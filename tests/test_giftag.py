import csv
import json

import pytest

from giftag import GiftagError
from giftag.annotate import annotate
from giftag.build import _fetch_sub_families, build, dbcan_cluster_key, dbcan_sub_layout, natural_key
from giftag.cli import main
from giftag.database import Database
from giftag.fetch import iter_hmm_records
from giftag.genes import genome_id
from giftag.search import Call, filter_overlaps, search_dbcan


def read_tsv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def observed(rows):
    return {(r["gene_id"], r["namespace"], r["accession"]) for r in rows}


@pytest.fixture(scope="module")
def protein_run(world, db, tmp_path_factory):
    out = tmp_path_factory.mktemp("run")
    fasta = world.protein_fasta(out / "genomeA.faa")
    annotate([fasta], out, db)
    return out, read_tsv(out / "giftag_markers.tsv")


# -- build ---------------------------------------------------------------------

def test_every_marker_gets_a_status(db):
    status = {(m["namespace"], m["accession"]): (m["status"], m["detail"])
              for m in read_tsv(db / "markers.tsv")}
    assert status[("KO", "K00001")][0] == "searchable"
    assert status[("KO", "K00003")][0] == "no_threshold"
    assert status[("KO", "K99999")][0] == "missing_from_source"
    assert status[("EC", "1.1.1.1")][0] == "unsupported_namespace"
    assert status[("CAZY", "GH77")][0] == "missing_from_source"
    # An NCBIfam version is part of the accession: a different build of the
    # profile is not silently substituted.
    assert status[("NCBIFAM", "NF000002.1")] == ("missing_from_source", "release has NF000002.2, not NF000002.1")


def test_unversioned_tigrfam_maps_to_the_release_version(db):
    rows = {(p["namespace"], p["accession"]): p for p in read_tsv(db / "profiles.tsv")}
    assert rows[("TIGRFAM", "TIGR00001")]["profile"] == "TIGR00001.3"
    assert rows[("NCBIFAM", "TIGR00001.3")]["profile"] == "TIGR00001.3"


def test_manifest_records_dbcan_z_from_the_full_libraries(db):
    manifest = json.loads((db / "giftag.json").read_text())
    dbcan = manifest["sources"]["dbcan"]
    # Z counts every upstream profile, including those not kept.
    assert dbcan["z_family"] == 3 and dbcan["z_sub"] == 4
    # Every cluster of a requested cluster's family is kept, requested or not.
    names = [name.split(".hmm")[0] for name, _ in iter_hmm_records(open(db / "dbcan_sub.hmm", "rb"))]
    assert names == ["GH5_e1", "GH5_e2", "GH9_e1", "GH9_e2"]
    assert manifest["gifter"]["markers"] == 17


def test_existing_database_needs_force(world, db):
    with pytest.raises(GiftagError, match="--force"):
        build(db, world.markers, kofam_dir=world.kofam)


def test_truncated_hmm_file_is_an_error(world):
    data = (world.dbcan / "dbCAN.hmm").read_bytes()
    with pytest.raises(GiftagError, match="truncated"):
        list(iter_hmm_records(iter(data[: len(data) - 1000].splitlines(keepends=True))))


# -- acceptance rules ----------------------------------------------------------

def test_kofam_adaptive_thresholds(protein_run):
    found = observed(protein_run[1])
    assert ("p_k1", "KO", "K00001") in found
    assert ("p_k2", "KO", "K00002") not in found          # threshold out of reach
    assert ("p_k4", "KO", "K00004") in found              # best-domain score
    assert not any(acc == "K00003" for _, _, acc in found)  # no threshold, never assigned


def test_trusted_and_gathering_cutoffs(protein_run):
    found = observed(protein_run[1])
    assert ("p_nf1", "NCBIFAM", "NF000001.2") in found
    assert ("p_nf3", "NCBIFAM", "NF000003.1") not in found   # below trusted cutoff
    assert ("p_tigr", "TIGRFAM", "TIGR00001") in found        # one hit, both labels
    assert ("p_tigr", "NCBIFAM", "TIGR00001.3") in found
    assert ("p_pf", "PFAM", "PF00001") in found


def test_dbcan_family_and_subfamily(protein_run):
    found = observed(protein_run[1])
    assert ("p_gh5_e1", "CAZY", "GH5") in found
    assert ("p_gh5_e1", "CAZY", "GH5_e1") in found


def test_a_stronger_unrequested_cluster_wins_its_region(protein_run):
    # p_gh5_e2 also matches GH5_e1, but GH5_e2 owns the region in run_dbcan's
    # overlap resolution; reporting GH5_e1 would be a false subfamily call.
    found = observed(protein_run[1])
    assert ("p_gh5_e2", "CAZY", "GH5") in found
    assert ("p_gh5_e2", "CAZY", "GH5_e1") not in found


def test_subfamily_gate_is_opt_in(world, db):
    database = Database(db)
    block = database.digitize([("p", world.proteins["p_gh5_e1"])])
    gh9_only = [h for h in database._load("dbcan") if h.name == "GH9.hmm"]
    # Without a GH5 family hit: run_dbcan still reports the cluster, and so
    # does giftag by default; the gate drops it.
    default = search_dbcan(gh9_only, database._sub_by_family(), block, 0, 3, 4)
    assert [c.profile for c in default] == ["GH5_e1"]
    assert search_dbcan(gh9_only, database._sub_by_family(), block, 0, 3, 4, gate=True) == []


def test_official_subfamily_model_counts_as_its_family(protein_run):
    # p_gh9_sub hits GH9_2 only. run_dbcan would report just "GH9_2"; the
    # protein is still a GH9, so the family marker and the GH9 clusters apply.
    rows = {(r["gene_id"], r["accession"]): r for r in protein_run[1]}
    assert rows[("p_gh9_sub", "GH9")]["profile"] == "GH9_2"
    assert ("p_gh9_sub", "GH9_e2") in rows


def test_unrelated_protein_has_no_markers(protein_run):
    assert not any(r["gene_id"] == "p_random" for r in protein_run[1])


def test_rows_carry_their_evidence(protein_run):
    row = next(r for r in protein_run[1] if r["accession"] == "K00001")
    assert row["genome_id"] == "genomeA"
    assert row["source"] == "kofam" and row["rule"] == "kofam_full"
    assert float(row["score"]) >= float(row["threshold"]) == 50.0


def test_filter_overlaps_matches_run_dbcan():
    def call(profile, start, end, evalue):
        return Call("g", "dbcan", profile, "dbcan_family", 0, evalue, 1e-15, 0.9, start, end)
    kept = filter_overlaps([call("A", 1, 100, 1e-20), call("B", 10, 110, 1e-30),
                            call("C", 200, 300, 1e-20), call("D", 290, 400, 1e-40)])
    # B replaces A (90% overlap, better E-value); C and D overlap by 10%, so both stay.
    assert [c.profile for c in kept] == ["B", "C", "D"]


# -- inputs and outputs ----------------------------------------------------------

def test_nucleotide_genome_is_gene_called(world, db, tmp_path):
    fasta = world.genome_fasta(tmp_path / "mag_001.fna", world.proteins["p_k1"])
    annotate([fasta], tmp_path / "out", db)
    rows = read_tsv(tmp_path / "out" / "giftag_markers.tsv")
    assert any(r["accession"] == "K00001" and r["gene_id"].startswith("contig1_") for r in rows)
    genome = read_tsv(tmp_path / "out" / "giftag_genomes.tsv")[0]
    assert genome["genome_id"] == "mag_001" and genome["gene_calling"] == "pyrodigal meta"
    assert (tmp_path / "out" / "proteins" / "mag_001.faa").exists()


def test_duplicate_genome_ids_are_refused(world, db, tmp_path):
    a = world.protein_fasta(tmp_path / "x.faa")
    (tmp_path / "sub").mkdir()
    b = world.protein_fasta(tmp_path / "sub" / "x.faa.gz".replace(".gz", ""))
    with pytest.raises(GiftagError, match="share genome IDs"):
        annotate([a, b], tmp_path / "out", db)


def test_genome_id_strips_fasta_and_gzip_suffixes():
    assert genome_id("dir/GCA_000001.1.fa.gz") == "GCA_000001.1"
    assert genome_id("bin.12.fna") == "bin.12"
    assert genome_id("proteins.faa") == "proteins"


def test_bare_input_flag_runs_annotate(world, db, tmp_path):
    fasta = world.protein_fasta(tmp_path / "g.faa", ["p_k1"])
    assert main(["-i", str(fasta), "-o", str(tmp_path / "out"), "-d", str(db)]) == 0
    assert observed(read_tsv(tmp_path / "out" / "giftag_markers.tsv")) == {("p_k1", "KO", "K00001")}


def test_missing_database_is_a_clean_error(tmp_path, capsys):
    assert main(["info", "-d", str(tmp_path / "nothing")]) == 1
    assert "giftag build" in capsys.readouterr().err


# -- dbCAN-sub by byte range -------------------------------------------------------

def test_natural_key_follows_dbcan_sub_order():
    names = ["AA1_e0", "AA3_e25", "AA3_e186", "AA11_e65", "CBM6_e70", "GH1_e3", "GH5_e12", "GH10_e0"]
    assert sorted(names, key=natural_key) == names


def test_ranged_fetch_returns_exactly_one_familys_profiles(world):
    from giftag.fetch import RangeFile
    records, _ = _fetch_sub_families(RangeFile(world.dbcan / "dbCAN_sub.hmm"), {"GH5"},
                                     {"GH5": 2, "GH9": 2})
    assert [name for name, _ in records] == ["GH5_e1.hmm|GH5:12|3.2.1.4:5", "GH5_e2.hmm|GH5:7"]
    streamed = (world.dbcan / "dbCAN_sub.hmm").read_bytes()
    assert b"".join(record for _, record in records) == streamed[: len(b"".join(r for _, r in records))]


def test_ranged_fetch_refuses_a_block_that_disagrees_with_the_layout(world):
    from giftag.fetch import RangeFile
    with pytest.raises(GiftagError, match="3 were expected"):
        _fetch_sub_families(RangeFile(world.dbcan / "dbCAN_sub.hmm"), {"GH5"}, {"GH5": 3, "GH9": 2})


def test_pinned_dbcan_release_has_a_layout():
    from giftag import sources
    layout = dbcan_sub_layout(sources.DBCAN_DEFAULT_RELEASE)
    assert sum(layout.values()) == 53411 and layout["GH13"] > 0
    assert dbcan_sub_layout("no_such_release") is None


def test_ranged_fetch_finds_every_family_in_a_large_library(world, tmp_path, monkeypatch):
    import random
    from giftag import build as B
    from giftag.fetch import RangeFile
    from conftest import make_hmm
    rng = random.Random(3)
    families = ["AA1", "AA3", "AA11", "CBM6", "GH1", "GH5", "GH10", "GH13", "PL9"]
    layout, chunks, expected = {}, [], {}
    template = make_hmm(rng, "X", "M" + "A" * 40)
    for family in sorted(families, key=natural_key):
        n = rng.randint(1, 12)
        layout[family] = n
        names = sorted((f"{family}_e{i}" for i in range(n)), key=natural_key)
        for name in names:
            # Profiles of uneven size, as in the real library.
            record = template.replace(b"NAME  X", f"NAME  {name}.hmm|{family}:3".encode())
            record = record.replace(b"\nLENG", b"\nDESC  " + b"x" * rng.randint(0, 3000) + b"\nLENG", 1)
            chunks.append(record)
            expected.setdefault(family, []).append(name)
    path = tmp_path / "dbCAN_sub.hmm"
    path.write_bytes(b"".join(chunks))
    monkeypatch.setattr(B, "_PROBE_SPAN", 512)
    monkeypatch.setattr(B, "_SCAN_SPAN", 2048)
    monkeypatch.setattr(B, "_COARSE_PROBES", 4)
    wanted = {"AA3", "GH1", "GH10", "PL9"}
    records, _ = _fetch_sub_families(RangeFile(path), wanted, layout)
    got = [dbcan_cluster_key(name) for name, _ in records]
    assert got == [n for f in sorted(wanted, key=natural_key) for n in expected[f]]


# -- command line --------------------------------------------------------------------

def test_normalize_rewrites_the_shorthand_and_repeated_inputs():
    from giftag.cli import normalize
    assert normalize(["-i", "a.fa", "b.fa", "-o", "out"]) == [
        "annotate", "-i", "a.fa", "-i", "b.fa", "-o", "out"]
    # Logging options stay in front of the command they precede.
    assert normalize(["-q", "--log-file", "run.log", "-i", "a.fa", "-o", "out"]) == [
        "-q", "--log-file", "run.log", "annotate", "-i", "a.fa", "-o", "out"]
    assert normalize(["build", "--force"]) == ["build", "--force"]
    assert normalize(["annotate", "--input", "a.fa", "b.fa"]) == [
        "annotate", "--input", "a.fa", "--input", "b.fa"]


def test_directory_input_and_log_file(world, db, tmp_path, capsys):
    genomes = tmp_path / "genomes"
    genomes.mkdir()
    world.protein_fasta(genomes / "a.faa", ["p_k1"])
    world.protein_fasta(genomes / "b.faa", ["p_pf"])
    (genomes / "notes.txt").write_text("not a genome")
    log = tmp_path / "logs" / "run.log"
    code = main(["--log-file", str(log), "-i", str(genomes), "-o", str(tmp_path / "out"), "-d", str(db)])
    assert code == 0
    ids = [r["genome_id"] for r in read_tsv(tmp_path / "out" / "giftag_genomes.tsv")]
    assert ids == ["a", "b"]
    text = log.read_text()
    # The file log is plain: UTC stamps, no Rich markup.
    assert "INFO    1/2 a · 1 proteins · 1 markers" in text and "[bold]" not in text
    assert "2 genomes" in capsys.readouterr().err


def test_quiet_run_prints_nothing_but_still_writes(world, db, tmp_path, capsys):
    fasta = world.protein_fasta(tmp_path / "g.faa", ["p_k1"])
    # The fixture database has unsearchable markers, so the one warning stays.
    assert main(["-q", "-i", str(fasta), "-o", str(tmp_path / "out"), "-d", str(db)]) == 0
    err = capsys.readouterr().err
    assert "cannot be searched" in err and "marker rows" not in err and "g ·" not in err
    assert (tmp_path / "out" / "giftag_run.json").exists()


def test_help_and_usage_errors(capsys):
    assert main([]) == 0
    assert "giftag build" in capsys.readouterr().out
    assert main(["annotate", "--help"]) == 0
    assert "starting" not in capsys.readouterr().err
    assert main(["annotate", "-o", "x"]) == 2
    assert main(["version"]) == 0
