import csv
import json

import pytest

from giftag import GiftagError
from giftag.annotate import annotate
from giftag.build import build
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
    assert dbcan["z_family"] == 2 and dbcan["z_sub"] == 3
    # Competitor clusters of a requested family are kept; other families are not.
    names = [name for name, _ in iter_hmm_records(open(db / "dbcan_sub.hmm", "rb"))]
    assert names == ["GH5_e1.hmm|GH5:12|3.2.1.4:5", "GH5_e2.hmm|GH5:7"]
    assert manifest["gifter"]["markers"] == 15


def test_existing_database_needs_force(world, db):
    with pytest.raises(GiftagError, match="--force"):
        build(db, world.markers, kofam_dir=world.kofam)


def test_truncated_hmm_file_is_an_error(world):
    data = (world.dbcan / "dbCAN.hmm").read_bytes()
    with pytest.raises(GiftagError, match="truncated"):
        list(iter_hmm_records(iter(data[: len(data) // 3].splitlines(keepends=True))))


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


def test_subfamily_needs_its_family(world, db):
    database = Database(db)
    block = database.digitize([("p", world.proteins["p_gh5_e1"])])
    gh9_only = [h for h in database._load("dbcan") if h.name == "GH9.hmm"]
    calls = search_dbcan(gh9_only, database._sub_by_family(), block, 0, 2, 3)
    assert calls == []


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
