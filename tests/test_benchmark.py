"""Check benchmark normalization and agreement on deliberately discordant data."""

import importlib.util
from pathlib import Path


spec = importlib.util.spec_from_file_location(
    "giftag_benchmark", Path(__file__).resolve().parents[1] / "benchmarks" / "compare.py")
bench = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bench)


def test_reference_normalizers(tmp_path):
    markers = {("KO", "K00001"), ("PFAM", "PF00001"), ("NCBIFAM", "NF000001.2"),
               ("NCBIFAM", "TIGR00001.3"),
               ("TIGRFAM", "TIGR00001"), ("CAZY", "GH43"), ("CAZY", "GH43_e2")}
    mapper = tmp_path / "mapper.tsv"
    mapper.write_text("g1\tK00001\ng2\t\ng3\tK99999\n")
    assert bench.normalize(mapper, "kofam-mapper", "genome", markers) == {
        ("genome", "g1", "KO", "K00001")}

    interpro = tmp_path / "interpro.tsv"
    interpro.write_text("g1\tmd5\t100\tPfam\tPF00001.7\tname\n"
                        "g2\tmd5\t100\tNCBIFAM\tNF000001.2\tname\n"
                        "g3\tmd5\t100\tNCBIFAM\tTIGR00001.3\tname\n")
    assert bench.normalize(interpro, "interproscan", "genome", markers) == {
        ("genome", "g1", "PFAM", "PF00001"),
        ("genome", "g2", "NCBIFAM", "NF000001.2"),
        ("genome", "g3", "NCBIFAM", "TIGR00001.3"),
        ("genome", "g3", "TIGRFAM", "TIGR00001")}

    family = tmp_path / "family.tsv"
    family.write_text("HMM Name\tTarget Name\nGH43_18.hmm\tg4\n")
    assert bench.normalize(family, "dbcan-family", "genome", markers) == {
        ("genome", "g4", "CAZY", "GH43")}
    sub = tmp_path / "sub.tsv"
    sub.write_text("Subfam Name\tTarget Name\nGH43_e2\tg5\n")
    assert bench.normalize(sub, "dbcan-sub", "genome", markers) == {
        ("genome", "g5", "CAZY", "GH43_e2")}

    domtblout = tmp_path / "hits.domtblout"
    domtblout.write_text("# target accession len query accession len ...\n"
                         "g1 - 100 Pfam_one PF00001.7 100 1e-30 100 0 1 1 "
                         "1e-30 1e-30 100 0 1 100 1 100 1 100 0.9 description\n")
    assert bench.normalize(domtblout, "hmmsearch-pfam", "genome", markers) == {
        ("genome", "g1", "PFAM", "PF00001")}


def test_comparison_counts_only_searchable_source_pairs(tmp_path):
    markers = {("KO", "K00001")}
    a = tmp_path / "giftag.tsv"
    b = tmp_path / "reference.tsv"
    a.write_text("genome_id\tgene_id\tnamespace\taccession\n"
                 "g\tp1\tKO\tK00001\n"
                 "g\tp2\tKO\tK00001\n"
                 "g\tp3\tCAZY\tGH43\n")
    b.write_text("genome_id\tgene_id\tnamespace\taccession\n"
                 "g\tp1\tKO\tK00001\n"
                 "g\tp4\tKO\tK00001\n")
    summary, differences, genomes, per_marker = bench.compare(a, b, markers, "kofam")
    assert (summary["shared_pairs"], summary["giftag_only_pairs"],
            summary["reference_only_pairs"], summary["jaccard"]) == (1, 1, 1, 1 / 3)
    assert {d["status"] for d in differences} == {"giftag_only", "reference_only"}
    assert genomes[0]["shared"] == 1
    assert per_marker[0]["accession"] == "K00001"
