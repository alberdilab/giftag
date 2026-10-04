"""Synthetic sources in each upstream's own layout, so the whole build and
annotate path runs offline against profiles whose answers are known."""

import io
import random
import tarfile

import pyhmmer
import pytest

from giftag.build import build

ABC = pyhmmer.easel.Alphabet.amino()
AA = "ACDEFGHIKLMNPQRSTVWY"
CODON = {
    "A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT", "G": "GGT", "H": "CAT",
    "I": "ATT", "K": "AAA", "L": "CTG", "M": "ATG", "N": "AAT", "P": "CCG", "Q": "CAG",
    "R": "CGT", "S": "TCT", "T": "ACC", "V": "GTT", "W": "TGG", "Y": "TAT",
}


def mutate(rng, seq, rate):
    return "".join(rng.choice(AA) if rng.random() < rate else c for c in seq)


def make_hmm(rng, name, seed, accession=None, trusted=None, gathering=None):
    members = [pyhmmer.easel.TextSequence(name=f"s{i}", sequence=mutate(rng, seed, 0.1))
               for i in range(8)]
    msa = pyhmmer.easel.TextMSA(name=name, sequences=members).digitize(ABC)
    builder = pyhmmer.plan7.Builder(ABC)
    hmm, _, _ = builder.build_msa(msa, pyhmmer.plan7.Background(ABC))
    hmm.name = name
    if accession:
        hmm.accession = accession
    if trusted:
        hmm.cutoffs.trusted = trusted
    if gathering:
        hmm.cutoffs.gathering = gathering
    buffer = io.BytesIO()
    hmm.write(buffer)
    return buffer.getvalue()


class World:
    """Seeds, proteins and source directories shared by the tests."""

    def __init__(self, root):
        rng = random.Random(7)
        seed = {k: "M" + "".join(rng.choice(AA) for _ in range(199)) for k in
                ("K00001", "K00002", "K00003", "K00004", "NF1", "NF3", "TIGR", "PF", "GH5", "GH9", "GH9s")}
        seed["GH5_e1"] = mutate(rng, seed["GH5"], 0.3)
        seed["GH5_e2"] = mutate(rng, seed["GH5"], 0.3)
        self.seed = seed
        self.proteins = {
            "p_k1": mutate(rng, seed["K00001"], 0.05),
            "p_k2": mutate(rng, seed["K00002"], 0.05),
            "p_k4": "".join(rng.choice(AA) for _ in range(150)) + mutate(rng, seed["K00004"], 0.05),
            "p_nf1": mutate(rng, seed["NF1"], 0.05),
            "p_nf3": mutate(rng, seed["NF3"], 0.05),
            "p_tigr": mutate(rng, seed["TIGR"], 0.05),
            "p_pf": mutate(rng, seed["PF"], 0.05),
            "p_gh5_e1": mutate(rng, seed["GH5_e1"], 0.05),
            "p_gh5_e2": mutate(rng, seed["GH5_e2"], 0.05),
            # Matches only the official subfamily model GH9_2, never bare GH9.
            "p_gh9_sub": mutate(rng, seed["GH9s"], 0.05),
            "p_random": "M" + "".join(rng.choice(AA) for _ in range(250)),
        }
        self.rng = rng
        self.root = root
        self.kofam = self._kofam(root / "kofam")
        self.ncbifam = self._ncbifam(root / "ncbifam")
        self.pfam = self._pfam(root / "pfam")
        self.dbcan = self._dbcan(root / "dbcan")
        self.markers = root / "markers.tsv"
        self.markers.write_text("namespace\taccession\tname\n" + "".join(
            f"{ns}\t{acc}\t\n" for ns, acc in [
                ("KO", "K00001"), ("KO", "K00002"), ("KO", "K00003"), ("KO", "K00004"),
                ("KO", "K99999"), ("NCBIFAM", "NF000001.2"), ("NCBIFAM", "NF000002.1"),
                ("NCBIFAM", "NF000003.1"), ("NCBIFAM", "TIGR00001.3"), ("TIGRFAM", "TIGR00001"),
                ("PFAM", "PF00001"), ("CAZY", "GH5"), ("CAZY", "GH5_e1"), ("CAZY", "GH77"),
                ("CAZY", "GH9"), ("CAZY", "GH9_e2"),
                ("EC", "1.1.1.1"),
            ]))

    def _kofam(self, d):
        d.mkdir()
        (d / "ko_list").write_text(
            "knum\tthreshold\tscore_type\tprofile_type\tF-measure\tnseq\tnseq_used\talen\tmlen\teff_nseq\tre/pos\tdefinition\n"
            "K00001\t50.0\tfull\tall\t1\t8\t8\t200\t200\t1\t1\tone\n"
            "K00002\t100000.0\tfull\tall\t1\t8\t8\t200\t200\t1\t1\tunreachable\n"
            "K00003\t-\t-\t-\t-\t8\t8\t200\t200\t1\t1\tno threshold\n"
            "K00004\t50.0\tdomain\tall\t1\t8\t8\t200\t200\t1\t1\tdomain scored\n")
        with tarfile.open(d / "profiles.tar.gz", "w:gz") as tar:
            for knum in ("K00001", "K00002", "K00003", "K00004"):
                data = make_hmm(self.rng, knum, self.seed[knum])
                info = tarfile.TarInfo(f"profiles/{knum}.hmm")
                info.size = len(data)
                tar.addfile(info, io.BytesIO(data))
        return d

    def _ncbifam(self, d):
        (d / "hmm_PGAP.HMM").mkdir(parents=True)
        (d / "hmm_PGAP.tsv").write_text(
            "#ncbi_accession\tsource_identifier\tlabel\n"
            "NF000001.2\t\tone\nNF000002.2\t\ttwo\nNF000003.1\t\tthree\n"
            "TIGR00001.2\t\told\nTIGR00001.3\t\ttigr\n")
        profiles = {"NF000001.2": ("NF1", (40.0, 40.0)), "NF000002.2": ("NF1", (40.0, 40.0)),
                    "NF000003.1": ("NF3", (5000.0, 5000.0)), "TIGR00001.3": ("TIGR", (40.0, 40.0))}
        for accession, (key, trusted) in profiles.items():
            (d / "hmm_PGAP.HMM" / f"{accession}.HMM").write_bytes(
                make_hmm(self.rng, f"n_{accession}", self.seed[key], accession, trusted=trusted))
        return d

    def _pfam(self, d):
        d.mkdir()
        (d / "PF00001.hmm").write_bytes(
            make_hmm(self.rng, "Pfam_one", self.seed["PF"], "PF00001.7", gathering=(30.0, 30.0)))
        return d

    def _dbcan(self, d):
        d.mkdir()
        (d / "dbCAN.hmm").write_bytes(
            make_hmm(self.rng, "GH5.hmm", self.seed["GH5"]) + make_hmm(self.rng, "GH9.hmm", self.seed["GH9"])
            + make_hmm(self.rng, "GH9_2.hmm", self.seed["GH9s"]))
        (d / "dbCAN_sub.hmm").write_bytes(
            make_hmm(self.rng, "GH5_e1.hmm|GH5:12|3.2.1.4:5", self.seed["GH5_e1"])
            + make_hmm(self.rng, "GH5_e2.hmm|GH5:7", self.seed["GH5_e2"])
            + make_hmm(self.rng, "GH9_e1.hmm|GH9:3", self.seed["GH9"])
            + make_hmm(self.rng, "GH9_e2.hmm|GH9:4", self.seed["GH9s"]))
        return d

    def protein_fasta(self, path, names=None):
        names = names or list(self.proteins)
        path.write_text("".join(f">{n}\n{self.proteins[n]}\n" for n in names))
        return path

    def genome_fasta(self, path, protein):
        """A contig carrying one ORF for `protein` in random flanking DNA."""
        rng = random.Random(11)
        flank = lambda n: "".join(rng.choice("ACGT") for _ in range(n))  # noqa: E731
        orf = "".join(CODON[a] for a in protein) + "TAA"
        path.write_text(f">contig1\n{flank(3000)}{orf}{flank(3000)}\n")
        return path


@pytest.fixture(scope="session")
def world(tmp_path_factory):
    return World(tmp_path_factory.mktemp("sources"))


@pytest.fixture(scope="session")
def db(world, tmp_path_factory):
    path = tmp_path_factory.mktemp("db") / "giftag"
    build(path, world.markers, kofam_dir=world.kofam, ncbifam_dir=world.ncbifam,
          pfam_dir=world.pfam, dbcan_dir=world.dbcan)
    return path
