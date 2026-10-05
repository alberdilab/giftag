"""Inspect one dbCAN-sub discrepancy using giftag's exact built HMMs and rules."""

import argparse
import json

import pyhmmer

from giftag.database import Database
from giftag.search import _dbcan_domains, filter_overlaps


def cluster_key(name):
    if isinstance(name, bytes):
        name = name.decode()
    return name.split("|", 1)[0].split(".hmm", 1)[0]


def sequence_at(fasta, gene):
    found = False
    parts = []
    with open(fasta) as handle:
        for line in handle:
            if line.startswith(">"):
                if found:
                    break
                found = line[1:].split()[0] == gene
            elif found:
                parts.append(line.strip())
    if not parts:
        raise ValueError(f"{gene} missing from {fasta}")
    return "".join(parts)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", required=True)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--gene", required=True)
    parser.add_argument("--cluster", required=True)
    args = parser.parse_args()
    db = Database(args.db)
    expected_sequence = sequence_at(args.fasta, args.gene)
    with pyhmmer.easel.SequenceFile(args.fasta, digital=True,
                                   alphabet=db.alphabet) as reader:
        all_sequences = reader.read_block()
    target = next((sequence for sequence in all_sequences
                   if sequence.name.decode() == args.gene), None)
    if target is None or len(target) != len(expected_sequence):
        raise ValueError(f"{args.gene}: unexpected FASTA sequence")
    block = pyhmmer.easel.DigitalSequenceBlock(db.alphabet, [target])
    family = args.cluster.rsplit("_e", 1)[0]
    group = [hmm for hmm in db._load("dbcan_sub")
             if cluster_key(hmm.name).rsplit("_e", 1)[0] == family]
    z = db.manifest["sources"]["dbcan"]["z_sub"]
    raw = _dbcan_domains("dbcan_sub", group, block, 1, z,
                         cluster_key, "dbcan_sub")
    kept = filter_overlaps(raw)
    single = next(hmm for hmm in group if cluster_key(hmm.name) == args.cluster)
    single_calls = _dbcan_domains("dbcan_sub", [single], block, 1, z,
                                  cluster_key, "dbcan_sub")
    def simplify(calls):
        return [dict(profile=c.profile, score=c.score, evalue=c.evalue,
                     coverage=c.coverage, target_from=c.target_from,
                     target_to=c.target_to) for c in calls]
    print(json.dumps(dict(gene=args.gene, cluster=args.cluster, sequence_length=len(block[0]),
                          family_models=len(group), z_sub=z,
                          single=simplify(single_calls),
                          family_raw=simplify(raw), family_kept=simplify(kept)), indent=2))


if __name__ == "__main__":
    main()
