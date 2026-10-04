"""Search proteins against each source with that source's own acceptance rule.

The sources calibrate their thresholds differently, and giftag keeps each rule
as its source defines it rather than inventing a common score:

- KOfam: the per-KO adaptive bit-score threshold from `ko_list`, applied to the
  full-sequence score or to the best single domain as `score_type` says. This
  is KofamScan's assignment rule (its `*` mark).
- NCBIfam/TIGRFAM: the profile's trusted cutoffs (HMMER `--cut_tc`).
- Pfam: the profile's gathering cutoffs (HMMER `--cut_ga`).
- dbCAN and dbCAN-sub: run_dbcan's domain rule (independent E-value below
  1e-15 and HMM coverage above 0.35), with E-values computed at run_dbcan's
  fixed Z, the profile count of the *full* upstream library, then run_dbcan's
  overlap resolution.

dbCAN-sub is searched on every protein, as run_dbcan does. `--gate-subfamilies`
restricts each family's clusters to proteins that carry a domain of that family
(its family model or one of its official subfamily models, such as `GH43_18`).
That is faster and stricter: on Bacteroides thetaiotaomicron VPI-5482 it keeps
144 of run_dbcan's 157 subfamily calls.
"""

from collections import defaultdict
from dataclasses import dataclass
from typing import Optional

import pyhmmer

from giftag import sources as S
from giftag.build import cazy_family, dbcan_cluster_key, dbcan_family_key


@dataclass
class Call:
    gene_id: str
    source: str
    profile: str
    rule: str
    score: float
    evalue: float
    threshold: float
    coverage: Optional[float] = None
    target_from: Optional[int] = None
    target_to: Optional[int] = None


def search_kofam(hmms, block, rules, cpus):
    calls = []
    for hits in pyhmmer.hmmsearch(hmms, block, cpus=cpus, T=0.0, domT=0.0):
        knum = hits.query.name
        rule, threshold = rules[("kofam", knum)]
        for hit in hits:
            domain = hit.best_domain
            if rule == "kofam_domain":
                score, evalue = domain.score, domain.i_evalue
            else:
                score, evalue = hit.score, hit.evalue
            if score >= threshold:
                calls.append(Call(hit.name, "kofam", knum, rule, score, evalue, threshold,
                                  target_from=domain.alignment.target_from,
                                  target_to=domain.alignment.target_to))
    return calls


def search_cutoff(source, hmms, block, rules, cpus, cutoff):
    """NCBIfam (trusted) and Pfam (gathering) model-specific cutoffs."""
    calls = []
    for hits in pyhmmer.hmmsearch(hmms, block, cpus=cpus, bit_cutoffs=cutoff):
        accession = hits.query.accession
        rule, threshold = rules[(source, accession)]
        for hit in hits.included:
            domains = list(hit.domains.included)
            if not domains:
                continue
            best = max(domains, key=lambda d: d.score)
            calls.append(Call(hit.name, source, accession, rule, hit.score, hit.evalue, threshold,
                              target_from=best.alignment.target_from,
                              target_to=best.alignment.target_to))
    return calls


def _dbcan_domains(source, hmms, block, cpus, z, key, rule):
    calls = []
    for hits in pyhmmer.hmmsearch(hmms, block, cpus=cpus, Z=z, domE=S.DBCAN_EVALUE):
        for hit in hits:
            for domain in hit.domains.included:
                aln = domain.alignment
                coverage = (aln.hmm_to - aln.hmm_from + 1) / aln.hmm_length
                if domain.i_evalue < S.DBCAN_EVALUE and coverage > S.DBCAN_COVERAGE:
                    calls.append(Call(hit.name, source, key(aln.hmm_name), rule, domain.score,
                                      domain.i_evalue, S.DBCAN_EVALUE, coverage,
                                      aln.target_from, aln.target_to))
    return calls


def filter_overlaps(calls, ratio=S.DBCAN_OVERLAP):
    """run_dbcan's `filter_overlaps`: per protein, walk domains by position and,
    where two overlap by more than `ratio` of either one's length, keep the one
    with the lower independent E-value. Reproduced as written, including that
    each domain is compared only with the last one kept."""
    by_gene = defaultdict(list)
    for call in calls:
        by_gene[call.gene_id].append(call)
    kept = []
    for gene in sorted(by_gene):
        keep = []
        for call in sorted(by_gene[gene], key=lambda c: (c.target_from, c.target_to)):
            if not keep:
                keep.append(call)
                continue
            last = keep[-1]
            overlap = min(last.target_to, call.target_to) - max(last.target_from, call.target_from)
            if overlap > 0:
                r_last = overlap / max(1, last.target_to - last.target_from)
                r_call = overlap / max(1, call.target_to - call.target_from)
                if r_last > ratio or r_call > ratio:
                    if last.evalue > call.evalue:
                        keep[-1] = call
                    continue
            keep.append(call)
        kept.extend(keep)
    return kept


def search_dbcan(family_hmms, sub_by_family, block, cpus, z_family, z_sub, gate=False):
    """dbCAN family calls, then dbCAN-sub calls.

    By default dbCAN-sub is searched on every protein, as run_dbcan does. With
    `gate`, a family's clusters are searched only on proteins that carry a
    domain of that family: faster, and stricter than run_dbcan."""
    families = filter_overlaps(_dbcan_domains(
        "dbcan", family_hmms, block, cpus, z_family, dbcan_family_key, "dbcan_family"))
    if not sub_by_family:
        return families
    if not gate:
        hmms = [hmm for family in sorted(sub_by_family) for hmm in sub_by_family[family]]
        return families + filter_overlaps(_dbcan_domains(
            "dbcan_sub", hmms, block, cpus, z_sub, dbcan_cluster_key, "dbcan_sub"))
    genes_by_family = defaultdict(set)
    for call in families:
        genes_by_family[cazy_family(call.profile)].add(call.gene_id)
    sequences = {seq.name: seq for seq in block}
    raw = []
    for family in sorted(genes_by_family):
        hmms = sub_by_family.get(family)
        if not hmms:
            continue
        subset = pyhmmer.easel.DigitalSequenceBlock(
            block.alphabet, [sequences[g] for g in sorted(genes_by_family[family])])
        raw.extend(_dbcan_domains("dbcan_sub", hmms, subset, cpus, z_sub,
                                  dbcan_cluster_key, "dbcan_sub"))
    return families + filter_overlaps(raw)
