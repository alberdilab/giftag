# dbCAN-sub discordance adjudication

Across the eight-genome full-library `run_dbcan` comparison, the only
dbCAN-sub reference-only marker was `GH95_e26` on
`btheta_vpi5482_p004502`. The run_dbcan hit had independent E-value `0.0`
and HMM coverage `0.9863`. The same profile in giftag's built database hits
that exact protein with score `1023.16` and coverage `0.9863`; the profile is
present and searchable.

The saved re-searches used the **same 56 stored GH95 HMMs and same protein**
under two installed pyhmmer versions. With pyhmmer 0.11.0 (run_dbcan 5.2.9),
the E-values of both `GH95_e26` and the stronger `GH95_e1` hit underflow to
zero. The upstream overlap rule keeps the first of these tied overlapping
hits, `GH95_e26`. With pyhmmer 0.12.3 (giftag 0.2.0), their E-values are
`2.03e-307` and `4.07e-315`; the same overlap rule keeps `GH95_e1`.
The full raw evidence is in the adjacent JSON files. The reference-only
`GH95_e26` marker changes **zero complete GIFT calls** when substituted into
the pinned gifter 0.7.3 evaluation for this genome.

This is an implementation-version and numerical tie difference, not evidence
for biological ground truth in either tool. The final concordance and GIFT
tables retain the difference.
