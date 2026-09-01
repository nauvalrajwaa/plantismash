# plantiSMASH `plants_extended` — changes (W1–W8)

Branch `feature/plant-bgc-extended` adds 15 HMM profiles + 2 rules on top of stock `plants` (108 profiles / 15 rules). This doc records the W1–W8 engine deltas, flags, and provenance.

## Suppressors + overlap margin (W4)

- GT overlap group (line 1 of `plants_extended/filterhmmdetails.txt`) contains the 5 suppressor profiles: `pfam_CESA`, `pfam_SUSY`, `pfam_PME`, `dbcan_GT43`, `dbcan_GT47`. They act by overlap-group best-hit-wins; no rule references them.
- New `plants_extended/suppressors.tsv` (`profile  margin  enabled`): all 5 at `margin=0.15`, `enabled=yes`.
- Rule: when the group winner is a suppressor and the best non-suppressor also overlaps (>20 aa), the suppressor only wins if `winner_bits >= (1+margin)*runner_up_bits`. Otherwise the detection hit survives and the suppressor hit is dropped. Every decision is logged (W3 `overlap_decisions.tsv`).

## Coverage / bitscore rules (W1–W2)

- `--min-hmm-coverage` default `0.25` (was `0.35`); logged at startup as `HMM filters …`.
- Coverage filtering applies only to HMMs with model length ≥ 300 aa (shorter models skip it unless they have an explicit `coverage_override` in `hit_thresholds.tsv`). Model length is taken from the HSP/SearchIO `seq_len` first; else parsed from the `.hmm` file `LENG` line at load.
- `plants_extended/hit_thresholds.tsv` (`profile  min_bits  coverage_override  note`, blank override = global): per-profile bitscore floors. A hit below BOTH its `hmmdetails.txt` cutoff and its `hit_thresholds.tsv` floor (when listed) is dropped before rule pools. Provisional floors for the 15 new profiles were calibrated from sec_met bitscore distributions in the 4 GBK run dirs (`PA_default` 94, `test25_PA` 86, `run1_28clusters` 28, `run2_24clusters_test3` 24) as `round(p10)` clamped `25–60`; suppressors note `suppressor`. Use `evaluation/calibrate_cutoffs.py` to reproduce.
- `--min-hmm-bitscore` (default `35`, `0` disables): hits below it never count toward rule quorums.

## Flags

| flag | default | meaning |
|---|---|---|
| `--min-hmm-coverage` | `0.25` | global coverage floor |
| `--min-hmm-bitscore` | `35` | global quorum floor |
| `--attribution-log` | `on` | write `<output>/attribution/` |

## New files

- `plants_extended/hit_thresholds.tsv`, `plants_extended/suppressors.tsv`
- `generic_modules/hmm_detection/attribution.py` + wiring in `run_antismash.py` and `hmm_detection/__init__.py`
- `evaluation/calibrate_cutoffs.py` (stdlib only)
- `plants_extended/cluster_rules.txt`: `cazyme_rich` window raised `3/0.5 → 10/2` kb (least-invasive path per W5 option 2, so dynamic multiplier need not change)

## `cazyme_rich` change (W5)

- Window `10 kb / 2 kb` in `plants_extended/cluster_rules.txt`. Combined with dynamic-cutoff multiplier already in `apply_cluster_rules`, this lets dbcan-anchored clusters (e.g. melon r2c001, Euphorbia c002) reach quorum without inflating genome-wide FPs. Marginal-quorum flag (below) catches the edge cases.
- Marginal quorum: rule satisfied with exactly `N` distinct genes and ALL counted hits `< 50` bits → appends `Marginal quorum: <rule>` to the cluster `/note` and sets `marginal_quorum` in `cluster_summary.tsv`. Surfaced in the attribution layer.

## Attribution & observability (W3)

Writes into `<output>/attribution/` when `--attribution-log on` (default):

- `hits.tsv`: `contig  gene_id  profile  bits  evalue  aln_len  model_len  coverage  decision  reason` (`kept` / `coverage_fail` / `bits_fail` / `overlap_loss` / `suppressor_margin_fail`)
- `overlap_decisions.tsv`: every overlap-group resolution
- `cluster_summary.tsv`: per-cluster type, gene counts, marginal flag, `gene:profile:bits` detail

`sec_met` `Domains detected` line now carries `, aln/model=<aln>/<model>, cov=<cov>` AFTER the existing `profile (E-value ..., bitscore ..., seeds ...)` prefix so `profile (E-value` regex still matches byte-for-byte.

## Provenance & combo mode

Run provenance-aware analyses with `--models "plants,plants_extended"` to keep baseline and extended label pools separated while sharing cutoffs/filters. The plant type is rewritten to `putative` for provenance in the GBK writer.

## Case studies

- Melon `28 → 24` (default → `plants_extended` at 0.25): 6 verified-fake clusters lost (c005 suppressed CESA, remainder coverage-failed); no ≥100-bit core lost.
- Euphorbia hirta `94 → 86 → 72` (default → 0.25 → 0.35): `0.35` destroys 8/14 clusters containing ≥100-bit cores (up to 414 bits), confirming `0.25` as the calibrated default. `86 → 72` pure-filter losses are the calibration set for W1.

## Model fidelity

No HMM profiles or rules were added/removed (model stays 108 / 15).
