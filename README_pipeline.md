# plants_extended — plant BGC prediction pipeline

Extended fork of [plantiSMASH 2.0](https://github.com/plantismash/plantismash) (itself an
antiSMASH derivative) that detects biosynthetic gene clusters in **plant genomes** using a
combined HMM library:

| Source | Version | Profiles | Role |
|---|---|---|---|
| plantiSMASH curated | bundled | 93 | core + tailoring enzymes (terpene synthases, P450s, MTs, UGTs, STS/CHS, …) |
| dbCAN CAZyme HMMdb | `db_v5-2-9_5-5-2026` (S3: `https://dbcan.s3.us-west-2.amazonaws.com/db_v5-2-9_5-5-2026/dbCAN.hmm`) | 7 selected families | glycosylation / activation tailoring signatures |
| Pfam-A | release 37.0 (`https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz`) | 3 selected domains | conjugation tailoring (sulfotransferase, GST N/C) |

Detection rules, physical-distance clustering logic and the HMM engine are **stock
plantiSMASH** — the extension is entirely data-driven via a new model directory.

---

## What was added

### 1. New detection model: `plants_extended/`

Location: `antismash/generic_modules/hmm_detection/plants_extended/`

Model directories are auto-discovered — any subdirectory of `hmm_detection/` that contains
`hmmdetails.txt` + `cluster_rules.txt` registers as a selectable detection model. No engine
code changes were needed.

Contents:
- `*.hmm` — 106 profiles (93 copied from `plants/` + 13 new: 7 dbCAN, 6 Pfam, see below)
- `hmmdetails.txt` — profile id → description / cutoff / filename (4 tab-separated columns;
  `-1` = use the HMM's own trusted cutoffs; uncalibrated dbCAN models fall back to bit-score ≥ 20)
- `cluster_rules.txt` — 14 cluster rules (12 baseline + 2 new)
- `filterhmmdetails.txt` — overlap groups so redundant hits on the same gene collapse to the best one

New profiles added (prefixed by source to keep namespaces clean):

```
dbcan_GT1                    UDP-glycosyltransferase family (CAZy GT1)
dbcan_GT2_cellulose_synt     Cellulose synthase / GT2-family glycosyltransferase
dbcan_GH1                    Beta-glucosidase (CAZy GH1)
dbcan_GH3                    Glycoside hydrolase family 3
dbcan_AA3                    GMC-family oxidoreductase (CAZy AA3)
dbcan_AA7                    Glucooligosaccharide oxidase (CAZy AA7)
dbcan_CE1                    Acetylxylan esterase (CAZy CE1)
pfam_Sulfotransferase        PF00685
pfam_GST_C                   PF00043
pfam_GST_N                   PF02798
```

New cluster rules (baselines untouched for comparability with stock plantiSMASH):

| Rule | Logic | Rationale |
|---|---|---|
| `glycoside` | `minimum(4, [general pool + GT/GH], [GT-signatures, GH1/GH3-signatures])`, 5 kb cutoff / 1 kb ext | Requires **both** a glycosyltransferase-class gene *and* a β-glucosidase-class gene nearby — the classic glycosylated-metabolite pattern (saponins, cyanogenic/benzoxazinoid glucosides). Avoids false positives from ubiquitous single UGTs. |
| `cazyme_rich` | `minimum(4, [pool + all dbCAN], [AA7/AA3/CE1/SulTra/GST])` | Catches CAZyme-enriched oxidative/conjugation tailoring regions lacking canonical cores. Experimental — calibrate before trusting. |

Cross-source overlap groups added to `filterhmmdetails.txt` (a gene hitting several
equivalent profiles keeps only the best hit):

- GT group extended: `+ dbcan_GT1, dbcan_GT2_cellulose_synt`
- `GMC_oxred_N/C + dbcan_AA3`
- `Abhydrolase_3, COesterase + dbcan_CE1`
- new: `Glyco_hydro_1 + dbcan_GH1`
- new: `pfam_GST_C + pfam_GST_N` (N- and C-terminal domains of one protein)

### 2. Rule DSL semantics (this fork)

Important — this fork's parser differs from upstream antiSMASH:

```
minimum(N, [required pool],[essential groups])
```

- **Required pool**: comma-separated profile list; at least **N distinct genes** must carry hits from this pool.
- **Essential groups**: comma-separated groups; **every group** must be satisfied by ≥1 gene in the neighbourhood. Inside a group, `/` separates alternatives.
- Rules are joined by ` or `. Bare profile names are complete rules (e.g. `cyclopeptide = BURP`).
- Clustering distance: genes within `cutoff` kb of a core gene join the cluster; boundaries extend by `extension` kb.

### 3. Upstream bug fixes applied

Two dead references inherited from stock `plants/` were corrected in **both** model dirs:

- `cluster_rules.txt`: saccharide rule referenced `Glycos_transf_28`; profile is `Glyco_transf_28`
- `filterhmmdetails.txt`: group listed `2OG-Fell_Oxy` (lowercase L); profile is `2OG-FeII_Oxy` (capital I)

---

## Installation

Requires conda/mamba. From the repository root:

```bash
# Linux / Intel macOS
mamba env create -f environment.yml

# Apple Silicon (recommended pins have no arm64 builds -> run under Rosetta)
CONDA_SUBDIR=osx-64 mamba env create -f environment.yml

conda activate plantismash

# Optional: MOODS motif-scoring backend for TFBS detection (C extension,
# fails on some modern toolchains; tfbs_finder degrades gracefully without it)
pip install -e .[tfbs]
```

Notes:
- The environment installs **this fork** editable (`pip install -e .`). Do not install the
  PyPI `plantismash` package into the same env — it would shadow this code.
- Python is pinned `<3.9` and biopython `==1.76` (the code uses `Bio.Alphabet`, removed in
  biopython ≥ 1.78). Do not bump these casually.
- Tool databases (Pfam for fullhmmer/smcogs etc.) still need the stock downloader:
  `python download_databases.py`

---

## Usage

Input: an annotated plant genome — GFF3 + genome FASTA (+ protein FASTA), or a GenBank file.

**Same CLI as stock plantiSMASH** — the extended model is selected with one extra flag:

```bash
# GFF3 input (extended model)
python run_antismash.py \
    --taxon plants \
    --enabled-detection-models plants_extended \
    --cpus 8 \
    genome.gff3 proteins.faa

# GenBank input, drop-in replacement for your usual plantiSMASH invocation
plantismash \
    --taxon plants \
    --enabled-detection-models plants_extended \
    --limit -1 \
    --verbose \
    --clusterblast \
    --clusterblastdir "${DB_DIR}" \
    --disable_subgroup \
    --min-hmm-coverage 0.35 \
    --cpus "${THREADS}" \
    --outputfolder "${OUT_DIR}" \
    "${INPUT_GBK}"
```

Notes:
- Omitting `--enabled-detection-models` gives stock plantiSMASH behaviour (model `plants`).
- `--min-hmm-coverage 0.35` is already the default; pass it explicitly for clarity, or set `0` to disable the new filter entirely.
- Cluster types in outputs are namespaced `plants_extended/<rule>`; baseline rules keep their usual names.
- All standard antiSMASH output modules apply (HTML report, GenBank, cluster JSON).

### Detection parameters (inherited knobs — all apply to `plants_extended` rules)

| Knob | Default under `--taxon plants` | Effect |
|---|---|---|
| Rule cutoff (kb) | 5 kb (all rules incl. new ones) | Max distance from a core gene for another core to join/merge a cluster |
| Rule extension (kb) | 1 kb (20 kb for cyclopeptide) | Final boundary padding after last core gene |
| `--dynamic-cutoff` | **ON** | Scales the kb window per core gene by local intergenic density (radius = 10 neighbouring genes): sparse regions get wider windows, dense ones tighter |
| `--cutoff-multiplier` | 1.0 | Global tuning multiplier on every rule's window |
| `--gene-num-cutoff N` | 0 (= only immediate neighbours qualify) | **Max non-BGC intervening genes**: merges cores separated by ≤ N genes even beyond the kb window; combined with kb via OR unless `-only` is set |
| `--gene-num-cutoff-only` | off | Replaces kb logic entirely with intervening-gene counting |
| CD-HIT tandem filter | **ON**, identity 0.5 | Collapses tandem repeat arrays so one array ≈ one hit |
| `--min-domain-number` | **2** (forced) | Minimum distinct pool-profile hits required before `minimum()` rules are even evaluated |

Annotation-poor assemblies benefit from e.g. `--gene-num-cutoff 5 --gene-num-cutoff-only`;
leave defaults to stay closest to published plantiSMASH behaviour.


### Rebuilding / extending the model directory

`bash_scripts/build_plants_extended.py` regenerates the non-baseline parts of
`plants_extended/` from raw sources (downloads dbCAN/Pfam into `/tmp/bgc_meta/` first):

```bash
curl -L -o /tmp/bgc_meta/dbCAN.hmm "https://dbcan.s3.us-west-2.amazonaws.com/db_v5-2-9_5-5-2026/dbCAN.hmm"
curl -L -o /tmp/bgc_meta/Pfam-A.hmm.gz "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
python bash_scripts/build_plants_extended.py antismash/generic_modules/hmm_detection/plants_extended
```

The script is idempotent (re-runs don't duplicate entries). To add more families: extend
`DBCAN_SELECT` / `PFAM_SELECT` dicts in the script, re-run, then update
`hmmdetails.txt` / `filterhmmdetails.txt` / `cluster_rules.txt` accordingly.

To add any HMM ad hoc: drop a `.hmm` file into the model dir — missing entries are
auto-appended to `hmmdetails.txt` on first run (DESC line becomes the description,
trusted cutoffs used).

---

## Validation status

- ✅ Model dir cross-consistency checked (all hmmdetails filenames exist, internal HMM
  NAMEs match ids, every rule/filter reference resolves) — zero dangling references.
- ✅ Package loads in the target env; all 14 `plants_extended/*` cluster types discovered
  and parsed by the rule engine.
- ⚠️ End-to-end run on a real genome **not yet executed** (skipped by request).
  Suggested validation target: *Arabidopsis thaliana* chr5 (~17.15 Mb) contains the
  well-characterised thalianol BGC (THAS/THAH/THAD) — a correct run should recover it as a
  terpene cluster under both `plants` and `plants_extended`.

## Caveats

- dbCAN cutoff handling: `dbcan_GT2_cellulose_synt` ships trusted cutoffs (TC) and
  `pfam_*` models ship GA lines → used automatically. The other 6 dbCAN models have no
  calibration lines, so an explicit bit-score **25** is set in `hmmdetails.txt` column 3
  (was: engine fallback of 20). Provisional — recalibrate on real plant hit distributions.
- **Alignment-coverage filter** (`--min-hmm-coverage`, default **0.35**): a hit whose
  HMM-coordinate span covers < 35% of the profile is discarded before rule evaluation —
  rejects partial/low-complexity domain hits that can score well against long models.
  Set `--min-hmm-coverage 0` to restore legacy (no-filter) behaviour.
- `cazyme_rich` is experimental; expect it to be noisy until calibrated on known genomes.
- Pfam/dbCAN subsets here are deliberately small and curated for plant specialized
  metabolism. Scanning the full Pfam-A (20k models) through this engine scales poorly
  because profiles are scanned per-gene-set.

### Primary-metabolism suppression

`pfam_CESA` (PF03552), `pfam_SUSY`/GT1-group (PF00534) and `pfam_PME` (PF01095) are
**suppressor profiles**: they appear in no rule, but sit in the same overlap group as the
glycosyltransferase labels. On genuine cellulose-synthase / sucrose-synthase / pectin-
esterase genes, the primary-metabolism hit outcompetes the CAZyme label, so those genes
carry no BGC type and cannot satisfy `cazyme_rich`/`glycoside` signatures or pool quotas.
Zero engine changes — pure Layer-3 best-hit-wins behaviour.
