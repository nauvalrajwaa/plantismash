# recruitment_miner Reference Data Directory

This directory stores the reference database for `recruitment_miner` (Feature B: target-guided mining / primary metabolism paralog recruitment).

## Expected Database Files

When fully populated by `bash_scripts/build_recruitment_db.py`, this folder contains:
1. `essential_proteins.faa` — FASTA of essential and curated recruited enzyme protein sequences (from Arabidopsis / TAIR10).
2. `reference_metadata.tsv` — TSV file detailing locus metadata:
   - `agi_locus`: AGI identifier (e.g. `AT1G01010`)
   - `source`: Evidence source (`EMB`, `OGEE`, or `CURATED_FAMILY`)
   - `family`: Enzyme / gene family annotation (e.g. `OSC`, `CPS`, `KSL`, `HMGR`, `SPDS`, `PMT`, `SAMS`, `ACCase-CT`, `tubulin`, `EMB_essential`)
   - `evidence_note`: Description or phenotypic evidence
3. `essential_proteins.dmnd` — DIAMOND database binary built with `diamond makedb --in essential_proteins.faa -d essential_proteins.dmnd`.

## Building / Updating the Database

Run:
```bash
python bash_scripts/build_recruitment_db.py \
  --tair-pep /path/to/TAIR10_pep_20101214.fasta \
  --emb-loci /path/to/SeedGenes_Ath_Mutants.txt \
  --ogee /path/to/ogee_athaliana_essential.tsv \
  --families /path/to/curated_families.tsv \
  --out-dir antismash/generic_modules/recruitment_miner/data/
```

### Input File Formats:
- `--emb-loci`: Plain text with one AGI locus per line (e.g. `AT1G01040`).
- `--ogee`: TSV export from OGEE containing locus and essentiality status.
- `--families`: TSV formatted as `family<TAB>agi_id[;agi_id...]` or `agi_id<TAB>family`.
- `--tair-pep`: Standard TAIR10 peptide FASTA.
