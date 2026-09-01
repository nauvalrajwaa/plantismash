# recruitment_miner Reference Data Directory

This directory stores the reference database for `recruitment_miner` (Feature B: target-guided mining / primary metabolism paralog recruitment).

## Database Files

1. `essential_proteins.faa` — FASTA containing 519 curated essential and recruited enzyme protein sequences (from *Arabidopsis thaliana* / TAIR10 proteome).
2. `reference_metadata.tsv` — TSV file detailing locus metadata (522 unique loci):
   - `agi_locus`: AGI identifier (e.g. `AT1G01010`)
   - `source`: Evidence source (`EMB`, `OGEE`, or `CURATED_FAMILY`)
   - `family`: Enzyme / gene family annotation (`OSC`, `CPS`, `KSL`, `HMGR`, `SPDS/PMT`, `SAMS`, `ACCase-CT`, `tubulin`, `EMB_essential`)
   - `evidence_note`: Description or phenotypic evidence
3. `essential_proteins.dmnd` — DIAMOND database binary built with `diamond makedb --in essential_proteins.faa -d essential_proteins.dmnd`.

## Automated Build / Update

To rebuild or refresh the database automatically:
```bash
python bash_scripts/build_recruitment_db.py
```
This will automatically fetch TAIR10 proteome sequences from Ensembl Plants and essential loci from SeedGenes, combine with curated enzyme families, generate the FASTA/metadata files, and compile the DIAMOND index.

### Custom / Offline Build:
```bash
python bash_scripts/build_recruitment_db.py \
  --tair-pep /path/to/TAIR10_pep_20101214.fasta \
  --emb-loci /path/to/SeedGenes_Ath_Mutants.txt \
  --ogee /path/to/ogee_athaliana_essential.tsv \
  --families /path/to/curated_families.tsv \
  --out-dir antismash/generic_modules/recruitment_miner/data/
```

