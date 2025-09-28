# Metadata

See the file `full_metadata.csv`. 

# List of BED files and documentation on enrichment strategies

## `atacseq_bed`

These `bed` files were obtained from the publication [Integrating chromatin accessibility states in the design of targeted sequencing panels for liquid biopsy](https://www.nature.com/articles/s41598-022-14675-z). For short, we call them `bed` file number 7 and number 8 (accessible regions and inaccessible regions, respectively). For these regions, we first filter our `BAM` files; keep fragments which fall into the regions only and discard other. We then construct our set of fragmentomics features (`FLEN, EM, ND`) as usual. See the notebook `process_bed_files.ipynb` for pre-processing script generating `bed` file from the `.xls` files.

Save feature folder name: `highdepth_atacseq_20250815`

## `bin600`

Take the nucleosome map `rpr_map_EXP0779.sorted.bed` and expand each nucleosome position to 600bp (+-300bp up-and-downstream). 

Save feature folder name: `bin600`

## `panel_7_8_vs_TSS` 

Intersect the panel 7 and 8 from `atacseq_bed` above with the 100bp up-and-downstream flanking regions of TSS. TSS coordinates obtained from `biomart` database, accessed via `R` package. 

Save feature folder nane: `highdepth_atacseq_overlap_tss`

## `methylation_regions`

Save feature folder name: `highdepth_methylbed`

## `TCGA_atac_seq_beds`
Run the analysis on TCGA data from the repo https://github.com/hieunguyen4193/ATAC-seq-TCGA-data-analysis 

Save feature folder name: `highdepth_atacseq_TCGA`

## `panel_7_8_vs_nucleosomeMap`

Save feature folder name: `panel_7_8_vs_nucleosomeMap`

## `custom_genes_tss` 
Use scripts from the repo `k-mer-feature`: https://github.com/hieunguyen4193/k-mer-feature. Run function `generate_TSS_promoters_for_input_genes.R`.

Save feature folder name: `custom_genes_TSS`

# List of "feature counts" for TSS bed files
