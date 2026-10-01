# Transcriptional Response of Inositol Pyrophosphate Mutants Under Phosphate Starved Conditions

This repository contains an RNA-seq workflow for studying how inositol pyrophosphate pathway mutations affect the transcriptional response to phosphate availability in Candida albicans. The project compares wild-type (WT) cells with mutants lacking KCS1, VIP1, or both genes under phosphate-replete (+Pi) and phosphate-starved (−Pi) conditions.
The repository brings together the kcs1Δ and vip1Δ phosphate-response analyses. The kcs1Δvip1Δ double-mutant analysis will be added as the project progresses. Each trial is analyzed separately using the same general workflow, with its own sample assignments, statistical model, and output files.

## Phase 1: Upstream HPC preprocessing
The upstream workflow processes single-end RNA-seq reads on a high-performance computing cluster to produce gene-level count matrices:
- Quality control (FastQC): Checks read quality and adapter contamination. Trimming was skipped for the current datasets because the reads were high quality.
- Reference indexing and alignment (STAR): Builds an index for the C. albicans reference genome and aligns the reads.
- Alignment processing (SAMtools): Converts SAM files to BAM format and sorts alignments by genomic coordinates.
- Gene quantification (featureCounts): Counts reads assigned to genes and produces a count matrix for each trial.

## Phase 2: Downstream R analysis
The same general analysis is applied to each mutant trial:
- Count preparation and filtering: Imports the featureCounts matrix, sums rows with identical gene IDs, removes tRNA and rRNA features, and filters genes with low counts before modeling.
- Gene annotation: Matches systematic IDs to Candida Genome Database (CGD) gene names. Genes without a matching name retain their identifiers.
- Sample quality checks: Uses variance-stabilized counts to generate PCA plots, sample correlation heatmaps, and gene-expression heatmaps. Library-size plots summarize the counts retained for analysis.
- Differential expression (DESeq2): Fits the model ~ Genotype * Pi separately for each trial. Pairwise comparisons measure differences among the four conditions, while the genotype-by-phosphate interaction tests whether the response to phosphate starvation differs between WT and the mutant.
- Phosphate-response analysis: Compares −Pi versus +Pi within each genotype. Response plots, Venn diagrams, and gene tables describe shared and genotype-specific changes.
- Directional GO enrichment (clusterProfiler and GO.db): Analyzes genes increased or decreased during phosphate starvation separately. CGD annotations provide gene-to-GO associations, and GO.db supplies term names. Enrichment uses the analyzed genes that map to the relevant GO ontology as the background. Outputs summarize Biological Process, Molecular Function, and Cellular Component enrichment.

## Outputs
Each trial saves its results separately. Outputs include sample metadata, annotated differential expression tables, significant-gene lists, QC plots, volcano and MA plots, heatmaps, phosphate-response plots, Venn diagrams, gene tables, GO enrichment results, and saved analysis objects.
Keeping the trial outputs separate preserves the experimental design of each dataset while allowing the same analysis approach to be used across all three mutants. Comparisons across mutants should also account for differences between trials.
Refer to the individual R scripts for the complete downstream analyses. Upstream preprocessing is run separately before importing the featureCounts matrices into R.
