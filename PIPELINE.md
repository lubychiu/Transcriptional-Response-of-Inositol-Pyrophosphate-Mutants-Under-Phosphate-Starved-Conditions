# Transcriptional Differences Between Inositol Pyrophosphate Pathway Mutants in *Candida albicans* Under Phosphate Starved Conditions
## Visualizations
Done in RStudio.

### WT vs KCS1 analysis
'''bash
############################################################
# CANDIDA ALBICANS RNA-seq ANALYSIS
#
# EXPERIMENT:
#
# WT    +Pi     3 biological replicates
# WT    -Pi     3 biological replicates
# KCS1  +Pi     3 biological replicates
# KCS1  -Pi     3 biological replicates
#
# 12 samples total
#
# SAMPLE ORDER AFTER SORTING:
#
# 1-3     WT +Pi
# 4-6     WT -Pi
# 7-9     KCS1 +Pi
# 10-12   KCS1 -Pi
#
# ANALYSIS:
#
# - Remove tRNA/rRNA features
# - Map systematic IDs to CGD gene names
# - DESeq2 using all 12 biological replicates
# - PCA
# - sample correlation
# - all 6 possible pairwise comparisons
# - volcano plots for all comparisons
# - MA plots for all comparisons
# - heatmaps for all comparisons
# - condition-level VST heatmap
# - Genotype x Pi interaction analysis
# - phosphate-response analysis
#
############################################################


############################################################
# 0. LOAD REQUIRED PACKAGES
############################################################

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

library(DESeq2)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(dplyr)
library(tidyr)
library(tibble)
library(matrixStats)


############################################################
# 1. OUTPUT DIRECTORIES
############################################################

output_dir <- "RNAseq_DE_Results"

dirs <- c(
  "01_QC",
  "02_DE_results",
  "03_Volcano",
  "04_MA_plots",
  "05_Heatmaps",
  "06_PCA",
  "07_Sample_correlations",
  "08_Significant_gene_lists",
  "09_Response_analysis"
)

dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

for (d in dirs) {
  dir.create(
    file.path(output_dir, d),
    showWarnings = FALSE,
    recursive = TRUE
  )
}


############################################################
# 2. LOAD CGD GENE ANNOTATION
############################################################

print("======================================================")
print("LOADING CGD GENE ANNOTATION")
print("======================================================")

if (!file.exists("cgd_features.tab")) {
  
  stop(
    paste0(
      "\nERROR: cgd_features.tab was not found.\n",
      "Place cgd_features.tab in the same working directory ",
      "as this R script.\n"
    )
  )
  
}

cgd_raw <- read.table(
  "cgd_features.tab",
  header = FALSE,
  sep = "\t",
  quote = "",
  comment.char = "!",
  fill = TRUE,
  stringsAsFactors = FALSE
)


############################################################
# CGD ANNOTATION STRUCTURE
#
# Column 1 = systematic ID
# Column 2 = gene name
# Column 4 = feature type
############################################################

mapping_table <- data.frame(
  
  systematic_id = trimws(
    as.character(cgd_raw[, 1])
  ),
  
  gene_symbol = trimws(
    as.character(cgd_raw[, 2])
  ),
  
  feature_type = trimws(
    as.character(cgd_raw[, 4])
  ),
  
  stringsAsFactors = FALSE
  
)


############################################################
# CLEAN GENE NAMES
############################################################

mapping_table$gene_symbol[
  is.na(mapping_table$gene_symbol) |
    mapping_table$gene_symbol == "" |
    mapping_table$gene_symbol ==
    mapping_table$systematic_id
] <- NA


############################################################
# NORMALIZE IDS
############################################################

normalize_id <- function(x) {
  
  x <- trimws(
    as.character(x)
  )
  
  x <- toupper(
    gsub(
      "[_\\-]",
      "",
      x
    )
  )
  
  return(x)
  
}


mapping_table$id_norm <- normalize_id(
  mapping_table$systematic_id
)


############################################################
# REMOVE DUPLICATE CGD ENTRIES
#
# Prefer entries that actually have gene names.
############################################################

mapping_table <- mapping_table[
  order(
    is.na(mapping_table$gene_symbol)
  ),
  ,
  drop = FALSE
]

mapping_table <- mapping_table[
  !duplicated(mapping_table$id_norm),
  ,
  drop = FALSE
]


############################################################
# ANNOTATION SUMMARY
############################################################

print(
  paste(
    "Total CGD features:",
    nrow(mapping_table)
  )
)

print(
  paste(
    "CGD features with gene names:",
    sum(
      !is.na(mapping_table$gene_symbol)
    )
  )
)


############################################################
# 3. LOAD FEATURECOUNTS MATRIX
############################################################

print("======================================================")
print("LOADING FEATURECOUNTS MATRIX")
print("======================================================")

if (!file.exists("candida_counts_matrix.txt")) {
  
  stop(
    paste0(
      "\nERROR: candida_counts_matrix.txt was not found.\n",
      "Place it in the same working directory as this script.\n"
    )
  )
  
}


############################################################
# ROBUST FEATURECOUNTS READER
#
# Avoids:
# "more columns than column names"
############################################################

read_featurecounts <- function(file) {
  
  lines <- readLines(
    file,
    warn = FALSE
  )
  
  lines <- lines[
    nzchar(
      trimws(lines)
    )
  ]
  
  lines <- lines[
    !grepl(
      "^#",
      trimws(lines)
    )
  ]
  
  if (length(lines) < 2) {
    stop(
      "ERROR: FeatureCounts file does not contain enough data."
    )
  }
  
  
  ##########################################################
  # HEADER
  ##########################################################
  
  header <- strsplit(
    lines[1],
    "\t",
    fixed = TRUE
  )[[1]]
  
  header[1] <- sub(
    "^\ufeff",
    "",
    header[1]
  )
  
  
  if (length(header) != 18) {
    
    stop(
      paste0(
        "\nERROR: Expected 18 columns in featureCounts header.\n",
        "Found ",
        length(header),
        " columns.\n\n",
        "Header detected:\n",
        paste(
          header,
          collapse = "\n"
        )
      )
    )
    
  }
  
  
  ##########################################################
  # DATA
  ##########################################################
  
  data_lines <- lines[-1]
  
  split_lines <- strsplit(
    data_lines,
    "\t",
    fixed = TRUE
  )
  
  
  ##########################################################
  # REMOVE TRAILING EMPTY FIELDS
  ##########################################################
  
  split_lines <- lapply(
    split_lines,
    function(x) {
      
      while (
        length(x) > 0 &&
        trimws(
          x[length(x)]
        ) == ""
      ) {
        
        x <- x[
          -length(x)
        ]
        
      }
      
      x
      
    }
  )
  
  
  field_counts <- lengths(
    split_lines
  )
  
  
  bad_rows <- which(
    field_counts != 18
  )
  
  
  if (length(bad_rows) > 0) {
    
    first_bad <- bad_rows[1]
    
    stop(
      paste0(
        "\nERROR: FeatureCounts contains a row with ",
        field_counts[first_bad],
        " columns instead of 18.\n\n",
        "Problem occurs at data row ",
        first_bad,
        ".\n\n",
        "Line:\n",
        data_lines[first_bad]
      )
    )
    
  }
  
  
  matrix_data <- do.call(
    rbind,
    split_lines
  )
  
  
  counts_df <- as.data.frame(
    matrix_data,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  colnames(counts_df) <- header
  
  return(counts_df)
  
}


counts_raw <- read_featurecounts(
  "candida_counts_matrix.txt"
)


print(
  paste(
    "Total columns detected:",
    ncol(counts_raw)
  )
)

print(
  "Column names detected:"
)

print(
  colnames(counts_raw)
)


############################################################
# 4. IDENTIFY AND REORDER SAMPLE COLUMNS
############################################################

print("======================================================")
print("IDENTIFYING RNA-seq SAMPLE COLUMNS")
print("======================================================")


############################################################
# FEATURECOUNTS STRUCTURE
#
# 1 = Geneid
# 2 = Chr
# 3 = Start
# 4 = End
# 5 = Strand
# 6 = Length
# 7-18 = samples
############################################################

sample_columns_original <- 7:18

original_sample_names <- colnames(
  counts_raw
)[sample_columns_original]


print(
  "Original sample order in featureCounts file:"
)

print(
  original_sample_names
)


############################################################
# EXTRACT N44VRL SAMPLE NUMBERS
############################################################

sample_numbers <- suppressWarnings(
  as.numeric(
    sub(
      ".*N44VRL_([0-9]+)_.*",
      "\\1",
      original_sample_names
    )
  )
)


print(
  "Detected N44VRL sample numbers:"
)

print(
  sample_numbers
)


############################################################
# VERIFY SAMPLE NUMBERS
############################################################

if (
  length(sample_numbers) != 12 ||
  any(is.na(sample_numbers)) ||
  !setequal(
    sample_numbers,
    1:12
  )
) {
  
  stop(
    paste0(
      "\nERROR: Could not correctly identify N44VRL ",
      "samples 1-12.\n\n",
      "Detected sample names:\n",
      paste(
        original_sample_names,
        collapse = "\n"
      ),
      "\n\nDetected sample numbers:\n",
      paste(
        sample_numbers,
        collapse = ", "
      )
    )
  )
  
}


############################################################
# REORDER BY N44VRL SAMPLE NUMBER
############################################################

sample_columns <- sample_columns_original[
  order(
    sample_numbers
  )
]

original_sample_names <- colnames(
  counts_raw
)[sample_columns]


print("======================================================")
print("CORRECTED SAMPLE ORDER")
print("======================================================")


for (i in seq_along(original_sample_names)) {
  
  print(
    paste(
      i,
      ":",
      original_sample_names[i]
    )
  )
  
}


############################################################
# BIOLOGICAL SAMPLE ORDER
############################################################

expected_genotype <- c(
  "WT",
  "WT",
  "WT",
  "WT",
  "WT",
  "WT",
  "KCS1",
  "KCS1",
  "KCS1",
  "KCS1",
  "KCS1",
  "KCS1"
)

expected_pi <- c(
  "plus",
  "plus",
  "plus",
  "minus",
  "minus",
  "minus",
  "plus",
  "plus",
  "plus",
  "minus",
  "minus",
  "minus"
)


############################################################
# DETERMINE GENOTYPE FROM SAMPLE NAMES
############################################################

detected_genotype <- ifelse(
  
  grepl(
    "WT",
    original_sample_names,
    ignore.case = TRUE
  ),
  
  "WT",
  
  ifelse(
    
    grepl(
      "KCS1",
      original_sample_names,
      ignore.case = TRUE
    ),
    
    "KCS1",
    
    NA
    
  )
  
)


############################################################
# DETERMINE PI CONDITION
############################################################

detected_pi <- ifelse(
  
  grepl(
    "\\+Pi",
    original_sample_names,
    ignore.case = TRUE
  ),
  
  "plus",
  
  ifelse(
    
    grepl(
      "-Pi",
      original_sample_names,
      ignore.case = TRUE
    ),
    
    "minus",
    
    NA
    
  )
  
)


############################################################
# VERIFY SAMPLE INFORMATION
############################################################

if (
  any(is.na(detected_genotype)) ||
  any(is.na(detected_pi))
) {
  
  print(original_sample_names)
  print(detected_genotype)
  print(detected_pi)
  
  stop(
    "\nERROR: Could not determine genotype or Pi condition from sample names."
  )
  
}


############################################################
# FINAL BIOLOGICAL ORDER CHECK
############################################################

if (
  !identical(
    detected_genotype,
    expected_genotype
  ) ||
  !identical(
    detected_pi,
    expected_pi
  )
) {
  
  print(
    "Expected genotype:"
  )
  
  print(
    expected_genotype
  )
  
  print(
    "Detected genotype:"
  )
  
  print(
    detected_genotype
  )
  
  print(
    "Expected Pi:"
  )
  
  print(
    expected_pi
  )
  
  print(
    "Detected Pi:"
  )
  
  print(
    detected_pi
  )
  
  stop(
    "\nERROR: Sample order does not match the expected experimental design."
  )
  
}


print(
  "SUCCESS: Sample order is correct."
)


############################################################
# 5. EXTRACT GENE IDS AND COUNTS
############################################################

counts_clean <- counts_raw[
  ,
  c(1, sample_columns),
  drop = FALSE
]

colnames(counts_clean)[1] <- "Geneid"


############################################################
# CONVERT COUNTS TO NUMERIC
############################################################

for (i in 2:ncol(counts_clean)) {
  
  counts_clean[[i]] <- as.numeric(
    counts_clean[[i]]
  )
  
}


############################################################
# CHECK FOR NA VALUES
############################################################

if (
  any(
    is.na(
      counts_clean[, -1]
    )
  )
) {
  
  stop(
    "\nERROR: NA values detected in count matrix."
  )
  
}


############################################################
# COLLAPSE DUPLICATE GENE IDS
############################################################

print(
  "Collapsing duplicate gene IDs..."
)

counts_fixed <- aggregate(
  . ~ Geneid,
  data = counts_clean,
  FUN = sum
)


############################################################
# 6. REMOVE tRNA AND rRNA GENES
############################################################

print("======================================================")
print("REMOVING tRNA / rRNA FEATURES")
print("======================================================")


############################################################
# NORMALIZE COUNT IDS
############################################################

stripped_keys <- gsub(
  "^CAALFM_",
  "",
  counts_fixed$Geneid,
  ignore.case = TRUE
)

normalized_matrix_keys <- normalize_id(
  stripped_keys
)


############################################################
# MATCH FEATURE TYPES
############################################################

matched_features <- mapping_table$feature_type[
  match(
    normalized_matrix_keys,
    mapping_table$id_norm
  )
]


############################################################
# MATCH GENE NAMES
############################################################

matched_symbols <- mapping_table$gene_symbol[
  match(
    normalized_matrix_keys,
    mapping_table$id_norm
  )
]


############################################################
# IDENTIFY tRNA / rRNA
############################################################

is_trna_rrna <- (
  
  grepl(
    "tRNA|rRNA",
    counts_fixed$Geneid,
    ignore.case = TRUE
  ) |
    
    grepl(
      "tRNA|rRNA",
      matched_features,
      ignore.case = TRUE
    ) |
    
    grepl(
      "^tRNA|^rRNA",
      matched_symbols,
      ignore.case = TRUE
    )
  
)

is_trna_rrna[
  is.na(is_trna_rrna)
] <- FALSE


print(
  paste(
    "Features before tRNA/rRNA filtering:",
    nrow(counts_fixed)
  )
)

print(
  paste(
    "tRNA/rRNA features removed:",
    sum(is_trna_rrna)
  )
)


counts_fixed <- counts_fixed[
  !is_trna_rrna,
  ,
  drop = FALSE
]


print(
  paste(
    "Features remaining:",
    nrow(counts_fixed)
  )
)


############################################################
# CREATE COUNT MATRIX
############################################################

rownames(counts_fixed) <- counts_fixed$Geneid

counts <- counts_fixed[
  ,
  -1,
  drop = FALSE
]


############################################################
# 7. CREATE MASTER GENE ANNOTATION
############################################################

print("======================================================")
print("ANNOTATING GENES")
print("======================================================")


original_gene_ids <- rownames(
  counts
)


############################################################
# STRIP CAALFM PREFIX
############################################################

stripped_ids <- gsub(
  "^CAALFM_",
  "",
  original_gene_ids,
  ignore.case = TRUE
)


normalized_ids <- normalize_id(
  stripped_ids
)


############################################################
# MATCH CGD GENE NAMES
############################################################

matched_gene_names <- mapping_table$gene_symbol[
  match(
    normalized_ids,
    mapping_table$id_norm
  )
]


############################################################
# DISPLAY NAME
#
# If a CGD gene name exists, use it.
# Otherwise retain the systematic ID.
############################################################

display_gene_names <- ifelse(
  
  !is.na(matched_gene_names) &
    matched_gene_names != "",
  
  matched_gene_names,
  
  stripped_ids
  
)


############################################################
# MASTER ANNOTATION TABLE
############################################################

gene_annotation <- data.frame(
  
  Systematic_ID = original_gene_ids,
  
  Gene_Name = display_gene_names,
  
  stringsAsFactors = FALSE
  
)


############################################################
# MAPPING SUMMARY
############################################################

mapped_count <- sum(
  !is.na(matched_gene_names) &
    matched_gene_names != ""
)

unmapped_count <- sum(
  is.na(matched_gene_names) |
    matched_gene_names == ""
)


print(
  paste(
    "Genes mapped to CGD gene names:",
    mapped_count
  )
)

print(
  paste(
    "Genes without a CGD gene name:",
    unmapped_count
  )
)


############################################################
# SAVE ANNOTATION
############################################################

write.csv(
  
  gene_annotation,
  
  file.path(
    output_dir,
    "02_DE_results",
    "Gene_ID_to_Gene_Name_Annotation.csv"
  ),
  
  row.names = FALSE
  
)


############################################################
# HELPER FUNCTION FOR GENE NAME LOOKUP
############################################################

get_gene_names <- function(ids) {
  
  names_out <- gene_annotation$Gene_Name[
    match(
      ids,
      gene_annotation$Systematic_ID
    )
  ]
  
  missing <- (
    is.na(names_out) |
      names_out == ""
  )
  
  names_out[
    missing
  ] <- ids[
    missing
  ]
  
  return(
    names_out
  )
  
}


############################################################
# 8. SAMPLE METADATA
############################################################

short_sample_names <- c(
  
  "WT_Pi_1",
  "WT_Pi_2",
  "WT_Pi_3",
  
  "WT_minusPi_1",
  "WT_minusPi_2",
  "WT_minusPi_3",
  
  "KCS1_Pi_1",
  "KCS1_Pi_2",
  "KCS1_Pi_3",
  
  "KCS1_minusPi_1",
  "KCS1_minusPi_2",
  "KCS1_minusPi_3"
  
)


############################################################
# HUMAN-READABLE CONDITION LABELS
#
# THESE ARE THE LABELS THAT WILL APPEAR ON FIGURES.
############################################################

condition_labels <- c(
  
  "WT+Pi",
  "WT+Pi",
  "WT+Pi",
  
  "WT-Pi",
  "WT-Pi",
  "WT-Pi",
  
  "KCS1+Pi",
  "KCS1+Pi",
  "KCS1+Pi",
  
  "KCS1-Pi",
  "KCS1-Pi",
  "KCS1-Pi"
  
)


replicate_vector <- c(
  1, 2, 3,
  1, 2, 3,
  1, 2, 3,
  1, 2, 3
)


############################################################
# RENAME COUNT MATRIX
############################################################

colnames(counts) <- short_sample_names


############################################################
# METADATA
############################################################

metadata <- data.frame(
  
  Original_Sample_ID = original_sample_names,
  
  Sample_Name = short_sample_names,
  
  Display_Name = condition_labels,
  
  Genotype = factor(
    detected_genotype,
    levels = c(
      "WT",
      "KCS1"
    )
  ),
  
  Pi = factor(
    detected_pi,
    levels = c(
      "plus",
      "minus"
    )
  ),
  
  Replicate = replicate_vector,
  
  stringsAsFactors = FALSE
  
)


rownames(metadata) <- short_sample_names


############################################################
# INTERNAL CONDITION VARIABLE
############################################################

metadata$Condition <- factor(
  
  paste(
    metadata$Genotype,
    metadata$Pi,
    sep = "_"
  ),
  
  levels = c(
    "WT_plus",
    "WT_minus",
    "KCS1_plus",
    "KCS1_minus"
  )
  
)


############################################################
# CONDITION LABEL FACTOR
############################################################

metadata$Condition_Label <- factor(
  
  condition_labels,
  
  levels = c(
    "WT+Pi",
    "WT-Pi",
    "KCS1+Pi",
    "KCS1-Pi"
  )
  
)


print("======================================================")
print("SAMPLE METADATA")
print("======================================================")

print(
  metadata
)


############################################################
# VERIFY REPLICATES
############################################################

replicate_table <- table(
  metadata$Condition
)

print(
  replicate_table
)

if (!all(replicate_table == 3)) {
  
  stop(
    "ERROR: Every condition must contain exactly 3 replicates."
  )
  
}


write.csv(
  
  metadata,
  
  file.path(
    output_dir,
    "sample_metadata.csv"
  ),
  
  row.names = TRUE
  
)


############################################################
# 9. FILTER LOW-COUNT GENES
############################################################

print("======================================================")
print("FILTERING LOW-COUNT GENES")
print("======================================================")


keep_genes <- rowSums(
  counts
) >= 10


print(
  paste(
    "Genes before filtering:",
    nrow(counts)
  )
)

print(
  paste(
    "Genes passing count filter:",
    sum(keep_genes)
  )
)

print(
  paste(
    "Genes removed:",
    sum(!keep_genes)
  )
)


counts_filtered <- counts[
  keep_genes,
  ,
  drop = FALSE
]


############################################################
# 10. DESEQ2
############################################################

print("======================================================")
print("CREATING DESEQ2 OBJECT")
print("======================================================")


dds <- DESeqDataSetFromMatrix(
  
  countData = round(
    as.matrix(
      counts_filtered
    )
  ),
  
  colData = metadata,
  
  design = ~ Genotype * Pi
  
)


dds <- DESeq(
  dds
)


saveRDS(
  
  dds,
  
  file.path(
    output_dir,
    "02_DE_results",
    "DESeq2_object.rds"
  )
  
)


############################################################
# 11. VST
############################################################

print("======================================================")
print("VST TRANSFORMATION")
print("======================================================")


vsd <- vst(
  dds,
  blind = FALSE
)

vst_matrix <- assay(
  vsd
)


write.csv(
  
  vst_matrix,
  
  file.path(
    output_dir,
    "01_QC",
    "VST_expression_matrix.csv"
  ),
  
  row.names = TRUE
  
)


############################################################
# 12. LIBRARY SIZE QC
############################################################

library_sizes <- colSums(
  counts(dds)
)


library_size_df <- data.frame(
  
  Sample = names(library_sizes),
  
  Condition = metadata[
    names(library_sizes),
    "Display_Name"
  ],
  
  Library_Size = as.numeric(
    library_sizes
  ),
  
  stringsAsFactors = FALSE
  
)


write.csv(
  
  library_size_df,
  
  file.path(
    output_dir,
    "01_QC",
    "Library_sizes.csv"
  ),
  
  row.names = FALSE
  
)


p_library <- ggplot(
  
  library_size_df,
  
  aes(
    x = Sample,
    y = Library_Size
  )
  
) +
  
  geom_col() +
  
  theme_bw() +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  ) +
  
  labs(
    title = "RNA-seq library sizes",
    x = "Sample",
    y = "Total assigned reads"
  )


ggsave(
  
  file.path(
    output_dir,
    "01_QC",
    "Library_sizes.png"
  ),
  
  p_library,
  
  width = 10,
  height = 6,
  dpi = 300
  
)


############################################################
# 13. PCA
############################################################

print("======================================================")
print("CREATING PCA")
print("======================================================")


pca_data <- plotPCA(
  
  vsd,
  
  intgroup = c(
    "Genotype",
    "Pi"
  ),
  
  returnData = TRUE
  
)


percent_variance <- round(
  
  100 *
    attr(
      pca_data,
      "percentVar"
    )
  
)


p_pca <- ggplot(
  
  pca_data,
  
  aes(
    x = PC1,
    y = PC2,
    label = name
  )
  
) +
  
  geom_point(
    size = 4
  ) +
  
  geom_text_repel(
    size = 3
  ) +
  
  theme_bw() +
  
  labs(
    
    title = "PCA of RNA-seq samples",
    
    x = paste0(
      "PC1: ",
      percent_variance[1],
      "% variance"
    ),
    
    y = paste0(
      "PC2: ",
      percent_variance[2],
      "% variance"
    )
    
  )


ggsave(
  
  file.path(
    output_dir,
    "06_PCA",
    "PCA_PC1_PC2.png"
  ),
  
  p_pca,
  
  width = 9,
  height = 7,
  dpi = 300
  
)


############################################################
# 14. SAMPLE CORRELATION
############################################################

cor_matrix <- cor(
  
  vst_matrix,
  
  method = "pearson"
  
)


write.csv(
  
  cor_matrix,
  
  file.path(
    output_dir,
    "07_Sample_correlations",
    "Sample_Pearson_correlations.csv"
  )
  
)


png(
  
  file.path(
    output_dir,
    "07_Sample_correlations",
    "Sample_Pearson_correlations.png"
  ),
  
  width = 2200,
  height = 2000,
  res = 300
  
)


cor_annotation <- data.frame(
  
  Condition = metadata[
    colnames(vst_matrix),
    "Display_Name"
  ]
  
)

rownames(cor_annotation) <- colnames(
  vst_matrix
)


pheatmap(
  
  cor_matrix,
  
  annotation_col = cor_annotation,
  
  annotation_row = cor_annotation,
  
  clustering_distance_rows = "correlation",
  
  clustering_distance_cols = "correlation",
  
  main = "Sample Pearson correlations",
  
  fontsize = 8
  
)


dev.off()


############################################################
# 15. DESEQ2 COEFFICIENTS
############################################################

print("======================================================")
print("DESEQ2 COEFFICIENTS")
print("======================================================")


results_names <- resultsNames(
  dds
)

print(
  results_names
)


############################################################
# GENOTYPE COEFFICIENT
############################################################

genotype_coef <- results_names[
  results_names == "Genotype_KCS1_vs_WT"
]


############################################################
# PI COEFFICIENT
############################################################

pi_coef <- results_names[
  results_names == "Pi_minus_vs_plus"
]


############################################################
# INTERACTION COEFFICIENT
############################################################

interaction_coef <- results_names[
  grepl(
    "Genotype.*Pi|Pi.*Genotype",
    results_names
  )
]


############################################################
# VERIFY COEFFICIENTS
############################################################

if (length(genotype_coef) != 1) {
  
  stop(
    paste0(
      "Could not identify genotype coefficient.\n",
      "Available coefficients:\n",
      paste(
        results_names,
        collapse = "\n"
      )
    )
  )
  
}


if (length(pi_coef) != 1) {
  
  stop(
    paste0(
      "Could not identify Pi coefficient.\n",
      "Available coefficients:\n",
      paste(
        results_names,
        collapse = "\n"
      )
    )
  )
  
}


if (length(interaction_coef) != 1) {
  
  stop(
    paste0(
      "Could not identify interaction coefficient.\n",
      "Available coefficients:\n",
      paste(
        results_names,
        collapse = "\n"
      )
    )
  )
  
}


############################################################
# PRINT COEFFICIENTS
############################################################

print(
  paste(
    "Genotype coefficient:",
    genotype_coef
  )
)

print(
  paste(
    "Pi coefficient:",
    pi_coef
  )
)

print(
  paste(
    "Interaction coefficient:",
    interaction_coef
  )
)
  
  ############################################################
  # 16. ALL SIX PAIRWISE COMPARISONS
  ############################################################
  #
  # Reference condition:
  #
  # WT+Pi = 0
  #
  # WT-Pi = Pi
  #
  # KCS1+Pi = Genotype
  #
  # KCS1-Pi =
  # Genotype + Pi + Interaction
  #
  ############################################################
  
  
  ############################################################
  # 1. WT+Pi vs WT-Pi
  #
  # Result is WT-Pi - WT+Pi
  ############################################################
  
  res_WT_minus_vs_WT_plus <- results(
    
    dds,
    
    name = pi_coef
    
  )
  
  
  ############################################################
  # 2. WT+Pi vs KCS1+Pi
  #
  # Result is KCS1+Pi - WT+Pi
  ############################################################
  
  res_KCS1_plus_vs_WT_plus <- results(
    
    dds,
    
    name = genotype_coef
    
  )
  
  
  ############################################################
  # 3. WT+Pi vs KCS1-Pi
  #
  # Result is KCS1-Pi - WT+Pi
  #
  # = Genotype + Pi + Interaction
  ############################################################
  
  res_KCS1_minus_vs_WT_plus <- results(
    
    dds,
    
    contrast = list(
      c(
        genotype_coef,
        pi_coef,
        interaction_coef
      )
    )
    
  )
  
  
  ############################################################
  # 4. WT-Pi vs KCS1+Pi
  #
  # Result is KCS1+Pi - WT-Pi
  #
  # = Genotype - Pi
  ############################################################
  
  res_KCS1_plus_vs_WT_minus <- results(
    
    dds,
    
    contrast = list(
      c(
        genotype_coef
      ),
      c(
        pi_coef
      )
    )
    
  )
  
  
  ############################################################
  # 5. WT-Pi vs KCS1-Pi
  #
  # Result is KCS1-Pi - WT-Pi
  #
  # = Genotype + Interaction
  ############################################################
  
  res_KCS1_minus_vs_WT_minus <- results(
    
    dds,
    
    contrast = list(
      c(
        genotype_coef,
        interaction_coef
      )
    )
    
  )
  
  
  ############################################################
  # 6. KCS1+Pi vs KCS1-Pi
  #
  # Result is KCS1-Pi - KCS1+Pi
  #
  # = Pi + Interaction
  ############################################################
  
  res_KCS1_minus_vs_KCS1_plus <- results(
    
    dds,
    
    contrast = list(
      c(
        pi_coef,
        interaction_coef
      )
    )
    
  )
  
  
  ############################################################
  # 7. FORMAL GENOTYPE x PI INTERACTION
  ############################################################
  
  res_interaction <- results(
    
    dds,
    
    name = interaction_coef
    
  )
  
  
  ############################################################
  # 17. ANNOTATE ALL DE RESULTS
  ############################################################
  
  annotate_results <- function(
    result,
    comparison_label
  ) {
    
    df <- as.data.frame(
      result
    )
    
    df$Systematic_ID <- rownames(
      df
    )
    
    df$Gene_Name <- get_gene_names(
      df$Systematic_ID
    )
    
    df$Comparison <- comparison_label
    
    
    ##########################################################
    # Put gene information first
    ##########################################################
    
    df <- df[
      ,
      c(
        "Systematic_ID",
        "Gene_Name",
        "Comparison",
        setdiff(
          colnames(df),
          c(
            "Systematic_ID",
            "Gene_Name",
            "Comparison"
          )
        )
      ),
      drop = FALSE
    ]
    
    return(df)
    
  }
  
  
  res_WT_minus_vs_WT_plus_df <- annotate_results(
    res_WT_minus_vs_WT_plus,
    "WT-Pi vs WT+Pi"
  )
  
  
  res_KCS1_plus_vs_WT_plus_df <- annotate_results(
    res_KCS1_plus_vs_WT_plus,
    "KCS1+Pi vs WT+Pi"
  )
  
  
  res_KCS1_minus_vs_WT_plus_df <- annotate_results(
    res_KCS1_minus_vs_WT_plus,
    "KCS1-Pi vs WT+Pi"
  )
  
  
  res_KCS1_plus_vs_WT_minus_df <- annotate_results(
    res_KCS1_plus_vs_WT_minus,
    "KCS1+Pi vs WT-Pi"
  )
  
  
  res_KCS1_minus_vs_WT_minus_df <- annotate_results(
    res_KCS1_minus_vs_WT_minus,
    "KCS1-Pi vs WT-Pi"
  )
  
  
  res_KCS1_minus_vs_KCS1_plus_df <- annotate_results(
    res_KCS1_minus_vs_KCS1_plus,
    "KCS1-Pi vs KCS1+Pi"
  )
  
  
  res_interaction_df <- annotate_results(
    res_interaction,
    "Genotype × Pi interaction"
  )
  
  
  ############################################################
  # CREATE NAMED RESULT LIST
  ############################################################
  
  all_results <- list(
    
    "WT-Pi_vs_WT+Pi" =
      res_WT_minus_vs_WT_plus_df,
    
    "KCS1+Pi_vs_WT+Pi" =
      res_KCS1_plus_vs_WT_plus_df,
    
    "KCS1-Pi_vs_WT+Pi" =
      res_KCS1_minus_vs_WT_plus_df,
    
    "KCS1+Pi_vs_WT-Pi" =
      res_KCS1_plus_vs_WT_minus_df,
    
    "KCS1-Pi_vs_WT-Pi" =
      res_KCS1_minus_vs_WT_minus_df,
    
    "KCS1-Pi_vs_KCS1+Pi" =
      res_KCS1_minus_vs_KCS1_plus_df,
    
    "Genotype_x_Pi_interaction" =
      res_interaction_df
    
  )
  
  
  ############################################################
  # 18. SAVE ALL DE RESULTS
  ############################################################
  
  for (comparison_name in names(all_results)) {
    
    write.csv(
      
      all_results[[comparison_name]],
      
      file.path(
        output_dir,
        "02_DE_results",
        paste0(
          "DE_",
          comparison_name,
          ".csv"
        )
      ),
      
      row.names = FALSE
      
    )
    
  }
  
  
  ############################################################
  # 19. SIGNIFICANT GENE LISTS
  ############################################################
  
  get_significant_genes <- function(
    result_df,
    padj_cutoff = 0.05
  ) {
    
    result_df %>%
      
      filter(
        !is.na(padj),
        padj < padj_cutoff
      ) %>%
      
      arrange(
        padj
      )
    
  }
  
  
  for (comparison_name in names(all_results)) {
    
    sig_genes <- get_significant_genes(
      all_results[[comparison_name]]
    )
    
    
    write.csv(
      
      sig_genes,
      
      file.path(
        output_dir,
        "08_Significant_gene_lists",
        paste0(
          "Significant_",
          comparison_name,
          ".csv"
        )
      ),
      
      row.names = FALSE
      
    )
    
  }
  
  
  ############################################################
  # 20. VOLCANO PLOT FUNCTION
  ############################################################
  
  make_volcano <- function(
    result_df,
    plot_title,
    output_file,
    fc_cutoff = 1,
    padj_cutoff = 0.05
  ) {
    
    
    plot_df <- result_df %>%
      
      filter(
        !is.na(log2FoldChange),
        !is.na(padj)
      )
    
    
    if (nrow(plot_df) == 0) {
      return(NULL)
    }
    
    
    plot_df$Significance <- "Not significant"
    
    
    plot_df$Significance[
      plot_df$padj < padj_cutoff &
        plot_df$log2FoldChange >= fc_cutoff
    ] <- "Upregulated"
    
    
    plot_df$Significance[
      plot_df$padj < padj_cutoff &
        plot_df$log2FoldChange <= -fc_cutoff
    ] <- "Downregulated"
    
    
    plot_df$Label <- NA_character_
    
    
    label_genes <- plot_df %>%
      
      filter(
        padj < padj_cutoff,
        abs(log2FoldChange) >= fc_cutoff
      ) %>%
      
      arrange(
        padj
      ) %>%
      
      head(15)
    
    
    if (nrow(label_genes) > 0) {
      
      plot_df$Label[
        match(
          label_genes$Systematic_ID,
          plot_df$Systematic_ID
        )
      ] <- label_genes$Gene_Name
      
    }
    
    
    p <- ggplot(
      
      plot_df,
      
      aes(
        x = log2FoldChange,
        y = -log10(padj),
        shape = Significance
      )
      
    ) +
      
      geom_point(
        alpha = 0.6,
        size = 1.5
      ) +
      
      geom_vline(
        xintercept = c(
          -fc_cutoff,
          fc_cutoff
        ),
        linetype = "dashed"
      ) +
      
      geom_hline(
        yintercept = -log10(
          padj_cutoff
        ),
        linetype = "dashed"
      ) +
      
      geom_text_repel(
        aes(
          label = Label
        ),
        na.rm = TRUE,
        size = 3,
        max.overlaps = 20
      ) +
      
      theme_bw() +
      
      labs(
        title = plot_title,
        x = "log2 fold change",
        y = "-log10 adjusted p-value"
      )
    
    
    ggsave(
      
      output_file,
      
      p,
      
      width = 9,
      height = 7,
      dpi = 300
      
    )
    
    
    return(p)
    
  }
  
  
  ############################################################
  # 21. CREATE VOLCANO PLOTS FOR ALL SIX COMPARISONS
  ############################################################
  
  volcano_titles <- c(
    
    "WT-Pi vs WT+Pi",
    "KCS1+Pi vs WT+Pi",
    "KCS1-Pi vs WT+Pi",
    "KCS1+Pi vs WT-Pi",
    "KCS1-Pi vs WT-Pi",
    "KCS1-Pi vs KCS1+Pi"
    
  )
  
  
  volcano_files <- c(
    
    "Volcano_WT-Pi_vs_WT+Pi.png",
    "Volcano_KCS1+Pi_vs_WT+Pi.png",
    "Volcano_KCS1-Pi_vs_WT+Pi.png",
    "Volcano_KCS1+Pi_vs_WT-Pi.png",
    "Volcano_KCS1-Pi_vs_WT-Pi.png",
    "Volcano_KCS1-Pi_vs_KCS1+Pi.png"
    
  )
  
  
  pairwise_results <- all_results[
    1:6
  ]
  
  
  for (i in seq_along(pairwise_results)) {
    
    make_volcano(
      
      pairwise_results[[i]],
      
      volcano_titles[i],
      
      file.path(
        output_dir,
        "03_Volcano",
        volcano_files[i]
      )
      
    )
    
  }
  
  
  ############################################################
  # INTERACTION VOLCANO
  ############################################################
  
  make_volcano(
    
    res_interaction_df,
    
    "Genotype × Pi interaction",
    
    file.path(
      output_dir,
      "03_Volcano",
      "Volcano_Genotype_x_Pi_interaction.png"
    )
    
  )
  
  
  ############################################################
  # 22. MA PLOT FUNCTION
  ############################################################
  
  make_ma_plot <- function(
    result_df,
    plot_title,
    output_file
  ) {
    
    
    plot_df <- result_df %>%
      
      filter(
        !is.na(baseMean),
        !is.na(log2FoldChange)
      )
    
    
    p <- ggplot(
      
      plot_df,
      
      aes(
        x = log10(baseMean + 1),
        y = log2FoldChange
      )
      
    ) +
      
      geom_point(
        alpha = 0.4,
        size = 1
      ) +
      
      geom_hline(
        yintercept = 0,
        linetype = "dashed"
      ) +
      
      theme_bw() +
      
      labs(
        title = plot_title,
        x = "log10 mean normalized expression",
        y = "log2 fold change"
      )
    
    
    ggsave(
      
      output_file,
      
      p,
      
      width = 9,
      height = 7,
      dpi = 300
      
    )
    
  }
  
  
  ############################################################
  # 23. MA PLOTS FOR ALL COMPARISONS
  ############################################################
  
  ma_files <- c(
    
    "MA_WT-Pi_vs_WT+Pi.png",
    "MA_KCS1+Pi_vs_WT+Pi.png",
    "MA_KCS1-Pi_vs_WT+Pi.png",
    "MA_KCS1+Pi_vs_WT-Pi.png",
    "MA_KCS1-Pi_vs_WT-Pi.png",
    "MA_KCS1-Pi_vs_KCS1+Pi.png"
    
  )
  
  
  for (i in seq_along(pairwise_results)) {
    
    make_ma_plot(
      
      pairwise_results[[i]],
      
      volcano_titles[i],
      
      file.path(
        output_dir,
        "04_MA_plots",
        ma_files[i]
      )
      
    )
    
  }
  
  
  make_ma_plot(
    
    res_interaction_df,
    
    "Genotype × Pi interaction",
    
    file.path(
      output_dir,
      "04_MA_plots",
      "MA_Genotype_x_Pi_interaction.png"
    )
    
  )
  
  
  ############################################################
  # 24. TOP VARIABLE GENE HEATMAP
  ############################################################
  
  print("======================================================")
  print("CREATING TOP VARIABLE GENE HEATMAP")
  print("======================================================")
  
  
  gene_variances <- rowVars(
    vst_matrix
  )
  
  
  top_n <- min(
    50,
    length(gene_variances)
  )
  
  
  top_variable_genes <- names(
    sort(
      gene_variances,
      decreasing = TRUE
    )
  )[1:top_n]
  
  
  heatmap_matrix <- vst_matrix[
    top_variable_genes,
    ,
    drop = FALSE
  ]
  
  
  heatmap_gene_names <- get_gene_names(
    rownames(
      heatmap_matrix
    )
  )
  
  
  rownames(
    heatmap_matrix
  ) <- make.unique(
    heatmap_gene_names
  )
  
  
  heatmap_scaled <- t(
    scale(
      t(
        heatmap_matrix
      )
    )
  )
  
  
  annotation_col <- data.frame(
    
    Condition = metadata[
      colnames(
        heatmap_scaled
      ),
      "Display_Name"
    ]
    
  )
  
  rownames(annotation_col) <- colnames(
    heatmap_scaled
  )
  
  
  png(
    
    file.path(
      output_dir,
      "05_Heatmaps",
      "Heatmap_Top_50_Variable_Genes.png"
    ),
    
    width = 2400,
    height = 2400,
    res = 300
    
  )
  
  
  pheatmap(
    
    heatmap_scaled,
    
    scale = "none",
    
    annotation_col = annotation_col,
    
    cluster_rows = TRUE,
    
    cluster_cols = TRUE,
    
    show_rownames = TRUE,
    
    fontsize_row = 6,
    
    fontsize_col = 9,
    
    main = "Top 50 most variable genes"
    
  )
  
  
  dev.off()
  
  
  ############################################################
  # 25. DE HEATMAP FUNCTION
  ############################################################
  
  make_de_heatmap <- function(
    result_df,
    plot_title,
    output_file,
    n_genes = 50
  ) {
    
    
    significant_genes <- result_df %>%
      
      filter(
        !is.na(padj),
        padj < 0.05
      ) %>%
      
      arrange(
        padj
      )
    
    
    if (nrow(significant_genes) < 2) {
      
      print(
        paste(
          "Too few significant genes for:",
          plot_title
        )
      )
      
      return(NULL)
      
    }
    
    
    selected_genes <- head(
      
      significant_genes$Systematic_ID,
      
      n_genes
      
    )
    
    
    selected_genes <- intersect(
      
      selected_genes,
      
      rownames(
        vst_matrix
      )
      
    )
    
    
    if (length(selected_genes) < 2) {
      return(NULL)
    }
    
    
    heatmap_data <- vst_matrix[
      
      selected_genes,
      ,
      drop = FALSE
      
    ]
    
    
    gene_names <- get_gene_names(
      
      rownames(
        heatmap_data
      )
      
    )
    
    
    rownames(
      heatmap_data
    ) <- make.unique(
      gene_names
    )
    
    
    heatmap_scaled <- t(
      
      scale(
        t(
          heatmap_data
        )
      )
      
    )
    
    
    ##########################################################
    # CONDITION LABELS
    ##########################################################
    
    annotation_col <- data.frame(
      
      Condition = metadata[
        colnames(
          heatmap_scaled
        ),
        "Display_Name"
      ]
      
    )
    
    
    rownames(annotation_col) <- colnames(
      heatmap_scaled
    )
    
    
    ##########################################################
    # FIXED-SIZE PNG
    ##########################################################
    
    png(
      
      output_file,
      
      width = 2400,
      height = 2600,
      res = 300
      
    )
    
    
    pheatmap(
      
      heatmap_scaled,
      
      scale = "none",
      
      annotation_col = annotation_col,
      
      cluster_rows = TRUE,
      
      cluster_cols = TRUE,
      
      show_rownames = TRUE,
      
      fontsize_row = 7,
      
      fontsize_col = 9,
      
      main = plot_title
      
    )
    
    
    dev.off()
    
  }
  
  
  ############################################################
  # 26. HEATMAPS FOR ALL SIX PAIRWISE COMPARISONS
  ############################################################
  
  heatmap_files <- c(
    
    "Heatmap_WT-Pi_vs_WT+Pi.png",
    "Heatmap_KCS1+Pi_vs_WT+Pi.png",
    "Heatmap_KCS1-Pi_vs_WT+Pi.png",
    "Heatmap_KCS1+Pi_vs_WT-Pi.png",
    "Heatmap_KCS1-Pi_vs_WT-Pi.png",
    "Heatmap_KCS1-Pi_vs_KCS1+Pi.png"
    
  )
  
  
  for (i in seq_along(pairwise_results)) {
    
    make_de_heatmap(
      
      pairwise_results[[i]],
      
      paste(
        "Differential expression:",
        volcano_titles[i]
      ),
      
      file.path(
        output_dir,
        "05_Heatmaps",
        heatmap_files[i]
      )
      
    )
    
  }
  
  
  ############################################################
  # INTERACTION HEATMAP
  ############################################################
  
  make_de_heatmap(
    
    res_interaction_df,
    
    "Genotype × Pi interaction",
    
    file.path(
      output_dir,
      "05_Heatmaps",
      "Heatmap_Genotype_x_Pi_interaction.png"
    )
    
  )
  
  ############################################################
  # 27. CONDITION-LEVEL VST HEATMAP
  ############################################################
  #
  # Replicates are averaged ONLY for visualization.
  #
  # Conditions:
  #
  # WT+Pi
  # WT-Pi
  # KCS1+Pi
  # KCS1-Pi
  #
  ############################################################
  
  print("======================================================")
  print("CREATING CONDITION-LEVEL HEATMAP")
  print("======================================================")
  
  
  ############################################################
  # CREATE CONDITION MEAN MATRIX
  ############################################################
  
  condition_means <- matrix(
    
    NA_real_,
    
    nrow = nrow(vst_matrix),
    
    ncol = 4
    
  )
  
  
  ############################################################
  # PRESERVE GENE IDs
  ############################################################
  
  rownames(condition_means) <- rownames(
    vst_matrix
  )
  
  
  ############################################################
  # CONDITION NAMES
  ############################################################
  
  colnames(condition_means) <- c(
    
    "WT+Pi",
    "WT-Pi",
    "KCS1+Pi",
    "KCS1-Pi"
    
  )
  
  
  ############################################################
  # IDENTIFY SAMPLE INDICES FOR EACH CONDITION
  ############################################################
  
  condition_indices <- list(
    
    "WT+Pi" = which(
      metadata$Condition == "WT_plus"
    ),
    
    "WT-Pi" = which(
      metadata$Condition == "WT_minus"
    ),
    
    "KCS1+Pi" = which(
      metadata$Condition == "KCS1_plus"
    ),
    
    "KCS1-Pi" = which(
      metadata$Condition == "KCS1_minus"
    )
    
  )
  
  
  ############################################################
  # VERIFY CONDITION INDICES
  ############################################################
  
  print("Condition sample indices:")
  
  for (condition_name in names(condition_indices)) {
    
    print(
      paste(
        condition_name,
        ":",
        paste(
          condition_indices[[condition_name]],
          collapse = ", "
        )
      )
    )
    
  }
  
  
  ############################################################
  # VERIFY THAT EACH CONDITION HAS 3 REPLICATES
  ############################################################
  
  condition_replicate_counts <- sapply(
    
    condition_indices,
    
    length
    
  )
  
  
  if (
    any(
      condition_replicate_counts != 3
    )
  ) {
    
    stop(
      paste0(
        "\nERROR: Condition-level heatmap does not have ",
        "exactly 3 replicates per condition.\n\n",
        paste(
          names(condition_replicate_counts),
          condition_replicate_counts,
          sep = ": ",
          collapse = "\n"
        )
      )
    )
    
  }
  
  
  ############################################################
  # CALCULATE MEAN VST EXPRESSION
  #
  # Replicates are averaged ONLY here for visualization.
  # DESeq2 statistics continue to use all biological replicates.
  ############################################################
  
  for (i in seq_along(condition_indices)) {
    
    condition_means[, i] <- rowMeans(
      
      vst_matrix[
        ,
        condition_indices[[i]],
        drop = FALSE
      ],
      
      na.rm = TRUE
      
    )
    
  }
  
  
  ############################################################
  # REMOVE GENES WITH NON-FINITE VALUES
  ############################################################
  
  finite_genes <- apply(
    
    condition_means,
    
    1,
    
    function(x) {
      
      all(
        is.finite(x)
      )
      
    }
    
  )
  
  
  print(
    paste(
      "Genes with finite values across all four conditions:",
      sum(finite_genes)
    )
  )
  
  
  condition_means <- condition_means[
    
    finite_genes,
    
    ,
    
    drop = FALSE
    
  ]
  
  
  ############################################################
  # CALCULATE VARIANCE ACROSS THE FOUR CONDITIONS
  ############################################################
  
  condition_variances <- matrixStats::rowVars(
    
    condition_means
    
  )
  
  
  ############################################################
  # REMOVE GENES WITH ZERO VARIANCE
  #
  # These genes would become NaN during row scaling.
  ############################################################
  
  variable_genes <- is.finite(
    
    condition_variances
    
  ) &
    
    condition_variances > 0
  
  
  print(
    paste(
      "Genes with variation across conditions:",
      sum(variable_genes)
    )
  )
  
  
  condition_variances <- condition_variances[
    
    variable_genes
    
  ]
  
  
  condition_means <- condition_means[
    
    variable_genes,
    
    ,
    
    drop = FALSE
    
  ]
  
  
  ############################################################
  # SELECT TOP 50 VARIABLE GENES
  #
  # IMPORTANT:
  # Use numeric row indices rather than names().
  # matrixStats::rowVars() does not necessarily preserve
  # row names as names of the returned variance vector.
  ############################################################
  
  n_condition_genes <- min(
    
    50,
    
    length(condition_variances)
    
  )
  
  
  if (
    n_condition_genes < 2
  ) {
    
    stop(
      paste0(
        "\nERROR: Fewer than 2 variable genes are available ",
        "for the condition-level heatmap."
      )
    )
    
  }
  
  
  top_condition_indices <- order(
    
    condition_variances,
    
    decreasing = TRUE
    
  )[
    seq_len(
      n_condition_genes
    )
  ]
  
  
  ############################################################
  # EXTRACT TOP VARIABLE GENES
  ############################################################
  
  condition_heatmap <- condition_means[
    
    top_condition_indices,
    
    ,
    
    drop = FALSE
    
  ]
  
  
  ############################################################
  # REPORT SELECTED GENES
  ############################################################
  
  print(
    paste(
      "Number of genes selected for condition heatmap:",
      nrow(condition_heatmap)
    )
  )
  
  
  ############################################################
  # MAP SYSTEMATIC IDs TO CGD GENE NAMES
  ############################################################
  
  condition_gene_names <- get_gene_names(
    
    rownames(
      condition_heatmap
    )
    
  )
  
  
  ############################################################
  # MAKE GENE NAMES UNIQUE
  ############################################################
  
  rownames(
    condition_heatmap
  ) <- make.unique(
    
    condition_gene_names
    
  )
  
  
  ############################################################
  # ROW-SCALE EXPRESSION
  #
  # Each gene is standardized across the four conditions.
  #
  # This shows relative expression patterns rather than
  # absolute VST expression.
  ############################################################
  
  condition_heatmap_scaled <- t(
    
    scale(
      
      t(
        condition_heatmap
      )
      
    )
    
  )
  
  
  ############################################################
  # REMOVE ANY ROWS THAT BECAME NON-FINITE AFTER SCALING
  ############################################################
  
  finite_scaled_genes <- apply(
    
    condition_heatmap_scaled,
    
    1,
    
    function(x) {
      
      all(
        is.finite(x)
      )
      
    }
    
  )
  
  
  condition_heatmap_scaled <- condition_heatmap_scaled[
    
    finite_scaled_genes,
    
    ,
    
    drop = FALSE
    
  ]
  
  
  ############################################################
  # VERIFY HEATMAP DIMENSIONS
  ############################################################
  
  print("Condition heatmap dimensions:")
  
  print(
    dim(
      condition_heatmap_scaled
    )
  )
  
  
  if (
    nrow(condition_heatmap_scaled) < 2
  ) {
    
    stop(
      "\nERROR: Fewer than 2 genes remain after scaling."
    )
    
  }
  
  
  if (
    ncol(condition_heatmap_scaled) != 4
  ) {
    
    stop(
      "\nERROR: Condition heatmap does not contain exactly 4 conditions."
    )
    
  }
  
  
  ############################################################
  # SAVE CONDITION MEAN MATRIX
  ############################################################
  
  write.csv(
    
    condition_means,
    
    file.path(
      output_dir,
      "05_Heatmaps",
      "Condition_Mean_VST_Expression.csv"
    ),
    
    row.names = TRUE
    
  )
  
  
  ############################################################
  # CREATE HEATMAP
  ############################################################
  
  png(
    
    file.path(
      output_dir,
      "05_Heatmaps",
      "Heatmap_Top_50_Variable_Condition_Means.png"
    ),
    
    width = 1800,
    height = 2400,
    res = 300
    
  )
  
  
  pheatmap(
    
    condition_heatmap_scaled,
    
    scale = "none",
    
    cluster_rows = TRUE,
    
    cluster_cols = FALSE,
    
    show_rownames = TRUE,
    
    show_colnames = TRUE,
    
    fontsize_row = 6,
    
    fontsize_col = 10,
    
    border_color = NA,
    
    main = "Top 50 variable genes: condition means"
    
  )
  
  
  dev.off()
  
  
  ############################################################
  # CONFIRM SUCCESS
  ############################################################
  
  print(
    "SUCCESS: Condition-level VST heatmap created."
  )
  
  print(
    paste(
      "Genes plotted:",
      nrow(condition_heatmap_scaled)
    )
  )
  
  print(
    "Conditions plotted:"
  )
  
  print(
    colnames(condition_heatmap_scaled)
  )
  
  ############################################################
  # 28. PHOSPHATE RESPONSE ANALYSIS
  ############################################################
  #
  # WT response:
  #
  # WT-Pi - WT+Pi
  #
  # KCS1 response:
  #
  # KCS1-Pi - KCS1+Pi
  #
  # The difference between these responses is tested formally
  # by the Genotype × Pi interaction.
  #
  ############################################################
  
  print("======================================================")
  print("PHOSPHATE RESPONSE ANALYSIS")
  print("======================================================")
  
  
  ############################################################
  # GET LOG2 FOLD CHANGES
  ############################################################
  
  WT_response <- res_WT_minus_vs_WT_plus_df$log2FoldChange
  
  KCS1_response <- res_KCS1_minus_vs_KCS1_plus_df$log2FoldChange
  
  
  response_table <- data.frame(
    
    Systematic_ID =
      res_WT_minus_vs_WT_plus_df$Systematic_ID,
    
    Gene_Name =
      res_WT_minus_vs_WT_plus_df$Gene_Name,
    
    WT_response =
      WT_response,
    
    WT_padj =
      res_WT_minus_vs_WT_plus_df$padj,
    
    KCS1_response =
      KCS1_response[
        match(
          res_WT_minus_vs_WT_plus_df$Systematic_ID,
          res_KCS1_minus_vs_KCS1_plus_df$Systematic_ID
        )
      ],
    
    KCS1_padj =
      res_KCS1_minus_vs_KCS1_plus_df$padj[
        match(
          res_WT_minus_vs_WT_plus_df$Systematic_ID,
          res_KCS1_minus_vs_KCS1_plus_df$Systematic_ID
        )
      ],
    
    Interaction_log2FC =
      res_interaction_df$log2FoldChange[
        match(
          res_WT_minus_vs_WT_plus_df$Systematic_ID,
          res_interaction_df$Systematic_ID
        )
      ],
    
    Interaction_padj =
      res_interaction_df$padj[
        match(
          res_WT_minus_vs_WT_plus_df$Systematic_ID,
          res_interaction_df$Systematic_ID
        )
      ],
    
    stringsAsFactors = FALSE
    
  )
  
  
  ############################################################
  # SAVE RESPONSE TABLE
  ############################################################
  
  write.csv(
    
    response_table,
    
    file.path(
      output_dir,
      "09_Response_analysis",
      "Phosphate_Response_Comparison.csv"
    ),
    
    row.names = FALSE
    
  )
  
  
  ############################################################
  # SIGNIFICANT INTERACTION GENES
  ############################################################
  
  significant_response_genes <- response_table %>%
    
    filter(
      !is.na(Interaction_padj),
      Interaction_padj < 0.05
    ) %>%
    
    arrange(
      Interaction_padj
    )
  
  
  write.csv(
    
    significant_response_genes,
    
    file.path(
      output_dir,
      "09_Response_analysis",
      "Significant_Phosphate_Response_Interaction_Genes.csv"
    ),
    
    row.names = FALSE
    
  )
  
  
  ############################################################
  # TOP 50 INTERACTION GENES
  ############################################################
  
  top_interaction_genes <- significant_response_genes %>%
    
    head(50)
  
  
  write.csv(
    
    top_interaction_genes,
    
    file.path(
      output_dir,
      "09_Response_analysis",
      "Top_50_Interaction_Genes.csv"
    ),
    
    row.names = FALSE
    
  )
  
  
  ############################################################
  # 29. PHOSPHATE RESPONSE HEATMAP
  ############################################################
  
  if (
    nrow(top_interaction_genes) >= 2
  ) {
    
    
    response_gene_ids <- intersect(
      
      top_interaction_genes$Systematic_ID,
      
      rownames(
        vst_matrix
      )
      
    )
    
    
    response_matrix <- condition_means[
      
      response_gene_ids,
      
      ,
      
      drop = FALSE
      
    ]
    
    
    rownames(
      response_matrix
    ) <- make.unique(
      
      get_gene_names(
        rownames(
          response_matrix
        )
      )
      
    )
    
    
    response_scaled <- t(
      
      scale(
        t(
          response_matrix
        )
      )
      
    )
    
    
    png(
      
      file.path(
        output_dir,
        "09_Response_analysis",
        "Heatmap_Top_Interaction_Genes.png"
      ),
      
      width = 1800,
      height = 2600,
      res = 300
      
    )
    
    
    pheatmap(
      
      response_scaled,
      
      scale = "none",
      
      cluster_rows = TRUE,
      
      cluster_cols = FALSE,
      
      show_rownames = TRUE,
      
      fontsize_row = 7,
      
      fontsize_col = 10,
      
      main =
        "Top genes with Genotype × Pi interaction"
      
    )
    
    
    dev.off()
    
  }
  
  
  ############################################################
  # 30. PHOSPHATE RESPONSE SCATTERPLOT
  ############################################################
  
  response_plot_df <- response_table %>%
    
    filter(
      !is.na(WT_response),
      !is.na(KCS1_response)
    )
  
  
  if (
    nrow(response_plot_df) > 0
  ) {
    
    
    p_response_scatter <- ggplot(
      
      response_plot_df,
      
      aes(
        x = WT_response,
        y = KCS1_response
      )
      
    ) +
      
      geom_point(
        alpha = 0.5,
        size = 2
      ) +
      
      geom_abline(
        slope = 1,
        intercept = 0,
        linetype = "dashed"
      ) +
      
      geom_vline(
        xintercept = 0,
        linetype = "dotted"
      ) +
      
      geom_hline(
        yintercept = 0,
        linetype = "dotted"
      ) +
      
      theme_bw() +
      
      labs(
        
        title =
          "WT vs KCS1 phosphate-starvation responses",
        
        x =
          "WT response: WT-Pi vs WT+Pi",
        
        y =
          "KCS1 response: KCS1-Pi vs KCS1+Pi"
        
      )
    
    
    ggsave(
      
      file.path(
        output_dir,
        "09_Response_analysis",
        "WT_vs_KCS1_Phosphate_Response.png"
      ),
      
      p_response_scatter,
      
      width = 9,
      height = 8,
      dpi = 300
      
    )
    
  }
  
  
  ############################################################
  # 31. PHOSPHATE RESPONSE PLOT
  #
  # Only significant interaction genes are shown.
  #
  ############################################################
  
  if (
    nrow(significant_response_genes) >= 1
  ) {
    
    
    response_long <- data.frame(
      
      Gene = rep(
        significant_response_genes$Gene_Name,
        2
      ),
      
      Genotype = rep(
        c(
          "WT",
          "KCS1"
        ),
        
        each =
          nrow(
            significant_response_genes
          )
      ),
      
      Response = c(
        
        significant_response_genes$WT_response,
        
        significant_response_genes$KCS1_response
        
      ),
      
      stringsAsFactors = FALSE
      
    )
    
    
    response_plot <- ggplot(
      
      response_long,
      
      aes(
        x = Response,
        y = reorder(
          Gene,
          Response
        ),
        shape = Genotype
      )
      
    ) +
      
      geom_point(
        size = 3
      ) +
      
      geom_vline(
        xintercept = 0,
        linetype = "dashed"
      ) +
      
      theme_bw() +
      
      labs(
        
        title =
          "Phosphate-response differences in significant interaction genes",
        
        x =
          "Change in VST expression: -Pi minus +Pi",
        
        y =
          "Gene"
        
      )
    
    
    ##########################################################
    # CAP HEIGHT TO PREVENT ggsave ERROR
    ##########################################################
    
    response_height <- min(
      
      max(
        6,
        0.25 *
          nrow(
            significant_response_genes
          )
      ),
      
      48
      
    )
    
    
    ggsave(
      
      file.path(
        output_dir,
        "09_Response_analysis",
        "Phosphate_response_plot.png"
      ),
      
      response_plot,
      
      width = 10,
      
      height = response_height,
      
      dpi = 300
      
    )
    
  }
  
  
  ############################################################
  # 32. SUMMARY TABLE
  ############################################################
  
  summary_table <- data.frame(
    
    Comparison = names(
      all_results
    ),
    
    Total_genes_tested = sapply(
      
      all_results,
      
      nrow
      
    ),
    
    Significant_padj_0.05 = sapply(
      
      all_results,
      
      function(df) {
        
        sum(
          !is.na(df$padj) &
            df$padj < 0.05
        )
        
      }
      
    ),
    
    Significant_padj_0.05_log2FC_1 = sapply(
      
      all_results,
      
      function(df) {
        
        sum(
          !is.na(df$padj) &
            df$padj < 0.05 &
            !is.na(df$log2FoldChange) &
            abs(df$log2FoldChange) >= 1
        )
        
      }
      
    ),
    
    stringsAsFactors = FALSE
    
  )
  
  
  write.csv(
    
    summary_table,
    
    file.path(
      output_dir,
      "02_DE_results",
      "DE_summary.csv"
    ),
    
    row.names = FALSE
    
  )
  
  
  ############################################################
  # 33. COMPARISON SUMMARY PRINT
  ############################################################
  
  print("======================================================")
  print("SIGNIFICANT GENE COUNTS")
  print("======================================================")
  
  
  print(
    summary_table
  )
  
  
  ############################################################
  # 34. SAVE VST OBJECT
  ############################################################
  
  saveRDS(
    
    vsd,
    
    file.path(
      output_dir,
      "01_QC",
      "VST_object.rds"
    )
    
  )
  
  
  ############################################################
  # 35. FINAL SUMMARY
  ############################################################
  
  print("")
  print("======================================================")
  print("RNA-seq ANALYSIS COMPLETE")
  print("======================================================")
  
  
  print(
    paste(
      "Total genes analyzed:",
      nrow(dds)
    )
  )
  
  
  print(
    paste(
      "WT+Pi replicates:",
      sum(
        metadata$Condition == "WT_plus"
      )
    )
  )
  
  
  print(
    paste(
      "WT-Pi replicates:",
      sum(
        metadata$Condition == "WT_minus"
      )
    )
  )
  
  
  print(
    paste(
      "KCS1+Pi replicates:",
      sum(
        metadata$Condition == "KCS1_plus"
      )
    )
  )
  
  
  print(
    paste(
      "KCS1-Pi replicates:",
      sum(
        metadata$Condition == "KCS1_minus"
      )
    )
  )
  
  
  print("")
  
  print(
    paste(
      "Genes with significant Genotype × Pi interaction:",
      sum(
        !is.na(
          res_interaction_df$padj
        ) &
          res_interaction_df$padj < 0.05
      )
    )
  )
  
  
  print("")
  
  print("Output directory:")
  
  print(
    normalizePath(
      output_dir
    )
  )
  
  
  print("")
  
  print("======================================================")
  print("PAIRWISE COMPARISONS GENERATED")
  print("======================================================")
  
  print(
    "1. WT-Pi vs WT+Pi"
  )
  
  print(
    "2. KCS1+Pi vs WT+Pi"
  )
  
  print(
    "3. KCS1-Pi vs WT+Pi"
  )
  
  print(
    "4. KCS1+Pi vs WT-Pi"
  )
  
  print(
    "5. KCS1-Pi vs WT-Pi"
  )
  
  print(
    "6. KCS1-Pi vs KCS1+Pi"
  )
  
  print(
    "7. Genotype × Pi interaction"
  )
  
  
  print("")
  print("======================================================")
  print("END OF ANALYSIS")
  print("======================================================")





  ############################################################
  # 27. CONDITION-LEVEL VST HEATMAP
  ############################################################
  #
  # Replicates are averaged ONLY for visualization.
  #
  # Conditions:
  #
  # WT+Pi
  # WT-Pi
  # KCS1+Pi
  # KCS1-Pi
  #
  ############################################################
  
  print("======================================================")
  print("CREATING CONDITION-LEVEL HEATMAP")
  print("======================================================")
  
  
  ############################################################
  # CREATE CONDITION MEAN MATRIX
  ############################################################
  
  condition_means <- matrix(
    
    NA_real_,
    
    nrow = nrow(vst_matrix),
    
    ncol = 4
    
  )
  
  
  ############################################################
  # PRESERVE GENE IDs
  ############################################################
  
  rownames(condition_means) <- rownames(
    vst_matrix
  )
  
  
  ############################################################
  # CONDITION NAMES
  ############################################################
  
  colnames(condition_means) <- c(
    
    "WT+Pi",
    "WT-Pi",
    "KCS1+Pi",
    "KCS1-Pi"
    
  )
  
  
  ############################################################
  # IDENTIFY SAMPLE INDICES FOR EACH CONDITION
  ############################################################
  
  condition_indices <- list(
    
    "WT+Pi" = which(
      metadata$Condition == "WT_plus"
    ),
    
    "WT-Pi" = which(
      metadata$Condition == "WT_minus"
    ),
    
    "KCS1+Pi" = which(
      metadata$Condition == "KCS1_plus"
    ),
    
    "KCS1-Pi" = which(
      metadata$Condition == "KCS1_minus"
    )
    
  )
  
  
  ############################################################
  # VERIFY CONDITION INDICES
  ############################################################
  
  print("Condition sample indices:")
  
  for (condition_name in names(condition_indices)) {
    
    print(
      paste(
        condition_name,
        ":",
        paste(
          condition_indices[[condition_name]],
          collapse = ", "
        )
      )
    )
    
  }
  
  
  ############################################################
  # VERIFY THAT EACH CONDITION HAS 3 REPLICATES
  ############################################################
  
  condition_replicate_counts <- sapply(
    
    condition_indices,
    
    length
    
  )
  
  
  if (
    any(
      condition_replicate_counts != 3
    )
  ) {
    
    stop(
      paste0(
        "\nERROR: Condition-level heatmap does not have ",
        "exactly 3 replicates per condition.\n\n",
        paste(
          names(condition_replicate_counts),
          condition_replicate_counts,
          sep = ": ",
          collapse = "\n"
        )
      )
    )
    
  }
  
  
  ############################################################
  # CALCULATE MEAN VST EXPRESSION
  #
  # Replicates are averaged ONLY here for visualization.
  # DESeq2 statistics continue to use all biological replicates.
  ############################################################
  
  for (i in seq_along(condition_indices)) {
    
    condition_means[, i] <- rowMeans(
      
      vst_matrix[
        ,
        condition_indices[[i]],
        drop = FALSE
      ],
      
      na.rm = TRUE
      
    )
    
  }
  
  
  ############################################################
  # REMOVE GENES WITH NON-FINITE VALUES
  ############################################################
  
  finite_genes <- apply(
    
    condition_means,
    
    1,
    
    function(x) {
      
      all(
        is.finite(x)
      )
      
    }
    
  )
  
  
  print(
    paste(
      "Genes with finite values across all four conditions:",
      sum(finite_genes)
    )
  )
  
  
  condition_means <- condition_means[
    
    finite_genes,
    
    ,
    
    drop = FALSE
    
  ]
  
  
  ############################################################
  # CALCULATE VARIANCE ACROSS THE FOUR CONDITIONS
  ############################################################
  
  condition_variances <- matrixStats::rowVars(
    
    condition_means
    
  )
  
  
  ############################################################
  # REMOVE GENES WITH ZERO VARIANCE
  #
  # These genes would become NaN during row scaling.
  ############################################################
  
  variable_genes <- is.finite(
    
    condition_variances
    
  ) &
    
    condition_variances > 0
  
  
  print(
    paste(
      "Genes with variation across conditions:",
      sum(variable_genes)
    )
  )
  
  
  condition_variances <- condition_variances[
    
    variable_genes
    
  ]
  
  
  condition_means <- condition_means[
    
    variable_genes,
    
    ,
    
    drop = FALSE
    
  ]
  
  
  ############################################################
  # SELECT TOP 50 VARIABLE GENES
  #
  # IMPORTANT:
  # Use numeric row indices rather than names().
  # matrixStats::rowVars() does not necessarily preserve
  # row names as names of the returned variance vector.
  ############################################################
  
  n_condition_genes <- min(
    
    50,
    
    length(condition_variances)
    
  )
  
  
  if (
    n_condition_genes < 2
  ) {
    
    stop(
      paste0(
        "\nERROR: Fewer than 2 variable genes are available ",
        "for the condition-level heatmap."
      )
    )
    
  }
  
  
  top_condition_indices <- order(
    
    condition_variances,
    
    decreasing = TRUE
    
  )[
    seq_len(
      n_condition_genes
    )
  ]
  
  
  ############################################################
  # EXTRACT TOP VARIABLE GENES
  ############################################################
  
  condition_heatmap <- condition_means[
    
    top_condition_indices,
    
    ,
    
    drop = FALSE
    
  ]
  
  
  ############################################################
  # REPORT SELECTED GENES
  ############################################################
  
  print(
    paste(
      "Number of genes selected for condition heatmap:",
      nrow(condition_heatmap)
    )
  )
  
  
  ############################################################
  # MAP SYSTEMATIC IDs TO CGD GENE NAMES
  ############################################################
  
  condition_gene_names <- get_gene_names(
    
    rownames(
      condition_heatmap
    )
    
  )
  
  
  ############################################################
  # MAKE GENE NAMES UNIQUE
  ############################################################
  
  rownames(
    condition_heatmap
  ) <- make.unique(
    
    condition_gene_names
    
  )
  
  
  ############################################################
  # ROW-SCALE EXPRESSION
  #
  # Each gene is standardized across the four conditions.
  #
  # This shows relative expression patterns rather than
  # absolute VST expression.
  ############################################################
  
  condition_heatmap_scaled <- t(
    
    scale(
      
      t(
        condition_heatmap
      )
      
    )
    
  )
  
  
  ############################################################
  # REMOVE ANY ROWS THAT BECAME NON-FINITE AFTER SCALING
  ############################################################
  
  finite_scaled_genes <- apply(
    
    condition_heatmap_scaled,
    
    1,
    
    function(x) {
      
      all(
        is.finite(x)
      )
      
    }
    
  )
  
  
  condition_heatmap_scaled <- condition_heatmap_scaled[
    
    finite_scaled_genes,
    
    ,
    
    drop = FALSE
    
  ]
  
  
  ############################################################
  # VERIFY HEATMAP DIMENSIONS
  ############################################################
  
  print("Condition heatmap dimensions:")
  
  print(
    dim(
      condition_heatmap_scaled
    )
  )
  
  
  if (
    nrow(condition_heatmap_scaled) < 2
  ) {
    
    stop(
      "\nERROR: Fewer than 2 genes remain after scaling."
    )
    
  }
  
  
  if (
    ncol(condition_heatmap_scaled) != 4
  ) {
    
    stop(
      "\nERROR: Condition heatmap does not contain exactly 4 conditions."
    )
    
  }
  
  
  ############################################################
  # SAVE CONDITION MEAN MATRIX
  ############################################################
  
  write.csv(
    
    condition_means,
    
    file.path(
      output_dir,
      "05_Heatmaps",
      "Condition_Mean_VST_Expression.csv"
    ),
    
    row.names = TRUE
    
  )
  
  
  ############################################################
  # CREATE HEATMAP
  ############################################################
  
  png(
    
    file.path(
      output_dir,
      "05_Heatmaps",
      "Heatmap_Top_50_Variable_Condition_Means.png"
    ),
    
    width = 1800,
    height = 2400,
    res = 300
    
  )
  
  
  pheatmap(
    
    condition_heatmap_scaled,
    
    scale = "none",
    
    cluster_rows = TRUE,
    
    cluster_cols = FALSE,
    
    show_rownames = TRUE,
    
    show_colnames = TRUE,
    
    fontsize_row = 6,
    
    fontsize_col = 10,
    
    border_color = NA,
    
    main = "Top 50 variable genes: condition means"
    
  )
  
  
  dev.off()
  
  
  ############################################################
  # CONFIRM SUCCESS
  ############################################################
  
  print(
    "SUCCESS: Condition-level VST heatmap created."
  )
  
  print(
    paste(
      "Genes plotted:",
      nrow(condition_heatmap_scaled)
    )
  )
  
  print(
    "Conditions plotted:"
  )
  
  print(
    colnames(condition_heatmap_scaled)
  )
  
  ############################################################
  # 28. PHOSPHATE RESPONSE ANALYSIS
  ############################################################
  #
  # WT response:
  #
  # WT-Pi - WT+Pi
  #
  # KCS1 response:
  #
  # KCS1-Pi - KCS1+Pi
  #
  # The difference between these responses is tested formally
  # by the Genotype × Pi interaction.
  #
  ############################################################
  
  print("======================================================")
  print("PHOSPHATE RESPONSE ANALYSIS")
  print("======================================================")
  
  
  ############################################################
  # GET LOG2 FOLD CHANGES
  ############################################################
  
  WT_response <- res_WT_minus_vs_WT_plus_df$log2FoldChange
  
  KCS1_response <- res_KCS1_minus_vs_KCS1_plus_df$log2FoldChange
  
  
  response_table <- data.frame(
    
    Systematic_ID =
      res_WT_minus_vs_WT_plus_df$Systematic_ID,
    
    Gene_Name =
      res_WT_minus_vs_WT_plus_df$Gene_Name,
    
    WT_response =
      WT_response,
    
    WT_padj =
      res_WT_minus_vs_WT_plus_df$padj,
    
    KCS1_response =
      KCS1_response[
        match(
          res_WT_minus_vs_WT_plus_df$Systematic_ID,
          res_KCS1_minus_vs_KCS1_plus_df$Systematic_ID
        )
      ],
    
    KCS1_padj =
      res_KCS1_minus_vs_KCS1_plus_df$padj[
        match(
          res_WT_minus_vs_WT_plus_df$Systematic_ID,
          res_KCS1_minus_vs_KCS1_plus_df$Systematic_ID
        )
      ],
    
    Interaction_log2FC =
      res_interaction_df$log2FoldChange[
        match(
          res_WT_minus_vs_WT_plus_df$Systematic_ID,
          res_interaction_df$Systematic_ID
        )
      ],
    
    Interaction_padj =
      res_interaction_df$padj[
        match(
          res_WT_minus_vs_WT_plus_df$Systematic_ID,
          res_interaction_df$Systematic_ID
        )
      ],
    
    stringsAsFactors = FALSE
    
  )
  
  
  ############################################################
  # SAVE RESPONSE TABLE
  ############################################################
  
  write.csv(
    
    response_table,
    
    file.path(
      output_dir,
      "09_Response_analysis",
      "Phosphate_Response_Comparison.csv"
    ),
    
    row.names = FALSE
    
  )
  
  
  ############################################################
  # SIGNIFICANT INTERACTION GENES
  ############################################################
  
  significant_response_genes <- response_table %>%
    
    filter(
      !is.na(Interaction_padj),
      Interaction_padj < 0.05
    ) %>%
    
    arrange(
      Interaction_padj
    )
  
  
  write.csv(
    
    significant_response_genes,
    
    file.path(
      output_dir,
      "09_Response_analysis",
      "Significant_Phosphate_Response_Interaction_Genes.csv"
    ),
    
    row.names = FALSE
    
  )
  
  
  ############################################################
  # TOP 50 INTERACTION GENES
  ############################################################
  
  top_interaction_genes <- significant_response_genes %>%
    
    head(50)
  
  
  write.csv(
    
    top_interaction_genes,
    
    file.path(
      output_dir,
      "09_Response_analysis",
      "Top_50_Interaction_Genes.csv"
    ),
    
    row.names = FALSE
    
  )
  
  
  ############################################################
  # 29. PHOSPHATE RESPONSE HEATMAP
  ############################################################
  
  if (
    nrow(top_interaction_genes) >= 2
  ) {
    
    
    response_gene_ids <- intersect(
      
      top_interaction_genes$Systematic_ID,
      
      rownames(
        vst_matrix
      )
      
    )
    
    
    response_matrix <- condition_means[
      
      response_gene_ids,
      
      ,
      
      drop = FALSE
      
    ]
    
    
    rownames(
      response_matrix
    ) <- make.unique(
      
      get_gene_names(
        rownames(
          response_matrix
        )
      )
      
    )
    
    
    response_scaled <- t(
      
      scale(
        t(
          response_matrix
        )
      )
      
    )
    
    
    png(
      
      file.path(
        output_dir,
        "09_Response_analysis",
        "Heatmap_Top_Interaction_Genes.png"
      ),
      
      width = 1800,
      height = 2600,
      res = 300
      
    )
    
    
    pheatmap(
      
      response_scaled,
      
      scale = "none",
      
      cluster_rows = TRUE,
      
      cluster_cols = FALSE,
      
      show_rownames = TRUE,
      
      fontsize_row = 7,
      
      fontsize_col = 10,
      
      main =
        "Top genes with Genotype × Pi interaction"
      
    )
    
    
    dev.off()
    
  }
  
  
  ############################################################
  # 30. PHOSPHATE RESPONSE SCATTERPLOT
  ############################################################
  
  response_plot_df <- response_table %>%
    
    filter(
      !is.na(WT_response),
      !is.na(KCS1_response)
    )
  
  
  if (
    nrow(response_plot_df) > 0
  ) {
    
    
    p_response_scatter <- ggplot(
      
      response_plot_df,
      
      aes(
        x = WT_response,
        y = KCS1_response
      )
      
    ) +
      
      geom_point(
        alpha = 0.5,
        size = 2
      ) +
      
      geom_abline(
        slope = 1,
        intercept = 0,
        linetype = "dashed"
      ) +
      
      geom_vline(
        xintercept = 0,
        linetype = "dotted"
      ) +
      
      geom_hline(
        yintercept = 0,
        linetype = "dotted"
      ) +
      
      theme_bw() +
      
      labs(
        
        title =
          "WT vs KCS1 phosphate-starvation responses",
        
        x =
          "WT response: WT-Pi vs WT+Pi",
        
        y =
          "KCS1 response: KCS1-Pi vs KCS1+Pi"
        
      )
    
    
    ggsave(
      
      file.path(
        output_dir,
        "09_Response_analysis",
        "WT_vs_KCS1_Phosphate_Response.png"
      ),
      
      p_response_scatter,
      
      width = 9,
      height = 8,
      dpi = 300
      
    )
    
  }
  
  
  ############################################################
  # 31. PHOSPHATE RESPONSE PLOT
  #
  # Only significant interaction genes are shown.
  #
  ############################################################
  
  if (
    nrow(significant_response_genes) >= 1
  ) {
    
    
    response_long <- data.frame(
      
      Gene = rep(
        significant_response_genes$Gene_Name,
        2
      ),
      
      Genotype = rep(
        c(
          "WT",
          "KCS1"
        ),
        
        each =
          nrow(
            significant_response_genes
          )
      ),
      
      Response = c(
        
        significant_response_genes$WT_response,
        
        significant_response_genes$KCS1_response
        
      ),
      
      stringsAsFactors = FALSE
      
    )
    
    
    response_plot <- ggplot(
      
      response_long,
      
      aes(
        x = Response,
        y = reorder(
          Gene,
          Response
        ),
        shape = Genotype
      )
      
    ) +
      
      geom_point(
        size = 3
      ) +
      
      geom_vline(
        xintercept = 0,
        linetype = "dashed"
      ) +
      
      theme_bw() +
      
      labs(
        
        title =
          "Phosphate-response differences in significant interaction genes",
        
        x =
          "Change in VST expression: -Pi minus +Pi",
        
        y =
          "Gene"
        
      )
    
    
    ##########################################################
    # CAP HEIGHT TO PREVENT ggsave ERROR
    ##########################################################
    
    response_height <- min(
      
      max(
        6,
        0.25 *
          nrow(
            significant_response_genes
          )
      ),
      
      48
      
    )
    
    
    ggsave(
      
      file.path(
        output_dir,
        "09_Response_analysis",
        "Phosphate_response_plot.png"
      ),
      
      response_plot,
      
      width = 10,
      
      height = response_height,
      
      dpi = 300
      
    )
    
  }
  
  
  ############################################################
  # 32. SUMMARY TABLE
  ############################################################
  
  summary_table <- data.frame(
    
    Comparison = names(
      all_results
    ),
    
    Total_genes_tested = sapply(
      
      all_results,
      
      nrow
      
    ),
    
    Significant_padj_0.05 = sapply(
      
      all_results,
      
      function(df) {
        
        sum(
          !is.na(df$padj) &
            df$padj < 0.05
        )
        
      }
      
    ),
    
    Significant_padj_0.05_log2FC_1 = sapply(
      
      all_results,
      
      function(df) {
        
        sum(
          !is.na(df$padj) &
            df$padj < 0.05 &
            !is.na(df$log2FoldChange) &
            abs(df$log2FoldChange) >= 1
        )
        
      }
      
    ),
    
    stringsAsFactors = FALSE
    
  )
  
  
  write.csv(
    
    summary_table,
    
    file.path(
      output_dir,
      "02_DE_results",
      "DE_summary.csv"
    ),
    
    row.names = FALSE
    
  )
  
  
  ############################################################
  # 33. COMPARISON SUMMARY PRINT
  ############################################################
  
  print("======================================================")
  print("SIGNIFICANT GENE COUNTS")
  print("======================================================")
  
  
  print(
    summary_table
  )
  
  
  ############################################################
  # 34. SAVE VST OBJECT
  ############################################################
  
  saveRDS(
    
    vsd,
    
    file.path(
      output_dir,
      "01_QC",
      "VST_object.rds"
    )
    
  )
  
  
  ############################################################
  # 35. FINAL SUMMARY
  ############################################################
  
  print("")
  print("======================================================")
  print("RNA-seq ANALYSIS COMPLETE")
  print("======================================================")
  
  
  print(
    paste(
      "Total genes analyzed:",
      nrow(dds)
    )
  )
  
  
  print(
    paste(
      "WT+Pi replicates:",
      sum(
        metadata$Condition == "WT_plus"
      )
    )
  )
  
  
  print(
    paste(
      "WT-Pi replicates:",
      sum(
        metadata$Condition == "WT_minus"
      )
    )
  )
  
  
  print(
    paste(
      "KCS1+Pi replicates:",
      sum(
        metadata$Condition == "KCS1_plus"
      )
    )
  )
  
  
  print(
    paste(
      "KCS1-Pi replicates:",
      sum(
        metadata$Condition == "KCS1_minus"
      )
    )
  )
  
  
  print("")
  
  print(
    paste(
      "Genes with significant Genotype × Pi interaction:",
      sum(
        !is.na(
          res_interaction_df$padj
        ) &
          res_interaction_df$padj < 0.05
      )
    )
  )
  
  
  print("")
  
  print("Output directory:")
  
  print(
    normalizePath(
      output_dir
    )
  )
  
  
  print("")
  
  print("======================================================")
  print("PAIRWISE COMPARISONS GENERATED")
  print("======================================================")
  
  print(
    "1. WT-Pi vs WT+Pi"
  )
  
  print(
    "2. KCS1+Pi vs WT+Pi"
  )
  
  print(
    "3. KCS1-Pi vs WT+Pi"
  )
  
  print(
    "4. KCS1+Pi vs WT-Pi"
  )
  
  print(
    "5. KCS1-Pi vs WT-Pi"
  )
  
  print(
    "6. KCS1-Pi vs KCS1+Pi"
  )
  
  print(
    "7. Genotype × Pi interaction"
  )
  
  
  print("")
  print("======================================================")
  print("END OF ANALYSIS")
  print("======================================================")


To double check p-values
 summary(res_df$padj)
  min(res_df$padj, na.rm = TRUE)











full script:
############################################################
# CANDIDA ALBICANS RNA-seq ANALYSIS
#
# EXPERIMENT:
#
# WT    +Pi     3 biological replicates
# WT    -Pi     3 biological replicates
# VIP1  +Pi     3 biological replicates
# VIP1  -Pi     3 biological replicates
#
# 12 samples total
#
# SAMPLE ORDER AFTER SORTING:
#
# 1-3     WT +Pi
# 4-6     WT -Pi
# 7-9     KCS1 +Pi
# 10-12   KCS1 -Pi
#
# ANALYSIS:
#
# - Remove tRNA/rRNA features
# - Map systematic IDs to CGD gene names
# - DESeq2 using all 12 biological replicates
# - PCA
# - sample correlation
# - all 6 possible pairwise comparisons
# - volcano plots for all comparisons
# - MA plots for all comparisons
# - heatmaps for all comparisons
# - condition-level VST heatmap
# - Genotype x Pi interaction analysis
# - phosphate-response analysis
#
############################################################


############################################################
# 0. LOAD REQUIRED PACKAGES
############################################################

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

library(DESeq2)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(dplyr)
library(tidyr)
library(tibble)
library(matrixStats)


############################################################
# 1. OUTPUT DIRECTORIES
############################################################

output_dir <- "RNAseq_DE_Results"

dirs <- c(
  "01_QC",
  "02_DE_results",
  "03_Volcano",
  "04_MA_plots",
  "05_Heatmaps",
  "06_PCA",
  "07_Sample_correlations",
  "08_Significant_gene_lists",
  "09_Response_analysis"
)

dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

for (d in dirs) {
  dir.create(
    file.path(output_dir, d),
    showWarnings = FALSE,
    recursive = TRUE
  )
}


############################################################
# 2. LOAD CGD GENE ANNOTATION
############################################################

print("======================================================")
print("LOADING CGD GENE ANNOTATION")
print("======================================================")

if (!file.exists("cgd_features.tab")) {
  
  stop(
    paste0(
      "\nERROR: cgd_features.tab was not found.\n",
      "Place cgd_features.tab in the same working directory ",
      "as this R script.\n"
    )
  )
  
}

cgd_raw <- read.table(
  "cgd_features.tab",
  header = FALSE,
  sep = "\t",
  quote = "",
  comment.char = "!",
  fill = TRUE,
  stringsAsFactors = FALSE
)


############################################################
# CGD ANNOTATION STRUCTURE
#
# Column 1 = systematic ID
# Column 2 = gene name
# Column 4 = feature type
############################################################

mapping_table <- data.frame(
  
  systematic_id = trimws(
    as.character(cgd_raw[, 1])
  ),
  
  gene_symbol = trimws(
    as.character(cgd_raw[, 2])
  ),
  
  feature_type = trimws(
    as.character(cgd_raw[, 4])
  ),
  
  stringsAsFactors = FALSE
  
)


############################################################
# CLEAN GENE NAMES
############################################################

mapping_table$gene_symbol[
  is.na(mapping_table$gene_symbol) |
    mapping_table$gene_symbol == "" |
    mapping_table$gene_symbol ==
    mapping_table$systematic_id
] <- NA


############################################################
# NORMALIZE IDS
############################################################

normalize_id <- function(x) {
  
  x <- trimws(
    as.character(x)
  )
  
  x <- toupper(
    gsub(
      "[_\\-]",
      "",
      x
    )
  )
  
  return(x)
  
}


mapping_table$id_norm <- normalize_id(
  mapping_table$systematic_id
)


############################################################
# REMOVE DUPLICATE CGD ENTRIES
#
# Prefer entries that actually have gene names.
############################################################

mapping_table <- mapping_table[
  order(
    is.na(mapping_table$gene_symbol)
  ),
  ,
  drop = FALSE
]

mapping_table <- mapping_table[
  !duplicated(mapping_table$id_norm),
  ,
  drop = FALSE
]


############################################################
# ANNOTATION SUMMARY
############################################################

print(
  paste(
    "Total CGD features:",
    nrow(mapping_table)
  )
)

print(
  paste(
    "CGD features with gene names:",
    sum(
      !is.na(mapping_table$gene_symbol)
    )
  )
)


############################################################
# 3. LOAD FEATURECOUNTS MATRIX
############################################################

print("======================================================")
print("LOADING FEATURECOUNTS MATRIX")
print("======================================================")

if (!file.exists("candida_counts_matrix.txt")) {
  
  stop(
    paste0(
      "\nERROR: candida_counts_matrix.txt was not found.\n",
      "Place it in the same working directory as this script.\n"
    )
  )
  
}


############################################################
# ROBUST FEATURECOUNTS READER
#
# Avoids:
# "more columns than column names"
############################################################

read_featurecounts <- function(file) {
  
  lines <- readLines(
    file,
    warn = FALSE
  )
  
  lines <- lines[
    nzchar(
      trimws(lines)
    )
  ]
  
  lines <- lines[
    !grepl(
      "^#",
      trimws(lines)
    )
  ]
  
  if (length(lines) < 2) {
    stop(
      "ERROR: FeatureCounts file does not contain enough data."
    )
  }
  
  
  ##########################################################
  # HEADER
  ##########################################################
  
  header <- strsplit(
    lines[1],
    "\t",
    fixed = TRUE
  )[[1]]
  
  header[1] <- sub(
    "^\ufeff",
    "",
    header[1]
  )
  
  
  if (length(header) != 18) {
    
    stop(
      paste0(
        "\nERROR: Expected 18 columns in featureCounts header.\n",
        "Found ",
        length(header),
        " columns.\n\n",
        "Header detected:\n",
        paste(
          header,
          collapse = "\n"
        )
      )
    )
    
  }
  
  
  ##########################################################
  # DATA
  ##########################################################
  
  data_lines <- lines[-1]
  
  split_lines <- strsplit(
    data_lines,
    "\t",
    fixed = TRUE
  )
  
  
  ##########################################################
  # REMOVE TRAILING EMPTY FIELDS
  ##########################################################
  
  split_lines <- lapply(
    split_lines,
    function(x) {
      
      while (
        length(x) > 0 &&
        trimws(
          x[length(x)]
        ) == ""
      ) {
        
        x <- x[
          -length(x)
        ]
        
      }
      
      x
      
    }
  )
  
  
  field_counts <- lengths(
    split_lines
  )
  
  
  bad_rows <- which(
    field_counts != 18
  )
  
  
  if (length(bad_rows) > 0) {
    
    first_bad <- bad_rows[1]
    
    stop(
      paste0(
        "\nERROR: FeatureCounts contains a row with ",
        field_counts[first_bad],
        " columns instead of 18.\n\n",
        "Problem occurs at data row ",
        first_bad,
        ".\n\n",
        "Line:\n",
        data_lines[first_bad]
      )
    )
    
  }
  
  
  matrix_data <- do.call(
    rbind,
    split_lines
  )
  
  
  counts_df <- as.data.frame(
    matrix_data,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  colnames(counts_df) <- header
  
  return(counts_df)
  
}


counts_raw <- read_featurecounts(
  "candida_counts_matrix.txt"
)


print(
  paste(
    "Total columns detected:",
    ncol(counts_raw)
  )
)

print(
  "Column names detected:"
)

print(
  colnames(counts_raw)
)


############################################################
# 4. IDENTIFY AND REORDER SAMPLE COLUMNS
############################################################

print("======================================================")
print("IDENTIFYING RNA-seq SAMPLE COLUMNS")
print("======================================================")


############################################################
# FEATURECOUNTS STRUCTURE
#
# 1 = Geneid
# 2 = Chr
# 3 = Start
# 4 = End
# 5 = Strand
# 6 = Length
# 7-18 = samples
############################################################

sample_columns_original <- 7:18

original_sample_names <- colnames(
  counts_raw
)[sample_columns_original]


print(
  "Original sample order in featureCounts file:"
)

print(
  original_sample_names
)


############################################################
# EXTRACT N44VRL SAMPLE NUMBERS
############################################################

sample_numbers <- suppressWarnings(
  as.numeric(
    sub(
      ".*N44VRL_([0-9]+)_.*",
      "\\1",
      original_sample_names
    )
  )
)


print(
  "Detected N44VRL sample numbers:"
)

print(
  sample_numbers
)


############################################################
# VERIFY SAMPLE NUMBERS
############################################################

if (
  length(sample_numbers) != 12 ||
  any(is.na(sample_numbers)) ||
  !setequal(
    sample_numbers,
    1:12
  )
) {
  
  stop(
    paste0(
      "\nERROR: Could not correctly identify N44VRL ",
      "samples 1-12.\n\n",
      "Detected sample names:\n",
      paste(
        original_sample_names,
        collapse = "\n"
      ),
      "\n\nDetected sample numbers:\n",
      paste(
        sample_numbers,
        collapse = ", "
      )
    )
  )
  
}


############################################################
# REORDER BY N44VRL SAMPLE NUMBER
############################################################

sample_columns <- sample_columns_original[
  order(
    sample_numbers
  )
]

original_sample_names <- colnames(
  counts_raw
)[sample_columns]


print("======================================================")
print("CORRECTED SAMPLE ORDER")
print("======================================================")


for (i in seq_along(original_sample_names)) {
  
  print(
    paste(
      i,
      ":",
      original_sample_names[i]
    )
  )
  
}


############################################################
# BIOLOGICAL SAMPLE ORDER
############################################################

expected_genotype <- c(
  "WT",
  "WT",
  "WT",
  "WT",
  "WT",
  "WT",
  "KCS1",
  "KCS1",
  "KCS1",
  "KCS1",
  "KCS1",
  "KCS1"
)

expected_pi <- c(
  "plus",
  "plus",
  "plus",
  "minus",
  "minus",
  "minus",
  "plus",
  "plus",
  "plus",
  "minus",
  "minus",
  "minus"
)


############################################################
# DETERMINE GENOTYPE FROM SAMPLE NAMES
############################################################

detected_genotype <- ifelse(
  
  grepl(
    "WT",
    original_sample_names,
    ignore.case = TRUE
  ),
  
  "WT",
  
  ifelse(
    
    grepl(
      "KCS1",
      original_sample_names,
      ignore.case = TRUE
    ),
    
    "KCS1",
    
    NA
    
  )
  
)


############################################################
# DETERMINE PI CONDITION
############################################################

detected_pi <- ifelse(
  
  grepl(
    "\\+Pi",
    original_sample_names,
    ignore.case = TRUE
  ),
  
  "plus",
  
  ifelse(
    
    grepl(
      "-Pi",
      original_sample_names,
      ignore.case = TRUE
    ),
    
    "minus",
    
    NA
    
  )
  
)


############################################################
# VERIFY SAMPLE INFORMATION
############################################################

if (
  any(is.na(detected_genotype)) ||
  any(is.na(detected_pi))
) {
  
  print(original_sample_names)
  print(detected_genotype)
  print(detected_pi)
  
  stop(
    "\nERROR: Could not determine genotype or Pi condition from sample names."
  )
  
}


############################################################
# FINAL BIOLOGICAL ORDER CHECK
############################################################

if (
  !identical(
    detected_genotype,
    expected_genotype
  ) ||
  !identical(
    detected_pi,
    expected_pi
  )
) {
  
  print(
    "Expected genotype:"
  )
  
  print(
    expected_genotype
  )
  
  print(
    "Detected genotype:"
  )
  
  print(
    detected_genotype
  )
  
  print(
    "Expected Pi:"
  )
  
  print(
    expected_pi
  )
  
  print(
    "Detected Pi:"
  )
  
  print(
    detected_pi
  )
  
  stop(
    "\nERROR: Sample order does not match the expected experimental design."
  )
  
}


print(
  "SUCCESS: Sample order is correct."
)


############################################################
# 5. EXTRACT GENE IDS AND COUNTS
############################################################

counts_clean <- counts_raw[
  ,
  c(1, sample_columns),
  drop = FALSE
]

colnames(counts_clean)[1] <- "Geneid"


############################################################
# CONVERT COUNTS TO NUMERIC
############################################################

for (i in 2:ncol(counts_clean)) {
  
  counts_clean[[i]] <- as.numeric(
    counts_clean[[i]]
  )
  
}


############################################################
# CHECK FOR NA VALUES
############################################################

if (
  any(
    is.na(
      counts_clean[, -1]
    )
  )
) {
  
  stop(
    "\nERROR: NA values detected in count matrix."
  )
  
}


############################################################
# COLLAPSE DUPLICATE GENE IDS
############################################################

print(
  "Collapsing duplicate gene IDs..."
)

counts_fixed <- aggregate(
  . ~ Geneid,
  data = counts_clean,
  FUN = sum
)


############################################################
# 6. REMOVE tRNA AND rRNA GENES
############################################################

print("======================================================")
print("REMOVING tRNA / rRNA FEATURES")
print("======================================================")


############################################################
# NORMALIZE COUNT IDS
############################################################

stripped_keys <- gsub(
  "^CAALFM_",
  "",
  counts_fixed$Geneid,
  ignore.case = TRUE
)

normalized_matrix_keys <- normalize_id(
  stripped_keys
)


############################################################
# MATCH FEATURE TYPES
############################################################

matched_features <- mapping_table$feature_type[
  match(
    normalized_matrix_keys,
    mapping_table$id_norm
  )
]


############################################################
# MATCH GENE NAMES
############################################################

matched_symbols <- mapping_table$gene_symbol[
  match(
    normalized_matrix_keys,
    mapping_table$id_norm
  )
]


############################################################
# IDENTIFY tRNA / rRNA
############################################################

is_trna_rrna <- (
  
  grepl(
    "tRNA|rRNA",
    counts_fixed$Geneid,
    ignore.case = TRUE
  ) |
    
    grepl(
      "tRNA|rRNA",
      matched_features,
      ignore.case = TRUE
    ) |
    
    grepl(
      "^tRNA|^rRNA",
      matched_symbols,
      ignore.case = TRUE
    )
  
)

is_trna_rrna[
  is.na(is_trna_rrna)
] <- FALSE


print(
  paste(
    "Features before tRNA/rRNA filtering:",
    nrow(counts_fixed)
  )
)

print(
  paste(
    "tRNA/rRNA features removed:",
    sum(is_trna_rrna)
  )
)


counts_fixed <- counts_fixed[
  !is_trna_rrna,
  ,
  drop = FALSE
]


print(
  paste(
    "Features remaining:",
    nrow(counts_fixed)
  )
)


############################################################
# CREATE COUNT MATRIX
############################################################

rownames(counts_fixed) <- counts_fixed$Geneid

counts <- counts_fixed[
  ,
  -1,
  drop = FALSE
]


############################################################
# 7. CREATE MASTER GENE ANNOTATION
############################################################

print("======================================================")
print("ANNOTATING GENES")
print("======================================================")


original_gene_ids <- rownames(
  counts
)


############################################################
# STRIP CAALFM PREFIX
############################################################

stripped_ids <- gsub(
  "^CAALFM_",
  "",
  original_gene_ids,
  ignore.case = TRUE
)


normalized_ids <- normalize_id(
  stripped_ids
)


############################################################
# MATCH CGD GENE NAMES
############################################################

matched_gene_names <- mapping_table$gene_symbol[
  match(
    normalized_ids,
    mapping_table$id_norm
  )
]


############################################################
# DISPLAY NAME
#
# If a CGD gene name exists, use it.
# Otherwise retain the systematic ID.
############################################################

display_gene_names <- ifelse(
  
  !is.na(matched_gene_names) &
    matched_gene_names != "",
  
  matched_gene_names,
  
  stripped_ids
  
)


############################################################
# MASTER ANNOTATION TABLE
############################################################

gene_annotation <- data.frame(
  
  Systematic_ID = original_gene_ids,
  
  Gene_Name = display_gene_names,
  
  stringsAsFactors = FALSE
  
)


############################################################
# MAPPING SUMMARY
############################################################

mapped_count <- sum(
  !is.na(matched_gene_names) &
    matched_gene_names != ""
)

unmapped_count <- sum(
  is.na(matched_gene_names) |
    matched_gene_names == ""
)


print(
  paste(
    "Genes mapped to CGD gene names:",
    mapped_count
  )
)

print(
  paste(
    "Genes without a CGD gene name:",
    unmapped_count
  )
)


############################################################
# SAVE ANNOTATION
############################################################

write.csv(
  
  gene_annotation,
  
  file.path(
    output_dir,
    "02_DE_results",
    "Gene_ID_to_Gene_Name_Annotation.csv"
  ),
  
  row.names = FALSE
  
)


############################################################
# HELPER FUNCTION FOR GENE NAME LOOKUP
############################################################

get_gene_names <- function(ids) {
  
  names_out <- gene_annotation$Gene_Name[
    match(
      ids,
      gene_annotation$Systematic_ID
    )
  ]
  
  missing <- (
    is.na(names_out) |
      names_out == ""
  )
  
  names_out[
    missing
  ] <- ids[
    missing
  ]
  
  return(
    names_out
  )
  
}


############################################################
# 8. SAMPLE METADATA
############################################################

short_sample_names <- c(
  
  "WT_Pi_1",
  "WT_Pi_2",
  "WT_Pi_3",
  
  "WT_minusPi_1",
  "WT_minusPi_2",
  "WT_minusPi_3",
  
  "KCS1_Pi_1",
  "KCS1_Pi_2",
  "KCS1_Pi_3",
  
  "KCS1_minusPi_1",
  "KCS1_minusPi_2",
  "KCS1_minusPi_3"
  
)


############################################################
# HUMAN-READABLE CONDITION LABELS
#
# THESE ARE THE LABELS THAT WILL APPEAR ON FIGURES.
############################################################

condition_labels <- c(
  
  "WT+Pi",
  "WT+Pi",
  "WT+Pi",
  
  "WT-Pi",
  "WT-Pi",
  "WT-Pi",
  
  "KCS1+Pi",
  "KCS1+Pi",
  "KCS1+Pi",
  
  "KCS1-Pi",
  "KCS1-Pi",
  "KCS1-Pi"
  
)


replicate_vector <- c(
  1, 2, 3,
  1, 2, 3,
  1, 2, 3,
  1, 2, 3
)


############################################################
# RENAME COUNT MATRIX
############################################################

colnames(counts) <- short_sample_names


############################################################
# METADATA
############################################################

metadata <- data.frame(
  
  Original_Sample_ID = original_sample_names,
  
  Sample_Name = short_sample_names,
  
  Display_Name = condition_labels,
  
  Genotype = factor(
    detected_genotype,
    levels = c(
      "WT",
      "KCS1"
    )
  ),
  
  Pi = factor(
    detected_pi,
    levels = c(
      "plus",
      "minus"
    )
  ),
  
  Replicate = replicate_vector,
  
  stringsAsFactors = FALSE
  
)


rownames(metadata) <- short_sample_names


############################################################
# INTERNAL CONDITION VARIABLE
############################################################

metadata$Condition <- factor(
  
  paste(
    metadata$Genotype,
    metadata$Pi,
    sep = "_"
  ),
  
  levels = c(
    "WT_plus",
    "WT_minus",
    "KCS1_plus",
    "KCS1_minus"
  )
  
)


############################################################
# CONDITION LABEL FACTOR
############################################################

metadata$Condition_Label <- factor(
  
  condition_labels,
  
  levels = c(
    "WT+Pi",
    "WT-Pi",
    "KCS1+Pi",
    "KCS1-Pi"
  )
  
)


print("======================================================")
print("SAMPLE METADATA")
print("======================================================")

print(
  metadata
)


############################################################
# VERIFY REPLICATES
############################################################

replicate_table <- table(
  metadata$Condition
)

print(
  replicate_table
)

if (!all(replicate_table == 3)) {
  
  stop(
    "ERROR: Every condition must contain exactly 3 replicates."
  )
  
}


write.csv(
  
  metadata,
  
  file.path(
    output_dir,
    "sample_metadata.csv"
  ),
  
  row.names = TRUE
  
)


############################################################
# 9. FILTER LOW-COUNT GENES
############################################################

print("======================================================")
print("FILTERING LOW-COUNT GENES")
print("======================================================")


keep_genes <- rowSums(
  counts
) >= 10


print(
  paste(
    "Genes before filtering:",
    nrow(counts)
  )
)

print(
  paste(
    "Genes passing count filter:",
    sum(keep_genes)
  )
)

print(
  paste(
    "Genes removed:",
    sum(!keep_genes)
  )
)


counts_filtered <- counts[
  keep_genes,
  ,
  drop = FALSE
]


############################################################
# 10. DESEQ2
############################################################

print("======================================================")
print("CREATING DESEQ2 OBJECT")
print("======================================================")


dds <- DESeqDataSetFromMatrix(
  
  countData = round(
    as.matrix(
      counts_filtered
    )
  ),
  
  colData = metadata,
  
  design = ~ Genotype * Pi
  
)


dds <- DESeq(
  dds
)


saveRDS(
  
  dds,
  
  file.path(
    output_dir,
    "02_DE_results",
    "DESeq2_object.rds"
  )
  
)


############################################################
# 11. VST
############################################################

print("======================================================")
print("VST TRANSFORMATION")
print("======================================================")


vsd <- vst(
  dds,
  blind = FALSE
)

vst_matrix <- assay(
  vsd
)


write.csv(
  
  vst_matrix,
  
  file.path(
    output_dir,
    "01_QC",
    "VST_expression_matrix.csv"
  ),
  
  row.names = TRUE
  
)


############################################################
# 12. LIBRARY SIZE QC
############################################################

library_sizes <- colSums(
  counts(dds)
)


library_size_df <- data.frame(
  
  Sample = names(library_sizes),
  
  Condition = metadata[
    names(library_sizes),
    "Display_Name"
  ],
  
  Library_Size = as.numeric(
    library_sizes
  ),
  
  stringsAsFactors = FALSE
  
)


write.csv(
  
  library_size_df,
  
  file.path(
    output_dir,
    "01_QC",
    "Library_sizes.csv"
  ),
  
  row.names = FALSE
  
)


p_library <- ggplot(
  
  library_size_df,
  
  aes(
    x = Sample,
    y = Library_Size
  )
  
) +
  
  geom_col() +
  
  theme_bw() +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  ) +
  
  labs(
    title = "RNA-seq library sizes",
    x = "Sample",
    y = "Total assigned reads"
  )


ggsave(
  
  file.path(
    output_dir,
    "01_QC",
    "Library_sizes.png"
  ),
  
  p_library,
  
  width = 10,
  height = 6,
  dpi = 300
  
)


############################################################
# 13. PCA
############################################################

print("======================================================")
print("CREATING PCA")
print("======================================================")


pca_data <- plotPCA(
  
  vsd,
  
  intgroup = c(
    "Genotype",
    "Pi"
  ),
  
  returnData = TRUE
  
)


percent_variance <- round(
  
  100 *
    attr(
      pca_data,
      "percentVar"
    )
  
)


p_pca <- ggplot(
  
  pca_data,
  
  aes(
    x = PC1,
    y = PC2,
    label = name
  )
  
) +
  
  geom_point(
    size = 4
  ) +
  
  geom_text_repel(
    size = 3
  ) +
  
  theme_bw() +
  
  labs(
    
    title = "PCA of RNA-seq samples",
    
    x = paste0(
      "PC1: ",
      percent_variance[1],
      "% variance"
    ),
    
    y = paste0(
      "PC2: ",
      percent_variance[2],
      "% variance"
    )
    
  )


ggsave(
  
  file.path(
    output_dir,
    "06_PCA",
    "PCA_PC1_PC2.png"
  ),
  
  p_pca,
  
  width = 9,
  height = 7,
  dpi = 300
  
)


############################################################
# 14. SAMPLE CORRELATION
############################################################

cor_matrix <- cor(
  
  vst_matrix,
  
  method = "pearson"
  
)


write.csv(
  
  cor_matrix,
  
  file.path(
    output_dir,
    "07_Sample_correlations",
    "Sample_Pearson_correlations.csv"
  )
  
)


png(
  
  file.path(
    output_dir,
    "07_Sample_correlations",
    "Sample_Pearson_correlations.png"
  ),
  
  width = 2200,
  height = 2000,
  res = 300
  
)


cor_annotation <- data.frame(
  
  Condition = metadata[
    colnames(vst_matrix),
    "Display_Name"
  ]
  
)

rownames(cor_annotation) <- colnames(
  vst_matrix
)


pheatmap(
  
  cor_matrix,
  
  annotation_col = cor_annotation,
  
  annotation_row = cor_annotation,
  
  clustering_distance_rows = "correlation",
  
  clustering_distance_cols = "correlation",
  
  main = "Sample Pearson correlations",
  
  fontsize = 8
  
)


dev.off()


############################################################
# 15. DESEQ2 COEFFICIENTS
############################################################

print("======================================================")
print("DESEQ2 COEFFICIENTS")
print("======================================================")


results_names <- resultsNames(
  dds
)

print(
  results_names
)


############################################################
# GENOTYPE COEFFICIENT
############################################################

genotype_coef <- results_names[
  results_names == "Genotype_KCS1_vs_WT"
]


############################################################
# PI COEFFICIENT
############################################################

pi_coef <- results_names[
  results_names == "Pi_minus_vs_plus"
]


############################################################
# INTERACTION COEFFICIENT
############################################################

interaction_coef <- results_names[
  grepl(
    "Genotype.*Pi|Pi.*Genotype",
    results_names
  )
]


############################################################
# VERIFY COEFFICIENTS
############################################################

if (length(genotype_coef) != 1) {
  
  stop(
    paste0(
      "Could not identify genotype coefficient.\n",
      "Available coefficients:\n",
      paste(
        results_names,
        collapse = "\n"
      )
    )
  )
  
}


if (length(pi_coef) != 1) {
  
  stop(
    paste0(
      "Could not identify Pi coefficient.\n",
      "Available coefficients:\n",
      paste(
        results_names,
        collapse = "\n"
      )
    )
  )
  
}


if (length(interaction_coef) != 1) {
  
  stop(
    paste0(
      "Could not identify interaction coefficient.\n",
      "Available coefficients:\n",
      paste(
        results_names,
        collapse = "\n"
      )
    )
  )
  
}


############################################################
# PRINT COEFFICIENTS
############################################################

print(
  paste(
    "Genotype coefficient:",
    genotype_coef
  )
)

print(
  paste(
    "Pi coefficient:",
    pi_coef
  )
)

print(
  paste(
    "Interaction coefficient:",
    interaction_coef
  )
)
  
  ############################################################
  # 16. ALL SIX PAIRWISE COMPARISONS
  ############################################################
  #
  # Reference condition:
  #
  # WT+Pi = 0
  #
  # WT-Pi = Pi
  #
  # KCS1+Pi = Genotype
  #
  # KCS1-Pi =
  # Genotype + Pi + Interaction
  #
  ############################################################
  
  
  ############################################################
  # 1. WT+Pi vs WT-Pi
  #
  # Result is WT-Pi - WT+Pi
  ############################################################
  
  res_WT_minus_vs_WT_plus <- results(
    
    dds,
    
    name = pi_coef
    
  )
  
  
  ############################################################
  # 2. WT+Pi vs KCS1+Pi
  #
  # Result is KCS1+Pi - WT+Pi
  ############################################################
  
  res_KCS1_plus_vs_WT_plus <- results(
    
    dds,
    
    name = genotype_coef
    
  )
  
  
  ############################################################
  # 3. WT+Pi vs KCS1-Pi
  #
  # Result is KCS1-Pi - WT+Pi
  #
  # = Genotype + Pi + Interaction
  ############################################################
  
  res_KCS1_minus_vs_WT_plus <- results(
    
    dds,
    
    contrast = list(
      c(
        genotype_coef,
        pi_coef,
        interaction_coef
      )
    )
    
  )
  
  
  ############################################################
  # 4. WT-Pi vs KCS1+Pi
  #
  # Result is KCS1+Pi - WT-Pi
  #
  # = Genotype - Pi
  ############################################################
  
  res_KCS1_plus_vs_WT_minus <- results(
    
    dds,
    
    contrast = list(
      c(
        genotype_coef
      ),
      c(
        pi_coef
      )
    )
    
  )
  
  
  ############################################################
  # 5. WT-Pi vs KCS1-Pi
  #
  # Result is KCS1-Pi - WT-Pi
  #
  # = Genotype + Interaction
  ############################################################
  
  res_KCS1_minus_vs_WT_minus <- results(
    
    dds,
    
    contrast = list(
      c(
        genotype_coef,
        interaction_coef
      )
    )
    
  )
  
  
  ############################################################
  # 6. KCS1+Pi vs KCS1-Pi
  #
  # Result is KCS1-Pi - KCS1+Pi
  #
  # = Pi + Interaction
  ############################################################
  
  res_KCS1_minus_vs_KCS1_plus <- results(
    
    dds,
    
    contrast = list(
      c(
        pi_coef,
        interaction_coef
      )
    )
    
  )
  
  
  ############################################################
  # 7. FORMAL GENOTYPE x PI INTERACTION
  ############################################################
  
  res_interaction <- results(
    
    dds,
    
    name = interaction_coef
    
  )
  
  
  ############################################################
  # 17. ANNOTATE ALL DE RESULTS
  ############################################################
  
  annotate_results <- function(
    result,
    comparison_label
  ) {
    
    df <- as.data.frame(
      result
    )
    
    df$Systematic_ID <- rownames(
      df
    )
    
    df$Gene_Name <- get_gene_names(
      df$Systematic_ID
    )
    
    df$Comparison <- comparison_label
    
    
    ##########################################################
    # Put gene information first
    ##########################################################
    
    df <- df[
      ,
      c(
        "Systematic_ID",
        "Gene_Name",
        "Comparison",
        setdiff(
          colnames(df),
          c(
            "Systematic_ID",
            "Gene_Name",
            "Comparison"
          )
        )
      ),
      drop = FALSE
    ]
    
    return(df)
    
  }
  
  
  res_WT_minus_vs_WT_plus_df <- annotate_results(
    res_WT_minus_vs_WT_plus,
    "WT-Pi vs WT+Pi"
  )
  
  
  res_KCS1_plus_vs_WT_plus_df <- annotate_results(
    res_KCS1_plus_vs_WT_plus,
    "KCS1+Pi vs WT+Pi"
  )
  
  
  res_KCS1_minus_vs_WT_plus_df <- annotate_results(
    res_KCS1_minus_vs_WT_plus,
    "KCS1-Pi vs WT+Pi"
  )
  
  
  res_KCS1_plus_vs_WT_minus_df <- annotate_results(
    res_KCS1_plus_vs_WT_minus,
    "KCS1+Pi vs WT-Pi"
  )
  
  
  res_KCS1_minus_vs_WT_minus_df <- annotate_results(
    res_KCS1_minus_vs_WT_minus,
    "KCS1-Pi vs WT-Pi"
  )
  
  
  res_KCS1_minus_vs_KCS1_plus_df <- annotate_results(
    res_KCS1_minus_vs_KCS1_plus,
    "KCS1-Pi vs KCS1+Pi"
  )
  
  
  res_interaction_df <- annotate_results(
    res_interaction,
    "Genotype × Pi interaction"
  )
  
  
  ############################################################
  # CREATE NAMED RESULT LIST
  ############################################################
  
  all_results <- list(
    
    "WT-Pi_vs_WT+Pi" =
      res_WT_minus_vs_WT_plus_df,
    
    "KCS1+Pi_vs_WT+Pi" =
      res_KCS1_plus_vs_WT_plus_df,
    
    "KCS1-Pi_vs_WT+Pi" =
      res_KCS1_minus_vs_WT_plus_df,
    
    "KCS1+Pi_vs_WT-Pi" =
      res_KCS1_plus_vs_WT_minus_df,
    
    "KCS1-Pi_vs_WT-Pi" =
      res_KCS1_minus_vs_WT_minus_df,
    
    "KCS1-Pi_vs_KCS1+Pi" =
      res_KCS1_minus_vs_KCS1_plus_df,
    
    "Genotype_x_Pi_interaction" =
      res_interaction_df
    
  )
  
  
  ############################################################
  # 18. SAVE ALL DE RESULTS
  ############################################################
  
  for (comparison_name in names(all_results)) {
    
    write.csv(
      
      all_results[[comparison_name]],
      
      file.path(
        output_dir,
        "02_DE_results",
        paste0(
          "DE_",
          comparison_name,
          ".csv"
        )
      ),
      
      row.names = FALSE
      
    )
    
  }
  
  
  ############################################################
  # 19. SIGNIFICANT GENE LISTS
  ############################################################
  
  get_significant_genes <- function(
    result_df,
    padj_cutoff = 0.05
  ) {
    
    result_df %>%
      
      filter(
        !is.na(padj),
        padj < padj_cutoff
      ) %>%
      
      arrange(
        padj
      )
    
  }
  
  
  for (comparison_name in names(all_results)) {
    
    sig_genes <- get_significant_genes(
      all_results[[comparison_name]]
    )
    
    
    write.csv(
      
      sig_genes,
      
      file.path(
        output_dir,
        "08_Significant_gene_lists",
        paste0(
          "Significant_",
          comparison_name,
          ".csv"
        )
      ),
      
      row.names = FALSE
      
    )
    
  }
  
 
  
  
  ############################################################
  # 20. VOLCANO PLOT FUNCTION
  ############################################################
  
  make_volcano <- function(
    result_df,
    plot_title,
    output_file,
    fc_cutoff = 1,
    padj_cutoff = 0.05,
    y_max = 50
  ) {
    
    ##########################################################
    # Check that required columns exist
    ##########################################################
    
    required_cols <- c(
      "log2FoldChange",
      "padj"
    )
    
    missing_cols <- setdiff(
      required_cols,
      colnames(result_df)
    )
    
    if (length(missing_cols) > 0) {
      
      stop(
        paste(
          "Missing required columns:",
          paste(missing_cols, collapse = ", ")
        )
      )
      
    }
    
    
    ##########################################################
    # Keep genes with valid fold changes and adjusted p-values
    ##########################################################
    
    plot_df <- result_df %>%
      filter(
        !is.na(log2FoldChange),
        !is.na(padj),
        padj > 0
      )
    
    
    if (nrow(plot_df) == 0) {
      
      warning(
        paste(
          "No valid genes available for:",
          plot_title
        )
      )
      
      return(NULL)
      
    }
    
    
    ##########################################################
    # Calculate -log10 adjusted p-value
    ##########################################################
    
    plot_df$neg_log10_padj <- -log10(
      plot_df$padj
    )
    
    
    ##########################################################
    # Assign significance categories
    ##########################################################
    
    plot_df$Significance <- "Not significant"
    
    
    plot_df$Significance[
      plot_df$padj < padj_cutoff &
        plot_df$log2FoldChange >= fc_cutoff
    ] <- "Upregulated"
    
    
    plot_df$Significance[
      plot_df$padj < padj_cutoff &
        plot_df$log2FoldChange <= -fc_cutoff
    ] <- "Downregulated"
    
    
    ##########################################################
    # Identify genes above the displayed y-axis limit
    ##########################################################
    
    plot_df$Above_Y_Max <- (
      plot_df$neg_log10_padj > y_max
    )
    
    
    n_above <- sum(
      plot_df$Above_Y_Max,
      na.rm = TRUE
    )
    
    
    ##########################################################
    # Create display value
    #
    # The TRUE -log10(padj) is retained in
    # plot_df$neg_log10_padj.
    #
    # Only the value used for plotting is capped.
    ##########################################################
    
    plot_df$neg_log10_padj_plot <- pmin(
      plot_df$neg_log10_padj,
      y_max
    )
    
    
    ##########################################################
    # Select genes to label
    ##########################################################
    
    plot_df$Label <- NA_character_
    
    
    if (
      "Systematic_ID" %in% colnames(plot_df) &&
      "Gene_Name" %in% colnames(plot_df)
    ) {
      
      label_genes <- plot_df %>%
        filter(
          padj < padj_cutoff,
          abs(log2FoldChange) >= fc_cutoff,
          !is.na(Gene_Name),
          Gene_Name != ""
        ) %>%
        arrange(padj) %>%
        head(40)
      
      
      if (nrow(label_genes) > 0) {
        
        label_indices <- match(
          label_genes$Systematic_ID,
          plot_df$Systematic_ID
        )
        
        plot_df$Label[label_indices] <-
          label_genes$Gene_Name
        
      }
      
    }
    
    
    ##########################################################
    # Create volcano plot
    ##########################################################
    
    p <- ggplot(
      plot_df,
      aes(
        x = log2FoldChange,
        y = neg_log10_padj_plot,
        shape = Significance
      )
    ) +
      
      ########################################################
    # All genes
    ########################################################
    
    geom_point(
      alpha = 0.6,
      size = 1.5
    ) +
      
      ########################################################
    # Fold-change thresholds
    ########################################################
    
    geom_vline(
      xintercept = c(
        -fc_cutoff,
        fc_cutoff
      ),
      linetype = "dashed"
    ) +
      
      ########################################################
    # Adjusted p-value threshold
    ########################################################
    
    geom_hline(
      yintercept = -log10(padj_cutoff),
      linetype = "dashed"
    ) +
      
      ########################################################
    # Mark genes above y-axis display limit
    #
    # Triangle indicates that the true value is > y_max.
    ########################################################
    
    geom_point(
      data = plot_df %>%
        filter(Above_Y_Max),
      
      aes(
        x = log2FoldChange,
        y = y_max
      ),
      
      shape = 24,
      size = 2.5,
      
      inherit.aes = FALSE
    ) +
      
      ########################################################
    # Gene labels
    ########################################################
    
    geom_text_repel(
      aes(
        label = Label
      ),
      
      na.rm = TRUE,
      size = 3,
      max.overlaps = 20
    ) +
      
      ########################################################
    # Display y-axis from 0 to 50
    #
    # This does NOT change the actual statistical results.
    ########################################################
    
    coord_cartesian(
      ylim = c(
        0,
        y_max
      ),
      clip = "off"
    ) +
      
      ########################################################
    # Theme
    ########################################################
    
    theme_bw() +
      
      labs(
        title = plot_title,
        x = "log2 fold change",
        y = expression(-log[10]("adjusted p-value")),
        shape = "Significance"
      )
    
    
    ##########################################################
    # Add annotation for genes above y-axis limit
    ##########################################################
    
    if (n_above > 0) {
      
      p <- p +
        annotate(
          "text",
          x = Inf,
          y = y_max,
          label = paste0(
            "▲ ",
            n_above,
            " genes > ",
            y_max
          ),
          hjust = 1.05,
          vjust = -0.5,
          size = 3.5
        )
      
    }
    
    
    ##########################################################
    # Save plot
    ##########################################################
    
    ggsave(
      filename = output_file,
      plot = p,
      width = 8,
      height = 6,
      dpi = 300
    )
    
    
    ##########################################################
    # Return plot
    ##########################################################
    
    return(p)
    
  }
  
  
  ############################################################
  # 21. CREATE VOLCANO PLOTS FOR ALL SIX PAIRWISE COMPARISONS
  ############################################################
  
  print("======================================================")
  print("CREATING VOLCANO PLOTS")
  print("======================================================")
  
  volcano_titles <- c(
    "WT-Pi vs WT+Pi",
    "KCS1+Pi vs WT+Pi",
    "KCS1-Pi vs WT+Pi",
    "KCS1+Pi vs WT-Pi",
    "KCS1-Pi vs WT-Pi",
    "KCS1-Pi vs KCS1+Pi"
  )
  
  volcano_files <- c(
    "Volcano_WT-Pi_vs_WT+Pi.png",
    "Volcano_KCS1+Pi_vs_WT+Pi.png",
    "Volcano_KCS1-Pi_vs_WT+Pi.png",
    "Volcano_KCS1+Pi_vs_WT-Pi.png",
    "Volcano_KCS1-Pi_vs_WT-Pi.png",
    "Volcano_KCS1-Pi_vs_KCS1+Pi.png"
  )
  
  ############################################################
  # USE ANNOTATED DATA FRAMES
  ############################################################
  
  pairwise_results <- list(
    res_WT_minus_vs_WT_plus_df,
    res_KCS1_plus_vs_WT_plus_df,
    res_KCS1_minus_vs_WT_plus_df,
    res_KCS1_plus_vs_WT_minus_df,
    res_KCS1_minus_vs_WT_minus_df,
    res_KCS1_minus_vs_KCS1_plus_df
  )
  
  ############################################################
  # CHECK THAT GENE NAMES ARE PRESENT
  ############################################################
  
  print("Checking gene-name annotations:")
  
  for (i in seq_along(pairwise_results)) {
    
    print(
      paste(
        volcano_titles[i],
        "columns:",
        paste(
          c("Systematic_ID", "Gene_Name") %in%
            colnames(pairwise_results[[i]]),
          collapse = ", "
        )
      )
    )
    
    print(
      paste(
        "Number of gene names:",
        sum(
          !is.na(pairwise_results[[i]]$Gene_Name) &
            pairwise_results[[i]]$Gene_Name != ""
        )
      )
    )
  }
  
  ############################################################
  # CREATE PAIRWISE VOLCANO PLOTS
  ############################################################
  
  for (i in seq_along(pairwise_results)) {
    
    make_volcano(
      result_df = pairwise_results[[i]],
      plot_title = volcano_titles[i],
      output_file = file.path(
        output_dir,
        "03_Volcano",
        volcano_files[i]
      )
    )
  }
  
  ############################################################
  # INTERACTION VOLCANO PLOT
  ############################################################
  
  make_volcano(
    result_df = res_interaction_df,
    plot_title = "Genotype × Pi interaction",
    output_file = file.path(
      output_dir,
      "03_Volcano",
      "Volcano_Genotype_x_Pi_interaction.png"
    )
  )
  
  
  ############################################################
  # 22. MA PLOT FUNCTION
  ############################################################
  
  make_ma_plot <- function(
    result_df,
    plot_title,
    output_file
  ) {
    
    
    plot_df <- result_df %>%
      
      filter(
        !is.na(baseMean),
        !is.na(log2FoldChange)
      )
    
    
    p <- ggplot(
      
      plot_df,
      
      aes(
        x = log10(baseMean + 1),
        y = log2FoldChange
      )
      
    ) +
      
      geom_point(
        alpha = 0.4,
        size = 1
      ) +
      
      geom_hline(
        yintercept = 0,
        linetype = "dashed"
      ) +
      
      theme_bw() +
      
      labs(
        title = plot_title,
        x = "log10 mean normalized expression",
        y = "log2 fold change"
      )
    
    
    ggsave(
      
      output_file,
      
      p,
      
      width = 9,
      height = 7,
      dpi = 300
      
    )
    
  }
  
  
  ############################################################
  # 23. MA PLOTS FOR ALL COMPARISONS
  ############################################################
  
  ma_files <- c(
    
    "MA_WT-Pi_vs_WT+Pi.png",
    "MA_KCS1+Pi_vs_WT+Pi.png",
    "MA_KCS1-Pi_vs_WT+Pi.png",
    "MA_KCS1+Pi_vs_WT-Pi.png",
    "MA_KCS1-Pi_vs_WT-Pi.png",
    "MA_KCS1-Pi_vs_KCS1+Pi.png"
    
  )
  
  
  for (i in seq_along(pairwise_results)) {
    
    make_ma_plot(
      
      pairwise_results[[i]],
      
      volcano_titles[i],
      
      file.path(
        output_dir,
        "04_MA_plots",
        ma_files[i]
      )
      
    )
    
  }
  
  
  make_ma_plot(
    
    res_interaction_df,
    
    "Genotype × Pi interaction",
    
    file.path(
      output_dir,
      "04_MA_plots",
      "MA_Genotype_x_Pi_interaction.png"
    )
    
  )
  
  
  ############################################################
  # 24. TOP VARIABLE GENE HEATMAP
  ############################################################
  
  print("======================================================")
  print("CREATING TOP VARIABLE GENE HEATMAP")
  print("======================================================")
  
  
  gene_variances <- rowVars(
    vst_matrix
  )
  
  
  top_n <- min(
    50,
    length(gene_variances)
  )
  
  
  top_variable_genes <- names(
    sort(
      gene_variances,
      decreasing = TRUE
    )
  )[1:top_n]
  
  
  heatmap_matrix <- vst_matrix[
    top_variable_genes,
    ,
    drop = FALSE
  ]
  
  
  heatmap_gene_names <- get_gene_names(
    rownames(
      heatmap_matrix
    )
  )
  
  
  rownames(
    heatmap_matrix
  ) <- make.unique(
    heatmap_gene_names
  )
  
  
  heatmap_scaled <- t(
    scale(
      t(
        heatmap_matrix
      )
    )
  )
  
  
  annotation_col <- data.frame(
    
    Condition = metadata[
      colnames(
        heatmap_scaled
      ),
      "Display_Name"
    ]
    
  )
  
  rownames(annotation_col) <- colnames(
    heatmap_scaled
  )
  
  
  png(
    
    file.path(
      output_dir,
      "05_Heatmaps",
      "Heatmap_Top_50_Variable_Genes.png"
    ),
    
    width = 2400,
    height = 2400,
    res = 300
    
  )
  
  
  pheatmap(
    
    heatmap_scaled,
    
    scale = "none",
    
    annotation_col = annotation_col,
    
    cluster_rows = TRUE,
    
    cluster_cols = TRUE,
    
    show_rownames = TRUE,
    
    fontsize_row = 6,
    
    fontsize_col = 9,
    
    main = "Top 50 most variable genes"
    
  )
  
  
  dev.off()
  
  
  ############################################################
  # 25. DE HEATMAP FUNCTION
  ############################################################
  
  make_de_heatmap <- function(
    result_df,
    plot_title,
    output_file,
    n_genes = 50
  ) {
    
    
    significant_genes <- result_df %>%
      
      filter(
        !is.na(padj),
        padj < 0.05
      ) %>%
      
      arrange(
        padj
      )
    
    
    if (nrow(significant_genes) < 2) {
      
      print(
        paste(
          "Too few significant genes for:",
          plot_title
        )
      )
      
      return(NULL)
      
    }
    
    
    selected_genes <- head(
      
      significant_genes$Systematic_ID,
      
      n_genes
      
    )
    
    
    selected_genes <- intersect(
      
      selected_genes,
      
      rownames(
        vst_matrix
      )
      
    )
    
    
    if (length(selected_genes) < 2) {
      return(NULL)
    }
    
    
    heatmap_data <- vst_matrix[
      
      selected_genes,
      ,
      drop = FALSE
      
    ]
    
    
    gene_names <- get_gene_names(
      
      rownames(
        heatmap_data
      )
      
    )
    
    
    rownames(
      heatmap_data
    ) <- make.unique(
      gene_names
    )
    
    
    heatmap_scaled <- t(
      
      scale(
        t(
          heatmap_data
        )
      )
      
    )
    
    
    ##########################################################
    # CONDITION LABELS
    ##########################################################
    
    annotation_col <- data.frame(
      
      Condition = metadata[
        colnames(
          heatmap_scaled
        ),
        "Display_Name"
      ]
      
    )
    
    
    rownames(annotation_col) <- colnames(
      heatmap_scaled
    )
    
    
    ##########################################################
    # FIXED-SIZE PNG
    ##########################################################
    
    png(
      
      output_file,
      
      width = 2400,
      height = 2600,
      res = 300
      
    )
    
    
    pheatmap(
      
      heatmap_scaled,
      
      scale = "none",
      
      annotation_col = annotation_col,
      
      cluster_rows = TRUE,
      
      cluster_cols = TRUE,
      
      show_rownames = TRUE,
      
      fontsize_row = 7,
      
      fontsize_col = 9,
      
      main = plot_title
      
    )
    
    
    dev.off()
    
  }
  
  
  ############################################################
  # 26. HEATMAPS FOR ALL SIX PAIRWISE COMPARISONS
  ############################################################
  
  heatmap_files <- c(
    
    "Heatmap_WT-Pi_vs_WT+Pi.png",
    "Heatmap_KCS1+Pi_vs_WT+Pi.png",
    "Heatmap_KCS1-Pi_vs_WT+Pi.png",
    "Heatmap_KCS1+Pi_vs_WT-Pi.png",
    "Heatmap_KCS1-Pi_vs_WT-Pi.png",
    "Heatmap_KCS1-Pi_vs_KCS1+Pi.png"
    
  )
  
  
  for (i in seq_along(pairwise_results)) {
    
    make_de_heatmap(
      
      pairwise_results[[i]],
      
      paste(
        "Differential expression:",
        volcano_titles[i]
      ),
      
      file.path(
        output_dir,
        "05_Heatmaps",
        heatmap_files[i]
      )
      
    )
    
  }
  
  
  ############################################################
  # INTERACTION HEATMAP
  ############################################################
  
  make_de_heatmap(
    
    res_interaction_df,
    
    "Genotype × Pi interaction",
    
    file.path(
      output_dir,
      "05_Heatmaps",
      "Heatmap_Genotype_x_Pi_interaction.png"
    )
    
  )
  
  ############################################################
  # 27. CONDITION-LEVEL VST HEATMAP
  ############################################################
  #
  # Replicates are averaged ONLY for visualization.
  #
  # Conditions:
  #
  # WT+Pi
  # WT-Pi
  # KCS1+Pi
  # KCS1-Pi
  #
  ############################################################
  
  print("======================================================")
  print("CREATING CONDITION-LEVEL HEATMAP")
  print("======================================================")
  
  
  ############################################################
  # CREATE CONDITION MEAN MATRIX
  ############################################################
  
  condition_means <- matrix(
    
    NA_real_,
    
    nrow = nrow(vst_matrix),
    
    ncol = 4
    
  )
  
  
  ############################################################
  # PRESERVE GENE IDs
  ############################################################
  
  rownames(condition_means) <- rownames(
    vst_matrix
  )
  
  
  ############################################################
  # CONDITION NAMES
  ############################################################
  
  colnames(condition_means) <- c(
    
    "WT+Pi",
    "WT-Pi",
    "KCS1+Pi",
    "KCS1-Pi"
    
  )
  
  
  ############################################################
  # IDENTIFY SAMPLE INDICES FOR EACH CONDITION
  ############################################################
  
  condition_indices <- list(
    
    "WT+Pi" = which(
      metadata$Condition == "WT_plus"
    ),
    
    "WT-Pi" = which(
      metadata$Condition == "WT_minus"
    ),
    
    "KCS1+Pi" = which(
      metadata$Condition == "KCS1_plus"
    ),
    
    "KCS1-Pi" = which(
      metadata$Condition == "KCS1_minus"
    )
    
  )
  
  
  ############################################################
  # VERIFY CONDITION INDICES
  ############################################################
  
  print("Condition sample indices:")
  
  for (condition_name in names(condition_indices)) {
    
    print(
      paste(
        condition_name,
        ":",
        paste(
          condition_indices[[condition_name]],
          collapse = ", "
        )
      )
    )
    
  }
  
  
  ############################################################
  # VERIFY THAT EACH CONDITION HAS 3 REPLICATES
  ############################################################
  
  condition_replicate_counts <- sapply(
    
    condition_indices,
    
    length
    
  )
  
  
  if (
    any(
      condition_replicate_counts != 3
    )
  ) {
    
    stop(
      paste0(
        "\nERROR: Condition-level heatmap does not have ",
        "exactly 3 replicates per condition.\n\n",
        paste(
          names(condition_replicate_counts),
          condition_replicate_counts,
          sep = ": ",
          collapse = "\n"
        )
      )
    )
    
  }
  
  
  ############################################################
  # CALCULATE MEAN VST EXPRESSION
  #
  # Replicates are averaged ONLY here for visualization.
  # DESeq2 statistics continue to use all biological replicates.
  ############################################################
  
  for (i in seq_along(condition_indices)) {
    
    condition_means[, i] <- rowMeans(
      
      vst_matrix[
        ,
        condition_indices[[i]],
        drop = FALSE
      ],
      
      na.rm = TRUE
      
    )
    
  }
  
  
  ############################################################
  # REMOVE GENES WITH NON-FINITE VALUES
  ############################################################
  
  finite_genes <- apply(
    
    condition_means,
    
    1,
    
    function(x) {
      
      all(
        is.finite(x)
      )
      
    }
    
  )
  
  
  print(
    paste(
      "Genes with finite values across all four conditions:",
      sum(finite_genes)
    )
  )
  
  
  condition_means <- condition_means[
    
    finite_genes,
    
    ,
    
    drop = FALSE
    
  ]
  
  
  ############################################################
  # CALCULATE VARIANCE ACROSS THE FOUR CONDITIONS
  ############################################################
  
  condition_variances <- matrixStats::rowVars(
    
    condition_means
    
  )
  
  
  ############################################################
  # REMOVE GENES WITH ZERO VARIANCE
  #
  # These genes would become NaN during row scaling.
  ############################################################
  
  variable_genes <- is.finite(
    
    condition_variances
    
  ) &
    
    condition_variances > 0
  
  
  print(
    paste(
      "Genes with variation across conditions:",
      sum(variable_genes)
    )
  )
  
  
  condition_variances <- condition_variances[
    
    variable_genes
    
  ]
  
  
  condition_means <- condition_means[
    
    variable_genes,
    
    ,
    
    drop = FALSE
    
  ]
  
  
  ############################################################
  # SELECT TOP 50 VARIABLE GENES
  #
  # IMPORTANT:
  # Use numeric row indices rather than names().
  # matrixStats::rowVars() does not necessarily preserve
  # row names as names of the returned variance vector.
  ############################################################
  
  n_condition_genes <- min(
    
    50,
    
    length(condition_variances)
    
  )
  
  
  if (
    n_condition_genes < 2
  ) {
    
    stop(
      paste0(
        "\nERROR: Fewer than 2 variable genes are available ",
        "for the condition-level heatmap."
      )
    )
    
  }
  
  
  top_condition_indices <- order(
    
    condition_variances,
    
    decreasing = TRUE
    
  )[
    seq_len(
      n_condition_genes
    )
  ]
  
  
  ############################################################
  # EXTRACT TOP VARIABLE GENES
  ############################################################
  
  condition_heatmap <- condition_means[
    
    top_condition_indices,
    
    ,
    
    drop = FALSE
    
  ]
  
  
  ############################################################
  # REPORT SELECTED GENES
  ############################################################
  
  print(
    paste(
      "Number of genes selected for condition heatmap:",
      nrow(condition_heatmap)
    )
  )
  
  
  ############################################################
  # MAP SYSTEMATIC IDs TO CGD GENE NAMES
  ############################################################
  
  condition_gene_names <- get_gene_names(
    
    rownames(
      condition_heatmap
    )
    
  )
  
  
  ############################################################
  # MAKE GENE NAMES UNIQUE
  ############################################################
  
  rownames(
    condition_heatmap
  ) <- make.unique(
    
    condition_gene_names
    
  )
  
  
  ############################################################
  # ROW-SCALE EXPRESSION
  #
  # Each gene is standardized across the four conditions.
  #
  # This shows relative expression patterns rather than
  # absolute VST expression.
  ############################################################
  
  condition_heatmap_scaled <- t(
    
    scale(
      
      t(
        condition_heatmap
      )
      
    )
    
  )
  
  
  ############################################################
  # REMOVE ANY ROWS THAT BECAME NON-FINITE AFTER SCALING
  ############################################################
  
  finite_scaled_genes <- apply(
    
    condition_heatmap_scaled,
    
    1,
    
    function(x) {
      
      all(
        is.finite(x)
      )
      
    }
    
  )
  
  
  condition_heatmap_scaled <- condition_heatmap_scaled[
    
    finite_scaled_genes,
    
    ,
    
    drop = FALSE
    
  ]
  
  
  ############################################################
  # VERIFY HEATMAP DIMENSIONS
  ############################################################
  
  print("Condition heatmap dimensions:")
  
  print(
    dim(
      condition_heatmap_scaled
    )
  )
  
  
  if (
    nrow(condition_heatmap_scaled) < 2
  ) {
    
    stop(
      "\nERROR: Fewer than 2 genes remain after scaling."
    )
    
  }
  
  
  if (
    ncol(condition_heatmap_scaled) != 4
  ) {
    
    stop(
      "\nERROR: Condition heatmap does not contain exactly 4 conditions."
    )
    
  }
  
  
  ############################################################
  # SAVE CONDITION MEAN MATRIX
  ############################################################
  
  write.csv(
    
    condition_means,
    
    file.path(
      output_dir,
      "05_Heatmaps",
      "Condition_Mean_VST_Expression.csv"
    ),
    
    row.names = TRUE
    
  )
  
  
  ############################################################
  # CREATE HEATMAP
  ############################################################
  
  png(
    
    file.path(
      output_dir,
      "05_Heatmaps",
      "Heatmap_Top_50_Variable_Condition_Means.png"
    ),
    
    width = 1800,
    height = 2400,
    res = 300
    
  )
  
  
  pheatmap(
    
    condition_heatmap_scaled,
    
    scale = "none",
    
    cluster_rows = TRUE,
    
    cluster_cols = FALSE,
    
    show_rownames = TRUE,
    
    show_colnames = TRUE,
    
    fontsize_row = 6,
    
    fontsize_col = 10,
    
    border_color = NA,
    
    main = "Top 50 variable genes: condition means"
    
  )
  
  
  dev.off()
  
  
  ############################################################
  # CONFIRM SUCCESS
  ############################################################
  
  print(
    "SUCCESS: Condition-level VST heatmap created."
  )
  
  print(
    paste(
      "Genes plotted:",
      nrow(condition_heatmap_scaled)
    )
  )
  
  print(
    "Conditions plotted:"
  )
  
  print(
    colnames(condition_heatmap_scaled)
  )
  
  ############################################################
  # 28. PHOSPHATE RESPONSE ANALYSIS
  ############################################################
  #
  # WT response:
  #
  # WT-Pi - WT+Pi
  #
  # KCS1 response:
  #
  # KCS1-Pi - KCS1+Pi
  #
  # The difference between these responses is tested formally
  # by the Genotype × Pi interaction.
  #
  ############################################################
  
  print("======================================================")
  print("PHOSPHATE RESPONSE ANALYSIS")
  print("======================================================")
  
  
  ############################################################
  # GET LOG2 FOLD CHANGES
  ############################################################
  
  WT_response <- res_WT_minus_vs_WT_plus_df$log2FoldChange
  
  KCS1_response <- res_KCS1_minus_vs_KCS1_plus_df$log2FoldChange
  
  
  response_table <- data.frame(
    
    Systematic_ID =
      res_WT_minus_vs_WT_plus_df$Systematic_ID,
    
    Gene_Name =
      res_WT_minus_vs_WT_plus_df$Gene_Name,
    
    WT_response =
      WT_response,
    
    WT_padj =
      res_WT_minus_vs_WT_plus_df$padj,
    
    KCS1_response =
      KCS1_response[
        match(
          res_WT_minus_vs_WT_plus_df$Systematic_ID,
          res_KCS1_minus_vs_KCS1_plus_df$Systematic_ID
        )
      ],
    
    KCS1_padj =
      res_KCS1_minus_vs_KCS1_plus_df$padj[
        match(
          res_WT_minus_vs_WT_plus_df$Systematic_ID,
          res_KCS1_minus_vs_KCS1_plus_df$Systematic_ID
        )
      ],
    
    Interaction_log2FC =
      res_interaction_df$log2FoldChange[
        match(
          res_WT_minus_vs_WT_plus_df$Systematic_ID,
          res_interaction_df$Systematic_ID
        )
      ],
    
    Interaction_padj =
      res_interaction_df$padj[
        match(
          res_WT_minus_vs_WT_plus_df$Systematic_ID,
          res_interaction_df$Systematic_ID
        )
      ],
    
    stringsAsFactors = FALSE
    
  )
  
  
  ############################################################
  # SAVE RESPONSE TABLE
  ############################################################
  
  write.csv(
    
    response_table,
    
    file.path(
      output_dir,
      "09_Response_analysis",
      "Phosphate_Response_Comparison.csv"
    ),
    
    row.names = FALSE
    
  )
  
  
  ############################################################
  # SIGNIFICANT INTERACTION GENES
  ############################################################
  
  significant_response_genes <- response_table %>%
    
    filter(
      !is.na(Interaction_padj),
      Interaction_padj < 0.05
    ) %>%
    
    arrange(
      Interaction_padj
    )
  
  
  write.csv(
    
    significant_response_genes,
    
    file.path(
      output_dir,
      "09_Response_analysis",
      "Significant_Phosphate_Response_Interaction_Genes.csv"
    ),
    
    row.names = FALSE
    
  )
  
  
  ############################################################
  # TOP 50 INTERACTION GENES
  ############################################################
  
  top_interaction_genes <- significant_response_genes %>%
    
    head(50)
  
  
  write.csv(
    
    top_interaction_genes,
    
    file.path(
      output_dir,
      "09_Response_analysis",
      "Top_50_Interaction_Genes.csv"
    ),
    
    row.names = FALSE
    
  )
  
  
  ############################################################
  # 29. PHOSPHATE RESPONSE HEATMAP
  ############################################################
  
  if (
    nrow(top_interaction_genes) >= 2
  ) {
    
    
    response_gene_ids <- intersect(
      
      top_interaction_genes$Systematic_ID,
      
      rownames(
        vst_matrix
      )
      
    )
    
    
    response_matrix <- condition_means[
      
      response_gene_ids,
      
      ,
      
      drop = FALSE
      
    ]
    
    
    rownames(
      response_matrix
    ) <- make.unique(
      
      get_gene_names(
        rownames(
          response_matrix
        )
      )
      
    )
    
    
    response_scaled <- t(
      
      scale(
        t(
          response_matrix
        )
      )
      
    )
    
    
    png(
      
      file.path(
        output_dir,
        "09_Response_analysis",
        "Heatmap_Top_Interaction_Genes.png"
      ),
      
      width = 1800,
      height = 2600,
      res = 300
      
    )
    
    
    pheatmap(
      
      response_scaled,
      
      scale = "none",
      
      cluster_rows = TRUE,
      
      cluster_cols = FALSE,
      
      show_rownames = TRUE,
      
      fontsize_row = 7,
      
      fontsize_col = 10,
      
      main =
        "Top genes with Genotype × Pi interaction"
      
    )
    
    
    dev.off()
    
  }
  
  
  ############################################################
  # 30. PHOSPHATE RESPONSE SCATTERPLOT
  ############################################################
  
  response_plot_df <- response_table %>%
    
    filter(
      !is.na(WT_response),
      !is.na(KCS1_response)
    )
  
  
  if (
    nrow(response_plot_df) > 0
  ) {
    
    
    p_response_scatter <- ggplot(
      
      response_plot_df,
      
      aes(
        x = WT_response,
        y = KCS1_response
      )
      
    ) +
      
      geom_point(
        alpha = 0.5,
        size = 2
      ) +
      
      geom_abline(
        slope = 1,
        intercept = 0,
        linetype = "dashed"
      ) +
      
      geom_vline(
        xintercept = 0,
        linetype = "dotted"
      ) +
      
      geom_hline(
        yintercept = 0,
        linetype = "dotted"
      ) +
      
      theme_bw() +
      
      labs(
        
        title =
          "WT vs KCS1 phosphate-starvation responses",
        
        x =
          "WT response: WT-Pi vs WT+Pi",
        
        y =
          "KCS1 response: KCS1-Pi vs KCS1+Pi"
        
      )
    
    
    ggsave(
      
      file.path(
        output_dir,
        "09_Response_analysis",
        "WT_vs_KCS1_Phosphate_Response.png"
      ),
      
      p_response_scatter,
      
      width = 9,
      height = 8,
      dpi = 300
      
    )
    
  }
  
  
  ############################################################
  # 31. PHOSPHATE RESPONSE PLOT
  #
  # Only significant interaction genes are shown.
  #
  ############################################################
  
  if (
    nrow(significant_response_genes) >= 1
  ) {
    
    
    response_long <- data.frame(
      
      Gene = rep(
        significant_response_genes$Gene_Name,
        2
      ),
      
      Genotype = rep(
        c(
          "WT",
          "KCS1"
        ),
        
        each =
          nrow(
            significant_response_genes
          )
      ),
      
      Response = c(
        
        significant_response_genes$WT_response,
        
        significant_response_genes$KCS1_response
        )
      ),
      
      stringsAsFactors = FALSE
      
    )
    
    
    response_plot <- ggplot(
      
      response_long,
      
      aes(
        x = Response,
        y = reorder(
          Gene,
          Response
        ),
        shape = Genotype
      )
      
    ) +
      
      geom_point(
        size = 3
      ) +
      
      geom_vline(
        xintercept = 0,
        linetype = "dashed"
      ) +
      
      theme_bw() +
      
      labs(
        
        title =
          "Phosphate-response differences in significant interaction genes",
        
        x =
          "Change in VST expression: -Pi minus +Pi",
        
        y =
          "Gene"
        
      )
    
    
    ##########################################################
    # CAP HEIGHT TO PREVENT ggsave ERROR
    ##########################################################
    
    response_height <- min(
      
      max(
        6,
        0.25 *
          nrow(
            significant_response_genes
          )
      ),
      
      48
      
    )
    
    
    ggsave(
      
      file.path(
        output_dir,
        "09_Response_analysis",
        "Phosphate_response_plot.png"
      ),
      
      response_plot,
      
      width = 10,
      
      height = response_height,
      
      dpi = 300
      
    )
    
  }
  
  
  ############################################################
  # 32. SUMMARY TABLE
  ############################################################
  
  summary_table <- data.frame(
    
    Comparison = names(
      all_results
    ),
    
    Total_genes_tested = sapply(
      
      all_results,
      
      nrow
      
    ),
    
    Significant_padj_0.05 = sapply(
      
      all_results,
      
      function(df) {
        
        sum(
          !is.na(df$padj) &
            df$padj < 0.05
        )
        
      }
      
    ),
    
    Significant_padj_0.05_log2FC_1 = sapply(
      
      all_results,
      
      function(df) {
        
        sum(
          !is.na(df$padj) &
            df$padj < 0.05 &
            !is.na(df$log2FoldChange) &
            abs(df$log2FoldChange) >= 1
        )
        
      }
      
    ),
    
    stringsAsFactors = FALSE
    
  )
  
  
  write.csv(
    
    summary_table,
    
    file.path(
      output_dir,
      "02_DE_results",
      "DE_summary.csv"
    ),
    
    row.names = FALSE
    
  )
  
  
  ############################################################
  # 33. COMPARISON SUMMARY PRINT
  ############################################################
  
  print("======================================================")
  print("SIGNIFICANT GENE COUNTS")
  print("======================================================")
  
  
  print(
    summary_table
  )
  
  
  ############################################################
  # 34. SAVE VST OBJECT
  ############################################################
  
  saveRDS(
    
    vsd,
    
    file.path(
      output_dir,
      "01_QC",
      "VST_object.rds"
    )
    
  )
  
  ############################################################
  # 22. VENN DIAGRAMS: PHOSPHATE-STARVATION RESPONSE
  #
  # All comparisons are:
  #       -Pi vs +Pi
  #
  # +Pi is the reference.
  #
  # Positive log2FC = UPREGULATED under -Pi
  # Negative log2FC = DOWNREGULATED under -Pi
  #
  # We compare:
  #
  # UPREGULATED:
  #   WT-only vs KCS1-only vs Common
  #
  # DOWNREGULATED:
  #   WT-only vs KCS1-only vs Common
  #
  # No cross-direction comparisons are made.
  ############################################################
  
  library(VennDiagram)
  library(grid)
  
  ############################################################
  # CREATE OUTPUT DIRECTORY
  ############################################################
  
  venn_dir <- file.path(
    output_dir,
    "04_Venn"
  )
  
  if (!dir.exists(venn_dir)) {
    dir.create(
      venn_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )
  }
  
  ############################################################
  # USE -Pi VS +Pi RESULTS
  ############################################################
  
  wt_response <- res_WT_minus_vs_WT_plus_df %>%
    filter(
      !is.na(padj),
      !is.na(log2FoldChange)
    )
  
  kcs1_response <- res_KCS1_minus_vs_KCS1_plus_df %>%
    filter(
      !is.na(padj),
      !is.na(log2FoldChange)
    )
  
  ############################################################
  # DEFINE UPREGULATED GENES UNDER -Pi
  ############################################################
  
  wt_up <- wt_response %>%
    filter(
      padj < 0.05,
      log2FoldChange >= 1
    ) %>%
    pull(Systematic_ID)
  
  kcs1_up <- kcs1_response %>%
    filter(
      padj < 0.05,
      log2FoldChange >= 1
    ) %>%
    pull(Systematic_ID)
  
  ############################################################
  # DEFINE DOWNREGULATED GENES UNDER -Pi
  ############################################################
  
  wt_down <- wt_response %>%
    filter(
      padj < 0.05,
      log2FoldChange <= -1
    ) %>%
    pull(Systematic_ID)
  
  kcs1_down <- kcs1_response %>%
    filter(
      padj < 0.05,
      log2FoldChange <= -1
    ) %>%
    pull(Systematic_ID)
  
  ############################################################
  # FIND COMMON AND STRAIN-SPECIFIC GENES
  ############################################################
  
  # UPREGULATED
  common_up <- intersect(
    wt_up,
    kcs1_up
  )
  
  wt_only_up <- setdiff(
    wt_up,
    kcs1_up
  )
  
  kcs1_only_up <- setdiff(
    kcs1_up,
    wt_up
  )
  
  # DOWNREGULATED
  common_down <- intersect(
    wt_down,
    kcs1_down
  )
  
  wt_only_down <- setdiff(
    wt_down,
    kcs1_down
  )
  
  kcs1_only_down <- setdiff(
    kcs1_down,
    wt_down
  )
  
  ############################################################
  # PRINT COUNTS
  ############################################################
  
  cat("\n======================================================\n")
  cat("PHOSPHATE-STARVATION RESPONSE: -Pi VS +Pi\n")
  cat("======================================================\n")
  
  cat("\nUPREGULATED UNDER -Pi\n")
  cat("----------------------\n")
  cat("WT total:       ", length(wt_up), "\n")
  cat("KCS1 total:     ", length(kcs1_up), "\n")
  cat("Common:         ", length(common_up), "\n")
  cat("WT-only:        ", length(wt_only_up), "\n")
  cat("KCS1-only:      ", length(kcs1_only_up), "\n")
  
  cat("\nDOWNREGULATED UNDER -Pi\n")
  cat("------------------------\n")
  cat("WT total:       ", length(wt_down), "\n")
  cat("KCS1 total:     ", length(kcs1_down), "\n")
  cat("Common:         ", length(common_down), "\n")
  cat("WT-only:        ", length(wt_only_down), "\n")
  cat("KCS1-only:      ", length(kcs1_only_down), "\n")
  
  ############################################################
  # VENN FUNCTION
  ############################################################
  
  make_venn <- function(
    set1,
    set2,
    name1,
    name2,
    plot_title,
    output_file
  ) {
    
    venn <- venn.diagram(
      x = list(
        set1 = set1,
        set2 = set2
      ),
      
      filename = NULL,
      
      width = 3200,
      height = 2800,
      resolution = 300,
      
      fill = c(
        "grey80",
        "grey60"
      ),
      
      alpha = 0.5,
      
      cex = 1.8,
      fontface = "bold",
      
      category.names = c(
        name1,
        name2
      ),
      
      cat.cex = 1.2,
      cat.fontface = "bold",
      
      cat.dist = c(
        0.08,
        0.08
      ),
      
      cat.pos = c(
        -15,
        15
      ),
      
      margin = 0.25,
      
      main = plot_title,
      main.cex = 1.5,
      main.fontface = "bold"
    )
    
    png(
      filename = output_file,
      width = 3200,
      height = 2800,
      res = 300
    )
    
    grid.newpage()
    
    pushViewport(
      viewport(
        x = 0.5,
        y = 0.47,
        width = 0.78,
        height = 0.70
      )
    )
    
    grid.draw(venn)
    
    popViewport()
    
    dev.off()
    
    cat(
      "Saved:",
      output_file,
      "\n"
    )
  }
  
  ############################################################
  # VENN 1: UPREGULATED UNDER -Pi
  ############################################################
  
  make_venn(
    set1 = wt_up,
    set2 = kcs1_up,
    
    name1 = "WT Upregulated",
    name2 = "KCS1 Upregulated",
    
    plot_title =
      "Genes Upregulated Under -Pi",
    
    output_file =
      file.path(
        venn_dir,
        "Venn_Common_Upregulated_WT_vs_KCS1.png"
      )
  )
  
  ############################################################
  # VENN 2: DOWNREGULATED UNDER -Pi
  ############################################################
  
  make_venn(
    set1 = wt_down,
    set2 = kcs1_down,
    
    name1 = "WT Downregulated",
    name2 = "KCS1 Downregulated",
    
    plot_title =
      "Genes Downregulated Under -Pi",
    
    output_file =
      file.path(
        venn_dir,
        "Venn_Common_Downregulated_WT_vs_KCS1.png"
      )
  )
  
  ############################################################
  # SAVE COMMON / STRAIN-SPECIFIC GENE LISTS
  ############################################################
  
  write.csv(
    data.frame(
      Systematic_ID = common_up,
      Gene_Name = get_gene_names(common_up)
    ),
    file.path(
      venn_dir,
      "Common_Upregulated_WT_KCS1.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(
    data.frame(
      Systematic_ID = wt_only_up,
      Gene_Name = get_gene_names(wt_only_up)
    ),
    file.path(
      venn_dir,
      "WT_Only_Upregulated.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(
    data.frame(
      Systematic_ID = kcs1_only_up,
      Gene_Name = get_gene_names(kcs1_only_up)
    ),
    file.path(
      venn_dir,
      "KCS1_Only_Upregulated.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(
    data.frame(
      Systematic_ID = common_down,
      Gene_Name = get_gene_names(common_down)
    ),
    file.path(
      venn_dir,
      "Common_Downregulated_WT_KCS1.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(
    data.frame(
      Systematic_ID = wt_only_down,
      Gene_Name = get_gene_names(wt_only_down)
    ),
    file.path(
      venn_dir,
      "WT_Only_Downregulated.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(
    data.frame(
      Systematic_ID = kcs1_only_down,
      Gene_Name = get_gene_names(kcs1_only_down)
    ),
    file.path(
      venn_dir,
      "KCS1_Only_Downregulated.csv"
    ),
    row.names = FALSE
  )
  
  ############################################################
  # FINISHED
  ############################################################
  
  cat("\n======================================================\n")
  cat("SECTION 22 COMPLETE\n")
  cat("======================================================\n")
  cat("Venn diagrams saved in:\n")
  cat(venn_dir, "\n")
  
  ############################################################
  # 23. TOP 50 GENES FOR WT VS KCS1 UNDER -Pi
  #
  # +Pi is the reference condition.
  #
  # Positive log2FC = higher under -Pi
  # Negative log2FC = lower under -Pi
  ############################################################
  
  ############################################################
  # DEFINE RESPONSE TABLES
  ############################################################
  
  wt_response <- res_WT_minus_vs_WT_plus_df %>%
    filter(
      !is.na(padj),
      !is.na(log2FoldChange)
    )
  
  kcs1_response <- res_KCS1_minus_vs_KCS1_plus_df %>%
    filter(
      !is.na(padj),
      !is.na(log2FoldChange)
    )
  
  ############################################################
  # DEFINE SIGNIFICANT GENES
  ############################################################
  
  # WT genes increased under -Pi
  wt_up <- wt_response %>%
    filter(
      padj < 0.05,
      log2FoldChange >= 1
    ) %>%
    pull(Systematic_ID)
  
  # WT genes decreased under -Pi
  wt_down <- wt_response %>%
    filter(
      padj < 0.05,
      log2FoldChange <= -1
    ) %>%
    pull(Systematic_ID)
  
  # KCS1 genes increased under -Pi
  kcs1_up <- kcs1_response %>%
    filter(
      padj < 0.05,
      log2FoldChange >= 1
    ) %>%
    pull(Systematic_ID)
  
  # KCS1 genes decreased under -Pi
  kcs1_down <- kcs1_response %>%
    filter(
      padj < 0.05,
      log2FoldChange <= -1
    ) %>%
    pull(Systematic_ID)
  
  ############################################################
  # FIND WT/KCS1 OVERLAPS
  ############################################################
  
  # Upregulated under -Pi
  common_up <- intersect(
    wt_up,
    kcs1_up
  )
  
  wt_only_up <- setdiff(
    wt_up,
    kcs1_up
  )
  
  kcs1_only_up <- setdiff(
    kcs1_up,
    wt_up
  )
  
  # Downregulated under -Pi
  common_down <- intersect(
    wt_down,
    kcs1_down
  )
  
  wt_only_down <- setdiff(
    wt_down,
    kcs1_down
  )
  
  kcs1_only_down <- setdiff(
    kcs1_down,
    wt_down
  )
  
  ############################################################
  # FUNCTION FOR WT/KCS1-ONLY GENES
  ############################################################
  
  get_top50 <- function(
    result_df,
    gene_ids,
    category
  ) {
    
    result_df %>%
      filter(
        Systematic_ID %in% gene_ids
      ) %>%
      arrange(
        padj,
        desc(abs(log2FoldChange))
      ) %>%
      mutate(
        Category = category
      ) %>%
      select(
        Category,
        Systematic_ID,
        Gene_Name,
        log2FoldChange,
        padj
      ) %>%
      slice_head(n = 50)
  }
  
  ############################################################
  # TOP 50 WT-ONLY UPREGULATED
  ############################################################
  
  top50_wt_only_up <- get_top50(
    wt_response,
    wt_only_up,
    "WT-only Upregulated under -Pi"
  )
  
  ############################################################
  # TOP 50 KCS1-ONLY UPREGULATED
  ############################################################
  
  top50_kcs1_only_up <- get_top50(
    kcs1_response,
    kcs1_only_up,
    "KCS1-only Upregulated under -Pi"
  )
  
  ############################################################
  # TOP 50 COMMON UPREGULATED
  #
  # Rank by WT padj, then show BOTH WT and KCS1 statistics.
  ############################################################
  
  top50_common_up <- wt_response %>%
    filter(
      Systematic_ID %in% common_up
    ) %>%
    select(
      Systematic_ID,
      Gene_Name,
      WT_log2FoldChange = log2FoldChange,
      WT_padj = padj
    ) %>%
    inner_join(
      kcs1_response %>%
        filter(
          Systematic_ID %in% common_up
        ) %>%
        select(
          Systematic_ID,
          KCS1_log2FoldChange = log2FoldChange,
          KCS1_padj = padj
        ),
      by = "Systematic_ID"
    ) %>%
    arrange(
      WT_padj,
      desc(abs(WT_log2FoldChange))
    ) %>%
    slice_head(n = 50)
  
  ############################################################
  # TOP 50 WT-ONLY DOWNREGULATED
  ############################################################
  
  top50_wt_only_down <- get_top50(
    wt_response,
    wt_only_down,
    "WT-only Downregulated under -Pi"
  )
  
  ############################################################
  # TOP 50 KCS1-ONLY DOWNREGULATED
  ############################################################
  
  top50_kcs1_only_down <- get_top50(
    kcs1_response,
    kcs1_only_down,
    "KCS1-only Downregulated under -Pi"
  )
  
  ############################################################
  # TOP 50 COMMON DOWNREGULATED
  #
  # Rank by WT padj, then show BOTH WT and KCS1 statistics.
  ############################################################
  
  top50_common_down <- wt_response %>%
    filter(
      Systematic_ID %in% common_down
    ) %>%
    select(
      Systematic_ID,
      Gene_Name,
      WT_log2FoldChange = log2FoldChange,
      WT_padj = padj
    ) %>%
    inner_join(
      kcs1_response %>%
        filter(
          Systematic_ID %in% common_down
        ) %>%
        select(
          Systematic_ID,
          KCS1_log2FoldChange = log2FoldChange,
          KCS1_padj = padj
        ),
      by = "Systematic_ID"
    ) %>%
    arrange(
      WT_padj,
      desc(abs(WT_log2FoldChange))
    ) %>%
    slice_head(n = 50)
  
  ############################################################
  # SAVE CSV FILES
  ############################################################
  
  write.csv(
    top50_wt_only_up,
    file.path(
      venn_dir,
      "Top50_WT_Only_Upregulated_under_minusPi.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(
    top50_kcs1_only_up,
    file.path(
      venn_dir,
      "Top50_KCS1_Only_Upregulated_under_minusPi.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(
    top50_common_up,
    file.path(
      venn_dir,
      "Top50_Common_Upregulated_under_minusPi.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(
    top50_wt_only_down,
    file.path(
      venn_dir,
      "Top50_WT_Only_Downregulated_under_minusPi.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(
    top50_kcs1_only_down,
    file.path(
      venn_dir,
      "Top50_KCS1_Only_Downregulated_under_minusPi.csv"
    ),
    row.names = FALSE
  )
  
  write.csv(
    top50_common_down,
    file.path(
      venn_dir,
      "Top50_Common_Downregulated_under_minusPi.csv"
    ),
    row.names = FALSE
  )
  
  ############################################################
  # PRINT COUNTS
  ############################################################
  
  cat("\n======================================================\n")
  cat("TOP 50 GENES: -Pi VS +Pi\n")
  cat("======================================================\n")
  
  cat("\nUPREGULATED UNDER -Pi\n")
  cat(
    "WT-only:   ",
    nrow(top50_wt_only_up),
    "\n"
  )
  cat(
    "KCS1-only: ",
    nrow(top50_kcs1_only_up),
    "\n"
  )
  cat(
    "Common:    ",
    nrow(top50_common_up),
    "\n"
  )
  
  cat("\nDOWNREGULATED UNDER -Pi\n")
  cat(
    "WT-only:   ",
    nrow(top50_wt_only_down),
    "\n"
  )
  cat(
    "KCS1-only: ",
    nrow(top50_kcs1_only_down),
    "\n"
  )
  cat(
    "Common:    ",
    nrow(top50_common_down),
    "\n"
  )
  
  cat("\nFiles saved to:\n")
  cat(venn_dir, "\n")
  
  ############################################################
  # 24. CREATE EXCEL FILE OF TOP 50 VENN GENES
  ############################################################
  
  library(openxlsx)
  
  excel_file <- file.path(
    venn_dir,
    "Top50_WT_vs_KCS1_under_minusPi.xlsx"
  )
  
  ############################################################
  # CREATE WORKBOOK
  ############################################################
  
  wb <- createWorkbook()
  
  ############################################################
  # 1. WT-ONLY UPREGULATED
  ############################################################
  
  addWorksheet(
    wb,
    "WT-only Up"
  )
  
  writeData(
    wb,
    "WT-only Up",
    top50_wt_only_up
  )
  
  ############################################################
  # 2. KCS1-ONLY UPREGULATED
  ############################################################
  
  addWorksheet(
    wb,
    "KCS1-only Up"
  )
  
  writeData(
    wb,
    "KCS1-only Up",
    top50_kcs1_only_up
  )
  
  ############################################################
  # 3. COMMON UPREGULATED
  ############################################################
  
  addWorksheet(
    wb,
    "Common Up"
  )
  
  writeData(
    wb,
    "Common Up",
    top50_common_up
  )
  
  ############################################################
  # 4. WT-ONLY DOWNREGULATED
  ############################################################
  
  addWorksheet(
    wb,
    "WT-only Down"
  )
  
  writeData(
    wb,
    "WT-only Down",
    top50_wt_only_down
  )
  
  ############################################################
  # 5. KCS1-ONLY DOWNREGULATED
  ############################################################
  
  addWorksheet(
    wb,
    "KCS1-only Down"
  )
  
  writeData(
    wb,
    "KCS1-only Down",
    top50_kcs1_only_down
  )
  
  ############################################################
  # 6. COMMON DOWNREGULATED
  ############################################################
  
  addWorksheet(
    wb,
    "Common Down"
  )
  
  writeData(
    wb,
    "Common Down",
    top50_common_down
  )
  
  ############################################################
  # FORMAT ALL SHEETS
  ############################################################
  
  for (sheet in names(wb)) {
    
    # Bold header
    addStyle(
      wb,
      sheet = sheet,
      style = createStyle(
        textDecoration = "bold"
      ),
      rows = 1,
      cols = 1:10,
      gridExpand = TRUE
    )
    
    # Freeze header row
    freezePane(
      wb,
      sheet = sheet,
      firstRow = TRUE
    )
    
    # Auto-width columns
    setColWidths(
      wb,
      sheet = sheet,
      cols = 1:10,
      widths = "auto"
    )
  }
  
  ############################################################
  # SAVE EXCEL FILE
  ############################################################
  
  saveWorkbook(
    wb,
    excel_file,
    overwrite = TRUE
  )
  
  ############################################################
  # CONFIRM
  ############################################################
  
  cat("\n======================================================\n")
  cat("EXCEL FILE CREATED\n")
  cat("======================================================\n")
  cat("File:\n")
  cat(excel_file, "\n")
  
  ############################################################
  # 36. FINAL SUMMARY
  ############################################################
  
  print("")
  print("======================================================")
  print("RNA-seq ANALYSIS COMPLETE")
  print("======================================================")
  
  
  print(
    paste(
      "Total genes analyzed:",
      nrow(dds)
    )
  )
  
  
  print(
    paste(
      "WT+Pi replicates:",
      sum(
        metadata$Condition == "WT_plus"
      )
    )
  )
  
  
  print(
    paste(
      "WT-Pi replicates:",
      sum(
        metadata$Condition == "WT_minus"
      )
    )
  )
  
  
  print(
    paste(
      "KCS1+Pi replicates:",
      sum(
        metadata$Condition == "KCS1_plus"
      )
    )
  )
  
  
  print(
    paste(
      "KCS1-Pi replicates:",
      sum(
        metadata$Condition == "KCS1_minus"
      )
    )
  )
  
  
  print("")
  
  print(
    paste(
      "Genes with significant Genotype × Pi interaction:",
      sum(
        !is.na(
          res_interaction_df$padj
        ) &
          res_interaction_df$padj < 0.05
      )
    )
  )
  
  
  print("")
  
  print("Output directory:")
  
  print(
    normalizePath(
      output_dir
    )
  )
  
  
  print("")
  
  print("======================================================")
  print("PAIRWISE COMPARISONS GENERATED")
  print("======================================================")
  
  print(
    "1. WT-Pi vs WT+Pi"
  )
  
  print(
    "2. KCS1+Pi vs WT+Pi"
  )
  
  print(
    "3. KCS1-Pi vs WT+Pi"
  )
  
  print(
    "4. KCS1+Pi vs WT-Pi"
  )
  
  print(
    "5. KCS1-Pi vs WT-Pi"
  )
  
  print(
    "6. KCS1-Pi vs KCS1+Pi"
  )
  
  print(
    "7. Genotype × Pi interaction"
  )
  
  
  print("")
  print("======================================================")
  print("END OF ANALYSIS")
  print("======================================================")
  
  
  summary(res_df$padj)
  min(res_df$padj, na.rm = TRUE)
  
  
  # Read CGD GO annotation file
  gene_association <- read.delim(
    "gene_association.cgd",
    header = FALSE,
    stringsAsFactors = FALSE,
    sep = "\t",
    quote = "",
    comment.char = "!"
  )
  
  # Read CGD feature annotation file
  cgd_features <- read.delim(
    "cgd_features.tab",
    header = FALSE,
    stringsAsFactors = FALSE,
    sep = "\t",
    quote = "",
    comment.char = ""
  )
  
  # Inspect both files
  head(gene_association)
  head(cgd_features)
  
  # Check dimensions
  dim(gene_association)
  dim(cgd_features)
  
  # Look at the actual non-comment lines of cgd_features.tab
  cgd_lines <- readLines("cgd_features.tab")
  
  # Remove comment/header lines
  cgd_lines <- cgd_lines[!grepl("^!", cgd_lines)]
  
  # Look at the first few actual data lines
  head(cgd_lines)
  
  # Check the IDs used in your RNA-seq results
  head(rownames(counts_filtered))
  head(res_WT_minus_vs_WT_plus_df$Systematic_ID)
  
  # Check the systematic IDs contained in the GO annotation file
  head(gene_association$V11)
  
  
  
  ############################################################
  # CHECK GO RESULTS
  ############################################################
  
  cat("\n")
  cat("============================================\n")
  cat("GO RESULT CHECK\n")
  cat("============================================\n")
  
  for (result_name in names(go_results)) {
    
    result <- go_results[[result_name]]
    
    if (is.null(result)) {
      
      cat(result_name, ": NULL\n")
      
    } else {
      
      result_df <- as.data.frame(result)
      
      cat(
        result_name,
        ":",
        nrow(result_df),
        "GO terms\n"
      )
    }
  }
  
  ############################################################
  # RE-RUN GO ENRICHMENT
  ############################################################
  
  go_results <- list()
  
  for (gene_set_name in names(gene_sets_cgd)) {
    
    cat("\n")
    cat("============================================\n")
    cat("PROCESSING:", gene_set_name, "\n")
    cat("============================================\n")
    
    # Biological Process
    bp_name <- paste0(gene_set_name, "_BP")
    
    go_results[[bp_name]] <- run_GO(
      genes = gene_sets_cgd[[gene_set_name]],
      gene_set_name = gene_set_name,
      ontology = "P"
    )
    
    # Molecular Function
    mf_name <- paste0(gene_set_name, "_MF")
    
    go_results[[mf_name]] <- run_GO(
      genes = gene_sets_cgd[[gene_set_name]],
      gene_set_name = gene_set_name,
      ontology = "F"
    )
    
    # Cellular Component
    cc_name <- paste0(gene_set_name, "_CC")
    
    go_results[[cc_name]] <- run_GO(
      genes = gene_sets_cgd[[gene_set_name]],
      gene_set_name = gene_set_name,
      ontology = "C"
    )
  }
  
  cat("\n")
  cat("============================================\n")
  cat("GO ENRICHMENT COMPLETE\n")
  cat("============================================\n")
  
  cat(
    "Number of GO result objects:",
    length(go_results),
    "\n"
  )
  
  ############################################################
  # REBUILD GENE ASSOCIATION TABLE
  ############################################################
  
  gene_association <- read.delim(
    text = paste(
      go_data_lines,
      collapse = "\n"
    ),
    header = FALSE,
    stringsAsFactors = FALSE,
    sep = "\t",
    quote = "",
    fill = TRUE
  )
  
  cat(
    "Gene association dimensions:",
    nrow(gene_association),
    "rows x",
    ncol(gene_association),
    "columns\n"
  )
  
  
  ############################################################
  # REBUILD GO ANNOTATIONS USING CORRECT GAF COLUMNS
  ############################################################
  
  go_annotations <- gene_association %>%
    dplyr::transmute(
      CGD_ID = trimws(as.character(V2)),
      Gene_Name = trimws(as.character(V3)),
      GO_ID = trimws(as.character(V5)),
      Evidence = trimws(as.character(V7)),
      Ontology = trimws(as.character(V9))
    ) %>%
    dplyr::filter(
      !is.na(CGD_ID),
      CGD_ID != "",
      !is.na(GO_ID),
      GO_ID != "",
      grepl("^GO:", GO_ID),
      Ontology %in% c("P", "F", "C")
    ) %>%
    dplyr::distinct()
  
  ############################################################
  # CHECK GO ANNOTATIONS
  ############################################################
  
  cat("GO annotations:", nrow(go_annotations), "\n")
  
  cat(
    "Unique CGD IDs:",
    length(unique(go_annotations$CGD_ID)),
    "\n"
  )
  
  cat(
    "Unique GO terms:",
    length(unique(go_annotations$GO_ID)),
    "\n"
  )
  
  cat("\nOntology counts:\n")
  print(table(go_annotations$Ontology, useNA = "ifany"))
  
  ############################################################
  # CHECK OVERLAP WITH RNA-SEQ BACKGROUND
  ############################################################
  
  cat(
    "\nRNA-seq background genes:",
    length(background_cgd),
    "\n"
  )
  
  cat(
    "Background genes with GO annotations:",
    sum(background_cgd %in% go_annotations$CGD_ID),
    "\n"
  )
  
  
  ############################################################
  # REBUILD GO TERM MAPPINGS
  ############################################################
  
  TERM2GENE_BP <- go_annotations %>%
    dplyr::filter(Ontology == "P") %>%
    
    
    ############################################################
  # VERIFY GO MAPPING TABLES
  ############################################################
  
  cat("BP mappings:", nrow(TERM2GENE_BP), "\n")
  cat("MF mappings:", nrow(TERM2GENE_MF), "\n")
  cat("CC mappings:", nrow(TERM2GENE_CC), "\n")
  cat("GO term names:", nrow(TERM2NAME), "\n")
  
  ############################################################
  # GO ENRICHMENT FUNCTION
  ############################################################
  
  run_GO <- function(
    genes,
    gene_set_name,
    ontology
  ) {
    
    ##########################################################
    # SELECT ONTOLOGY
    ##########################################################
    
    if (ontology == "P") {
      term2gene <- TERM2GENE_BP
    } else if (ontology == "F") {
      term2gene <- TERM2GENE_MF
    } else if (ontology == "C") {
      term2gene <- TERM2GENE_CC
    } else {
      stop("Ontology must be P, F, or C.")
    }
    
    ##########################################################
    # CLEAN MAPPING TABLE
    ##########################################################
    
    term2gene <- as.data.frame(
      term2gene,
      stringsAsFactors = FALSE
    )
    
    colnames(term2gene) <- c(
      "GO_ID",
      "CGD_ID"
    )
    
    ##########################################################
    # MATCH GENES TO GO ANNOTATIONS
    ##########################################################
    
    genes <- intersect(
      genes,
      unique(term2gene$CGD_ID)
    )
    
    universe <- intersect(
      background_cgd,
      unique(term2gene$CGD_ID)
    )
    
    cat(
      gene_set_name,
      "-",
      ontology,
      ":",
      length(genes),
      "GO-mapped genes;",
      length(universe),
      "background genes\n"
    )
    
    ##########################################################
    # STOP IF NOT ENOUGH GENES
    ##########################################################
    
    if (
      length(genes) == 0 ||
      length(universe) == 0
    ) {
      return(NULL)
    }
    
    ##########################################################
    # RUN ENRICHMENT
    ##########################################################
    
    ego <- clusterProfiler::enricher(
      gene = genes,
      universe = universe,
      TERM2GENE = term2gene,
      TERM2NAME = TERM2NAME,
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      minGSSize = 5,
      maxGSSize = 500
    )
    
    ##########################################################
    # ADD METADATA
    ##########################################################
    
    if (!is.null(ego)) {
      
      if (nrow(as.data.frame(ego)) > 0) {
        
        ego@result$Gene_Set <- gene_set_name
        ego@result$Ontology <- ontology
        
      }
    }
    
    return(ego)
  }
  
  ############################################################
  # RUN ALL GO ENRICHMENT ANALYSES
  ############################################################
  
  go_results <- list()
  
  for (gene_set_name in names(gene_sets_cgd)) {
    
    ##########################################################
    # BIOLOGICAL PROCESS
    ##########################################################
    
    bp_name <- paste0(
      gene_set_name,
      "_BP"
    )
    
    go_results[[bp_name]] <- run_GO(
      genes = gene_sets_cgd[[gene_set_name]],
      gene_set_name = gene_set_name,
      ontology = "P"
    )
    
    ##########################################################
    # MOLECULAR FUNCTION
    ##########################################################
    
    mf_name <- paste0(
      gene_set_name,
      "_MF"
    )
    
    go_results[[mf_name]] <- run_GO(
      genes = gene_sets_cgd[[gene_set_name]],
      gene_set_name = gene_set_name,
      ontology = "F"
    )
    
    ##########################################################
    # CELLULAR COMPONENT
    ##########################################################
    
    cc_name <- paste0(
      gene_set_name,
      "_CC"
    )
    
    go_results[[cc_name]] <- run_GO(
      genes = gene_sets_cgd[[gene_set_name]],
      gene_set_name = gene_set_name,
      ontology = "C"
    )
  }
  
  ############################################################
  # CHECK NUMBER OF RESULTS
  ############################################################
  
  cat(
    "\nNumber of GO result objects:",
    length(go_results),
    "\n"
  )
  
  ############################################################
  # CHECK GO ENRICHMENT RESULTS
  ############################################################
  
  for (result_name in names(go_results)) {
    
    result <- go_results[[result_name]]
    
    if (is.null(result)) {
      cat(result_name, ": NULL\n")
      next
    }
    
    result_df <- as.data.frame(result)
    
    sig <- result_df %>%
      dplyr::filter(
        !is.na(p.adjust),
        p.adjust < 0.05
      )
    
    cat(
      result_name,
      ":",
      nrow(result_df),
      "total terms;",
      nrow(sig),
      "significant terms\n"
    )
  }
  
  
  ############################################################
  # PRINT ALL SIGNIFICANT GO TERMS
  ############################################################
  
  for (result_name in names(go_results)) {
    
    result <- go_results[[result_name]]
    
    if (is.null(result)) {
      next
    }
    
    result_df <- as.data.frame(result)
    
    sig <- result_df %>%
      dplyr::filter(
        !is.na(p.adjust),
        p.adjust < 0.05
      ) %>%
      dplyr::arrange(p.adjust)
    
    if (nrow(sig) == 0) {
      next
    }
    
    cat("\n")
    cat("============================================================\n")
    cat(result_name, "\n")
    cat("============================================================\n")
    
    print(
      sig %>%
        dplyr::select(
          ID,
          Description,
          GeneRatio,
          BgRatio,
          Count,
          pvalue,
          p.adjust
        )
    )
  }
  
  
  ############################################################
  # CREATE COMBINED GO RESULTS TABLE
  ############################################################
  
  combined_GO <- dplyr::bind_rows(
    lapply(
      go_results,
      function(x) {
        
        if (is.null(x)) {
          return(NULL)
        }
        
        df <- as.data.frame(x)
        
        if (nrow(df) == 0) {
          return(NULL)
        }
        
        df
      }
    )
  )
  
  ############################################################
  # KEEP SIGNIFICANT RESULTS
  ############################################################
  
  combined_GO_sig <- combined_GO %>%
    dplyr::filter(
      !is.na(p.adjust),
      p.adjust < 0.05
    ) %>%
    dplyr::arrange(
      Ontology,
      Gene_Set,
      p.adjust
    )
  
  cat(
    "Total significant GO results:",
    nrow(combined_GO_sig),
    "\n"
  )
  
  ############################################################
  # CREATE EXCEL WORKBOOK
  ############################################################
  
  go_output_dir <- "05_GO_Analysis"
  
  if (!dir.exists(go_output_dir)) {
    dir.create(
      go_output_dir,
      recursive = TRUE
    )
  }
  
  wb <- openxlsx::createWorkbook()
  
  ############################################################
  # ALL SIGNIFICANT GO TERMS
  ############################################################
  
  openxlsx::addWorksheet(
    wb,
    "All_Significant_GO"
  )
  
  openxlsx::writeData(
    wb,
    "All_Significant_GO",
    combined_GO_sig
  )
  
  ############################################################
  # BIOLOGICAL PROCESS
  ############################################################
  
  bp_results <- combined_GO_sig %>%
    dplyr::filter(
      Ontology == "P"
    )
  
  openxlsx::addWorksheet(
    wb,
    "Biological_Process"
  )
  
  openxlsx::writeData(
    wb,
    "Biological_Process",
    bp_results
  )
  
  ############################################################
  # MOLECULAR FUNCTION
  ############################################################
  
  mf_results <- combined_GO_sig %>%
    dplyr::filter(
      Ontology == "F"
    )
  
  openxlsx::addWorksheet(
    wb,
    "Molecular_Function"
  )
  
  openxlsx::writeData(
    wb,
    "Molecular_Function",
    mf_results
  )
  
  ############################################################
  # CELLULAR COMPONENT
  ############################################################
  
  cc_results <- combined_GO_sig %>%
    dplyr::filter(
      Ontology == "C"
    )
  
  openxlsx::addWorksheet(
    wb,
    "Cellular_Component"
  )
  
  openxlsx::writeData(
    wb,
    "Cellular_Component",
    cc_results
  )
  
  ############################################################
  # ONE SHEET PER GENE SET
  ############################################################
  
  for (gene_set_name in names(gene_sets_cgd)) {
    
    gene_set_results <- combined_GO_sig %>%
      dplyr::filter(
        Gene_Set == gene_set_name
      )
    
    if (nrow(gene_set_results) == 0) {
      next
    }
    
    sheet_name <- substr(
      gene_set_name,
      1,
      31
    )
    
    openxlsx::addWorksheet(
      wb,
      sheet_name
    )
    
    openxlsx::writeData(
      wb,
      sheet_name,
      gene_set_results
    )
  }
  
  ############################################################
  # SAVE WORKBOOK
  ############################################################
  
  excel_file <- file.path(
    go_output_dir,
    "GO_Enrichment_Results.xlsx"
  )
  
  openxlsx::saveWorkbook(
    wb,
    excel_file,
    overwrite = TRUE
  )
  
  cat(
    "Excel workbook saved to:",
    excel_file,
    "\n"
  )
  
  ############################################################
  # CREATE GO DOTPLOTS
  ############################################################
  
  go_plot_dir <- file.path(
    go_output_dir,
    "Plots"
  )
  
  if (!dir.exists(go_plot_dir)) {
    dir.create(
      go_plot_dir,
      recursive = TRUE
    )
  }
  
  ############################################################
  # GENERATE ONE PLOT PER GENE SET
  ############################################################
  
  for (gene_set_name in names(gene_sets_cgd)) {
    
    ##########################################################
    # GET RESULTS FOR THIS GENE SET
    ##########################################################
    
    result_df <- combined_GO_sig %>%
      dplyr::filter(
        Gene_Set == gene_set_name
      ) %>%
      dplyr::arrange(
        p.adjust
      )
    
    if (nrow(result_df) == 0) {
      cat(
        "No significant GO terms:",
        gene_set_name,
        "\n"
      )
      next
    }
    
    ##########################################################
    # KEEP TOP 15 TERMS
    ##########################################################
    
    result_df <- result_df %>%
      dplyr::slice_head(
        n = 15
      )
    
    ##########################################################
    # CREATE LABEL
    ##########################################################
    
    result_df <- result_df %>%
      dplyr::mutate(
        Description = factor(
          Description,
          levels = rev(Description)
        ),
        neg_log10_padj = -log10(p.adjust)
      )
    
    ##########################################################
    # PLOT
    ##########################################################
    
    p <- ggplot2::ggplot(
      result_df,
      ggplot2::aes(
        x = neg_log10_padj,
        y = Description,
        size = Count,
        color = Ontology
      )
    ) +
      ggplot2::geom_point() +
      ggplot2::labs(
        title = paste0(
          "GO Enrichment: ",
          gene_set_name
        ),
        x = "-log10(adjusted p-value)",
        y = "GO Term",
        size = "Gene Count",
        color = "Ontology"
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          hjust = 0.5,
          face = "bold"
        ),
        axis.text.y = ggplot2::element_text(
          size = 9
        ),
        legend.position = "right"
      )
    
    ##########################################################
    # SAVE
    ##########################################################
    
    plot_name <- gsub(
      "[^A-Za-z0-9_-]",
      "_",
      gene_set_name
    )
    
    output_file <- file.path(
      go_plot_dir,
      paste0(
        plot_name,
        "_GO_dotplot.png"
      )
    )
    
    ggplot2::ggsave(
      filename = output_file,
      plot = p,
      width = 11,
      height = 8,
      dpi = 300
    )
    
    cat(
      "Saved:",
      output_file,
      "\n"
    )
  }
  
  ############################################################
  # GO ANALYSIS: PHOSPHATE STARVATION RESPONSE
  #
  # WT -Pi vs WT +Pi
  # KCS1 -Pi vs KCS1 +Pi
  #
  # +Pi is the baseline
  ############################################################
  
  ############################################################
  # SIGNIFICANCE THRESHOLDS
  ############################################################
  
  padj_cutoff <- 0.05
  log2fc_cutoff <- 1
  
  ############################################################
  # WT PHOSPHATE STARVATION
  ############################################################
  
  WT_starvation <- res_WT_minus_vs_WT_plus_df %>%
    dplyr::filter(
      !is.na(padj),
      !is.na(log2FoldChange)
    )
  
  WT_starvation_up <- WT_starvation %>%
    dplyr::filter(
      padj < padj_cutoff,
      log2FoldChange >= log2fc_cutoff
    ) %>%
    dplyr::pull(Systematic_ID) %>%
    unique()
  
  WT_starvation_down <- WT_starvation %>%
    dplyr::filter(
      padj < padj_cutoff,
      log2FoldChange <= -log2fc_cutoff
    ) %>%
    dplyr::pull(Systematic_ID) %>%
    unique()
  
  ############################################################
  # KCS1 PHOSPHATE STARVATION
  ############################################################
  
  KCS1_starvation <- res_KCS1_minus_vs_KCS1_plus_df %>%
    dplyr::filter(
      !is.na(padj),
      !is.na(log2FoldChange)
    )
  
  KCS1_starvation_up <- KCS1_starvation %>%
    dplyr::filter(
      padj < padj_cutoff,
      log2FoldChange >= log2fc_cutoff
    ) %>%
    dplyr::pull(Systematic_ID) %>%
    unique()
  
  KCS1_starvation_down <- KCS1_starvation %>%
    dplyr::filter(
      padj < padj_cutoff,
      log2FoldChange <= -log2fc_cutoff
    ) %>%
    dplyr::pull(Systematic_ID) %>%
    unique()
  
  ############################################################
  # CHECK RNA-SEQ GENE COUNTS
  ############################################################
  
  cat(
    "WT starvation Up:",
    length(WT_starvation_up),
    "\n"
  )
  
  cat(
    "WT starvation Down:",
    length(WT_starvation_down),
    "\n"
  )
  
  cat(
    "KCS1 starvation Up:",
    length(KCS1_starvation_up),
    "\n"
  )
  
  cat(
    "KCS1 starvation Down:",
    length(KCS1_starvation_down),
    "\n"
  )
  
  
  ############################################################
  # MAP PHOSPHATE-STARVATION GENE SETS TO CGD IDs
  ############################################################
  
  convert_to_cgd <- function(genes) {
    
    rna_mapping %>%
      dplyr::filter(
        RNAseq_ID %in% genes,
        !is.na(CGD_ID),
        CGD_ID != ""
      ) %>%
      dplyr::pull(CGD_ID) %>%
      unique()
  }
  
  ############################################################
  # CONVERT ALL FOUR GENE SETS
  ############################################################
  
  WT_starvation_up_cgd <- convert_to_cgd(
    WT_starvation_up
  )
  
  WT_starvation_down_cgd <- convert_to_cgd(
    WT_starvation_down
  )
  
  KCS1_starvation_up_cgd <- convert_to_cgd(
    KCS1_starvation_up
  )
  
  KCS1_starvation_down_cgd <- convert_to_cgd(
    KCS1_starvation_down
  )
  
  ############################################################
  # CHECK MAPPING
  ############################################################
  
  cat(
    "WT starvation Up:",
    length(WT_starvation_up),
    "RNA genes ->",
    length(WT_starvation_up_cgd),
    "CGD genes\n"
  )
  
  cat(
    "WT starvation Down:",
    length(WT_starvation_down),
    "RNA genes ->",
    length(WT_starvation_down_cgd),
    "CGD genes\n"
  )
  
  cat(
    "KCS1 starvation Up:",
    length(KCS1_starvation_up),
    "RNA genes ->",
    length(KCS1_starvation_up_cgd),
    "CGD genes\n"
  )
  
  cat(
    "KCS1 starvation Down:",
    length(KCS1_starvation_down),
    "RNA genes ->",
    length(KCS1_starvation_down_cgd),
    "CGD genes\n"
  )
  
  
  ############################################################
  # DEFINE PHOSPHATE-STARVATION GO GENE SETS
  ############################################################
  
  starvation_gene_sets <- list(
    
    WT_Starvation_Up = WT_starvation_up_cgd,
    
    WT_Starvation_Down = WT_starvation_down_cgd,
    
    KCS1_Starvation_Up = KCS1_starvation_up_cgd,
    
    KCS1_Starvation_Down = KCS1_starvation_down_cgd
  )
  
  ############################################################
  # CHECK GO ANNOTATION COVERAGE
  ############################################################
  
  for (gene_set_name in names(starvation_gene_sets)) {
    
    genes <- starvation_gene_sets[[gene_set_name]]
    
    mapped <- intersect(
      genes,
      unique(go_annotations$CGD_ID)
    )
    
    cat(
      gene_set_name,
      ":",
      length(genes),
      "genes;",
      length(mapped),
      "with GO annotations\n"
    )
  }
  
  
  ############################################################
  # GO ENRICHMENT: PHOSPHATE STARVATION RESPONSE
  ############################################################
  
  starvation_go_results <- list()
  
  for (gene_set_name in names(starvation_gene_sets)) {
    
    ##########################################################
    # Biological Process
    ##########################################################
    
    bp_name <- paste0(gene_set_name, "_BP")
    
    starvation_go_results[[bp_name]] <- run_GO(
      genes = starvation_gene_sets[[gene_set_name]],
      gene_set_name = gene_set_name,
      ontology = "P"
    )
    
    ##########################################################
    # Molecular Function
    ##########################################################
    
    mf_name <- paste0(gene_set_name, "_MF")
    
    starvation_go_results[[mf_name]] <- run_GO(
      genes = starvation_gene_sets[[gene_set_name]],
      gene_set_name = gene_set_name,
      ontology = "F"
    )
    
    ##########################################################
    # Cellular Component
    ##########################################################
    
    cc_name <- paste0(gene_set_name, "_CC")
    
    starvation_go_results[[cc_name]] <- run_GO(
      genes = starvation_gene_sets[[gene_set_name]],
      gene_set_name = gene_set_name,
      ontology = "C"
    )
  }
  
  
  ############################################################
  # PRINT SIGNIFICANT GO TERMS
  ############################################################
  
  for (result_name in names(starvation_go_results)) {
    
    result <- starvation_go_results[[result_name]]
    
    if (is.null(result)) {
      next
    }
    
    result_df <- as.data.frame(result)
    
    sig <- result_df %>%
      dplyr::filter(
        !is.na(p.adjust),
        p.adjust < 0.05
      ) %>%
      dplyr::arrange(p.adjust)
    
    cat("\n============================================================\n")
    cat(result_name, "\n")
    cat("============================================================\n")
    
    if (nrow(sig) == 0) {
      
      cat("No significant GO terms.\n")
      
    } else {
      
      print(
        sig %>%
          dplyr::select(
            ID,
            Description,
            GeneRatio,
            BgRatio,
            Count,
            pvalue,
            p.adjust
          )
      )
    }
  }
  
  ############################################################
  # EXCEL EXPORT: PHOSPHATE STARVATION GO ANALYSIS
  ############################################################
  
  library(openxlsx)
  library(dplyr)
  
  ############################################################
  # CREATE OUTPUT DIRECTORY
  ############################################################
  
  go_output_dir <- "05_GO_Analysis"
  
  if (!dir.exists(go_output_dir)) {
    dir.create(go_output_dir, recursive = TRUE)
  }
  
  ############################################################
  # EXCEL WORKBOOK
  ############################################################
  
  wb <- openxlsx::createWorkbook()
  
  ############################################################
  # SUMMARY TABLE
  ############################################################
  
  summary_df <- data.frame(
    Gene_Set = c(
      "WT_Starvation_Up",
      "WT_Starvation_Down",
      "KCS1_Starvation_Up",
      "KCS1_Starvation_Down"
    ),
    Significant_Genes = c(
      length(WT_starvation_up_cgd),
      length(WT_starvation_down_cgd),
      length(KCS1_starvation_up_cgd),
      length(KCS1_starvation_down_cgd)
    ),
    GO_Annotated_Genes = c(
      length(intersect(
        WT_starvation_up_cgd,
        unique(go_annotations$CGD_ID)
      )),
      length(intersect(
        WT_starvation_down_cgd,
        unique(go_annotations$CGD_ID)
      )),
      length(intersect(
        KCS1_starvation_up_cgd,
        unique(go_annotations$CGD_ID)
      )),
      length(intersect(
        KCS1_starvation_down_cgd,
        unique(go_annotations$CGD_ID)
      ))
    ),
    stringsAsFactors = FALSE
  )
  
  openxlsx::addWorksheet(
    wb,
    "Summary"
  )
  
  openxlsx::writeData(
    wb,
    "Summary",
    summary_df
  )
  
  ############################################################
  # FUNCTION TO PREPARE GO RESULTS
  ############################################################
  
  prepare_GO_table <- function(result) {
    
    if (is.null(result)) {
      return(
        data.frame(
          Message = "No GO enrichment result",
          stringsAsFactors = FALSE
        )
      )
    }
    
    df <- as.data.frame(result)
    
    if (nrow(df) == 0) {
      return(
        data.frame(
          Message = "No significant GO terms",
          stringsAsFactors = FALSE
        )
      )
    }
    
    df <- df %>%
      dplyr::mutate(
        GeneRatio = as.character(GeneRatio),
        BgRatio = as.character(BgRatio),
        pvalue = as.numeric(pvalue),
        p.adjust = as.numeric(p.adjust),
        qvalue = as.numeric(qvalue)
      ) %>%
      dplyr::arrange(p.adjust)
    
    return(df)
  }
  
  ############################################################
  # WRITE ALL GO RESULTS
  ############################################################
  
  for (result_name in names(starvation_go_results)) {
    
    result_df <- prepare_GO_table(
      starvation_go_results[[result_name]]
    )
    
    # Make Excel sheet name <= 31 characters
    sheet_name <- result_name
    
    if (nchar(sheet_name) > 31) {
      sheet_name <- substr(sheet_name, 1, 31)
    }
    
    openxlsx::addWorksheet(
      wb,
      sheet_name
    )
    
    openxlsx::writeData(
      wb,
      sheet_name,
      result_df,
      withFilter = TRUE
    )
    
    openxlsx::setColWidths(
      wb,
      sheet_name,
      cols = 1:ncol(result_df),
      widths = "auto"
    )
  }
  
  ############################################################
  # FORMAT SUMMARY
  ############################################################
  
  header_style <- openxlsx::createStyle(
    textDecoration = "bold",
    halign = "center",
    border = "Bottom"
  )
  
  openxlsx::addStyle(
    wb,
    "Summary",
    header_style,
    rows = 1,
    cols = 1:ncol(summary_df),
    gridExpand = TRUE
  )
  
  openxlsx::setColWidths(
    wb,
    "Summary",
    cols = 1:ncol(summary_df),
    widths = "auto"
  )
  
  ############################################################
  # SAVE WORKBOOK
  ############################################################
  
  excel_file <- file.path(
    go_output_dir,
    "Phosphate_Starvation_GO_Results.xlsx"
  )
  
  openxlsx::saveWorkbook(
    wb,
    excel_file,
    overwrite = TRUE
  )
  
  cat("\n============================================================\n")
  cat("GO EXCEL WORKBOOK SAVED\n")
  cat("============================================================\n")
  cat(excel_file, "\n")
  
  
  ############################################################
  # GO BAR PLOTS: PHOSPHATE STARVATION
  ############################################################
  
  library(ggplot2)
  library(dplyr)
  
  ############################################################
  # PLOT OUTPUT DIRECTORY
  ############################################################
  
  plot_dir <- file.path(
    go_output_dir,
    "GO_Plots"
  )
  
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
  }
  
  ############################################################
  # FUNCTION TO PREPARE PLOT DATA
  ############################################################
  
  prepare_GO_plot <- function(
    result_names,
    ontology_labels,
    top_n = 10
  ) {
    
    plot_list <- list()
    
    for (i in seq_along(result_names)) {
      
      result_name <- result_names[i]
      ontology_label <- ontology_labels[i]
      
      result <- starvation_go_results[[result_name]]
      
      if (is.null(result)) {
        next
      }
      
      df <- as.data.frame(result)
      
      if (nrow(df) == 0) {
        next
      }
      
      df <- df %>%
        dplyr::filter(
          !is.na(p.adjust),
          p.adjust < 0.05
        ) %>%
        dplyr::mutate(
          p.adjust = as.numeric(p.adjust),
          Ontology = ontology_label
        ) %>%
        dplyr::arrange(p.adjust) %>%
        dplyr::slice_head(n = top_n)
      
      if (nrow(df) > 0) {
        plot_list[[length(plot_list) + 1]] <- df
      }
    }
    
    if (length(plot_list) == 0) {
      return(NULL)
    }
    
    combined <- dplyr::bind_rows(plot_list)
    
    combined <- combined %>%
      dplyr::mutate(
        neg_log10_FDR = -log10(p.adjust),
        Description = factor(
          Description,
          levels = rev(unique(Description))
        )
      )
    
    return(combined)
  }
  
  ############################################################
  # FUNCTION TO CREATE A PLOT
  ############################################################
  
  make_GO_barplot <- function(
    result_names,
    ontology_labels,
    plot_title,
    output_file,
    top_n = 10
  ) {
    
    plot_df <- prepare_GO_plot(
      result_names = result_names,
      ontology_labels = ontology_labels,
      top_n = top_n
    )
    
    if (is.null(plot_df)) {
      
      cat(
        "\nNo significant GO terms for:",
        plot_title,
        "\n"
      )
      
      return(NULL)
    }
    
    p <- ggplot(
      plot_df,
      aes(
        x = neg_log10_FDR,
        y = Description
      )
    ) +
      geom_col() +
      facet_wrap(
        ~Ontology,
        scales = "free_y"
      ) +
      labs(
        title = plot_title,
        x = expression(-log[10]("FDR")),
        y = "GO term"
      ) +
      theme_bw() +
      theme(
        plot.title = element_text(
          face = "bold",
          size = 14
        ),
        axis.text.y = element_text(
          size = 9
        ),
        strip.text = element_text(
          face = "bold"
        )
      )
    
    ggsave(
      filename = output_file,
      plot = p,
      width = 12,
      height = 8,
      dpi = 300
    )
    
    return(p)
  }
  
  ############################################################
  # WT STARVATION UP
  ############################################################
  
  WT_up_plot <- make_GO_barplot(
    result_names = c(
      "WT_Starvation_Up_BP",
      "WT_Starvation_Up_MF",
      "WT_Starvation_Up_CC"
    ),
    ontology_labels = c(
      "Biological Process",
      "Molecular Function",
      "Cellular Component"
    ),
    plot_title = "WT response to phosphate starvation: upregulated genes",
    output_file = file.path(
      plot_dir,
      "WT_Starvation_Up_GO.png"
    ),
    top_n = 10
  )
  
  ############################################################
  # WT STARVATION DOWN
  ############################################################
  
  WT_down_plot <- make_GO_barplot(
    result_names = c(
      "WT_Starvation_Down_BP",
      "WT_Starvation_Down_MF",
      "WT_Starvation_Down_CC"
    ),
    ontology_labels = c(
      "Biological Process",
      "Molecular Function",
      "Cellular Component"
    ),
    plot_title = "WT response to phosphate starvation: downregulated genes",
    output_file = file.path(
      plot_dir,
      "WT_Starvation_Down_GO.png"
    ),
    top_n = 10
  )
  
  ############################################################
  # KCS1 STARVATION UP
  ############################################################
  
  KCS1_up_plot <- make_GO_barplot(
    result_names = c(
      "KCS1_Starvation_Up_BP",
      "KCS1_Starvation_Up_MF",
      "KCS1_Starvation_Up_CC"
    ),
    ontology_labels = c(
      "Biological Process",
      "Molecular Function",
      "Cellular Component"
    ),
    plot_title = "KCS1 response to phosphate starvation: upregulated genes",
    output_file = file.path(
      plot_dir,
      "KCS1_Starvation_Up_GO.png"
    ),
    top_n = 10
  )
  
  ############################################################
  # KCS1 STARVATION DOWN
  ############################################################
  
  KCS1_down_plot <- make_GO_barplot(
    result_names = c(
      "KCS1_Starvation_Down_BP",
      "KCS1_Starvation_Down_MF",
      "KCS1_Starvation_Down_CC"
    ),
    ontology_labels = c(
      "Biological Process",
      "Molecular Function",
      "Cellular Component"
    ),
    plot_title = "KCS1 response to phosphate starvation: downregulated genes",
    output_file = file.path(
      plot_dir,
      "KCS1_Starvation_Down_GO.png"
    ),
    top_n = 10
  )
  
  
  ############################################################
  # COMPARISON PLOT
  # WT VS KCS1 STARVATION RESPONSE
  ############################################################
  
  comparison_df <- bind_rows(
    
    prepare_GO_plot(
      c(
        "WT_Starvation_Up_BP",
        "WT_Starvation_Up_MF",
        "WT_Starvation_Up_CC"
      ),
      c(
        "Biological Process",
        "Molecular Function",
        "Cellular Component"
      ),
      top_n = 10
    ) %>%
      mutate(
        Strain = "WT",
        Direction = "Up"
      ),
    
    prepare_GO_plot(
      c(
        "KCS1_Starvation_Up_BP",
        "KCS1_Starvation_Up_MF",
        "KCS1_Starvation_Up_CC"
      ),
      c(
        "Biological Process",
        "Molecular Function",
        "Cellular Component"
      ),
      top_n = 10
    ) %>%
      mutate(
        Strain = "KCS1",
        Direction = "Up"
      ),
    
    prepare_GO_plot(
      c(
        "WT_Starvation_Down_BP",
        "WT_Starvation_Down_MF",
        "WT_Starvation_Down_CC"
      ),
      c(
        "Biological Process",
        "Molecular Function",
        "Cellular Component"
      ),
      top_n = 10
    ) %>%
      mutate(
        Strain = "WT",
        Direction = "Down"
      ),
    
    prepare_GO_plot(
      c(
        "KCS1_Starvation_Down_BP",
        "KCS1_Starvation_Down_MF",
        "KCS1_Starvation_Down_CC"
      ),
      c(
        "Biological Process",
        "Molecular Function",
        "Cellular Component"
      ),
      top_n = 10
    ) %>%
      mutate(
        Strain = "KCS1",
        Direction = "Down"
      )
  )
  
  comparison_df <- comparison_df %>%
    mutate(
      Strain_Direction = paste(
        Strain,
        Direction
      )
    )
  
  comparison_plot <- ggplot(
    comparison_df,
    aes(
      x = neg_log10_FDR,
      y = Description
    )
  ) +
    geom_col() +
    facet_grid(
      Strain_Direction ~ Ontology,
      scales = "free_y"
    ) +
    labs(
      title = "GO enrichment during phosphate starvation",
      x = expression(-log[10]("FDR")),
      y = "GO term"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 14
      ),
      axis.text.y = element_text(
        size = 8
      ),
      strip.text = element_text(
        face = "bold"
      )
    )
  
  ggsave(
    filename = file.path(
      plot_dir,
      "WT_vs_KCS1_Starvation_GO_Comparison.png"
    ),
    plot = comparison_plot,
    width = 14,
    height = 14,
    dpi = 300
  )
  
  print(comparison_plot)
