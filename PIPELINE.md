# Transcriptional Differences Between Inositol Pyrophosphate Pathway Mutants in *Candida albicans* Under Phosphate Starved Conditions
## Visualizations
Done in RStudio.

### WT vs KCS1 analysis

```bash
# Candida albicans RNA-seq analysis
#
# 12 samples: WT and KCS1, each grown in +Pi and -Pi, with 3 biological replicates
# per condition. Samples are ordered as WT+Pi (1-3), WT-Pi (4-6),
# KCS1+Pi (7-9), and KCS1-Pi (10-12).
#
# The script filters tRNA/rRNA and low-count genes, runs DESeq2 with a
# Genotype * Pi design, generates QC/DE plots, and tests the genotype-by-Pi
# interaction to compare phosphate responses between WT and KCS1.

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

# 1. Load Required Packages

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

# 2. Output Directories

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

# 3. Load Cgd Gene Annotation

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

# Cgd annotation structure
#
# Column 1 = systematic ID
# Column 2 = gene name
# Column 4 = feature type

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

mapping_table$gene_symbol[
is.na(mapping_table$gene_symbol) |
    mapping_table$gene_symbol == "" |
    mapping_table$gene_symbol ==
    mapping_table$systematic_id
] <- NA

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

#
# Prefer entries that actually have gene names.

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

# 4. Load Featurecounts Matrix

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

# Robust featurecounts reader
#
# Avoids:
# "more columns than column names"

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

data_lines <- lines[-1]

split_lines <- strsplit(
    data_lines,
    "\t",
    fixed = TRUE
  )

# Remove trailing empty fields

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

# 5. Identify And Reorder Sample Columns

print("======================================================")
print("IDENTIFYING RNA-seq SAMPLE COLUMNS")
print("======================================================")

#
# 1 = Geneid
# 2 = Chr
# 3 = Start
# 4 = End
# 5 = Strand
# 6 = Length
# 7-18 = samples

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

# Extract n44vrl sample numbers

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

# Reorder by n44vrl sample number

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

# 6. Extract Gene Ids And Counts

counts_clean <- counts_raw[
  ,
c(1, sample_columns),
drop = FALSE
]

colnames(counts_clean)[1] <- "Geneid"

for (i in 2:ncol(counts_clean)) {

counts_clean[[i]] <- as.numeric(
    counts_clean[[i]]
  )

}

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

# Collapse duplicate gene IDs

print(
  "Collapsing duplicate gene IDs..."
)

counts_fixed <- aggregate(
  . ~ Geneid,
data = counts_clean,
FUN = sum
)

# 7. Remove Trna And Rrna Genes

print("======================================================")
print("REMOVING tRNA / rRNA FEATURES")
print("======================================================")

stripped_keys <- gsub(
  "^CAALFM_",
  "",
counts_fixed$Geneid,
ignore.case = TRUE
)

normalized_matrix_keys <- normalize_id(
stripped_keys
)

matched_features <- mapping_table$feature_type[
match(
    normalized_matrix_keys,
    mapping_table$id_norm
  )
]

matched_symbols <- mapping_table$gene_symbol[
match(
    normalized_matrix_keys,
    mapping_table$id_norm
  )
]

# IDENTIFY tRNA / rRNA

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

rownames(counts_fixed) <- counts_fixed$Geneid

counts <- counts_fixed[
  ,
  -1,
drop = FALSE
]

# 8. Create Master Gene Annotation

print("======================================================")
print("ANNOTATING GENES")
print("======================================================")

original_gene_ids <- rownames(
counts
)

stripped_ids <- gsub(
  "^CAALFM_",
  "",
original_gene_ids,
ignore.case = TRUE
)

normalized_ids <- normalize_id(
stripped_ids
)

# Match CGD gene names

matched_gene_names <- mapping_table$gene_symbol[
match(
    normalized_ids,
    mapping_table$id_norm
  )
]

#
# If a CGD gene name exists, use it.
# Otherwise retain the systematic ID.

display_gene_names <- ifelse(

  !is.na(matched_gene_names) &
    matched_gene_names != "",

matched_gene_names,

stripped_ids

)

gene_annotation <- data.frame(

Systematic_ID = original_gene_ids,

Gene_Name = display_gene_names,

stringsAsFactors = FALSE

)

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

write.csv(

gene_annotation,

file.path(
    output_dir,
    "02_DE_results",
    "Gene_ID_to_Gene_Name_Annotation.csv"
  ),

row.names = FALSE

)

# Helper function for gene name lookup

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

# 9. Sample Metadata

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

# Human-readable condition labels
#
# THESE ARE THE LABELS THAT WILL APPEAR ON FIGURES.

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

colnames(counts) <- short_sample_names

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

# Internal condition variable

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

# Condition label factor

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

# Verify replicates

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

# 10. Filter Low-Count Genes

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

# 11. Deseq2

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

# 12. Vst

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

# 13. Library Size Qc

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

# 14. Pca

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

# 15. Sample Correlation

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

# 16. Deseq2 Coefficients

print("======================================================")
print("DESEQ2 COEFFICIENTS")
print("======================================================")

results_names <- resultsNames(
dds
)

print(
results_names
)

# Genotype coefficient

genotype_coef <- results_names[
results_names == "Genotype_KCS1_vs_WT"
]

# Pi coefficient

pi_coef <- results_names[
results_names == "Pi_minus_vs_plus"
]

# Interaction coefficient

interaction_coef <- results_names[
grepl(
    "Genotype.*Pi|Pi.*Genotype",
    results_names
  )
]

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

# 17. All Six Pairwise Comparisons
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

# 17.1 WT+Pi vs WT-Pi
#
# Result is WT-Pi - WT+Pi

res_WT_minus_vs_WT_plus <- results(

    dds,

    name = pi_coef

  )

# 17.2 WT+Pi vs KCS1+Pi
#
# Result is KCS1+Pi - WT+Pi

res_KCS1_plus_vs_WT_plus <- results(

    dds,

    name = genotype_coef

  )

# 17.3 WT+Pi vs KCS1-Pi
#
# Result is KCS1-Pi - WT+Pi
#
# = Genotype + Pi + Interaction

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

# 17.4 WT-Pi vs KCS1+Pi
#
# Result is KCS1+Pi - WT-Pi
#
# = Genotype - Pi

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

# 17.5 WT-Pi vs KCS1-Pi
#
# Result is KCS1-Pi - WT-Pi
#
# = Genotype + Interaction

res_KCS1_minus_vs_WT_minus <- results(

    dds,

    contrast = list(
      c(
        genotype_coef,
        interaction_coef
      )
    )

  )

# 17.6 KCS1+Pi vs KCS1-Pi
#
# Result is KCS1-Pi - KCS1+Pi
#
# = Pi + Interaction

res_KCS1_minus_vs_KCS1_plus <- results(

    dds,

    contrast = list(
      c(
        pi_coef,
        interaction_coef
      )
    )

  )

# 17.7 Formal Genotype X Pi Interaction

res_interaction <- results(

    dds,

    name = interaction_coef

  )

# 18. Annotate All De Results

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

    # Put gene information first

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

# Create named result list

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

# 19. Save All De Results

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

# 20. Significant Gene Lists

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

# 21. Volcano Plot Function

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

# 22. Create Volcano Plots For All Six Comparisons

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

# Interaction volcano

make_volcano(

    res_interaction_df,

    "Genotype × Pi interaction",

    file.path(
      output_dir,
      "03_Volcano",
      "Volcano_Genotype_x_Pi_interaction.png"
    )

  )

# 23. Ma Plot Function

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

# 24. Ma Plots For All Comparisons

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

# 25. Top Variable Gene Heatmap

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

# 26. De Heatmap Function

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

# 27. Heatmaps For All Six Pairwise Comparisons

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

# Interaction heatmap

make_de_heatmap(

    res_interaction_df,

    "Genotype × Pi interaction",

    file.path(
      output_dir,
      "05_Heatmaps",
      "Heatmap_Genotype_x_Pi_interaction.png"
    )

  )

# 28. Condition-Level Vst Heatmap
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

print("======================================================")
print("CREATING CONDITION-LEVEL HEATMAP")
print("======================================================")

# Create condition mean matrix

condition_means <- matrix(

    NA_real_,

    nrow = nrow(vst_matrix),

    ncol = 4

  )

rownames(condition_means) <- rownames(
    vst_matrix
  )

colnames(condition_means) <- c(

    "WT+Pi",
    "WT-Pi",
    "KCS1+Pi",
    "KCS1-Pi"

  )

# Identify sample indices for each condition

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

# Verify that each condition has 3 replicates

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

# Calculate mean VST expression
#
# Replicates are averaged ONLY here for visualization.
# DESeq2 statistics continue to use all biological replicates.

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

# Remove genes with non-finite values

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

# Calculate variance across the four conditions

condition_variances <- matrixStats::rowVars(

    condition_means

  )

# Remove genes with zero variance
#
# These genes would become NaN during row scaling.

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

# Select top 50 variable genes
#
# IMPORTANT:
# Use numeric row indices rather than names().
# matrixStats::rowVars() does not necessarily preserve
# row names as names of the returned variance vector.

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

condition_heatmap <- condition_means[

    top_condition_indices,

    ,

    drop = FALSE

  ]

print(
    paste(
      "Number of genes selected for condition heatmap:",
      nrow(condition_heatmap)
    )
  )

# MAP SYSTEMATIC IDs TO CGD GENE NAMES

condition_gene_names <- get_gene_names(

    rownames(
      condition_heatmap
    )

  )

rownames(
    condition_heatmap
  ) <- make.unique(

    condition_gene_names

  )

# Row-scale expression
#
# Each gene is standardized across the four conditions.
#
# This shows relative expression patterns rather than
# absolute VST expression.

condition_heatmap_scaled <- t(

    scale(

      t(
        condition_heatmap
      )

    )

  )

# Remove any rows that became non-finite after scaling

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

# Save condition mean matrix

write.csv(

    condition_means,

    file.path(
      output_dir,
      "05_Heatmaps",
      "Condition_Mean_VST_Expression.csv"
    ),

    row.names = TRUE

  )

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

# 29. Phosphate Response Analysis
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

print("======================================================")
print("PHOSPHATE RESPONSE ANALYSIS")
print("======================================================")

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

write.csv(

    response_table,

    file.path(
      output_dir,
      "09_Response_analysis",
      "Phosphate_Response_Comparison.csv"
    ),

    row.names = FALSE

  )

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

# 30. Phosphate Response Heatmap

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

# 31. Phosphate Response Scatterplot

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

# 32. Phosphate Response Plot
#
# Only significant interaction genes are shown.
#

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

    # CAP HEIGHT TO PREVENT ggsave ERROR

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

# 33. Summary Table

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

# 34. Comparison Summary Print

print("======================================================")
print("SIGNIFICANT GENE COUNTS")
print("======================================================")

print(
    summary_table
  )

# 35. Save Vst Object

saveRDS(

    vsd,

    file.path(
      output_dir,
      "01_QC",
      "VST_object.rds"
    )

  )

# 36. Final Summary

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
```

### WT vs VIP analysis
```bash
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
 # Candida albicans: VIP1 phosphate-response trial
# Three biological replicates per condition, twelve samples total.
# After sorting N44VRL sample numbers:
#   1–3 WT +Pi; 4–6 WT -Pi; 7–9 VIP1 +Pi; 10–12 VIP1 -Pi.
# Run from the VIP1 trial folder containing candida_counts_matrix.txt
# and cgd_features.tab. This script writes to RNAseq_VIP1_DE_Results.
# GO analysis also requires gene_association.CGD, clusterProfiler, and GO.db.

# 1. Load required packages ---------------

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

# 2. Output directories ---------------

output_dir <- "RNAseq_VIP1_DE_Results"
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
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
for (d in dirs) {
  dir.create(file.path(output_dir, d), showWarnings = FALSE, recursive = TRUE)
}

# 3. Load CGD gene annotation ---------------

message("Loading CGD gene annotations...")
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
# Column 1 = systematic ID
# Column 2 = gene name
# Column 4 = feature type
mapping_table <- data.frame(
  systematic_id = trimws(as.character(cgd_raw[, 1])),
  gene_symbol = trimws(as.character(cgd_raw[, 2])),
  feature_type = trimws(as.character(cgd_raw[, 4])),
  stringsAsFactors = FALSE
)
mapping_table$gene_symbol[
  is.na(mapping_table$gene_symbol) |
  mapping_table$gene_symbol == "" |
  mapping_table$gene_symbol ==
  mapping_table$systematic_id
] <- NA
normalize_id <- function(x) {
  x <- trimws(as.character(x))
  x <- toupper(gsub("[_\\-]", "", x))
  return(x)
}
mapping_table$id_norm <- normalize_id(mapping_table$systematic_id)
# Prefer entries that actually have gene names.
mapping_table <- mapping_table[
  order(is.na(mapping_table$gene_symbol)),
  ,
  drop = FALSE
]
mapping_table <- mapping_table[
  !duplicated(mapping_table$id_norm),
  ,
  drop = FALSE
]
print(paste("Total CGD features:", nrow(mapping_table)))
print(paste("CGD features with gene names:", sum(!is.na(mapping_table$gene_symbol))))

# 4. Load featureCounts matrix ---------------

message("Loading the featureCounts matrix...")
if (!file.exists("candida_counts_matrix.txt")) {
  stop(
    paste0(
      "\nERROR: candida_counts_matrix.txt was not found.\n",
      "Place it in the same working directory as this script.\n"
    )
  )
}
# Avoids:
# "more columns than column names"
read_featurecounts <- function(file) {
  lines <- readLines(file, warn = FALSE)
  lines <- lines[
    nzchar(trimws(lines))
  ]
  lines <- lines[
    !grepl(
      "^#",
      trimws(lines)
    )
  ]
  if (length(lines) < 2) {
    stop("ERROR: FeatureCounts file does not contain enough data.")
  }
  header <- strsplit(lines[1], "\t", fixed = TRUE)[[1]]
  header[1] <- sub("^\ufeff", "", header[1])
  if (length(header) != 18) {
    stop(
      paste0(
        "\nERROR: Expected 18 columns in featureCounts header.\n",
        "Found ",
        length(header),
        " columns.\n\n",
        "Header detected:\n",
        paste(header, collapse = "\n")
      )
    )
  }
  data_lines <- lines[-1]
  split_lines <- strsplit(data_lines, "\t", fixed = TRUE)
  split_lines <- lapply(
    split_lines,
    function(x) {
      while (length(x) > 0 && trimws(x[length(x)]) == "") {
        x <- x[
          -length(x)
        ]
      }
      x
    }
  )
  field_counts <- lengths(split_lines)
  bad_rows <- which(field_counts != 18)
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
  matrix_data <- do.call(rbind, split_lines)
  counts_df <- as.data.frame(matrix_data, stringsAsFactors = FALSE, check.names = FALSE)
  colnames(counts_df) <- header
  return(counts_df)
}
counts_raw <- read_featurecounts("candida_counts_matrix.txt")
print(paste("Total columns detected:", ncol(counts_raw)))
print("Column names detected:")
print(colnames(counts_raw))

# 5. Identify and reorder sample columns ---------------

print("IDENTIFYING RNA-seq SAMPLE COLUMNS")
# 1 = Geneid
# 2 = Chr
# 3 = Start
# 4 = End
# 5 = Strand
# 6 = Length
# 7-18 = samples
sample_columns_original <- 7:18
original_sample_names <- colnames(counts_raw)[sample_columns_original]
print("Original sample order in featureCounts file:")
print(original_sample_names)
sample_numbers <- suppressWarnings(
  as.numeric(sub(".*N44VRL_([0-9]+)_.*", "\\1", original_sample_names))
)
print("Detected N44VRL sample numbers:")
print(sample_numbers)
if (length(sample_numbers) != 12 || any(is.na(sample_numbers)) || !setequal(sample_numbers, 1:12)) {
  stop(
    paste0(
      "\nERROR: Could not correctly identify N44VRL ",
      "samples 1-12.\n\n",
      "Detected sample names:\n",
      paste(original_sample_names, collapse = "\n"),
      "\n\nDetected sample numbers:\n",
      paste(sample_numbers, collapse = ", ")
    )
  )
}
sample_columns <- sample_columns_original[
  order(sample_numbers)
]
original_sample_names <- colnames(counts_raw)[sample_columns]
print("CORRECTED SAMPLE ORDER")
for (i in seq_along(original_sample_names)) {
  print(paste(i, ":", original_sample_names[i]))
}
expected_genotype <- rep(c("WT", "VIP1"), each = 6)
expected_pi <- rep(rep(c("plus", "minus"), each = 3), times = 2)
detected_genotype <- ifelse(
  grepl("WT", original_sample_names, ignore.case = TRUE),
  "WT",
  ifelse(grepl("VIP1", original_sample_names, ignore.case = TRUE), "VIP1", NA)
)
detected_pi <- ifelse(
  grepl("\\+Pi", original_sample_names, ignore.case = TRUE),
  "plus",
  ifelse(grepl("-Pi", original_sample_names, ignore.case = TRUE), "minus", NA)
)
if (any(is.na(detected_genotype)) || any(is.na(detected_pi))) {
  print(original_sample_names)
  print(detected_genotype)
  print(detected_pi)
  stop("\nERROR: Could not determine genotype or Pi condition from sample names.")
}
if (!identical(detected_genotype, expected_genotype) || !identical(detected_pi, expected_pi)) {
  print("Expected genotype:")
  print(expected_genotype)
  print("Detected genotype:")
  print(detected_genotype)
  print("Expected Pi:")
  print(expected_pi)
  print("Detected Pi:")
  print(detected_pi)
  stop("\nERROR: Sample order does not match the expected experimental design.")
}
print("SUCCESS: Sample order is correct.")

# 6. Extract gene ids and counts ---------------

counts_clean <- counts_raw[
  ,
  c(1, sample_columns),
  drop = FALSE
]
colnames(counts_clean)[1] <- "Geneid"
for (i in 2:ncol(counts_clean)) {
  counts_clean[[i]] <- as.numeric(counts_clean[[i]])
}
if (any(is.na(counts_clean[, -1]))) {
  stop("\nERROR: NA values detected in count matrix.")
}
print("Collapsing duplicate gene IDs...")
counts_fixed <- aggregate(. ~ Geneid, data = counts_clean, FUN = sum)

# 7. Remove tRNA and rRNA genes ---------------

print("REMOVING tRNA / rRNA FEATURES")
stripped_keys <- gsub("^CAALFM_", "", counts_fixed$Geneid, ignore.case = TRUE)
normalized_matrix_keys <- normalize_id(stripped_keys)
matched_features <- mapping_table$feature_type[
  match(normalized_matrix_keys, mapping_table$id_norm)
]
matched_symbols <- mapping_table$gene_symbol[
  match(normalized_matrix_keys, mapping_table$id_norm)
]
# IDENTIFY tRNA / rRNA
is_trna_rrna <- (
  grepl("tRNA|rRNA", counts_fixed$Geneid, ignore.case = TRUE) |
  grepl("tRNA|rRNA", matched_features, ignore.case = TRUE) |
  grepl("^tRNA|^rRNA", matched_symbols, ignore.case = TRUE)
)
is_trna_rrna[
  is.na(is_trna_rrna)
] <- FALSE
print(paste("Features before tRNA/rRNA filtering:", nrow(counts_fixed)))
print(paste("tRNA/rRNA features removed:", sum(is_trna_rrna)))
counts_fixed <- counts_fixed[
  !is_trna_rrna,
  ,
  drop = FALSE
]
print(paste("Features remaining:", nrow(counts_fixed)))
rownames(counts_fixed) <- counts_fixed$Geneid
counts <- counts_fixed[
  ,
  -1,
  drop = FALSE
]

# 8. Create master gene annotation ---------------

print("ANNOTATING GENES")
original_gene_ids <- rownames(counts)
stripped_ids <- gsub("^CAALFM_", "", original_gene_ids, ignore.case = TRUE)
normalized_ids <- normalize_id(stripped_ids)
matched_gene_names <- mapping_table$gene_symbol[
  match(normalized_ids, mapping_table$id_norm)
]
# If a CGD gene name exists, use it.
# Otherwise retain the systematic ID.
display_gene_names <- ifelse(
  !is.na(matched_gene_names) &
  matched_gene_names != "",
  matched_gene_names,
  stripped_ids
)
gene_annotation <- data.frame(
  Systematic_ID = original_gene_ids,
  Gene_Name = display_gene_names,
  stringsAsFactors = FALSE
)
mapped_count <- sum(!is.na(matched_gene_names) & matched_gene_names != "")
unmapped_count <- sum(is.na(matched_gene_names) | matched_gene_names == "")
print(paste("Genes mapped to CGD gene names:", mapped_count))
print(paste("Genes without a CGD gene name:", unmapped_count))
write.csv(
  gene_annotation,
  file.path(output_dir, "02_DE_results", "Gene_ID_to_Gene_Name_Annotation.csv"),
  row.names = FALSE
)
get_gene_names <- function(ids) {
  names_out <- gene_annotation$Gene_Name[
    match(ids, gene_annotation$Systematic_ID)
  ]
  missing <- (is.na(names_out) | names_out == "")
  names_out[
    missing
  ] <- ids[
    missing
  ]
  return(names_out)
}

# 9. Sample metadata ---------------

short_sample_names <- c(
  "WT_Pi_1",
  "WT_Pi_2",
  "WT_Pi_3",
  "WT_minusPi_1",
  "WT_minusPi_2",
  "WT_minusPi_3",
  "VIP1_Pi_1",
  "VIP1_Pi_2",
  "VIP1_Pi_3",
  "VIP1_minusPi_1",
  "VIP1_minusPi_2",
  "VIP1_minusPi_3"
)
condition_labels <- c(
  "WT+Pi",
  "WT+Pi",
  "WT+Pi",
  "WT-Pi",
  "WT-Pi",
  "WT-Pi",
  "VIP1+Pi",
  "VIP1+Pi",
  "VIP1+Pi",
  "VIP1-Pi",
  "VIP1-Pi",
  "VIP1-Pi"
)
replicate_vector <- c(1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3)
colnames(counts) <- short_sample_names
metadata <- data.frame(
  Original_Sample_ID = original_sample_names,
  Sample_Name = short_sample_names,
  Display_Name = condition_labels,
  Genotype = factor(detected_genotype, levels = c("WT", "VIP1")),
  Pi = factor(detected_pi, levels = c("plus", "minus")),
  Replicate = replicate_vector,
  stringsAsFactors = FALSE
)
rownames(metadata) <- short_sample_names
metadata$Condition <- factor(
  paste(metadata$Genotype, metadata$Pi, sep = "_"),
  levels = c("WT_plus", "WT_minus", "VIP1_plus", "VIP1_minus")
)
metadata$Condition_Label <- factor(
  condition_labels,
  levels = c("WT+Pi", "WT-Pi", "VIP1+Pi", "VIP1-Pi")
)
print("SAMPLE METADATA")
print(metadata)
replicate_table <- table(metadata$Condition)
print(replicate_table)
if (!all(replicate_table == 3)) {
  stop("ERROR: Every condition must contain exactly 3 replicates.")
}
write.csv(metadata, file.path(output_dir, "sample_metadata.csv"), row.names = TRUE)

# 10. Filter low-count genes ---------------

print("FILTERING LOW-COUNT GENES")
keep_genes <- rowSums(counts) >= 10
print(paste("Genes before filtering:", nrow(counts)))
print(paste("Genes passing count filter:", sum(keep_genes)))
print(paste("Genes removed:", sum(!keep_genes)))
counts_filtered <- counts[
  keep_genes,
  ,
  drop = FALSE
]

# 11. DESeq2 ---------------

print("CREATING DESEQ2 OBJECT")
dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts_filtered)),
  colData = metadata,
  design = ~ Genotype * Pi
)
dds <- DESeq(dds)
saveRDS(dds, file.path(output_dir, "02_DE_results", "DESeq2_object.rds"))

# 12. VST ---------------

print("VST TRANSFORMATION")
vsd <- vst(dds, blind = FALSE)
vst_matrix <- assay(vsd)
write.csv(
  vst_matrix,
  file.path(output_dir, "01_QC", "VST_expression_matrix.csv"),
  row.names = TRUE
)

# 13. Library size QC ---------------

library_sizes <- colSums(counts(dds))
library_size_df <- data.frame(
  Sample = names(library_sizes),
  Condition = metadata[
    names(library_sizes),
    "Display_Name"
  ],
  Library_Size = as.numeric(library_sizes),
  stringsAsFactors = FALSE
)
write.csv(library_size_df, file.path(output_dir, "01_QC", "Library_sizes.csv"), row.names = FALSE)
p_library <- ggplot(library_size_df, aes(x = Sample, y = Library_Size)) +
geom_col() +
theme_bw() +
theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
labs(title = "RNA-seq library sizes", x = "Sample", y = "Total assigned reads")
ggsave(
  file.path(output_dir, "01_QC", "Library_sizes.png"),
  p_library,
  width = 10,
  height = 6,
  dpi = 300
)

# 14. PCA ---------------

print("CREATING PCA")
pca_data <- plotPCA(vsd, intgroup = c("Genotype", "Pi"), returnData = TRUE)
percent_variance <- round(100 * attr(pca_data, "percentVar"))
p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, label = name)) +
geom_point(size = 4) +
geom_text_repel(size = 3) +
theme_bw() +
labs(
  title = "PCA of RNA-seq samples",
  x = paste0("PC1: ", percent_variance[1], "% variance"),
  y = paste0("PC2: ", percent_variance[2], "% variance")
)
ggsave(
  file.path(output_dir, "06_PCA", "PCA_PC1_PC2.png"),
  p_pca,
  width = 9,
  height = 7,
  dpi = 300
)

# 15. Sample correlation ---------------

cor_matrix <- cor(vst_matrix, method = "pearson")
write.csv(
  cor_matrix,
  file.path(output_dir, "07_Sample_correlations", "Sample_Pearson_correlations.csv")
)
png(
  file.path(output_dir, "07_Sample_correlations", "Sample_Pearson_correlations.png"),
  width = 2200,
  height = 2000,
  res = 300
)
cor_annotation <- data.frame(Condition = metadata[ colnames(vst_matrix), "Display_Name" ])
rownames(cor_annotation) <- colnames(vst_matrix)
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

# 16. DESeq2 coefficients ---------------

print("DESEQ2 COEFFICIENTS")
results_names <- resultsNames(dds)
print(results_names)
genotype_coef <- results_names[
  results_names == "Genotype_VIP1_vs_WT"
]
pi_coef <- results_names[
  results_names == "Pi_minus_vs_plus"
]
interaction_coef <- results_names[
  grepl("Genotype.*Pi|Pi.*Genotype", results_names)
]
if (length(genotype_coef) != 1) {
  stop(
    paste0(
      "Could not identify genotype coefficient.\n",
      "Available coefficients:\n",
      paste(results_names, collapse = "\n")
    )
  )
}
if (length(pi_coef) != 1) {
  stop(
    paste0(
      "Could not identify Pi coefficient.\n",
      "Available coefficients:\n",
      paste(results_names, collapse = "\n")
    )
  )
}
if (length(interaction_coef) != 1) {
  stop(
    paste0(
      "Could not identify interaction coefficient.\n",
      "Available coefficients:\n",
      paste(results_names, collapse = "\n")
    )
  )
}
print(paste("Genotype coefficient:", genotype_coef))
print(paste("Pi coefficient:", pi_coef))
print(paste("Interaction coefficient:", interaction_coef))

# 17. All six pairwise comparisons ---------------

# WT +Pi is the reference. VIP1 -Pi combines the genotype, Pi,
# and interaction coefficients.
# 17.1. WT+Pi vs WT-Pi
# Result is WT-Pi - WT+Pi
res_WT_minus_vs_WT_plus <- results(dds, name = pi_coef)
# 17.2. WT+Pi vs VIP1+Pi
# Result is VIP1+Pi - WT+Pi
res_VIP1_plus_vs_WT_plus <- results(dds, name = genotype_coef)
# 17.3. WT+Pi vs VIP1-Pi
# Result is VIP1-Pi - WT+Pi
# = Genotype + Pi + Interaction
res_VIP1_minus_vs_WT_plus <- results(
  dds,
  contrast = list(c(genotype_coef, pi_coef, interaction_coef))
)
# 17.4. WT-Pi vs VIP1+Pi
# Result is VIP1+Pi - WT-Pi
# = Genotype - Pi
res_VIP1_plus_vs_WT_minus <- results(dds, contrast = list(c(genotype_coef), c(pi_coef)))
# 17.5. WT-Pi vs VIP1-Pi
# Result is VIP1-Pi - WT-Pi
# = Genotype + Interaction
res_VIP1_minus_vs_WT_minus <- results(dds, contrast = list(c(genotype_coef, interaction_coef)))
# 17.6. VIP1+Pi vs VIP1-Pi
# Result is VIP1-Pi - VIP1+Pi
# = Pi + Interaction
res_VIP1_minus_vs_VIP1_plus <- results(dds, contrast = list(c(pi_coef, interaction_coef)))
# 17.7. Formal genotype x Pi interaction
res_interaction <- results(dds, name = interaction_coef)

# 18. Annotate all DE results ---------------

annotate_results <- function(result, comparison_label) {
  df <- as.data.frame(result)
  df$Systematic_ID <- rownames(df)
  df$Gene_Name <- get_gene_names(df$Systematic_ID)
  df$Comparison <- comparison_label
  # Put gene information first
  df <- df[
    ,
    c(
      "Systematic_ID",
      "Gene_Name",
      "Comparison",
      setdiff(colnames(df), c("Systematic_ID", "Gene_Name", "Comparison"))
    ),
    drop = FALSE
  ]
  return(df)
}
res_WT_minus_vs_WT_plus_df <- annotate_results(res_WT_minus_vs_WT_plus, "WT-Pi vs WT+Pi")
res_VIP1_plus_vs_WT_plus_df <- annotate_results(res_VIP1_plus_vs_WT_plus, "VIP1+Pi vs WT+Pi")
res_VIP1_minus_vs_WT_plus_df <- annotate_results(res_VIP1_minus_vs_WT_plus, "VIP1-Pi vs WT+Pi")
res_VIP1_plus_vs_WT_minus_df <- annotate_results(res_VIP1_plus_vs_WT_minus, "VIP1+Pi vs WT-Pi")
res_VIP1_minus_vs_WT_minus_df <- annotate_results(res_VIP1_minus_vs_WT_minus, "VIP1-Pi vs WT-Pi")
res_VIP1_minus_vs_VIP1_plus_df <- annotate_results(res_VIP1_minus_vs_VIP1_plus, "VIP1-Pi vs VIP1+Pi")
res_interaction_df <- annotate_results(res_interaction, "Genotype × Pi interaction")
all_results <- list(
  "WT-Pi_vs_WT+Pi" =
  res_WT_minus_vs_WT_plus_df,
  "VIP1+Pi_vs_WT+Pi" =
  res_VIP1_plus_vs_WT_plus_df,
  "VIP1-Pi_vs_WT+Pi" =
  res_VIP1_minus_vs_WT_plus_df,
  "VIP1+Pi_vs_WT-Pi" =
  res_VIP1_plus_vs_WT_minus_df,
  "VIP1-Pi_vs_WT-Pi" =
  res_VIP1_minus_vs_WT_minus_df,
  "VIP1-Pi_vs_VIP1+Pi" =
  res_VIP1_minus_vs_VIP1_plus_df,
  "Genotype_x_Pi_interaction" =
  res_interaction_df
)

# 19. Save all DE results ---------------

for (comparison_name in names(all_results)) {
  write.csv(
    all_results[[comparison_name]],
    file.path(output_dir, "02_DE_results", paste0("DE_", comparison_name, ".csv")),
    row.names = FALSE
  )
}

# 20. Significant gene lists ---------------

get_significant_genes <- function(result_df, padj_cutoff = 0.05) {
  result_df %>%
  filter(!is.na(padj), padj < padj_cutoff) %>%
  arrange(padj)
}
for (comparison_name in names(all_results)) {
  sig_genes <- get_significant_genes(all_results[[comparison_name]])
  write.csv(
    sig_genes,
    file.path(
      output_dir,
      "08_Significant_gene_lists",
      paste0("Significant_", comparison_name, ".csv")
    ),
    row.names = FALSE
  )
}

# 21. Volcano plot function ---------------

make_volcano <- function(
  result_df,
  plot_title,
  output_file,
  fc_cutoff = 1,
  padj_cutoff = 0.05,
  y_max = 50
) {
  # Check that required columns exist
  required_cols <- c("log2FoldChange", "padj")
  missing_cols <- setdiff(required_cols, colnames(result_df))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }
  # Keep genes with valid fold changes and adjusted p-values
  plot_df <- result_df %>%
  filter(!is.na(log2FoldChange), !is.na(padj), padj > 0)
  if (nrow(plot_df) == 0) {
    warning(paste("No valid genes available for:", plot_title))
    return(NULL)
  }
  # Calculate -log10 adjusted p-value
  plot_df$neg_log10_padj <- -log10(plot_df$padj)
  # Assign significance categories
  plot_df$Significance <- "Not significant"
  plot_df$Significance[
    plot_df$padj < padj_cutoff &
    plot_df$log2FoldChange >= fc_cutoff
  ] <- "Upregulated"
  plot_df$Significance[
    plot_df$padj < padj_cutoff &
    plot_df$log2FoldChange <= -fc_cutoff
  ] <- "Downregulated"
  # Identify genes above the displayed y-axis limit
  plot_df$Above_Y_Max <- (plot_df$neg_log10_padj > y_max)
  n_above <- sum(plot_df$Above_Y_Max, na.rm = TRUE)
  # Create display value
  # The TRUE -log10(padj) is retained in
  # plot_df$neg_log10_padj.
  # Only the value used for plotting is capped.
  plot_df$neg_log10_padj_plot <- pmin(plot_df$neg_log10_padj, y_max)
  # Select genes to label
  plot_df$Label <- NA_character_
  if ("Systematic_ID" %in% colnames(plot_df) && "Gene_Name" %in% colnames(plot_df)) {
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
      label_indices <- match(label_genes$Systematic_ID, plot_df$Systematic_ID)
      plot_df$Label[label_indices] <-
      label_genes$Gene_Name
    }
  }
  # Create volcano plot
  p <- ggplot(plot_df, aes(x = log2FoldChange, y = neg_log10_padj_plot, shape = Significance)) +
  # All genes
  geom_point(alpha = 0.6, size = 1.5) +
  # Fold-change thresholds
  geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed") +
  # Adjusted p-value threshold
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed") +
  # Mark genes above y-axis display limit
  # Triangle indicates that the true value is > y_max.
  geom_point(
    data = plot_df %>%
    filter(Above_Y_Max),
    aes(x = log2FoldChange, y = y_max),
    shape = 24,
    size = 2.5,
    inherit.aes = FALSE
  ) +
  # Gene labels
  geom_text_repel(aes(label = Label), na.rm = TRUE, size = 3, max.overlaps = 20) +
  # Display y-axis from 0 to 50
  # This does NOT change the actual statistical results.
  coord_cartesian(ylim = c(0, y_max), clip = "off") +
  # Theme
  theme_bw() +
  labs(
    title = plot_title,
    x = "log2 fold change",
    y = expression(-log[10]("adjusted p-value")),
    shape = "Significance"
  )
  # Add annotation for genes above y-axis limit
  if (n_above > 0) {
    p <- p +
    annotate(
      "text",
      x = Inf,
      y = y_max,
      label = paste0("▲ ", n_above, " genes > ", y_max),
      hjust = 1.05,
      vjust = -0.5,
      size = 3.5
    )
  }
  # Save plot
  ggsave(filename = output_file, plot = p, width = 8, height = 6, dpi = 300)
  # Return plot
  return(p)
}

# 22. Create volcano plots for all six pairwise comparisons ---------------

print("CREATING VOLCANO PLOTS")
volcano_titles <- c(
  "WT-Pi vs WT+Pi",
  "VIP1+Pi vs WT+Pi",
  "VIP1-Pi vs WT+Pi",
  "VIP1+Pi vs WT-Pi",
  "VIP1-Pi vs WT-Pi",
  "VIP1-Pi vs VIP1+Pi"
)
volcano_files <- c(
  "Volcano_WT-Pi_vs_WT+Pi.png",
  "Volcano_VIP1+Pi_vs_WT+Pi.png",
  "Volcano_VIP1-Pi_vs_WT+Pi.png",
  "Volcano_VIP1+Pi_vs_WT-Pi.png",
  "Volcano_VIP1-Pi_vs_WT-Pi.png",
  "Volcano_VIP1-Pi_vs_VIP1+Pi.png"
)
pairwise_results <- list(
  res_WT_minus_vs_WT_plus_df,
  res_VIP1_plus_vs_WT_plus_df,
  res_VIP1_minus_vs_WT_plus_df,
  res_VIP1_plus_vs_WT_minus_df,
  res_VIP1_minus_vs_WT_minus_df,
  res_VIP1_minus_vs_VIP1_plus_df
)
print("Checking gene-name annotations:")
for (i in seq_along(pairwise_results)) {
  print(
    paste(
      volcano_titles[i],
      "columns:",
      paste(c("Systematic_ID", "Gene_Name") %in% colnames(pairwise_results[[i]]), collapse = ", ")
    )
  )
  print(
    paste(
      "Number of gene names:",
      sum(!is.na(pairwise_results[[i]]$Gene_Name) & pairwise_results[[i]]$Gene_Name != "")
    )
  )
}
for (i in seq_along(pairwise_results)) {
  make_volcano(
    result_df = pairwise_results[[i]],
    plot_title = volcano_titles[i],
    output_file = file.path(output_dir, "03_Volcano", volcano_files[i])
  )
}
make_volcano(
  result_df = res_interaction_df,
  plot_title = "Genotype × Pi interaction",
  output_file = file.path(output_dir, "03_Volcano", "Volcano_Genotype_x_Pi_interaction.png")
)

# 23. MA plot function ---------------

make_ma_plot <- function(result_df, plot_title, output_file) {
  plot_df <- result_df %>%
  filter(!is.na(baseMean), !is.na(log2FoldChange))
  p <- ggplot(plot_df, aes(x = log10(baseMean + 1), y = log2FoldChange)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(title = plot_title, x = "log10 mean normalized expression", y = "log2 fold change")
  ggsave(output_file, p, width = 9, height = 7, dpi = 300)
}

# 24. MA plots for all comparisons ---------------

ma_files <- c(
  "MA_WT-Pi_vs_WT+Pi.png",
  "MA_VIP1+Pi_vs_WT+Pi.png",
  "MA_VIP1-Pi_vs_WT+Pi.png",
  "MA_VIP1+Pi_vs_WT-Pi.png",
  "MA_VIP1-Pi_vs_WT-Pi.png",
  "MA_VIP1-Pi_vs_VIP1+Pi.png"
)
for (i in seq_along(pairwise_results)) {
  make_ma_plot(
    pairwise_results[[i]],
    volcano_titles[i],
    file.path(output_dir, "04_MA_plots", ma_files[i])
  )
}
make_ma_plot(
  res_interaction_df,
  "Genotype × Pi interaction",
  file.path(output_dir, "04_MA_plots", "MA_Genotype_x_Pi_interaction.png")
)

# 25. Top variable gene heatmap ---------------

print("CREATING TOP VARIABLE GENE HEATMAP")
gene_variances <- rowVars(vst_matrix)
top_n <- min(50, length(gene_variances))
top_variable_genes <- names(sort(gene_variances, decreasing = TRUE))[1:top_n]
heatmap_matrix <- vst_matrix[
  top_variable_genes,
  ,
  drop = FALSE
]
heatmap_gene_names <- get_gene_names(rownames(heatmap_matrix))
rownames(heatmap_matrix) <- make.unique(heatmap_gene_names)
heatmap_scaled <- t(scale(t(heatmap_matrix)))
annotation_col <- data.frame(Condition = metadata[ colnames(heatmap_scaled), "Display_Name" ])
rownames(annotation_col) <- colnames(heatmap_scaled)
png(
  file.path(output_dir, "05_Heatmaps", "Heatmap_Top_50_Variable_Genes.png"),
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

# 26. DE heatmap function ---------------

make_de_heatmap <- function(result_df, plot_title, output_file, n_genes = 50) {
  significant_genes <- result_df %>%
  filter(!is.na(padj), padj < 0.05) %>%
  arrange(padj)
  if (nrow(significant_genes) < 2) {
    print(paste("Too few significant genes for:", plot_title))
    return(NULL)
  }
  selected_genes <- head(significant_genes$Systematic_ID, n_genes)
  selected_genes <- intersect(selected_genes, rownames(vst_matrix))
  if (length(selected_genes) < 2) {
    return(NULL)
  }
  heatmap_data <- vst_matrix[
    selected_genes,
    ,
    drop = FALSE
  ]
  gene_names <- get_gene_names(rownames(heatmap_data))
  rownames(heatmap_data) <- make.unique(gene_names)
  heatmap_scaled <- t(scale(t(heatmap_data)))
  annotation_col <- data.frame(Condition = metadata[ colnames(heatmap_scaled), "Display_Name" ])
  rownames(annotation_col) <- colnames(heatmap_scaled)
  png(output_file, width = 2400, height = 2600, res = 300)
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

# 27. Heatmaps for all six pairwise comparisons ---------------

heatmap_files <- c(
  "Heatmap_WT-Pi_vs_WT+Pi.png",
  "Heatmap_VIP1+Pi_vs_WT+Pi.png",
  "Heatmap_VIP1-Pi_vs_WT+Pi.png",
  "Heatmap_VIP1+Pi_vs_WT-Pi.png",
  "Heatmap_VIP1-Pi_vs_WT-Pi.png",
  "Heatmap_VIP1-Pi_vs_VIP1+Pi.png"
)
for (i in seq_along(pairwise_results)) {
  make_de_heatmap(
    pairwise_results[[i]],
    paste("Differential expression:", volcano_titles[i]),
    file.path(output_dir, "05_Heatmaps", heatmap_files[i])
  )
}
make_de_heatmap(
  res_interaction_df,
  "Genotype × Pi interaction",
  file.path(output_dir, "05_Heatmaps", "Heatmap_Genotype_x_Pi_interaction.png")
)

# 28. Condition-level VST heatmap ---------------

# Replicates are averaged ONLY for visualization.
# Conditions:
# WT+Pi
# WT-Pi
# VIP1+Pi
# VIP1-Pi
print("CREATING CONDITION-LEVEL HEATMAP")
condition_means <- matrix(NA_real_, nrow = nrow(vst_matrix), ncol = 4)
# PRESERVE GENE IDs
rownames(condition_means) <- rownames(vst_matrix)
colnames(condition_means) <- c("WT+Pi", "WT-Pi", "VIP1+Pi", "VIP1-Pi")
condition_indices <- list(
  "WT+Pi" = which(metadata$Condition == "WT_plus"),
  "WT-Pi" = which(metadata$Condition == "WT_minus"),
  "VIP1+Pi" = which(metadata$Condition == "VIP1_plus"),
  "VIP1-Pi" = which(metadata$Condition == "VIP1_minus")
)
print("Condition sample indices:")
for (condition_name in names(condition_indices)) {
  print(paste(condition_name, ":", paste(condition_indices[[condition_name]], collapse = ", ")))
}
condition_replicate_counts <- sapply(condition_indices, length)
if (any(condition_replicate_counts != 3)) {
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
# Replicates are averaged ONLY here for visualization.
# DESeq2 statistics continue to use all biological replicates.
for (i in seq_along(condition_indices)) {
  condition_means[, i] <- rowMeans(vst_matrix[ , condition_indices[[i]], drop = FALSE ], na.rm = TRUE)
}
finite_genes <- apply(condition_means, 1, function(x) { all(is.finite(x)) })
print(paste("Genes with finite values across all four conditions:", sum(finite_genes)))
condition_means <- condition_means[
  finite_genes,
  ,
  drop = FALSE
]
condition_variances <- matrixStats::rowVars(condition_means)
# These genes would become NaN during row scaling.
variable_genes <- is.finite(condition_variances) &
condition_variances > 0
print(paste("Genes with variation across conditions:", sum(variable_genes)))
condition_variances <- condition_variances[
  variable_genes
]
condition_means <- condition_means[
  variable_genes,
  ,
  drop = FALSE
]
# Use numeric row indices rather than names().
# matrixStats::rowVars() does not necessarily preserve
# row names as names of the returned variance vector.
n_condition_genes <- min(50, length(condition_variances))
if (n_condition_genes < 2) {
  stop(
    paste0(
      "\nERROR: Fewer than 2 variable genes are available ",
      "for the condition-level heatmap."
    )
  )
}
top_condition_indices <- order(condition_variances, decreasing = TRUE)[
  seq_len(n_condition_genes)
]
condition_heatmap <- condition_means[
  top_condition_indices,
  ,
  drop = FALSE
]
print(paste("Number of genes selected for condition heatmap:", nrow(condition_heatmap)))
# MAP SYSTEMATIC IDs TO CGD GENE NAMES
condition_gene_names <- get_gene_names(rownames(condition_heatmap))
rownames(condition_heatmap) <- make.unique(condition_gene_names)
# Each gene is standardized across the four conditions.
# This shows relative expression patterns rather than
# absolute VST expression.
condition_heatmap_scaled <- t(scale(t(condition_heatmap)))
finite_scaled_genes <- apply(condition_heatmap_scaled, 1, function(x) { all(is.finite(x)) })
condition_heatmap_scaled <- condition_heatmap_scaled[
  finite_scaled_genes,
  ,
  drop = FALSE
]
print("Condition heatmap dimensions:")
print(dim(condition_heatmap_scaled))
if (nrow(condition_heatmap_scaled) < 2) {
  stop("\nERROR: Fewer than 2 genes remain after scaling.")
}
if (ncol(condition_heatmap_scaled) != 4) {
  stop("\nERROR: Condition heatmap does not contain exactly 4 conditions.")
}
write.csv(
  condition_means,
  file.path(output_dir, "05_Heatmaps", "Condition_Mean_VST_Expression.csv"),
  row.names = TRUE
)
png(
  file.path(output_dir, "05_Heatmaps", "Heatmap_Top_50_Variable_Condition_Means.png"),
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
print("SUCCESS: Condition-level VST heatmap created.")
print(paste("Genes plotted:", nrow(condition_heatmap_scaled)))
print("Conditions plotted:")
print(colnames(condition_heatmap_scaled))

# 29. Phosphate response analysis ---------------

# WT response:
# WT-Pi - WT+Pi
# VIP1 response:
# VIP1-Pi - VIP1+Pi
# The difference between these responses is tested formally
# by the Genotype × Pi interaction.
print("PHOSPHATE RESPONSE ANALYSIS")
WT_response <- res_WT_minus_vs_WT_plus_df$log2FoldChange
VIP1_response <- res_VIP1_minus_vs_VIP1_plus_df$log2FoldChange
response_table <- data.frame(
  Systematic_ID = res_WT_minus_vs_WT_plus_df$Systematic_ID,
  Gene_Name = res_WT_minus_vs_WT_plus_df$Gene_Name,
  WT_response = WT_response,
  WT_padj = res_WT_minus_vs_WT_plus_df$padj,
  VIP1_response = VIP1_response[
    match(res_WT_minus_vs_WT_plus_df$Systematic_ID, res_VIP1_minus_vs_VIP1_plus_df$Systematic_ID)
  ],
  VIP1_padj = res_VIP1_minus_vs_VIP1_plus_df$padj[
    match(res_WT_minus_vs_WT_plus_df$Systematic_ID, res_VIP1_minus_vs_VIP1_plus_df$Systematic_ID)
  ],
  Interaction_log2FC = res_interaction_df$log2FoldChange[
    match(res_WT_minus_vs_WT_plus_df$Systematic_ID, res_interaction_df$Systematic_ID)
  ],
  Interaction_padj = res_interaction_df$padj[
    match(res_WT_minus_vs_WT_plus_df$Systematic_ID, res_interaction_df$Systematic_ID)
  ],
  stringsAsFactors = FALSE
)
write.csv(
  response_table,
  file.path(output_dir, "09_Response_analysis", "Phosphate_Response_Comparison.csv"),
  row.names = FALSE
)
significant_response_genes <- response_table %>%
filter(!is.na(Interaction_padj), Interaction_padj < 0.05) %>%
arrange(Interaction_padj)
write.csv(
  significant_response_genes,
  file.path(
    output_dir,
    "09_Response_analysis",
    "Significant_Phosphate_Response_Interaction_Genes.csv"
  ),
  row.names = FALSE
)
top_interaction_genes <- significant_response_genes %>%
head(50)
write.csv(
  top_interaction_genes,
  file.path(output_dir, "09_Response_analysis", "Top_50_Interaction_Genes.csv"),
  row.names = FALSE
)

# 30. Phosphate response heatmap ---------------

if (nrow(top_interaction_genes) >= 2) {
  response_gene_ids <- intersect(top_interaction_genes$Systematic_ID, rownames(vst_matrix))
  response_matrix <- condition_means[
    response_gene_ids,
    ,
    drop = FALSE
  ]
  rownames(response_matrix) <- make.unique(get_gene_names(rownames(response_matrix)))
  response_scaled <- t(scale(t(response_matrix)))
  png(
    file.path(output_dir, "09_Response_analysis", "Heatmap_Top_Interaction_Genes.png"),
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
    main = "Top genes with Genotype × Pi interaction"
  )
  dev.off()
}

# 31. Phosphate response scatterplot ---------------

response_plot_df <- response_table %>%
filter(!is.na(WT_response), !is.na(VIP1_response))
if (nrow(response_plot_df) > 0) {
  p_response_scatter <- ggplot(response_plot_df, aes(x = WT_response, y = VIP1_response)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  theme_bw() +
  labs(
    title = "WT vs VIP1 phosphate-starvation responses",
    x = "WT response: WT-Pi vs WT+Pi",
    y = "VIP1 response: VIP1-Pi vs VIP1+Pi"
  )
  ggsave(
    file.path(output_dir, "09_Response_analysis", "WT_vs_VIP1_Phosphate_Response.png"),
    p_response_scatter,
    width = 9,
    height = 8,
    dpi = 300
  )
}

# 32. Phosphate response plot ---------------

# Only significant interaction genes are shown.
if (nrow(significant_response_genes) >= 1) {
  response_long <- data.frame(
    Gene = rep(significant_response_genes$Gene_Name, 2),
    Genotype = rep(c("WT", "VIP1"), each = nrow(significant_response_genes)),
    Response = c(significant_response_genes$WT_response, significant_response_genes$VIP1_response),
    stringsAsFactors = FALSE
  )
  response_plot <- ggplot(
    response_long,
    aes(x = Response, y = reorder(Gene, Response), shape = Genotype)
  ) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(
    title = "Phosphate-response differences in significant interaction genes",
    x = "log2 fold change: -Pi versus +Pi",
    y = "Gene"
  )
  # Keep large gene sets within the ggsave height limit.
  response_height <- min(max(6, 0.25 * nrow(significant_response_genes)), 48)
  ggsave(
    file.path(output_dir, "09_Response_analysis", "Phosphate_response_plot.png"),
    response_plot,
    width = 10,
    height = response_height,
    dpi = 300
  )
}

# 33. Summary table ---------------

summary_table <- data.frame(
  Comparison = names(all_results),
  Total_genes_tested = sapply(all_results, nrow),
  Significant_padj_0.05 = sapply(
    all_results,
    function(df) {
      sum(!is.na(df$padj) & df$padj < 0.05)
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
  file.path(output_dir, "02_DE_results", "DE_summary.csv"),
  row.names = FALSE
)

# 34. Comparison summary print ---------------

print("SIGNIFICANT GENE COUNTS")
print(summary_table)

# 35. Save VST object ---------------

saveRDS(vsd, file.path(output_dir, "01_QC", "VST_object.rds"))

# 36. Venn diagrams: phosphate-starvation response ---------------

# All comparisons are:
# -Pi vs +Pi
# +Pi is the reference.
# Positive log2FC = UPREGULATED under -Pi
# Negative log2FC = DOWNREGULATED under -Pi
# We compare:
# WT-only vs VIP1-only vs Common
# WT-only vs VIP1-only vs Common
# No cross-direction comparisons are made.
library(VennDiagram)
library(grid)
venn_dir <- file.path(output_dir, "04_Venn")
if (!dir.exists(venn_dir)) {
  dir.create(venn_dir, recursive = TRUE, showWarnings = FALSE)
}
# USE -Pi VS +Pi RESULTS
wt_response <- res_WT_minus_vs_WT_plus_df %>%
filter(!is.na(padj), !is.na(log2FoldChange))
vip1_response <- res_VIP1_minus_vs_VIP1_plus_df %>%
filter(!is.na(padj), !is.na(log2FoldChange))
# DEFINE UPREGULATED GENES UNDER -Pi
wt_up <- wt_response %>%
filter(padj < 0.05, log2FoldChange >= 1) %>%
pull(Systematic_ID)
vip1_up <- vip1_response %>%
filter(padj < 0.05, log2FoldChange >= 1) %>%
pull(Systematic_ID)
# DEFINE DOWNREGULATED GENES UNDER -Pi
wt_down <- wt_response %>%
filter(padj < 0.05, log2FoldChange <= -1) %>%
pull(Systematic_ID)
vip1_down <- vip1_response %>%
filter(padj < 0.05, log2FoldChange <= -1) %>%
pull(Systematic_ID)
common_up <- intersect(wt_up, vip1_up)
wt_only_up <- setdiff(wt_up, vip1_up)
vip1_only_up <- setdiff(vip1_up, wt_up)
common_down <- intersect(wt_down, vip1_down)
wt_only_down <- setdiff(wt_down, vip1_down)
vip1_only_down <- setdiff(vip1_down, wt_down)
cat("PHOSPHATE-STARVATION RESPONSE: -Pi VS +Pi\n")
cat("\nUPREGULATED UNDER -Pi\n")
cat("----------------------\n")
cat("WT total:       ", length(wt_up), "\n")
cat("VIP1 total:     ", length(vip1_up), "\n")
cat("Common:         ", length(common_up), "\n")
cat("WT-only:        ", length(wt_only_up), "\n")
cat("VIP1-only:      ", length(vip1_only_up), "\n")
cat("\nDOWNREGULATED UNDER -Pi\n")
cat("------------------------\n")
cat("WT total:       ", length(wt_down), "\n")
cat("VIP1 total:     ", length(vip1_down), "\n")
cat("Common:         ", length(common_down), "\n")
cat("WT-only:        ", length(wt_only_down), "\n")
cat("VIP1-only:      ", length(vip1_only_down), "\n")
make_venn <- function(set1, set2, name1, name2, plot_title, output_file) {
  venn <- venn.diagram(
    x = list(set1 = set1, set2 = set2),
    filename = NULL,
    width = 3200,
    height = 2800,
    resolution = 300,
    fill = c("grey80", "grey60"),
    alpha = 0.5,
    cex = 1.8,
    fontface = "bold",
    category.names = c(name1, name2),
    cat.cex = 1.2,
    cat.fontface = "bold",
    cat.dist = c(0.08, 0.08),
    cat.pos = c(-15, 15),
    margin = 0.25,
    main = plot_title,
    main.cex = 1.5,
    main.fontface = "bold"
  )
  png(filename = output_file, width = 3200, height = 2800, res = 300)
  grid.newpage()
  pushViewport(viewport(x = 0.5, y = 0.47, width = 0.78, height = 0.70))
  grid.draw(venn)
  popViewport()
  dev.off()
  cat("Saved:", output_file, "\n")
}
# VENN 1: UPREGULATED UNDER -Pi
make_venn(
  set1 = wt_up,
  set2 = vip1_up,
  name1 = "WT Upregulated",
  name2 = "VIP1 Upregulated",
  plot_title = "Genes Upregulated Under -Pi",
  output_file = file.path(venn_dir, "Venn_Common_Upregulated_WT_vs_VIP1.png")
)
# VENN 2: DOWNREGULATED UNDER -Pi
make_venn(
  set1 = wt_down,
  set2 = vip1_down,
  name1 = "WT Downregulated",
  name2 = "VIP1 Downregulated",
  plot_title = "Genes Downregulated Under -Pi",
  output_file = file.path(venn_dir, "Venn_Common_Downregulated_WT_vs_VIP1.png")
)
write.csv(
  data.frame(Systematic_ID = common_up, Gene_Name = get_gene_names(common_up)),
  file.path(venn_dir, "Common_Upregulated_WT_VIP1.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(Systematic_ID = wt_only_up, Gene_Name = get_gene_names(wt_only_up)),
  file.path(venn_dir, "WT_Only_Upregulated.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(Systematic_ID = vip1_only_up, Gene_Name = get_gene_names(vip1_only_up)),
  file.path(venn_dir, "VIP1_Only_Upregulated.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(Systematic_ID = common_down, Gene_Name = get_gene_names(common_down)),
  file.path(venn_dir, "Common_Downregulated_WT_VIP1.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(Systematic_ID = wt_only_down, Gene_Name = get_gene_names(wt_only_down)),
  file.path(venn_dir, "WT_Only_Downregulated.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(Systematic_ID = vip1_only_down, Gene_Name = get_gene_names(vip1_only_down)),
  file.path(venn_dir, "VIP1_Only_Downregulated.csv"),
  row.names = FALSE
)
cat("SECTION 22 COMPLETE\n")
cat("Venn diagrams saved in:\n")
cat(venn_dir, "\n")

# 37. Top 50 genes for WT vs VIP1 under -Pi ---------------

# +Pi is the reference condition.
# Positive log2FC = higher under -Pi
# Negative log2FC = lower under -Pi
wt_response <- res_WT_minus_vs_WT_plus_df %>%
filter(!is.na(padj), !is.na(log2FoldChange))
vip1_response <- res_VIP1_minus_vs_VIP1_plus_df %>%
filter(!is.na(padj), !is.na(log2FoldChange))
# WT genes increased under -Pi
wt_up <- wt_response %>%
filter(padj < 0.05, log2FoldChange >= 1) %>%
pull(Systematic_ID)
# WT genes decreased under -Pi
wt_down <- wt_response %>%
filter(padj < 0.05, log2FoldChange <= -1) %>%
pull(Systematic_ID)
# VIP1 genes increased under -Pi
vip1_up <- vip1_response %>%
filter(padj < 0.05, log2FoldChange >= 1) %>%
pull(Systematic_ID)
# VIP1 genes decreased under -Pi
vip1_down <- vip1_response %>%
filter(padj < 0.05, log2FoldChange <= -1) %>%
pull(Systematic_ID)
# Upregulated under -Pi
common_up <- intersect(wt_up, vip1_up)
wt_only_up <- setdiff(wt_up, vip1_up)
vip1_only_up <- setdiff(vip1_up, wt_up)
# Downregulated under -Pi
common_down <- intersect(wt_down, vip1_down)
wt_only_down <- setdiff(wt_down, vip1_down)
vip1_only_down <- setdiff(vip1_down, wt_down)
get_top50 <- function(result_df, gene_ids, category) {
  result_df %>%
  filter(Systematic_ID %in% gene_ids) %>%
  arrange(padj, desc(abs(log2FoldChange))) %>%
  mutate(Category = category) %>%
  dplyr::select(Category, Systematic_ID, Gene_Name, log2FoldChange, padj) %>%
  slice_head(n = 50)
}
top50_wt_only_up <- get_top50(wt_response, wt_only_up, "WT-only Upregulated under -Pi")
top50_vip1_only_up <- get_top50(vip1_response, vip1_only_up, "VIP1-only Upregulated under -Pi")
# Rank by WT padj, then show BOTH WT and VIP1 statistics.
top50_common_up <- wt_response %>%
filter(Systematic_ID %in% common_up) %>%
dplyr::select(Systematic_ID, Gene_Name, WT_log2FoldChange = log2FoldChange, WT_padj = padj) %>%
inner_join(
  vip1_response %>%
  filter(Systematic_ID %in% common_up) %>%
  dplyr::select(Systematic_ID, VIP1_log2FoldChange = log2FoldChange, VIP1_padj = padj),
  by = "Systematic_ID"
) %>%
arrange(WT_padj, desc(abs(WT_log2FoldChange))) %>%
slice_head(n = 50)
top50_wt_only_down <- get_top50(wt_response, wt_only_down, "WT-only Downregulated under -Pi")
top50_vip1_only_down <- get_top50(vip1_response, vip1_only_down, "VIP1-only Downregulated under -Pi")
# Rank by WT padj, then show BOTH WT and VIP1 statistics.
top50_common_down <- wt_response %>%
filter(Systematic_ID %in% common_down) %>%
dplyr::select(Systematic_ID, Gene_Name, WT_log2FoldChange = log2FoldChange, WT_padj = padj) %>%
inner_join(
  vip1_response %>%
  filter(Systematic_ID %in% common_down) %>%
  dplyr::select(Systematic_ID, VIP1_log2FoldChange = log2FoldChange, VIP1_padj = padj),
  by = "Systematic_ID"
) %>%
arrange(WT_padj, desc(abs(WT_log2FoldChange))) %>%
slice_head(n = 50)
write.csv(
  top50_wt_only_up,
  file.path(venn_dir, "Top50_WT_Only_Upregulated_under_minusPi.csv"),
  row.names = FALSE
)
write.csv(
  top50_vip1_only_up,
  file.path(venn_dir, "Top50_VIP1_Only_Upregulated_under_minusPi.csv"),
  row.names = FALSE
)
write.csv(
  top50_common_up,
  file.path(venn_dir, "Top50_Common_Upregulated_under_minusPi.csv"),
  row.names = FALSE
)
write.csv(
  top50_wt_only_down,
  file.path(venn_dir, "Top50_WT_Only_Downregulated_under_minusPi.csv"),
  row.names = FALSE
)
write.csv(
  top50_vip1_only_down,
  file.path(venn_dir, "Top50_VIP1_Only_Downregulated_under_minusPi.csv"),
  row.names = FALSE
)
write.csv(
  top50_common_down,
  file.path(venn_dir, "Top50_Common_Downregulated_under_minusPi.csv"),
  row.names = FALSE
)
cat("TOP 50 GENES: -Pi VS +Pi\n")
cat("\nUPREGULATED UNDER -Pi\n")
cat("WT-only:   ", nrow(top50_wt_only_up), "\n")
cat("VIP1-only: ", nrow(top50_vip1_only_up), "\n")
cat("Common:    ", nrow(top50_common_up), "\n")
cat("\nDOWNREGULATED UNDER -Pi\n")
cat("WT-only:   ", nrow(top50_wt_only_down), "\n")
cat("VIP1-only: ", nrow(top50_vip1_only_down), "\n")
cat("Common:    ", nrow(top50_common_down), "\n")
cat("\nFiles saved to:\n")
cat(venn_dir, "\n")

# 38. Create excel file of top 50 venn genes ---------------

library(openxlsx)
excel_file <- file.path(venn_dir, "Top50_WT_vs_VIP1_under_minusPi.xlsx")
wb <- createWorkbook()
# WT-only upregulated
addWorksheet(wb, "WT-only Up")
writeData(wb, "WT-only Up", top50_wt_only_up)
# VIP1-only upregulated
addWorksheet(wb, "VIP1-only Up")
writeData(wb, "VIP1-only Up", top50_vip1_only_up)
# Common upregulated
addWorksheet(wb, "Common Up")
writeData(wb, "Common Up", top50_common_up)
# WT-only downregulated
addWorksheet(wb, "WT-only Down")
writeData(wb, "WT-only Down", top50_wt_only_down)
# VIP1-only downregulated
addWorksheet(wb, "VIP1-only Down")
writeData(wb, "VIP1-only Down", top50_vip1_only_down)
# Common downregulated
addWorksheet(wb, "Common Down")
writeData(wb, "Common Down", top50_common_down)
for (sheet in names(wb)) {
  # Bold header
  addStyle(
    wb,
    sheet = sheet,
    style = createStyle(textDecoration = "bold"),
    rows = 1,
    cols = 1:10,
    gridExpand = TRUE
  )
  # Freeze header row
  freezePane(wb, sheet = sheet, firstRow = TRUE)
  # Auto-width columns
  setColWidths(wb, sheet = sheet, cols = 1:10, widths = "auto")
}
saveWorkbook(wb, excel_file, overwrite = TRUE)
cat("EXCEL FILE CREATED\n")
cat("File:\n")
cat(excel_file, "\n")

# 39. GO annotations and identifier mapping -------------------------------

if (!file.exists("gene_association.cgd")) {
  stop("Main analysis saved. Add gene_association.cgd to run GO enrichment.")
}

library(clusterProfiler)
library(GO.db)

# GAF uses CGD database IDs; match those through names and systematic aliases.
gaf <- read.delim(
  "gene_association.cgd", header = FALSE, sep = "\t",
  quote = "", comment.char = "!", fill = TRUE,
  stringsAsFactors = FALSE
)
if (ncol(gaf) < 11) stop("Expected at least 11 columns in the CGD GAF file.")

go_annotations <- gaf %>%
transmute(
  CGD_ID = trimws(V2), Gene_Name = trimws(V3),
  GO_ID = trimws(V5), Ontology = trimws(V9), Qualifier = V4
) %>%
filter(
  CGD_ID != "", grepl("^GO:[0-9]+$", GO_ID),
  Ontology %in% c("P", "F", "C"),
  !grepl("(^|\\|)NOT($|\\|)", Qualifier)
) %>%
distinct()

alias_table <- bind_rows(lapply(seq_len(nrow(gaf)), function(i) {
      aliases <- unique(c(gaf$V2[i], gaf$V3[i], strsplit(gaf$V11[i], "|", fixed = TRUE)[[1]]))
      data.frame(CGD_ID = gaf$V2[i], Alias = aliases)
})) %>%
filter(!is.na(Alias), Alias != "") %>%
mutate(Key = normalize_id(Alias)) %>%
distinct(Key, CGD_ID)

# Ambiguous aliases are excluded so one RNA-seq gene cannot map arbitrarily.
alias_table <- alias_table %>%
group_by(Key) %>%
filter(n_distinct(CGD_ID) == 1) %>%
ungroup()

rna_ids <- rownames(counts_filtered)
rna_mapping <- data.frame(RNAseq_ID = rna_ids, CGD_ID = NA_character_)
for (i in seq_along(rna_ids)) {
  candidates <- c(rna_ids[i], sub("^CAALFM_", "", rna_ids[i]), get_gene_names(rna_ids[i]))
  matches <- unique(alias_table$CGD_ID[alias_table$Key %in% normalize_id(candidates)])
  if (length(matches) == 1) rna_mapping$CGD_ID[i] <- matches
}

convert_to_cgd <- function(genes) {
  unique(na.omit(rna_mapping$CGD_ID[rna_mapping$RNAseq_ID %in% genes]))
}
background_cgd <- convert_to_cgd(rna_ids)
go_output_dir <- file.path(output_dir, "10_GO_Analysis")
dir.create(go_output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(rna_mapping, file.path(go_output_dir, "RNAseq_to_CGD_mapping.csv"), row.names = FALSE)
message("GO mapping: ", sum(!is.na(rna_mapping$CGD_ID)), " of ", length(rna_ids), " genes.")

# 40. GO enrichment for phosphate starvation -------------------------------

get_response_genes <- function(result, direction) {
  result %>%
  filter(!is.na(padj), padj < 0.05, !is.na(log2FoldChange), direction * log2FoldChange >= 1) %>%
  pull(Systematic_ID) %>%
  unique() %>%
  convert_to_cgd()
}

starvation_gene_sets <- list(
  WT_Starvation_Up = get_response_genes(res_WT_minus_vs_WT_plus_df, 1),
  WT_Starvation_Down = get_response_genes(res_WT_minus_vs_WT_plus_df, -1),
  VIP1_Starvation_Up = get_response_genes(res_VIP1_minus_vs_VIP1_plus_df, 1),
  VIP1_Starvation_Down = get_response_genes(res_VIP1_minus_vs_VIP1_plus_df, -1)
)
ontology_codes <- c(BP = "P", MF = "F", CC = "C")
starvation_go_results <- list()

for (set_name in names(starvation_gene_sets)) {
  for (ontology_name in names(ontology_codes)) {
    term2gene <- go_annotations %>%
    filter(Ontology == ontology_codes[[ontology_name]]) %>%
    dplyr::select(GO_ID, CGD_ID) %>%
    distinct()
    universe <- intersect(background_cgd, term2gene$CGD_ID)
    genes <- intersect(starvation_gene_sets[[set_name]], universe)
    if (!length(genes) || !length(universe)) next

    term_ids <- unique(term2gene$GO_ID)
    term2name <- data.frame(
      GO_ID = term_ids,
      Description = vapply(term_ids, function(id) {
          term <- GO.db::GOTERM[[id]]
          if (is.null(term)) id else AnnotationDbi::Term(term)
        }, character(1))
    )
    result <- clusterProfiler::enricher(
      gene = genes, universe = universe,
      TERM2GENE = term2gene, TERM2NAME = term2name,
      pvalueCutoff = 0.05, qvalueCutoff = 0.05,
      pAdjustMethod = "BH", minGSSize = 5, maxGSSize = 500
    )
    if (is.null(result)) next
    df <- as.data.frame(result)
    if (!nrow(df)) next
    df$Gene_Set <- set_name
    df$Ontology <- ontology_name
    starvation_go_results[[paste(set_name, ontology_name, sep = "_")]] <- df
  }
}

# 41. Save GO tables and plots ---------------------------------------------

wb <- openxlsx::createWorkbook()
combined_go <- bind_rows(starvation_go_results)
openxlsx::addWorksheet(wb, "All_GO_results")
if (nrow(combined_go)) {
  openxlsx::writeData(wb, "All_GO_results", combined_go)
  write.csv(combined_go, file.path(go_output_dir, "Starvation_GO_results.csv"), row.names = FALSE)
} else {
  openxlsx::writeData(wb, "All_GO_results", "No GO terms passed the enrichment cutoffs.")
}
for (ontology_name in names(ontology_codes)) {
  openxlsx::addWorksheet(wb, ontology_name)
  if (nrow(combined_go)) {
    openxlsx::writeData(wb, ontology_name, filter(combined_go, Ontology == ontology_name))
  }
}
openxlsx::saveWorkbook(wb, file.path(go_output_dir, "Phosphate_Starvation_GO.xlsx"), overwrite = TRUE)

if (nrow(combined_go)) {
  plot_data <- combined_go %>%
  filter(!is.na(p.adjust), p.adjust < 0.05) %>%
  group_by(Gene_Set, Ontology) %>%
  slice_min(p.adjust, n = 10, with_ties = FALSE) %>%
  ungroup()
  for (set_name in unique(plot_data$Gene_Set)) {
    df <- filter(plot_data, Gene_Set == set_name)
    plot <- ggplot(df, aes(-log10(pmax(p.adjust, .Machine$double.xmin)),
        reorder(Description, -p.adjust), fill = Ontology)) +
    geom_col() + facet_wrap(~ Ontology, scales = "free_y", ncol = 1) +
    theme_bw() + labs(title = gsub("_", " ", set_name), x = "-log10(adjusted p-value)", y = NULL)
    ggsave(file.path(go_output_dir, paste0(set_name, "_GO.png")), plot,
      width = 10, height = max(6, 0.3 * nrow(df)), dpi = 300)
  }
}

# 42. Finish ---------------------------------------------------------------

message("VIP1 analysis complete. Results saved in: ", normalizePath(output_dir))
```
