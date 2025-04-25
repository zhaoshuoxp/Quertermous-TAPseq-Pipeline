#!/usr/bin/env Rscript

# --------------------
required_packages <- c("sceptre", "Seurat", "reshape2")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("Installing missing package: %s\n", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  library(pkg, character.only = TRUE)
}

# ----------------------------
args <- commandArgs(trailingOnly = TRUE)

# 
mtx <- ifelse(length(args) >= 1, args[1], stop("Error: missing cellranger matrix"))
guides_file <- ifelse(length(args) >= 2, args[2], stop("Error: missing file with guides id"))
cov_file <- ifelse(length(args) >= 3, args[3], stop("Error: missing covarates file"))
grna_target <- ifelse(length(args) >= 4, args[4], stop("Error: missing grna-target data frame"))
pos <- ifelse(length(args) >= 5, args[5], stop("Error: missing positive controls"))
discovery <- ifelse(length(args) >= 6, args[6], stop("Error: missing discovery pairs"))
cores <- ifelse(length(args) >= 7, as.numeric(args[7]), 20)  

# --------------------------------
cat("Loading data...\n")
counts <- Read10X(data.dir = mtx)
gene_exp <- counts$`Gene Expression`
crispr <- counts$`CRISPR Guide Capture`

guides <- read.table(guides_file, stringsAsFactors = FALSE)$V1
cov <- read.table(cov_file, header = TRUE, row.names = 1, sep = '\t')

crispr_clean <- crispr[rownames(crispr) %in% guides, colnames(crispr) %in% rownames(cov)]
gene_exp_clean <- gene_exp[, colnames(gene_exp) %in% rownames(cov)]

grna_target_data_frame <- read.table(grna_target, header = TRUE)
positive_control_pairs <- read.table(pos, header = TRUE)
discovery_pairs <- read.table(discovery, header = TRUE)

# ---------------------------
cat("Initializing sceptre object...\n")
sceptre_object <- import_data(
  response_matrix = gene_exp_clean,
  grna_matrix = crispr_clean,
  grna_target_data_frame = grna_target_data_frame,
  moi = 'high',
  extra_covariates = cov
)

cat("Setting analysis parameters...\n")
sceptre_object <- set_analysis_parameters(
  sceptre_object = sceptre_object,
  discovery_pairs = discovery_pairs,
  positive_control_pairs = positive_control_pairs,
  grna_integration_strategy = "union",
  side = 'both'
)

cat("Assigning gRNAs...\n")
sceptre_object <- assign_grnas(sceptre_object, method = "mixture", parallel = TRUE, n_processors = cores)

cat("Running quality control...\n")
sceptre_object <- run_qc(sceptre_object, p_mito_threshold = 0.075)

cat("Running calibration check...\n")
sceptre_object <- run_calibration_check(sceptre_object, parallel = TRUE, n_processors = cores)

cat("Running power check...\n")
sceptre_object <- run_power_check(sceptre_object, parallel = TRUE, n_processors = cores)

cat("Generating plots...\n")
plot(sceptre_object)
dev.off()

cat("Writing outputs...\n")
write_outputs_to_directory(
  sceptre_object = sceptre_object,
  directory = "./"
)

cat("Saving results...\n")
result <- get_result(sceptre_object, "run_discovery_analysis")

result$p.adj <- p.adjust(
  p = result$p_value,
  method = sceptre_object@multiple_testing_method
)

write.table(result, 'results.txt', sep = '\t', quote = FALSE, row.names = FALSE)

saveRDS(sceptre_object, 'sceptre.rds')

cat("Analysis complete!\n")
