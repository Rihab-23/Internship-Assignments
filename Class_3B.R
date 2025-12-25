# ===============================================================
#                    Task_BENSALEK Rihab
# ---------------------------------------------------------
# ---------------------------------------------------------
#     Module II: Introduction to Genomics Data Analysis
# ---------------------------------------------------------
# ---------------------------------------------------------
#         Assignment 4:   Microarray Data Analysis
# ===============================================================

# -------------------------------
# 0. Install and Load Packages
# -------------------------------
# Install Bioconductor packages
BiocManager::install(c("lumi", "arrayQualityMetrics", "illuminaio"), ask = FALSE)
BiocManager::install("ArrayExpress", force = TRUE)
BiocManager::install("lumi", force = TRUE)
BiocManager::install("arrayQualityMetrics", force = TRUE)
BiocManager::install("limma", force = TRUE)


# Install CRAN packages for data manipulation
install.packages("dplyr")
install.packages("ggplot2")

# Load Required Libraries
library(ArrayExpress) 
library(illuminaio)
library(Biobase)              # Download processed data matrix only 
library(lumi)                 # Illumina expression arrays
library(arrayQualityMetrics)  # QC reports
library(dplyr)                # data manipulation
library(limma)                # For normalization and differential expression
library(ggplot2)


# 1. Set Paths
# -------------------------------
idat_dir <- "C:/Users/hp/Documents/AI_Omics_Internship_2025/Module_II/E-MTAB-8148_(1)"
metadata_file1 <- "C:/Users/hp/Documents/AI_Omics_Internship_2025/Module_II/MAGE-TAB_Files/E-MTAB-8148.sdrf.txt"

# -------------------------------
# 2. List IDAT Files
# -------------------------------
idat_files <- list.files(idat_dir, pattern = "[iI][dD][aA][tT]$", full.names = TRUE)
cat("Total IDAT files found:", length(idat_files), "\n")


# -------------------------------
# 3. Load Metadata
# -------------------------------
phenotype_data <- read.delim(metadata_file1, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Optional: check column names and candidate columns for disease info
colnames(phenotype_data)
candidate_cols <- grep("Characteristics|Factor Value|Source Type", colnames(phenotype_data), value = TRUE)

# -------------------------------
# 4. Read IDAT Files in Batches
# -------------------------------
batch_size <- 50
intensity_matrix <- NULL
probe_names <- NULL
successful_files <- character(0)

for(i in seq(1, length(idat_files), by = batch_size)){
  batch_files <- idat_files[i:min(i+batch_size-1, length(idat_files))]
  cat("Reading batch", i, "to", min(i+batch_size-1, length(idat_files)), "\n")
  batch_list <- lapply(batch_files, function(f){
    ok <- try(readIDAT(f), silent = TRUE)
    if(inherits(ok, "try-error")) {
      warning("Failed to read: ", f)
      return(NULL)
    }
    ok
  })
  # remove NULLs (failed reads) and keep names of successful files
  valid_idx <- !sapply(batch_list, is.null)
  if(!any(valid_idx)) next
  batch_list <- batch_list[valid_idx]
  successful_files <- c(successful_files, batch_files[valid_idx])
  
  if(is.null(probe_names)) probe_names <- batch_list[[1]]$Manifest$Name
  batch_matrix <- sapply(batch_list, function(x){
    # adapt depending on readIDAT structure; try these options
    if(!is.null(x$Quants$Mean)) return(x$Quants$Mean)
    if(!is.null(x$Intensity)) return(x$Intensity)
    stop("Cannot find intensity vector in readIDAT output")
  })
  # ensure matrix orientation: probes × samples
  if(is.vector(batch_matrix)) batch_matrix <- matrix(batch_matrix, ncol = 1)
  if(is.null(intensity_matrix)) intensity_matrix <- batch_matrix else intensity_matrix <- cbind(intensity_matrix, batch_matrix)
  
  rm(batch_list, batch_matrix); gc()
}

# sanity checks
if(is.null(intensity_matrix)) stop("No IDATs were successfully read.")
rownames(intensity_matrix) <- probe_names
colnames(intensity_matrix) <- gsub("\\.idat$","", basename(successful_files))

cat("Final intensity matrix dims (probes x samples):", dim(intensity_matrix), "\n")
cat("Successful files read:", length(successful_files), "\n")

# 2) Align phenotype metadata to the actual samples we have
# Try to find column in phenotype_data that holds the raw filename (common names: Array.Data.File, Source.Name)
# We'll try Array.Data.File first, otherwise look for a column that contains sample names
md <- phenotype_data

# Detect candidate column that contains the file names (loosely)
cands <- c("Array.Data.File","Array.Data.File_name","Array.Data.File","Array.Data.File.URI","Source.Name")
found <- intersect(cands, colnames(md))
if(length(found) == 0){
  # fallback: look for any column which contains the first sample name (partial match)
  maybe <- sapply(md, function(col) any(grepl(gsub("\\.idat$","", basename(successful_files)[1]), col, ignore.case = TRUE)))
  if(any(maybe)){
    found <- names(maybe)[which(maybe)[1]]
    message("Auto-detected phenotype column: ", found)
  } else {
    message("No obvious phenotype column found. Showing head of metadata columns:")
    print(colnames(md))
    stop("Please identify which metadata column matches the sample names from IDAT filenames.")
  }
}

file_col <- found[1]
cat("Using phenotype column:", file_col, "\n")

# create matching vector: remove .idat suffixes from phenotype column values for matching
#phen_names <- gsub("\\.idat$","", as.character(md[[file_col]]))
#sample_names_actual <- colnames(intensity_matrix)

#ord <- match(sample_names_actual, phen_names)
if(any(is.na(ord))){
  missing_samples <- sample_names_actual[is.na(ord)]
  stop("The following samples are present in intensity matrix but not found in phenotype metadata column '", file_col, "':\n", paste(missing_samples, collapse = ", "))
}

rownames(intensity_matrix) <- probe_names
colnames(intensity_matrix) <- gsub("\\.idat$","", basename(successful_files))

cat("Intensity matrix dims (probes x samples):", dim(intensity_matrix), "\n")
cat("Successful files read:", length(successful_files), "\n")


# -------------------------------
# 5. Align Phenotype Metadata with Samples
# -------------------------------
# Match sample names between metadata and intensity matrix
phen_names <- gsub("\\.idat$","", phenotype_data$Array.Data.File)
ord <- match(colnames(intensity_matrix), phen_names)
md2 <- phenotype_data[ord, , drop = FALSE]
rownames(md2) <- colnames(intensity_matrix)

# -------------------------------
# 6. Build ExpressionSet
# -------------------------------
pheno_adf <- AnnotatedDataFrame(md2)
eset <- ExpressionSet(assayData = as.matrix(intensity_matrix), phenoData = pheno_adf)

head(colnames(intensity_matrix))
head(rownames(pheno_adf))
dim(intensity_matrix)
dim(pheno_adf)

# Quick checks to make sure data loaded successfully & correctly 
cat("ExpressionSet dimensions (probes x samples):", dim(exprs(eset)), "\n") # dim(exprs(eset)) → prints the number of probes × number of samples, tells how many features and samples are there
cat("pData rows:", nrow(pData(eset)), 
    "should equal ncol(exprs(eset)):", 
    ncol(exprs(eset)), "\n") # nrow(pData(eset)) vs ncol(exprs(eset)) → verifies that the number of rows in the phenotype data matches the number of samples in the expression matrix.

# -------------------------------
# 7. QC Before Normalization
# -------------------------------
arrayQualityMetrics(expressionset = eset,
                    outdir = "Results/QC_Raw_Data",
                    force = TRUE,
                    do.logtransform = TRUE)
# OUTPUT: HTML report with flagged arrays (outliers)

# -------------------------------
# 8. Normalization
# -------------------------------
exprs(eset) <- normalizeBetweenArrays(exprs(eset), method = "quantile")
eset_norm <- eset

# QC after normalization
arrayQualityMetrics(expressionset = eset_norm,
                    outdir = "Results/QC_Normalized_Data",
                    force = TRUE)
# OUTPUT: HTML report post-normalization

# -------------------------------
# 9. Filter Low-Intensity Probes
# -------------------------------
processed_data <- exprs(eset_norm)
row_median <- apply(processed_data, 1, median)

hist(row_median, breaks = 100, freq = FALSE, main = "Median Probe Intensity")
threshold <- 120  # chosen based on histogram inspection
abline(v = threshold, col = "blue", lwd = 2)

filtered_data <- processed_data[row_median > threshold, ]
cat("Probes remaining after filtering:", nrow(filtered_data), "\n")
summary(row_median) # check the distribution of probe medians

#Probes (features) present before filtering" -> nrow(exprs(eset_norm))
#Transcripts remained after filtering low-intensity probes -> nrow(filtered_data)

# -------------------------------
# 10. Define Experimental Groups
# -------------------------------
# Identify which metadata column contains disease info
groups <- factor(phenotype_data$Characteristics.disease.,
                 levels = c("Normal mucosa adjacent", "colorectal cancer"),
                 labels = c("normal", "colorectal cancer"))

# Assign groups to ExpressionSet
pData(eset_norm)$group <- groups
table(pData(eset_norm)$group)
# OUTPUT: counts of normal vs cancer samples

# -------------------------------
# 11. Summary for Assignment
# -------------------------------
nrow(exprs(eset_norm))          # Probes before filtering
nrow(filtered_data)             # Probes after filtering
length(colnames(eset_norm))     # Total samples
table(pData(eset_norm)$group)   # Sample counts by group

# -------------------------------
# 12.Saving the ExpressionSet after normalization
# -------------------------------
save(eset_norm, file = "eset_norm.RData")
# -------------------------------
# 14.Saving the entire R workspace
# -------------------------------
save.image(file = "BENSALEK_Rihab_Class3B_Assignment4.RData")

# -------------------------------
# 15.Loading Workspace
# -------------------------------
load("BENSALEK_Rihab_Class3B_Assignment4.RData")



























# -------------------------------
# 0. Install and Load Packages
# -------------------------------
BiocManager::install(c("lumi", "arrayQualityMetrics", "illuminaio"), ask = FALSE)
BiocManager::install("limma", force = TRUE)

install.packages("dplyr")
install.packages("ggplot2")

library(ArrayExpress)
library(illuminaio)
library(Biobase)              
library(lumi)                 
library(arrayQualityMetrics)  
library(dplyr)                
library(limma)                
library(ggplot2)

# -------------------------------
# 1. Set path to manually downloaded data
# -------------------------------
# You should manually download the raw data files (.idat or processed text files)
# from ArrayExpress website: https://www.ebi.ac.uk/arrayexpress/
# Save all files under a folder, e.g., "C:/Users/hp/Documents/AI_Omics_Internship_2025/MTAB8148/RawFiles"

data_dir <- "C:/Users/hp/Documents/AI_Omics_Internship_2025/Raw_data"

# -------------------------------
# 2. Load raw Illumina data from local folder
# -------------------------------
# If you have IDAT files:
idat_files <- list.files(data_dir, pattern = "*.idat$", full.names = TRUE)
length(idat_files)

# Read all IDAT files (illuminaio::readIDAT)
raw_data_list <- lapply(idat_files, readIDAT)

# Convert to ExpressionSet (manually, if needed):
# Example: If you have one-color Illumina BeadChips
# Extract signal intensity matrix
expr_matrix <- sapply(raw_data_list, function(x) x$Quants$AvgSignal)
rownames(expr_matrix) <- raw_data_list[[1]]$Quants$ProbeID

# Build ExpressionSet
pheno_file <- "C:/Users/hp/Documents/AI_Omics_Internship_2025/MTAB8148/E-MTAB-8148.sdrf.txt"
pheno_data <- read.delim(pheno_file, stringsAsFactors = TRUE)
pheno_data$group <- factor(pheno_data$Source.Name,
                           levels = c("gastric mucosa", "gastric adenocarcinoma"),
                           labels = c("normal", "cancer"))

library(Biobase)
raw_data <- ExpressionSet(assayData = expr_matrix,
                          phenoData = AnnotatedDataFrame(pheno_data))

# -------------------------------
# 3. Preprocessing
# -------------------------------
lumi_data <- lumiExpresso(raw_data)

# -------------------------------
# 4. Quality Control
# -------------------------------
dir.create("Results/QC_Raw_Data", showWarnings = FALSE)
dir.create("Results/QC_Normalized_Data", showWarnings = FALSE)

arrayQualityMetrics(
  expressionset = raw_data,
  outdir = "Results/QC_Raw_Data",
  force = TRUE
)

arrayQualityMetrics(
  expressionset = lumi_data,
  outdir = "Results/QC_Normalized_Data",
  force = TRUE
)

# -------------------------------
# 5. Filtering low-variance probes
# -------------------------------
expr_data <- exprs(lumi_data)
row_median <- apply(expr_data, 1, median)
threshold <- quantile(row_median, 0.25)
filtered_data <- expr_data[row_median > threshold, ]

# -------------------------------
# 6. Probe Annotation (Optional)
# -------------------------------
# Make sure illuminaHumanv3.db is installed
BiocManager::install("illuminaHumanv3.db")
library(illuminaHumanv3.db)

probe_ids <- rownames(filtered_data)

gene_symbols <- mapIds(
  illuminaHumanv3.db,
  keys = probe_ids,
  keytype = "PROBEID",
  column = "SYMBOL",
  multiVals = "first"
)

# -------------------------------
# 7. Save
# -------------------------------
save(lumi_data, file = "lumi_data.RData")
save(raw_data, file = "raw_data.RData")












dest <- "C:/Users/hp/Documents/AI_Omics_Internship_2025/MTAB8148"

#Download and Load Data Directly from ArrayExpress
options(ArrayExpress.download.method = "curl")
raw_data <- ArrayExpress("E-MTAB-8148")
raw_data

#Preprocessing Illumina Data
lumi_data <- lumiExpresso(raw_data)

#Quality Control
arrayQualityMetrics(
  expressionset = raw_data,
  outdir = "Results/QC_Raw_Data",
  force = TRUE
)

arrayQualityMetrics(
  expressionset = lumi_data,
  outdir = "Results/QC_Normalized_Data",
  force = TRUE
)

#Filtering Low-Variance Probes
expr_data <- exprs(lumi_data)
row_median <- apply(expr_data, 1, median)
threshold <- quantile(row_median, 0.25)  # Keep upper 75%
filtered_data <- expr_data[row_median > threshold, ]

#Add Phenotype Data
#ArrayExpress data often includes SDRF metadata automatically:
phenotype_data <- pData(lumi_data)

phenotype_data <- read.delim("E-MTAB-8148.sdrf.txt")

groups <- factor(phenotype_data$Source.Name,
                 levels = c("gastric mucosa", "gastric adenocarcinoma"),
                 labels = c("normal", "cancer"))

pData(lumi_data)$group <- groups

#Probe Annotation
probe_ids <- rownames(filtered_data)

gene_symbols <- mapIds(
  illuminaHumanv3.db,
  keys = probe_ids,
  keytype = "PROBEID",
  column = "SYMBOL",
  multiVals = "first"
)

#Save
save(lumi_data, file = "lumi_data.RData")
