# ---------------------------------------------------------
#                    Task_BENSALEK Rihab
# ---------------------------------------------------------
# Assignment 2: Differential Gene Expression Classification
# ---------------------------------------------------------

# SETTING WORKING DIRECTORY#
setwd("C:/Users/hp/Documents/AI_Omics_Internship_2025")
getwd()

# 1. Define input and output directories
input_dir <- "Raw_Data"        # folder containing input CSV files
output_dir <- "Results"        # folder where processed files will be saved

input_dir
output_dir

# 2. Create Raw_Data & Results folder if don't exist
if(!dir.exists(input_dir)){
  dir.create(input_dir)
}
if(!dir.exists(output_dir)){
  dir.create(output_dir)
}

# 3. List of files to process
files_to_process <- c("DEGs_Data_1.csv", "DEGs_Data_2.csv")  # input files
result_list <- list()   # empty list to store processed datasets

files_to_process
result_list

# 4. Function to classify genes based on logFC and padj
classify_gene <- function(logFC, padj){
  if (logFC > 1 & padj < 0.05) {        
    return("Upregulated")               # Condition 1 "Upregulated"
  } else if (logFC < -1 & padj < 0.05) {
    return("Downregulated")             # Condition 2 "Downregulated"
  } else {
    return("NA_ Not Significant")       # Otherwise 
  }
}

classify_gene

# 5. Process each file in a loop
for(files in files_to_process){
  cat("\nProcessing :", files , "\n")   # print which file is being processed
  
  # Build the path to the csv file and read it
  file_path <- file.path(input_dir, files)
  data <- read.csv(file_path, header = TRUE)
  cat("Files imported, Checking for missing values.....\n")
  
  # 6. Handle missing values in padj column
  if("padj" %in% names(data)){
    missing_count <- sum(is.na(data$padj))              # count missing padj
    
    # 0 < p_value < 1 padj ranges from 0 to 1, small values (close to 0) are significant results 
    #while large values (close to 1) are not significant results 
    cat("Missing values in 'padj' column:", missing_count , "\n")
    data$padj[is.na(data$padj)] <- 1                    # replace NA (missing padj) with 1 
  }
  
  # 7. Handle missing values in logFC column
  if("logFC" %in% names(data)){
    missing_values <- sum(is.na(data$logFC))            # count missing logFC
    
    cat("Missing values in 'logFC':", missing_values , "\n")    # Replace missing logFC with column mean
    data$logFC[is.na(data$logFC)] <- mean(data$logFC, na.rm = TRUE)   # na.rm = TRUE -> ignores NA
  }                     
  
  # 8. Add new column 'status' using classify_gene function
  # Apply function to each row using mapply (maps two vectors: logFC and padj)
  data$status <- mapply(classify_gene, data$logFC, data$padj)
  cat("Status column added successfully!\n")
  
  # 9. Save processed file into Results folder
  output_file_path <- file.path(output_dir, paste0("Classification_", files))
  write.csv(data, output_file_path, row.names = FALSE)
  cat("Results saved to:", output_file_path, "\n")
  
  output_file_path
  
  # 10. Print summary counts
  gene_counts <- table(data$status)
  cat("Summary counts for", files, ":\n")
  print(gene_counts)
  
  # Store result in list for later use
  result_list[[files]] <- data
  
  # 11. Print summary counts of classifications
  gene_counts <- table(data$status)
  cat("Summary counts for", files, ":\n")
  print(gene_counts)
}

# 12. Store each processed dataset separately from the list
result_1 <- result_list[[1]]
result_2 <- result_list[[2]]

result_1
result_2

# 13. Save R workspace for reproducibility
save.image(file = "BENSALEK_Rihab_Class2_Assignment.RData")
