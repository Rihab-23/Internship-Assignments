# ===================================================================
# -------------------------------------------------------------------
#                       Module I: Class 1C
# -------------------------------------------------------------------
# ===================================================================
# Retrieve the value in the console
gene <- "TP53"
gene 
print(gene)
# Store numeric values in one variable (vector)
expression_levels <- c(2.3, 4.6, 3.6, 7.2, 4.7)
expression_levels

# Import a CSV file as a variable
raw_data <- read.csv("C:/Users/hp/Documents/AI_Omics_Internship_2025/patient_info.csv")
raw_data
View(raw_data)

# Create a copy of raw_data
data <- raw_data
data
View(data)

# Remove columns
# Remove patient_id column
raw_data$patient_id <- NULL

# Create a new dataset without the first column
clean_data <- data[ ,-1]
clean_data
View(clean_data)

# Sort " age " from largest to smallest
sorted_age <- sort(raw_data$age, decreasing = TRUE)
sorted_age

# Sort 'age' from smallest to largest
sorted_age1 <- sort(raw_data$age, decreasing = FALSE)
sorted_age1

# Practice exercises
# Logical conditions: if & else
# 1 Check Cholesterol level
cholesterol <- 230
if (cholesterol < 190) {
  print("Low Cholesterol")
}

# 2 Blood pressure status
Systolic_bp <- 130
if (Systolic_bp < 120) {
  print("Normal systolic blood pressure")
} else {
  print("Higher systolic blood pressure")
}

# 3 Convert character columns to factors
# Loop For
clean_data <- raw_data
str(clean_data)
for (i in 1:ncol(clean_data)) {
  if (is.character(clean_data[[i]])) {
    clean_data[[i]] <- as.factor(clean_data[[i]])
  }
}  

# 4 Convert Factors to Numeric Codes
raw_data$smoker_numeric <- numeric(nrow(raw_data))

for (i in 1:nrow(raw_data)) {
  if (raw_data$smoker[i] == "Yes") {
    raw_data$smoker_numeric[i] <- 1
  } else if (raw_data$smoker[i] == "No") {
    raw_data$smoker_numeric[i] <- 0
  }
}

# Check structure of raw_data
str(data)
clean_data <- raw_data
clean_data
View(clean_data)
str(clean_data)
