#Task_BENSALEK Rihab#
#SETTING WORKING DIRECTORY#
setwd("C:/Users/hp/Documents/AI_Omics_Internship_2025")
getwd()
                          
#Creating Project Folder#
##Creating Project "Module_I##
dir.create("C:/Users/hp/Documents/AI_Omics_Internship_2025/Module_I")

#Creating subfolders#
dir.create("C:/Users/hp/Documents/AI_Omics_Internship_2025/Module_I/raw_data")
dir.create("C:/Users/hp/Documents/AI_Omics_Internship_2025/Module_I/clean_data")
dir.create("C:/Users/hp/Documents/AI_Omics_Internship_2025/Module_I/scripts")
dir.create("C:/Users/hp/Documents/AI_Omics_Internship_2025/Module_I/results")
dir.create("C:/Users/hp/Documents/AI_Omics_Internship_2025/Module_I/plots")

#Loading "Patient_info.csv" dataset into my R environment#
data <- read.csv("C:/Users/hp/Documents/AI_Omics_Internship_2025/patient_info.csv")
data
View(data)

#Inspecting the structure of "Patient_info.csv" dataset#
str(data)
summary(data)

#Inspecting rows and columns#
nrow(data)
ncol(data)
colnames(data)

#Identifying variables with incorrect or inconsistent data types#
sapply(data, class)

#Converting variables into appropriate data types#
data$patient_id <- as.character(data$patient_id)
data$age <- as.integer(data$age)
data$gender<- as.factor(data$gender)
data$diagnosis<- as.factor(data$diagnosis)
data$bmi <- as.numeric(data$bmi)
data$smoker<- as.factor(data$smoker)
str(data)

#Creating new variables#
##Creating Binary variable##
data$smoker_status <- ifelse(data$smoker == "Yes", 1, 0)
class(data$smoker_status)

##Converting the binary variable into a factor##
data$smoker_status <- factor(data$smoker_status, levels = c(0, 1))
class(data$smoker_status)

#Saving the clean dataset#
write.csv(data, file = "C:/Users/hp/Documents/AI_Omics_Internship_2025/Module_I/clean_data/patient_info_clean.csv")

#Saving the entire R workspace"
save.image(file = "BENSALEK_Rihab_ClassIb_Assignment.RData")

#Loading Workspace#
load("BENSALEK_Rihab_ClassIb_Assignment.RData")
