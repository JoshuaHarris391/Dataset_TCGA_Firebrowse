
########################################################################################
#### RSEM TPM Pre-processing 
########################################################################################

# reading in file
RSEM.TPM <- data.table::fread(file = "DATA/RSEM_Normalized/RSEM_Normalized.txt", header = T, sep = "\t", na.strings = c("", NA))
# Removing first Row 
RSEM.TPM <- RSEM.TPM[2:length(rownames(RSEM.TPM)), ]
# Splitting Column 1
RSEM.TPM <- tidyr::separate(RSEM.TPM, colnames(RSEM.TPM)[1],c("Gene", "Entrez_ID"), sep = "\\|")
# Excluding Rows with no gene annotation 
RSEM.TPM <- RSEM.TPM[!(RSEM.TPM$Gene == "?"), ]
# Converting Entrez_ID to rownames 
rownames(RSEM.TPM) <- RSEM.TPM$Entrez_ID

### Removing normal solid tissue samples
# Removing last 12 characters from IDs (. == any character), ({12} == 12 times),  ($ == end of the string)
colnames(RSEM.TPM)[3:length(colnames(RSEM.TPM))] <- gsub(".{13}$", "", colnames(RSEM.TPM[, 3:length(colnames(RSEM.TPM))]))
# Finding Solid Normal tissue samples ending with "11" 
RSEM.cols <- colnames(RSEM.TPM)
norm_sam <- RSEM.cols[grepl("11$", RSEM.cols)]  
# Removing normal tissue samples 
keep_cols <- !(colnames(RSEM.TPM) %in% norm_sam)
RSEM.TPM <- subset(RSEM.TPM, select = keep_cols)
# Calculating percentage of Normal Solid tumours 
nst <- length(norm_sam)*100/(length(RSEM.cols)-2)
print(paste(format(nst, digits = 5), "%", sep = " "))  # 9.2% of the patient IDs are normal solid tissues

### Removing Metastases "06"
# Finding Solid Normal tissue samples ending with "06" 
met_sam <- RSEM.cols[grepl("06$", RSEM.cols)] 
# Removing samples from metastases 
keep_cols <- !(colnames(RSEM.TPM) %in% met_sam)
RSEM.TPM <- subset(RSEM.TPM, select = keep_cols)
# Caluclating number of metastatic specimens 
ms <- length(met_sam)*100/length(RSEM.cols)
print(paste(format(ms, digits = 2), "%", sep = " "))  # 0.58 % of the patient IDs are metastatic tissues

# Checking that all remaning samples are primary solid tumours only ending in "01"
table(grepl("01", colnames(RSEM.TPM[, 3:ncol(RSEM.TPM)])))

# Converting columns to numeric
RSEM.TPM.nmat <- as.data.frame(apply(RSEM.TPM[, 3:length(colnames(RSEM.TPM))], 2, as.numeric))
RSEM.TPM <- cbind(RSEM.TPM[, 1:2], RSEM.TPM.nmat)

# removing tissue identifier 
colnames(RSEM.TPM)[3:length(colnames(RSEM.TPM))] <- gsub(".{3}$", "", colnames(RSEM.TPM[, 3:length(colnames(RSEM.TPM))]))



########################################################################################
#### EXP.COUNTS Pre-processing 
########################################################################################

# Reading in file 
EXP.COUNTS <- data.table::fread(file = "DATA/Merge_Expression/Merge_Expression.txt", header = T, na.strings = c("", NA))
# Creating gene.entrez lookup 
gene.entrez <- as.data.frame(EXP.COUNTS[2:length(rownames(EXP.COUNTS)), 1])
colnames(gene.entrez) <- "gene.entrez"
# Splitting DF 
gene.entrez <- tidyr::separate(gene.entrez, gene.entrez,c("Gene", "Entrez_ID"), sep = "\\|")
# Subsetting Raw EXP.COUNTS 
EXP.COUNTS <- as.data.frame(t(EXP.COUNTS))
EXP.COUNTS <- EXP.COUNTS[(EXP.COUNTS[, 1] == "raw_counts"), ]
# Transposing and subsetting raw EXP.COUNTS 
EXP.COUNTS <- t(EXP.COUNTS[, 2:length(colnames(EXP.COUNTS))])
# adding back gene annotation 
EXP.COUNTS <- cbind(gene.entrez, EXP.COUNTS)
# Converting rownames to entrez IDs 
rownames(EXP.COUNTS) <- EXP.COUNTS$Entrez_ID
# removing unannotated genes 
EXP.COUNTS <- EXP.COUNTS[!(EXP.COUNTS$Gene == "?"), ]

### Removing normal solid tissue samples
# Removing last 12 characters from IDs (. == any character), ({12} == 12 times),  ($ == end of the string)
colnames(EXP.COUNTS)[3:length(colnames(EXP.COUNTS))] <- gsub(".{13}$", "", colnames(EXP.COUNTS[, 3:length(colnames(EXP.COUNTS))]))
# Finding Solid Normal tissue samples ending with "11" 
COUNTS.cols <- colnames(EXP.COUNTS)
norm_sam <- COUNTS.cols[grepl("11$", COUNTS.cols)]  
# Removing normal tissue samples 
keep_cols <- !(colnames(EXP.COUNTS) %in% norm_sam)
EXP.COUNTS <- subset(EXP.COUNTS, select = keep_cols)

### Removing Metastases "06"
# Finding Solid Normal tissue samples ending with "06" 
met_sam <- COUNTS.cols[grepl("06$", COUNTS.cols)] 
# Removing samples from metastases 
keep_cols <- !(colnames(EXP.COUNTS) %in% met_sam)
EXP.COUNTS <- subset(EXP.COUNTS, select = keep_cols)

# Checking that all remaning samples are primary solid tumours only ending in "01"
table(grepl("01", colnames(EXP.COUNTS[, 3:ncol(EXP.COUNTS)])))

# Converting columns to numeric
Count.nmat <- as.data.frame(apply(EXP.COUNTS[, 3:length(colnames(EXP.COUNTS))], 2, as.numeric))
EXP.COUNTS <- cbind(EXP.COUNTS[, 1:2], Count.nmat)

# Removing tissue identifier
colnames(EXP.COUNTS)[3:length(colnames(EXP.COUNTS))] <- gsub(".{3}$", "", colnames(EXP.COUNTS[, 3:length(colnames(EXP.COUNTS))]))


########################################################################################
#### cBioPortal TCGA Clinical Pre-Processing 
########################################################################################
# Reading in Data 
CLINICAL <- read.delim(file = "DATA/Clin_cBio.txt", na.strings = c("", NA))
# Converting "." to "-"
rownames(CLINICAL) <- gsub("\\.", "-", rownames(CLINICAL))
# Adding Patient IDs as a Variable 
CLINICAL$PATIENT_ID <- rownames(CLINICAL)


### Removing normal solid tissue samples
# Finding Solid Normal tissue samples ending with "11" 
clin_sam <- rownames(CLINICAL)
# Checking to see if there are any samples that are normal tissue
table(grepl("11$", clin_sam)) # There are no normal tissue samples


### Removing Metastases "06"
# Finding Solid Normal tissue samples ending with "06" 
met_sam <- clin_sam[grepl("06$", clin_sam)] 
# Removing samples from metastases 
keep_rows <- !(rownames(CLINICAL) %in% met_sam)
CLINICAL <- subset(CLINICAL, keep_rows)
# Rearranging Columns 
col_len <- length(colnames(CLINICAL))
CLINICAL <- CLINICAL[, c(col_len, 1:(col_len-1))]
# removing tissue type identifier 
rownames(CLINICAL) <- gsub(".{3}$", "", rownames(CLINICAL))
# Updating Patient ID column 
CLINICAL$PATIENT_ID <- rownames(CLINICAL)


########################################################################################
#### Nature paper TCGA Clinical Pre-Processing 
########################################################################################
# Reading in data 
NATURE_CLIN <- read.csv(file = "DATA/Clinical/Nature_TCGA_Clin.csv", header = T, na.strings = c("", NA))

# Changing ER receptor status to binary (0 == neg, 1 == pos)
table(NATURE_CLIN$ER_Status)
NATURE_CLIN$ER_Status <- ifelse(NATURE_CLIN$ER_Status == "Positive", "Positive", 
                                ifelse(NATURE_CLIN$ER_Status == "Negative", "Negative", NA))
NATURE_CLIN$ER_Status <- as.character(NATURE_CLIN$ER_Status)
NATURE_CLIN$ER_Status <- factor(NATURE_CLIN$ER_Status, c("Negative", "Positive"))

# Changing PR receptor status to binary 
table(NATURE_CLIN$PR_Status)
NATURE_CLIN$PR_Status <- ifelse(NATURE_CLIN$PR_Status == "Positive", "Positive", 
                                ifelse(NATURE_CLIN$PR_Status == "Negative", "Negative", NA))
NATURE_CLIN$PR_Status <- as.character(NATURE_CLIN$PR_Status)
NATURE_CLIN$PR_Status <- factor(NATURE_CLIN$PR_Status, c("Negative", "Positive"))


# Changing HER2 receptor status to discrete (note that equivocal will be consider NA)
table(NATURE_CLIN$HER2_Final_Status)
NATURE_CLIN$HER2_Final_Status <- ifelse(NATURE_CLIN$HER2_Final_Status == "Positive", "Positive", 
                                ifelse(NATURE_CLIN$HER2_Final_Status == "Negative", "Negative", NA))
NATURE_CLIN$HER2_Final_Status <- as.character(NATURE_CLIN$HER2_Final_Status)
NATURE_CLIN$HER2_Final_Status <- factor(NATURE_CLIN$HER2_Final_Status, c("Negative", "Positive"))


# Creating "Node_Coded" to discrete
table(NATURE_CLIN$Node_Coded)
NATURE_CLIN$Node_Coded <- ifelse(NATURE_CLIN$Node_Coded == "Positive", "Positive", "Negative")
NATURE_CLIN$Node_Coded <- as.character(NATURE_CLIN$Node_Coded)
NATURE_CLIN$Node_Coded <- factor(NATURE_CLIN$Node_Coded, c("Negative", "Positive"))

# Creating an integer stage column 
table(NATURE_CLIN$AJCC_Stage)
NATURE_CLIN$Integer_Stage <- ifelse(NATURE_CLIN$AJCC_Stage == "Stage I", 1, 
                                    ifelse(NATURE_CLIN$AJCC_Stage == "Stage IA", 1, 
                                           ifelse(NATURE_CLIN$AJCC_Stage == "Stage IB", 1, 
                                                  ifelse(NATURE_CLIN$AJCC_Stage == "Stage II", 2, 
                                                         ifelse(NATURE_CLIN$AJCC_Stage == "Stage IIA", 2, 
                                                                ifelse(NATURE_CLIN$AJCC_Stage == "Stage IIB", 2, 
                                                                       ifelse(NATURE_CLIN$AJCC_Stage == "Stage III", 3, 
                                                                              ifelse(NATURE_CLIN$AJCC_Stage == "Stage IIIA", 3, 
                                                                                     ifelse(NATURE_CLIN$AJCC_Stage == "Stage IIIB", 3, 
                                                                                            ifelse(NATURE_CLIN$AJCC_Stage == "Stage IIIC", 3, 
                                                                                                   ifelse(NATURE_CLIN$AJCC_Stage == "Stage IV", 4, NA)))))))))))
NATURE_CLIN$Integer_Stage <- as.integer(NATURE_CLIN$Integer_Stage)
NATURE_CLIN$Integer_Stage <- factor(NATURE_CLIN$Integer_Stage, c(1, 2, 3, 4))

# Converting OS time from days to years 
NATURE_CLIN$OS_Time <- (NATURE_CLIN$OS_Time)/365 
NATURE_CLIN$OS_Time <- as.numeric(NATURE_CLIN$OS_Time)


########################################################################################
#### Removing unlinked cases between RSEM.TPM and "CLINICAL"
########################################################################################

## Since we know that RSEM.TPM has more values than CLINICAL df, I will subtract unlinked cases from RSEM.TPM 
# Creating lookup vectors for RSEM and CLINICAL DFs 
RSEM.lookvec <- colnames(RSEM.TPM)
RSEM.Clin.Lookvec <- rownames(CLINICAL)
# Creating Match lookup 
RSEM_Match <- match(RSEM.Clin.Lookvec, RSEM.lookvec)
# Subsetting by matched cases
RSEM_DF <- t(RSEM.TPM)
RSEM_DF <- RSEM_DF[c(RSEM_Match), ]
RSEM_DF <- t(RSEM_DF)
RSEM_DF <- as.data.frame(RSEM_DF)
# Converting Factors to Numeric
RSEM_DF <- data.frame(lapply(RSEM_DF, as.character), stringsAsFactors=FALSE)
RSEM_DF <- data.frame(lapply(RSEM_DF, as.numeric), stringsAsFactors=FALSE)
# Converting "." to "-"
colnames(RSEM_DF) <- gsub("\\.", "-", colnames(RSEM_DF))
# Combining TPM values with Gene and Entrez_ID columns 
RSEM_DF <- cbind(RSEM.TPM[, 1:2], RSEM_DF, stringsAsFactors = FALSE)
# Checking DF Properties
# View(RSEM_DF[1:10, 1:10])
# View(CLINICAL[1:10, 1:10])
# 
# str(RSEM_DF[1:10, 1:10])
# dim(RSEM_DF)


########################################################################################
#### Removing unlinked cases between "EXP.COUNTS" and "CLINICAL"
########################################################################################

## Since we know that clinical.df has more values than the EXP.COUNTS df, I will subtract unlinked cases from clinical.df and call it "CLIN_COUNTS" 
# Creating lookup vectors for RSEM and CLINICAL DFs 
Counts.lookvec <- colnames(EXP.COUNTS)
Counts.Clin.Lookvec <- rownames(CLINICAL)
# Creating Match lookup 
Counts_Match <-  match(Counts.lookvec, Counts.Clin.Lookvec)
# Removing NA values from match vector
Counts_Match <- na.omit(Counts_Match)
# Subsetting by matched cases
CLIN_COUNTS <- CLINICAL[c(Counts_Match), ]
# Checking DF Properties
# View(CLIN_COUNTS[1:10, 1:10])
# View(EXP.COUNTS[1:10, 1:10])
# 
# str(CLIN_COUNTS[1:10, 1:10])
# dim(CLIN_COUNTS)


########################################################################################
#### Removing unlinked cases between RSEM.TPM and "NATURE_CLIN"
########################################################################################

## Since we know that RSEM.TPM has more values than NATURE_CLIN df, I will subtract unlinked cases from RSEM.TPM 
# Creating lookup vectors for RSEM and NATURE_CLIN DFs 
RSEM.lookvec <- colnames(RSEM.TPM)
RSEM.Clin.Lookvec <- NATURE_CLIN$Complete_TCGA_ID
# Creating Match lookup 
RSEM_Match <- match(RSEM.Clin.Lookvec, RSEM.lookvec)
# Subsetting by matched cases
RSEM_DF_NATURE <- t(RSEM.TPM)
RSEM_DF_NATURE <- RSEM_DF_NATURE[c(RSEM_Match), ]
RSEM_DF_NATURE <- t(RSEM_DF_NATURE)
RSEM_DF_NATURE <- as.data.frame(RSEM_DF_NATURE)
# Converting Factors to Numeric
RSEM_DF_NATURE <- data.frame(lapply(RSEM_DF_NATURE, as.character), stringsAsFactors=FALSE)
RSEM_DF_NATURE <- data.frame(lapply(RSEM_DF_NATURE, as.numeric), stringsAsFactors=FALSE)
# Converting "." to "-"
colnames(RSEM_DF_NATURE) <- gsub("\\.", "-", colnames(RSEM_DF_NATURE))
# Combining TPM values with Gene and Entrez_ID columns 
RSEM_DF_NATURE <- cbind(RSEM.TPM[, 1:2], RSEM_DF_NATURE, stringsAsFactors = FALSE)
# # Checking DF Properties
# View(RSEM_DF_NATURE[1:10, 1:10])
# View(NATURE_CLIN[1:10, 1:10])
# 
# str(RSEM_DF_NATURE)
# dim(RSEM_DF_NATURE)


########################################################################################
#### Removing unlinked cases between "EXP.COUNTS" and "NATURE_CLIN"
########################################################################################

## Since we know that NATURE_CLIN.df has more values than the EXP.COUNTS df, I will subtract unlinked cases from NATURE_CLIN.df and call it "CLIN_COUNTS" 
# Creating lookup vectors for RSEM and NATURE_CLIN DFs 
Counts.lookvec <- colnames(EXP.COUNTS)
Counts.Clin.Lookvec <- NATURE_CLIN$Complete_TCGA_ID
# Creating Match lookup 
Counts_Match <-  match(Counts.lookvec, Counts.Clin.Lookvec)
# Removing NA values from match vector
Counts_Match <- na.omit(Counts_Match)
# Subsetting by matched cases
CLIN_COUNTS_NATURE <- NATURE_CLIN[c(Counts_Match), ]
# # Checking DF Properties
# View(CLIN_COUNTS_NATURE[1:10, 1:10])
# View(EXP.COUNTS[1:10, 1:10])
# 
# str(CLIN_COUNTS_NATURE)
# dim(CLIN_COUNTS_NATURE)



########################################################################################
#### Creating Post Processed Clinical DF
########################################################################################
clin_factors <- c(
  "PATIENT_ID",
  "AJCC_NODES_PATHOLOGIC_PN",
  "AJCC_PATHOLOGIC_TUMOR_STAGE", 
  "DFS_MONTHS", 
  "DFS_STATUS", 
  "ER_STATUS_BY_IHC", 
  "HER2_IHC_SCORE",
  "IHC_HER2", 
  "OS_MONTHS", 
  "OS_STATUS", 
  "PR_STATUS_BY_IHC"
)
CLIN_POST <- subset(CLINICAL, select = clin_factors)


# Changing Node Status to integer 
# Converting Node Negative == 0, Node positive == 1, Node NX == NA
CLIN_POST$INTEGER_NODE_STATUS <- as.integer(ifelse(CLIN_POST$AJCC_NODES_PATHOLOGIC_PN == "NX", NA, 
                                                        ifelse(grepl("^N0", CLIN_POST$AJCC_NODES_PATHOLOGIC_PN, perl=TRUE), 0, 1)))
CLIN_POST$INTEGER_NODE_STATUS <- factor(CLIN_POST$INTEGER_NODE_STATUS, c(0,1))

# Changing Stage to Integers
CLIN_POST$INTEGER_STAGE <- ifelse(CLIN_POST$AJCC_PATHOLOGIC_TUMOR_STAGE == "Stage I", 1, 
                                              ifelse(CLIN_POST$AJCC_PATHOLOGIC_TUMOR_STAGE == "Stage IA", 1, 
                                                     ifelse(CLIN_POST$AJCC_PATHOLOGIC_TUMOR_STAGE == "Stage IB", 1, 
                                                            ifelse(CLIN_POST$AJCC_PATHOLOGIC_TUMOR_STAGE == "Stage II", 2, 
                                                                   ifelse(CLIN_POST$AJCC_PATHOLOGIC_TUMOR_STAGE == "Stage IIA", 2, 
                                                                          ifelse(CLIN_POST$AJCC_PATHOLOGIC_TUMOR_STAGE == "Stage IIB", 2, 
                                                                                 ifelse(CLIN_POST$AJCC_PATHOLOGIC_TUMOR_STAGE == "Stage III", 3, 
                                                                                        ifelse(CLIN_POST$AJCC_PATHOLOGIC_TUMOR_STAGE == "Stage IIIA", 3, 
                                                                                               ifelse(CLIN_POST$AJCC_PATHOLOGIC_TUMOR_STAGE == "Stage IIIB", 3, 
                                                                                                      ifelse(CLIN_POST$AJCC_PATHOLOGIC_TUMOR_STAGE == "Stage IIIC", 3, 
                                                                                                             ifelse(CLIN_POST$AJCC_PATHOLOGIC_TUMOR_STAGE == "Stage IV", 4, NA)))))))))))

CLIN_POST$INTEGER_STAGE <- as.integer(CLIN_POST$INTEGER_STAGE)
CLIN_POST$INTEGER_STAGE <- factor(CLIN_POST$INTEGER_STAGE, c(1,2,3,4))

# Converting HER2_IHC Score to factored characters 
CLIN_POST$HER2_IHC_SCORE <- as.integer(CLIN_POST$HER2_IHC_SCORE)
CLIN_POST$HER2_IHC_SCORE <- factor(CLIN_POST$HER2_IHC_SCORE, c("0", "1", "2", "3"))

# Converting Survival status to binary 
CLIN_POST$DFS_STATUS <- ifelse(CLIN_POST$DFS_STATUS == "DiseaseFree", 0, 1)
CLIN_POST$DFS_STATUS <- as.integer(CLIN_POST$DFS_STATUS)

CLIN_POST$OS_STATUS <- ifelse(CLIN_POST$OS_STATUS == "LIVING", 0, 1)
CLIN_POST$OS_STATUS <- as.integer(CLIN_POST$OS_STATUS)


# Converting ER Status to binary 
table(CLIN_POST$ER_STATUS_BY_IHC)
CLIN_POST$ER_STATUS_BY_IHC <- ifelse(CLIN_POST$ER_STATUS_BY_IHC == "Positive", "Positive", 
                                ifelse(CLIN_POST$ER_STATUS_BY_IHC == "Negative", "Negative", NA))
CLIN_POST$ER_STATUS_BY_IHC <- as.character(CLIN_POST$ER_STATUS_BY_IHC)
CLIN_POST$ER_STATUS_BY_IHC <- factor(CLIN_POST$ER_STATUS_BY_IHC, c("Negative", "Positive"))

# Converting PR Status to binary 
table(CLIN_POST$PR_STATUS_BY_IHC)
CLIN_POST$PR_STATUS_BY_IHC <- ifelse(CLIN_POST$PR_STATUS_BY_IHC == "Positive", "Positive", 
                                     ifelse(CLIN_POST$PR_STATUS_BY_IHC == "Negative", "Negative", NA))
CLIN_POST$PR_STATUS_BY_IHC <- as.character(CLIN_POST$PR_STATUS_BY_IHC)
CLIN_POST$PR_STATUS_BY_IHC <- factor(CLIN_POST$PR_STATUS_BY_IHC, c("Negative", "Positive"))

# Converting HER2 Status to binary (Equivocol values are regarded as NA)
table(CLIN_POST$IHC_HER2)
CLIN_POST$IHC_HER2 <- ifelse(CLIN_POST$IHC_HER2 == "Positive", "Positive", 
                                     ifelse(CLIN_POST$IHC_HER2 == "Negative", "Negative", NA))
CLIN_POST$IHC_HER2 <- as.character(CLIN_POST$IHC_HER2)
CLIN_POST$IHC_HER2 <- factor(CLIN_POST$IHC_HER2, c("Negative", "Positive"))

########################################################################################
#### Creating R Data file
########################################################################################

# Writing Feather files 
# install.packages("feather")
library(feather)
write_feather(RSEM_DF, path = "DATA/Processed/RSEM_DF.txt")
write_feather(EXP.COUNTS, path = "DATA/Processed/EXP_COUNTS.txt")
write_feather(CLINICAL, path  = "DATA/Processed/CLINICAL.txt")
write_feather(CLIN_COUNTS, path  = "DATA/Processed/CLIN_COUNTS.txt")
write_feather(CLIN_POST, path  = "DATA/Processed/CLIN_POST.txt")
write_feather(NATURE_CLIN, path  = "DATA/Processed/NATURE_CLIN.txt")
write_feather(RSEM_DF_NATURE, path  = "DATA/Processed/NATURE_RSEM_DF.txt")
write_feather(CLIN_COUNTS_NATURE, path  = "DATA/Processed/NATURE_CLIN_COUNTS.txt")


# Reading in files 
rm(list=ls())
RSEM_DF <- read_feather(path = "DATA/Processed/RSEM_DF.txt")
EXP_COUNTS <- read_feather(path = "DATA/Processed/EXP_COUNTS.txt")
CLINICAL <- read_feather(path  = "DATA/Processed/CLINICAL.txt")
CLIN_COUNTS <- read_feather(path  = "DATA/Processed/CLIN_COUNTS.txt")
CLIN_POST <- read_feather(path  = "DATA/Processed/CLIN_POST.txt")
NATURE_CLIN <- read_feather(path = "DATA/Processed/NATURE_CLIN.txt")
NATURE_RSEM_DF <- read_feather(path = "DATA/Processed/NATURE_RSEM_DF.txt")
NATURE_CLIN_COUNTS <- read_feather(path = "DATA/Processed/NATURE_CLIN_COUNTS.txt")
# Saving Global Env
save.image(file = "DATA/Processed/FB_BRCA_TCGA.RData")
rm(list=ls())
# Loading Glob Env
load(file = "DATA/Processed/FB_BRCA_TCGA.RData")


# Viewing files 
View(CLIN_COUNTS[1:10, 1:10])
View(CLIN_POST[1:10, 1:10])
View(CLINICAL[1:10, 1:10])
View(EXP_COUNTS[1:10, 1:10])
View(RSEM_DF[1:10, 1:10])
View(NATURE_CLIN[1:10, 1:10])
View(NATURE_RSEM_DF[1:10, 1:10])
View(NATURE_CLIN_COUNTS[1:10, 1:10])










