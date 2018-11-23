install.packages("gdata")
library("gdata")

# Readin in file 
nat_clin <- read.xls("/Users/joshua_harris/Dropbox/Research/PhD/Bioinformatics/Results/TCGA/FireBrowse_TCGA_BRCA/FB_BRCA_TCGA/DATA/Clinical/Nature_TCGA_Clin.xls")
View(nat_clin)

colnames(nat_clin) <- gsub("\\.", "_" ,colnames(nat_clin))
colnames(nat_clin) <- gsub("\\__", "_" ,colnames(nat_clin))
colnames(nat_clin) <- gsub("\\_$", "" ,colnames(nat_clin))

# Saving file as CSV 
write.csv(nat_clin, file = "DATA/Clinical/Nature_TCGA_Clin.csv", row.names = F)

# Checking Written file

nat_imp <- read.csv(file = "DATA/Clinical/Nature_TCGA_Clin.csv", header = T)
View(nat_imp)
