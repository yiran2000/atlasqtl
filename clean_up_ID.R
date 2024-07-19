rm(list=ls())
library(dplyr)
library(data.table)
library(tictoc)
library(peakRAM)
library(missMethods)
library(mice)
library(ggplot2)
library(corrplot)
library(cluster)


geno_path = "/rds/user/yl2021/hpc-work/myukbb"
pheno_path = "~/rds/rds-mrc-bsu-csoP2nj6Y6Y/biobank/proteomics"
save_path = "/rds/user/yl2021/hpc-work/myukbb/imputed_Y"

task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

#-------------------------------------------------------------------------------
#load data
#sample id
# sample_ls = fread(file.path(geno_path, "/sample_qc_pass_id.txt"))[[2]]

#proteomics
olink = fread(file.path(pheno_path, "/olink_data.txt.gz"))
olink_wide = dcast(olink, eid+ ins_index ~ protein_id, value.var = "result")
pheno = olink_wide[ins_index==0]

#select samples after QC
sample_qc_pass_id = fread(file.path(geno_path, "sample_qc_pass_id.txt")) %>% rename(FID = V1, IID = V2)
pheno = pheno %>% filter(eid %in% sample_qc_pass_id$IID)

rm(olink)
rm(olink_wide)

#sample sequence
sample = fread(file.path(geno_path, "/imputed_chr1_ld1.sample")) 
sample = sample[-1, ] %>% rename(FID = ID_1, IID = ID_2, SEX = sex)

pheno = pheno[match(sample$IID, pheno[["eid"]]), ] #now the sample sequence of sample and pheno are exactly the same

sum(sample$IID == pheno[["eid"]])


well_plate = fread(file.path(pheno_path, "ukb678244.csv")) %>% 
  filter(eid %in% pheno$eid)%>% 
  select("eid", "30901-0.0", "30902-0.0") %>% 
  rename(PlateID = "30901-0.0", WellID = "30902-0.0") #we have 4 people doesn't have corresponding info in here

pheno = pheno %>% filter(eid %in% well_plate$eid)

#--------------------------------------------------------------------------------
#remove rows with high level of missingness
mr_Y2 <- rowMeans(is.na(pheno))
selected_row <- which(mr_Y2 < 0.4)

#remove columns of missingness > 0.2
mr_Y <- colMeans(is.na(pheno[selected_row,,drop = FALSE]))
selected_col <- names(mr_Y)[mr_Y < 0.2]

Y <- pheno[selected_row,..selected_col, drop = FALSE]

ID = Y$eid
Y = Y[,3:ncol(Y)]
rownames(Y) = ID
colnames(Y) = paste0("protein_", colnames(Y))# 36636 samples remained

saveRDS(Y, file.path(geno_path, "Y_original.rds"))
saveRDS(ID, file.path(geno_path, "sample_ID.rds"))