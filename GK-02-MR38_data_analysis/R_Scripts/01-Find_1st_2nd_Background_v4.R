# Description: Take output from 01-Find_1st_2nd_Background and compile SNP report
# Details: For further analysis we need for every strain a list of the animals Alleles compared to the BG(s)
# Type: Script

options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx8192m"))
options(scipen=999)
library("dplyr")
library("tidyr")
library("ggplot2")
library("janitor")
library("xlsx")
library("stringr")

# Set working directory ####
setwd("S:/AG/AG-Knauf/Daten/AG Knauf/Gabriel Kirchgraber/GK-01-MUGA/")
# setwd("C:/Users/sweet/Dropbox/Arbeit/GK-01-MUGA/")

# Set File Paths ####
output_dir_STR = "./GK-02-MR38_data_analysis/R_Scripts/Output/"
plots_dir_STR = "./GK-02-MR38_data_analysis/R_Scripts/Plots/"
data_Neo_dir_STR = "./_Data/Neogene_Files/"
MR38_dir_STR = "./GK-02-MR38_data_analysis/"
data_dir = paste(MR38_dir_STR, "Find_Background_Data/Single_Animals/", sep = '')

final_report_STR = "Charite_Berlin_MURCOMV02_20191231_FinalReport.txt"
SNP_Map_STR = "SNP_Map_v3.csv"
consensus_STR = "MiniMUGAv2 Consensus Calls.txt"
report_STR = "BC_SNP_analysis.xlsx"

final_report_file_STR = paste(data_Neo_dir_STR, final_report_STR, sep = '')
SNP_map_file_STR = paste(data_Neo_dir_STR, SNP_Map_STR, sep = '')
consensus_calls_file_STR = paste(data_Neo_dir_STR, consensus_STR, sep = '')
report_file_STR   = paste(output_dir_STR, report_STR, sep = '')


# Load data ####
SNP_map_file_DF = read.csv(SNP_map_file_STR)
SNP_map_file_DF = select(SNP_map_file_DF, Name, Chromosome, Position, refsnp_id, ensembl_gene_stable_id, reg_feature_stable_id)
consensus_DF = read.delim(consensus_calls_file_STR)
consensus_DF = select(consensus_DF, -c(Chromosome, Position..b38.))
combined_DF = left_join(SNP_map_file_DF, consensus_DF, by = c('Name' = 'Marker'))


# Loop through all files and append info to SNP_map. 
# From Raw_Data Tab append Allele column and rename.
# Save 1st background name. Later, remove all backgrounds that did not occur
file_list_CV = list.files(data_dir, recursive = TRUE)
# file_list_CV = refile_list_CV[1]
# fle = "2968 DNA.xlsx"
col_to_keep_CV = c("Name", "Chromosome", "Position", "refsnp_id", "ensembl_gene_stable_id", "reg_feature_stable_id")
SNP_to_keep_CV = c()
for (fle in file_list_CV){
  ident_STR = str_replace(fle, ".xlsx", "")
  ident_STR = str_replace(ident_STR, " ", "_")
  Raw_Data_DF = read.xlsx(paste(data_dir,fle, sep = ""),sheetName = "Raw_Data", startRow = 3, colIndex = c(2,246))
  # Raw_Data_DF = select(Raw_Data_DF, c(SNP.Name, Allele))
  Raw_Data_DF = rename(Raw_Data_DF, !!ident_STR := Allele)
  combined_DF = left_join(combined_DF, Raw_Data_DF, by = c('Name' = 'SNP.Name'))
  primary_bg_DF = read.xlsx(paste(data_dir,fle, sep = ""), sheetName = "Count_Data_1st_BG", colIndex = 3)
  col_to_keep_CV = append(col_to_keep_CV, names(primary_bg_DF))
  col_to_keep_CV = append(col_to_keep_CV, ident_STR)
  SNP_to_keep_CV = unique(append(SNP_to_keep_CV, Raw_Data_DF$SNP.Name))
  print (paste(fle, " finished!"))
}

col_to_keep_CV = append(col_to_keep_CV, "X129S6.SvEvTac")
combined_DF = combined_DF[, names(combined_DF) %in% col_to_keep_CV]
combined_DF = combined_DF[combined_DF$Name %in% SNP_to_keep_CV, ]

## Set up report
wb = createWorkbook()
s1 = createSheet(wb, 'Overview')
s2 = createSheet(wb, 'Data')
s3 = createSheet(wb, 'Description')

s2_rows = createRow(s2, rowIndex = 1:1)
s2_cells = createCell(s2_rows, colIndex = 1)
setCellValue(s2_cells[[1,1]], "All SNPs from MR38 that were different from their primary background.")
addDataFrame(combined_DF, sheet = s2, row.names = FALSE, startRow = 2)

s3_rows = createRow(s3, rowIndex = 1:3)
s3_cells = createCell(s3_rows, colIndex = 1)
setCellValue(s3_cells[[1,1]], "# Description: A report on the differences per SNP of all Samples in MR38 to their primary background.")
setCellValue(s3_cells[[2,1]], "# Details:")
setCellValue(s3_cells[[3,1]], "# Type: Report")

saveWorkbook(wb, report_file_STR)









