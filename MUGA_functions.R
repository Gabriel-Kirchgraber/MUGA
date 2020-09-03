# Description: Functions for analysing the MUGA data
# Details: load_clean_data_FUN
# Type: Script

# library('dplyr')
# library("tidyr")
# # Set working directory ####
# setwd("S:/AG/AG-Knauf/Daten/AG Knauf/Gabriel Kirchgraber/GK-01-MUGA/")
# # setwd("C:/Users/sweet/Dropbox/Arbeit/GK-01-MUGA/")
# 
# # Set File Paths ####
# output_dir_STR = "./GK-02-MR38_data_analysis/R_Scripts/Output/"
# plots_dir_STR = "./GK-02-MR38_data_analysis/R_Scripts/Plots/"
# data_Neo_dir_STR = "./_Data/Neogene_Files/"
# 
# final_report_STR = "Charite_Berlin_MURCOMV02_20191231_FinalReport.txt"
# SNP_Map_STR = "SNP_Map_v3.csv"
# consensus_STR = "MiniMUGAv2 Consensus Calls.txt"
# fun_STR = "MUGA_functions.R"
# ideogram_STR = "3xBL6.jpg"
# 
# final_report_file_STR = paste(data_Neo_dir_STR, final_report_STR, sep = '')
# SNP_map_file_STR = paste(data_Neo_dir_STR, SNP_Map_STR, sep = '')
# consensus_calls_file_STR = paste(data_Neo_dir_STR, consensus_STR, sep = '')
# fun_file = paste('./', fun_STR, sep = '')
# ideogram_file_STR = paste(plots_dir_STR, ideogram_STR, sep='')
# 
# sample_to_keep_CV = c("C57BL/6N (Chr)","C57BL/6N (Chr)A","C57BL/6N (FEM)","C57BL/6N (FEM)A","C57BL/6N (Janvier)","C57BL/6N (Janvier)A")
# BG_to_keep_CV = c('C57BL.6NCrl')
# SNP_to_keep_CV = c()
# 
# sample_to_keep_CV = c()
# BG_to_keep_CV = c()
# SNP_to_keep_CV = c()

# load_clean_data_FUN: Function to load and clean data. Output is combined_DF and track_SNPs_TB ####
load_clean_data_FUN = function(final_report_file_STR,
                               consensus_calls_file_STR,
                               sample_to_keep_CV = NULL,
                               BG_to_keep_CV = NULL,
                               SNP_to_keep_CV = NULL,
                               SNP_map_file_STR = NULL) {
  ## read data into dataframe
  final_report_DF = read.delim(final_report_file_STR, skip = 9)
  consensus_DF = read.delim(consensus_calls_file_STR)
  
  ## Filter for samples
  if (!is.null(sample_to_keep_CV)) {
    final_report_DF = filter(final_report_DF, Sample.ID %in% sample_to_keep_CV)
  }
  
  ## rename columns
  final_report_DF = rename(final_report_DF, A1 = Allele1...Forward, A2 = Allele2...Forward)
  consensus_DF = rename(consensus_DF, Position = Position..b38.)
  
  ## filter for backgrounds
  if (!is.null(BG_to_keep_CV)) {
    BG_to_keep_CV = append(c('Marker', 'Chromosome', 'Position'), all_of(BG_to_keep_CV))
    consensus_DF = select(consensus_DF, BG_to_keep_CV)
  }
  
  ## drop unused columns in final_report
  final_report_DF = subset(final_report_DF, select = -c(X, Y, Theta, X.Raw, Y.Raw, R))
  
  ## Make all nucleotides uppercase in consensus calls
  last_column = dim(consensus_DF)[2]
  consensus_DF[, 3:last_column] = mutate_all(consensus_DF[, 3:last_column], .funs = toupper) #Start at Position. If only one background 4:4 gives error
  
  ## Merge consensus calls with final_report
  combined_DF = left_join(final_report_DF, consensus_DF, by = c('SNP.Name' = 'Marker'))
  
  ## Improve column names. Some have / inside, interferes with including them in file paths
  names(combined_DF) = make.names(names(combined_DF))
  
  ## Make position numeric
  combined_DF$Position = as.numeric(combined_DF$Position )
  
  ## Create dataframe to track dropped SNPs
  track_SNPs_DF = data.frame(
    Action = 'Initialize',
    Dropped_SNPs = 0,
    Remaining_SNPs = length(unique(combined_DF$SNP.Name))
  )
  
  if (is.null(SNP_to_keep_CV)) {
    ## Drop constructs
    SNPs_beforeINT = length(unique(combined_DF$SNP.Name))
    combined_DF = filter(combined_DF,!is.na(Chromosome))
    SNPsafterINT = length(unique(combined_DF$SNP.Name))
    track_SNPs_DF = add_row(
      track_SNPs_DF,
      Action = 'Drop constructs',
      Dropped_SNPs = SNPs_beforeINT - SNPsafterINT,
      Remaining_SNPs = SNPsafterINT
    )
    
    ## Replace '-' with N
    levels(combined_DF$A1) = c(levels(combined_DF$A1), 'N')
    levels(combined_DF$A2) = c(levels(combined_DF$A2), 'N')
    combined_DF = combined_DF %>% mutate(A1 = replace(A1, A1 == '-', 'N'))
    combined_DF = combined_DF %>% mutate(A2 = replace(A2, A2 == '-', 'N'))
    
    ## Drop duplicated markers
    SNPs_beforeINT = length(unique(combined_DF$SNP.Name))
    combined_DF = distinct(combined_DF, Sample.ID, Position, Chromosome, .keep_all = TRUE)
    SNPsafterINT = length(unique(combined_DF$SNP.Name))
    track_SNPs_DF = add_row(
      track_SNPs_DF,
      Action = 'Drop duplicates',
      Dropped_SNPs = SNPs_beforeINT - SNPsafterINT,
      Remaining_SNPs = SNPsafterINT
    )
    
    ## Drop SNPs where every Allele is 'N' (GC.Score == 0)
    SNPs_beforeINT = length(unique(combined_DF$SNP.Name))
    GC_TB = combined_DF %>% group_by(SNP.Name) %>% summarise(gcscore = sum(GC.Score))
    GC_TB = filter(GC_TB, gcscore > 0)
    combined_DF = combined_DF[combined_DF$SNP.Name %in% GC_TB$SNP.Name, ]
    SNPsafterINT = length(unique(combined_DF$SNP.Name))
    track_SNPs_DF = add_row(
      track_SNPs_DF,
      Action = 'Dropped SNPs with all N',
      Dropped_SNPs = SNPs_beforeINT - SNPsafterINT,
      Remaining_SNPs = SNPsafterINT
    )
    
  } else {
    SNPs_beforeINT = length(unique(combined_DF$SNP.Name))
    combined_DF = filter(combined_DF, SNP.Name %in% SNP_to_keep_CV)
    SNPsafterINT = length(unique(combined_DF$SNP.Name))
    track_SNPs_DF = add_row(
      track_SNPs_DF,
      Action = 'Dropped SNPs by file input',
      Dropped_SNPs = SNPs_beforeINT - SNPsafterINT,
      Remaining_SNPs = SNPsafterINT
    )
  }
  
  ## Add info for allele: G/G -> G, A/G -> H, N/G -> N [does not happen], N/N ->N
  combined_DF$Allele = NA
  combined_DF$Allele[combined_DF$A1 == combined_DF$A2] = filter(combined_DF, A1 == A2)$A1
  combined_DF$Allele[is.na(combined_DF$Allele)] = 'H'
  
  ## Add factors and levels in correct order for sorting
  combined_DF$Chromosome = factor(combined_DF$Chromosome, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "MT", "PAR", "X", "Y"))
  
  ## Add info from SNP_map
  if(!is.null(SNP_map_file_STR)){
    SNP_map_DF = read.csv(SNP_map_file_STR)
    SNP_map_DF = select(SNP_map_DF, Name, refsnp_id, ensembl_gene_stable_id, reg_feature_stable_id)
    combined_DF = left_join(combined_DF, SNP_map_DF, by = c('SNP.Name' = 'Name'))
  }
  
  ## Return combined_DF and track_SNPs_DF
  return(list(combined_DF, track_SNPs_DF))
  
}