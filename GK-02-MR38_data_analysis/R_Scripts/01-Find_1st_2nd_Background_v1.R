# Description: Compare MUGA data from a mouse to all the consensus backgrounds and determine primary and secondary background
# Details: To find secondary background, take the background where the least differences were found and compare those SNPs to all backgrounds.
# Type: Script

options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx8192m"))
options(scipen=999)
library("dplyr")
library("tidyr")
library("ggplot2")
library("janitor")
library("xlsx")

# Set working directory ####
# setwd("S:/AG/AG-Knauf/Daten/AG Knauf/Gabriel Kirchgraber/GK-01-MUGA/")
setwd("C:/Users/sweet/Dropbox/Arbeit/GK-01-MUGA/")

# Set File Paths ####
output_dir_STR = "./GK-02-MR38_data_analysis/R_Scripts/Output/"
plots_dir_STR = "./GK-02-MR38_data_analysis/R_Scripts/Plots/"
data_Neo_dir_STR = "./_Data/Neogene_Files/"

final_report_STR = "Charite_Berlin_MURCOMV02_20191231_FinalReport.txt"
SNP_Map_STR = "SNP_Map_v3.csv"
consensus_STR = "MiniMUGAv2 Consensus Calls.txt"
fun_STR = "MUGA_functions.R"
ideogram_STR = "Background_comparison_1st_2nd_background_ideogram.jpg"
report_STR = "Background_comparison_1st_2nd_background.xlsx"

final_report_file_STR = paste(data_Neo_dir_STR, final_report_STR, sep = '')
SNP_map_file_STR = paste(data_Neo_dir_STR, SNP_Map_STR, sep = '')
consensus_calls_file_STR = paste(data_Neo_dir_STR, consensus_STR, sep = '')
fun_file = paste('./', fun_STR, sep = '')
ideogram_file_STR = paste(plots_dir_STR, ideogram_STR, sep = '')
report_file_STR   = paste(output_dir_STR, report_STR, sep = '')

# Load MUGA functions and set parameters ####
source(fun_file)
sample_to_keep_CV = c()
BG_to_keep_CV = c()
SNP_to_keep_CV = c()

# Load data ####
temp = load_clean_data_FUN(
  final_report_file_STR = final_report_file_STR,
  consensus_calls_file_STR = consensus_calls_file_STR
)
combined_DF = temp[[1]]
track_SNPs_DF = temp[[2]]
rm(temp)

# Function to compare a set of SNPs from a group of mice to all background Allelels. ####
# First pool Allele information of mice in group. Throw out SNPs were not all mice have the same allele
# Then group them and count differences to the backgrounds

cmpr_to_background_FUN = function (combined_DF,
                                   track_SNPs_DF,
                                   SNP_to_keep_CV = NULL,
                                   BG_to_keep_CV = NULL,
                                   sample_to_keep_CV = NULL) {
  # Filter dataset
  if (!is.null(SNP_to_keep_CV)) {
    combined_DF = combined_DF %>% filter(SNP.Name %in% SNP_to_keep_CV)
  }
  if (!is.null(BG_to_keep_CV)) {
    BG_to_keep_CV = append(
      c(
        "SNP.Name",
        "Sample.ID",
        "A1",
        "A2",
        "GC.Score",
        "Chromosome",
        "Position",
        "Allele"
      ),
      BG_to_keep_CV
    )
    combined_DF = combined_DF %>% select(BG_to_keep_CV)
  }
  if (!is.null(sample_to_keep_CV)) {
    combined_DF = combined_DF %>% filter(Sample.ID %in% sample_to_keep_CV)
  }
  
  # Drop SNPs where alleles for mice of interest are not the same
  # Summarize dataframe so that only SNPs where the alleles are different for the samples have more than one line
  # Count number of appearances of the SNP.Name and remove all that appear more than once
  # -> List of SNPs where all Alleles are the same
  SNPs_beforeINT = length(unique(combined_DF$SNP.Name))
  keep_DF = combined_DF %>% group_by(SNP.Name, Chromosome, Position, Allele) %>% summarise()
  keep_VC = filter(count(keep_DF, SNP.Name), n == 1)$SNP.Name
  count_DF = combined_DF[combined_DF$SNP.Name %in% keep_VC,]
  SNPsafterINT = length(unique(count_DF$SNP.Name))
  track_SNPs_DF = add_row(
    track_SNPs_DF,
    Action = 'Drop SNPs that are not same in all samples',
    Dropped_SNPs = SNPs_beforeINT - SNPsafterINT,
    Remaining_SNPs = SNPsafterINT
  )
  
  
  # Only one row per SNP needed
  count_DF = distinct(select(count_DF,-c(Sample.ID, GC.Score, A1, A2)), .keep_all = TRUE)
  count_DF_copy = count_DF
  
  # Compare sample allele to background alleles
  last_column = dim(count_DF)[2]
  cols_VC = names(count_DF)[4:(last_column - 1)]
  for (nme in cols_VC) {
    count_DF[nme] = count_DF$Allele != count_DF[nme]
  }
  
  # Count differences per chromosome
  count_DF = count_DF %>% group_by(Chromosome) %>% summarise(across(all_of(cols_VC), sum))
  
  ## Add totals
  count_DF = adorn_totals(count_DF, name = 'Total')
  
  ## Sort by totals
  sort_order = 'Chromosome'
  last_row = dim(count_DF)[1]
  last_column = dim(count_DF)[2]
  sort_order = append(sort_order, names(sort(count_DF[last_row, 2:last_column])))
  count_DF = count_DF[sort_order]
  
  ## Add rank
  rank_row = c('Rank', seq(1, last_column - 1, 1))
  count_DF = rbind(count_DF, rank_row)
  
  ## Add SNP total
  SNP_total = track_SNPs_DF$Remaining_SNPs[length(track_SNPs_DF$Remaining_SNPs)]
  total_row = c('SNPs', rep(SNP_total, times = last_column - 1))
  count_DF = rbind(count_DF, total_row)
  
  # Get the SNPs that make up the difference to the background with the fewest differences
  diff_SNP_CV = count_DF_copy[count_DF_copy$Allele != count_DF_copy[names(count_DF)[2]], ]$SNP.Name
  
  # Return relevant data
  return(list(
    names(count_DF)[2],
    diff_SNP_CV,
    count_DF,
    count_DF_copy,
    track_SNPs_DF
  ))
  
}

# Start analysis ####

Find_Backgrounds_FUN = function(combined_DF,
                                track_SNPs_DF,
                                sample_to_keep_CV) {
  temp = cmpr_to_background_FUN(combined_DF, track_SNPs_DF, sample_to_keep_CV = sample_to_keep_CV)
  Frst_BG_STR = temp[[1]]
  Frst_BG_SNP_CV = temp[[2]]
  Frst_BG_count_DF = temp[[3]]
  Frst_BG_raw_DF = temp[[4]]
  Frst_BG_track_DF = temp[[5]]
  rm(temp)
  
  temp = cmpr_to_background_FUN(
    combined_DF = combined_DF,
    track_SNPs_DF = track_SNPs_DF,
    SNP_to_keep_CV = Frst_BG_SNP_CV,
    sample_to_keep_CV = sample_to_keep_CV
  )
  
  Scnd_BG_STR = temp[[1]]
  Scnd_BG_SNP_CV = temp[[2]]
  Scnd_BG_count_DF = temp[[3]]
  Scnd_BG_raw_DF = temp[[4]]
  Scnd_BG_track_DF = temp[[5]]
  rm(temp)
  return(
    list(
      Frst_BG_STR,
      Frst_BG_SNP_CV,
      Frst_BG_count_DF,
      Frst_BG_raw_DF,
      Frst_BG_track_DF,
      Scnd_BG_STR,
      Scnd_BG_SNP_CV,
      Scnd_BG_count_DF,
      Scnd_BG_raw_DF,
      Scnd_BG_track_DF
    )
  )
}

sample_to_keep_CV = c("Casp 1A", "Casp 1B", "Casp 1")
temp = Find_Backgrounds_FUN(combined_DF, track_SNPs_DF, sample_to_keep_CV)




# create_report_FUN: Export Report function. Creates a Excel report of the analysis. ####
create_report_FUN = function(combined_DF, track_SNPs_TB, df_1, file_ident_STR, bg_to_keep_CV=NULL){
  
  output_file_report_STR = paste(output_dir_STR, "Background_analysis_",file_ident_STR, "_report.xlsx", sep="")
  
  # Create Excel report
  ## Add totals
  df_1 = adorn_totals(df_1, name = 'Total')
  
  ## Sort by totals
  sort_order = 'Chromosome'
  last_row = dim(df_1)[1]
  last_column = dim(df_1)[2]
  sort_order = append(sort_order, names(sort(df_1[last_row,2:last_column])))
  df_1 = df_1[sort_order]
  
  ## Add rank
  rank_row = c('Rank', seq(1,last_column-1,1))
  df_1 = rbind(df_1, rank_row)
  
  ### Add identifier
  ident_row = c('Sample', rep(file_ident_STR, times = last_column-1))
  df_1 = rbind(df_1, ident_row)
  
  ## Add SNP total
  SNP_total = track_SNPs_TB$Remaining_SNPs[length(track_SNPs_TB$Remaining_SNPs)]
  total_row = c('SNPs', rep(SNP_total, times = last_column-1 ))
  df_1 = rbind(df_1, total_row)
  df_1_copy = df_1
  df_1 = df_1_copy
  
  ## If bg_to_keep is supplied drop unused columns
  if (!is.null(bg_to_keep_CV)){
    df_1 = select(df_1, bg_to_keep_CV)
  }
  
  ## Set up report
  wb = createWorkbook()
  s1 = createSheet(wb, 'Background_comparison')
  s2 = createSheet(wb, 'Dropped_SNPs')
  s3 = createSheet(wb, 'Raw_Data')
  
  addDataFrame(df_1, sheet = s1, row.names = TRUE)
  addDataFrame(track_SNPs_TB, sheet = s2, row.names = TRUE)
  addDataFrame(combined_DF, sheet = s3, row.names = FALSE)
  
  saveWorkbook(wb, output_file_report_STR)
  
  # write.csv2(combined_DF, output_file_raw_STR , na = '')
  # write.csv2(df, output_file_diff_per_chr_STR, na = '')
  # write.csv2(track_SNPs_TB, output_file_dropped_SNP_STR, na = '')
}
# Start Script ####
sample_to_keep_STR = 'Casp 1A'
load_clean_data_FUN(final_report_file_STR, consensus_calls_file_STR, sample_to_keep_STR)

# If creating more than one report use a copy so you don't have to run load_clean_data multiple times
combined_copy_DF = combined_DF
track_SNPs_copy_TB = track_SNPs_TB
analysis_FUN(combined_DF, track_SNPs_copy_TB)




