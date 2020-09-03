options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx8192m"))
library("dplyr")
library("tidyr")
library("ggplot2")
library("janitor")
library("xlsx")
options(scipen=999)
# Compare set of mice to all backgrounds

# Set working directory ####
setwd("S:/AG/AG-Knauf/Daten/AG Knauf/Gabriel Kirchgraber/GK-01.00-MUGA/")
# setwd("C:/Users/sweet/Dropbox/Arbeit/GK-01.00-MUGA/")

# Set File Paths ####
output_dir_STR = "./GK-01.20-MR38_data_analysis/R_Scripts/Output/"
plots_dir_STR = "./GK-01.20-MR38_data_analysis/R_Scripts/Plots/"
data_Neo_dir_STR = "./.Data/Neogene_Files/"

final_report_file_STR = paste(data_Neo_dir_STR, "Charite_Berlin_MURCOMV02_20191231_FinalReport.txt", sep ='')
consensus_calls_file_STR = paste(data_Neo_dir_STR, "MiniMUGAv2 Consensus Calls.txt", sep = '')

# load_clean_data_FUN: Function to load and clean data. Output is combined_DF and track_SNPs_TB ####
load_clean_data_FUN = function(final_report_file_STR, consensus_calls_file_STR){
  # Load data ####
  combined_DF = NA
  track_SNPs_TB = NA
  ## read data into dataframe
  final_report_DF = read.delim(final_report_file_STR, skip=9)
  consensus_DF = read.delim(consensus_calls_file_STR)
  
  ## rename columns
  final_report_DF = rename(final_report_DF, A1 = Allele1...Forward, A2 = Allele2...Forward)
  consensus_DF = rename(consensus_DF, Position = Position..b38.)
  
  ## drop unused columns
  final_report_DF = subset(final_report_DF, select = -c(X, Y, Theta, X.Raw, Y.Raw, R))
  
  ## Make all nucleotides uppercase in consensus calls
  last_column = dim(consensus_DF)[2]
  consensus_DF[,4:last_column] = mutate_all(consensus_DF[,4:last_column], .funs = toupper)
  
  ## Merge consensus calls with final_report
  combined_DF = left_join(final_report_DF, consensus_DF, by = c('SNP.Name' = 'Marker'))
  
  # Clean Data ####
  ## Create tibble to track dropped SNPs
  track_SNPs_TB = tibble(Action = 'Initialize', Dropped_SNPs = 0, Remaining_SNPs = length(unique(combined_DF$SNP.Name)))
  
  ## Drop constructs
  SNPs_beforeINT = length(unique(combined_DF$SNP.Name))
  combined_DF = filter(combined_DF, !is.na(Chromosome))
  SNPsafterINT = length(unique(combined_DF$SNP.Name))
  track_SNPs_TB = add_row(track_SNPs_TB, Action = 'Drop constructs', Dropped_SNPs = SNPs_beforeINT - SNPsafterINT, Remaining_SNPs = SNPsafterINT)
  
  ## Replace '-' with N
  levels(combined_DF$A1) = c(levels(combined_DF$A1), 'N')
  levels(combined_DF$A2) = c(levels(combined_DF$A2), 'N')
  combined_DF = combined_DF %>% mutate(A1 = replace(A1, A1 == '-', 'N'))
  combined_DF = combined_DF %>% mutate(A2 = replace(A2, A2 == '-', 'N'))
  
  ## Drop duplicated markers
  SNPs_beforeINT = length(unique(combined_DF$SNP.Name))
  ### Protect markers with NA for Chromosome and Position from distinct. 
  ### All SNPS where Chromosome is NA get a unique negative Number for position
  # combined_DF[is.na(combined_DF$Chromosome),]$Position = seq(-1, -dim(combined_DF[is.na(combined_DF$Chromosome),])[1], -1)
  ### Remove duplicated 
  combined_DF = distinct(combined_DF, Sample.ID, Position, Chromosome, .keep_all = TRUE)
  ### Change neg. Position values back to NA
  # combined_DF[combined_DF$Position<0,]$Position = NA
  SNPsafterINT = length(unique(combined_DF$SNP.Name))
  track_SNPs_TB = add_row(track_SNPs_TB, Action = 'Drop duplicates', Dropped_SNPs = SNPs_beforeINT - SNPsafterINT, Remaining_SNPs = SNPsafterINT)
  
  ## Drop SNPs where every Allele is 'N' (GC.Score == 0)
  SNPs_beforeINT = length(unique(combined_DF$SNP.Name))
  GC_TB = combined_DF %>% group_by(SNP.Name) %>% summarise(gcscore = sum(GC.Score))
  GC_TB = filter(GC_TB, gcscore > 0)
  combined_DF = combined_DF[combined_DF$SNP.Name %in% GC_TB$SNP.Name,]
  SNPsafterINT = length(unique(combined_DF$SNP.Name))
  track_SNPs_TB = add_row(track_SNPs_TB, Action = 'Drop SNPs with all N', Dropped_SNPs = SNPs_beforeINT - SNPsafterINT, Remaining_SNPs = SNPsafterINT)
  
  ## Drop SNPs where at least one SNP has N
  SNPs_beforeINT = length(unique(combined_DF$SNP.Name))
  ### 1. Create copy of combined_DF to hold list auf SNPs with at least one N call
  ###    Drop unused columns. A2 can be dropped because if A1 = N, A2 is always N
  N_DF = subset(combined_DF, select = c(SNP.Name, Sample.ID, A1))
  ### 2. Replace all not N calls with NA
  N_DF = N_DF %>% mutate(A1 = replace(A1, A1 != 'N', NA))
  ### 3. Group by SNP.Name and summarise by A1
  ###    All SNPs are grouped by SNP, SNPs with more than one N get multiple entries. SNPs with no N and only NAs get one entry
  N_DF = N_DF %>% group_by(SNP.Name) %>% summarise(A1)
  ### 4. Drop all SNPs that are not N
  N_DF = filter(N_DF, A1 == 'N')
  ### 5. Remove duplicates
  ###    List now only contains SNPs with at least one N call
  N_DF = unique(N_DF$SNP.Name)
  combined_DF = combined_DF[!(combined_DF$SNP.Name %in% N_DF),]
  SNPsafterINT = length(unique(combined_DF$SNP.Name))
  track_SNPs_TB = add_row(track_SNPs_TB, Action = 'Drop SNPs with at least one N', Dropped_SNPs = SNPs_beforeINT - SNPsafterINT, Remaining_SNPs = SNPsafterINT)
  combined_DF <<- combined_DF
  track_SNPs_TB <<- track_SNPs_TB
}
# analysis_FUN: Function to analyse data. Output is combined_DF and track_SNPs_TB ####
analysis_FUN = function(combined_DF, track_SNPs_TB){
  # Add information about group Allele. All animals in Group with same nucleotid, get G/A/T/C, when heterozygous get H #### 
  # only works if no more Ns in dataset
  combined_DF$Allele = NA
  combined_DF$Allele[combined_DF$A1 == combined_DF$A2] = filter(combined_DF, A1==A2)$A1
  combined_DF$Allele[is.na(combined_DF$Allele)] = 'H'
  
  # Drop unused columns
  combined_DF = subset(combined_DF, select = -c(A1, A2, GC.Score))
  
  # Drop SNPs where alleles for mice of interest are not the same ####
  # Summarize dataframe so that only SNPs where the alleles are different for the samples have more than one line
  # Count number of appearances of the SNP.Name and remove all that appear more than once 
  # -> List of SNPs where all Alleles are the same
  SNPs_beforeINT = length(unique(combined_DF$SNP.Name))
  keep_DF = combined_DF %>% group_by(SNP.Name, Chromosome, Position, Allele) %>% summarise()
  keep_VC = filter(count(keep_DF, SNP.Name),n==1)$SNP.Name
  combined_DF = combined_DF[combined_DF$SNP.Name %in% keep_VC, ]
  SNPsafterINT = length(unique(combined_DF$SNP.Name))
  track_SNPs_TB = add_row(track_SNPs_TB, Action = 'Drop SNPs that are not same in all samples', Dropped_SNPs = SNPs_beforeINT - SNPsafterINT, Remaining_SNPs = SNPsafterINT)
  
  ## List of SNPs to keep is used in BL6NCrl vs J Background comparison analysis
  #SNPs_to_keep_VC = unique(combined_DF$SNP.Name)
  #save(SNPs_to_keep_VC, file = 'C:/Users/sweet/Dropbox/Arbeit/MUGA/DATA/MR38_Background/Background_comparison_SNPs_to_keep.Rdata')
  
  # Only one row per SNP needed
  combined_DF$Chromosome = factor(combined_DF$Chromosome, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "MT", "PAR", "X", "Y"))
  combined_DF = distinct(select(combined_DF, !2), .keep_all = TRUE)
  
  # Compare sample Allele to background allele ####
  ## use bracket notation to add columns by variable name
  df = combined_DF
  last_column = dim(combined_DF)[2]
  cols_VC = names(combined_DF)[4:(last_column-1)]
  for (nme in cols_VC){
    df[nme] = df$Allele != df[nme]
  }
  
  df = df %>% group_by(Chromosome) %>% summarise(across(all_of(cols_VC), sum))
  
  df_1 <<- df
  combined_DF <<- combined_DF
  track_SNPs_TB <<- track_SNPs_TB
}

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
load_clean_data_FUN(final_report_file_STR, consensus_calls_file_STR)
# If creating more than one report use a copy so you don't have to run load_clean_data multiple times
combined_copy_DF = combined_DF
track_SNPs_copy_TB = track_SNPs_TB

# Remove Samples
# samples_to_keep_VC = c(194, 195, 196)
# samples_to_keep_VC = c("Mut439A","Mut439B","Mut439")
# samples_to_keep_VC = c('Casp 1A', 'Casp 1B', 'Casp 1')
# samples_to_keep_VC = c('IL1aA', 'IL1aAB', 'IL1aA')
# samples_to_keep_VC = c('IL1bA', 'IL1bAB', 'IL1bA')

#Remove columns
bg_to_keep_CV = c('Chromosome', 'C57BL.6J', 'C57BL.6NCrl','C57BL.6NJ')

# combined_DF = filter(combined_copy_DF, Sample.ID %in% samples_to_keep_VC)
# analysis_FUN(combined_DF, track_SNPs_copy_TB)
# 
# ## Set identifying part of File name
# file_ident_STR ="NLRP3"
# create_report_FUN(combined_DF, track_SNPs_TB, df_1, file_ident_STR, bg_to_keep_CV)

# samples_to_keep_VC = c('Casp 1A', 'Casp 1B', 'Casp 1')
# combined_DF = filter(combined_copy_DF, Sample.ID %in% samples_to_keep_VC)
# analysis_FUN(combined_DF, track_SNPs_copy_TB)
# file_ident_STR ="Casp1"
# create_report_FUN(combined_DF, track_SNPs_TB, df_1, file_ident_STR, bg_to_keep_CV)

samples_to_keep_VC = c('IL1aA', 'IL1aAB', 'IL1aA')
combined_DF = filter(combined_copy_DF, Sample.ID %in% samples_to_keep_VC)
analysis_FUN(combined_DF, track_SNPs_copy_TB)
file_ident_STR ="IL1a"
create_report_FUN(combined_DF, track_SNPs_TB, df_1, file_ident_STR, bg_to_keep_CV)

samples_to_keep_VC = c('IL1bA', 'IL1bAB', 'IL1bA')
combined_DF = filter(combined_copy_DF, Sample.ID %in% samples_to_keep_VC)
analysis_FUN(combined_DF, track_SNPs_copy_TB)
file_ident_STR ="IL1b"
create_report_FUN(combined_DF, track_SNPs_TB, df_1, file_ident_STR, bg_to_keep_CV)

final_report_file_STR = paste(data_Neo_dir_STR, "Charite_Berlin_MURCOMV02_20200608_FinalReport.txt", sep ='')
load_clean_data_FUN(final_report_file_STR, consensus_calls_file_STR)
samples_to_keep_VC = c(194, 195, 196)
combined_DF = filter(combined_DF, Sample.ID %in% samples_to_keep_VC)
analysis_FUN(combined_DF, track_SNPs_TB)
file_ident_STR ="protected"
create_report_FUN(combined_DF, track_SNPs_TB, df_1, file_ident_STR, bg_to_keep_CV)


