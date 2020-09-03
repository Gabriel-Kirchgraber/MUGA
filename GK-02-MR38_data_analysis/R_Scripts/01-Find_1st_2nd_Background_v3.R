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
setwd("S:/AG/AG-Knauf/Daten/AG Knauf/Gabriel Kirchgraber/GK-01-MUGA/")
# setwd("C:/Users/sweet/Dropbox/Arbeit/GK-01-MUGA/")

# Set File Paths ####
output_dir_STR = "./GK-02-MR38_data_analysis/R_Scripts/Output/"
plots_dir_STR = "./GK-02-MR38_data_analysis/R_Scripts/Plots/"
data_Neo_dir_STR = "./_Data/Neogene_Files/"
MR38_dir_STR = "./GK-02-MR38_data_analysis/"
data_dir = paste(MR38_dir_STR, "Find_Background_Data/", sep = '')

final_report_STR = "Charite_Berlin_MURCOMV02_20191231_FinalReport.txt"
SNP_Map_STR = "SNP_Map_v3.csv"
consensus_STR = "MiniMUGAv2 Consensus Calls.txt"
fun_STR = "MUGA_functions.R"
ideogram_STR = "Background_comparison_1st_2nd_background_ideogram.jpg"
report_STR = "Background_comparison_1st_2nd_background_pooled.xlsx"
Results_MiniMUGA_STR = "MR-Results MiniMUGA.csv"

final_report_file_STR = paste(data_Neo_dir_STR, final_report_STR, sep = '')
SNP_map_file_STR = paste(data_Neo_dir_STR, SNP_Map_STR, sep = '')
consensus_calls_file_STR = paste(data_Neo_dir_STR, consensus_STR, sep = '')
fun_file = paste('./', fun_STR, sep = '')
ideogram_file_STR = paste(plots_dir_STR, ideogram_STR, sep = '')
report_file_STR   = paste(output_dir_STR, report_STR, sep = '')
Results_MiniMUGA_file_STR = paste(MR38_dir_STR, Results_MiniMUGA_STR, sep = '')

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
  # This also removes all SNPs where at least one Allele is N. Samples where all SNPs are N stay.
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
  
  # remove SNPs with N Allele. These are SNPs that were not caught while loading data, where SNPs where all samples had N were removed.
  # They happen when only our samples of interest have all Ns.
  SNPs_beforeINT = length(unique(count_DF$SNP.Name))
  count_DF = filter(count_DF, Allele != 'N')
  count_DF_copy = count_DF
  SNPsafterINT = length(unique(count_DF$SNP.Name))
  track_SNPs_DF = add_row(
    track_SNPs_DF,
    Action = 'Drop SNPs where one sample had N',
    Dropped_SNPs = SNPs_beforeINT - SNPsafterINT,
    Remaining_SNPs = SNPsafterINT
  )
  count_DF_raw = count_DF
  
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
  count_DF_raw = filter(count_DF_raw, SNP.Name %in% diff_SNP_CV)
  
  # Return relevant data
  return(list(
    names(count_DF)[2],
    diff_SNP_CV,
    count_DF,
    count_DF_raw,
    track_SNPs_DF
  ))
  
}



# Start analysis ####

create_report_FUN = function(MR38_dir_STR, sample_to_keep_CV, temp){
  ## Set up report
  wb = createWorkbook()
  s1 = createSheet(wb, 'Raw_Data')
  s2 = createSheet(wb, 'Count_Data_1st_BG')
  s3 = createSheet(wb, 'SNP_List')
  s4 = createSheet(wb, 'Count_Data_2nd_BG')
  s5 = createSheet(wb, 'Dropped_SNPs')
  s6 = createSheet(wb, 'Description')
  
  sample_to_keep_CV = paste(sample_to_keep_CV, collapse =" ")
  
  s1_rows = createRow(s1, rowIndex = 1:2)
  s1_cells = createCell(s1_rows, colIndex = 1)
  setCellValue(s1_cells[[1,1]], "Raw Data from 01-Find_1st_2nd_Background_search_clusters_v1.R")
  setCellValue(s1_cells[[2,1]], paste("Sample is ", sample_to_keep_CV, " from MR38", sep =''))
  addDataFrame(temp[[4]], sheet = s1, row.names = TRUE, startRow = 3)
  
  s2_rows = createRow(s2, rowIndex = 1:2)
  s2_cells = createCell(s2_rows, colIndex = 1)
  setCellValue(s2_cells[[1,1]], "Count Data from 01-Find_1st_2nd_Background_search_clusters_v1.R")
  setCellValue(s2_cells[[2,1]], paste("Sample is ", sample_to_keep_CV, " from MR38. Differences between Sample and all backgrounds.", sep =''))
  addDataFrame(temp[[3]], sheet = s2, row.names = TRUE, startRow = 3)
  
  s3_rows = createRow(s3, rowIndex = 1:2)
  s3_cells = createCell(s3_rows, colIndex = 1)
  setCellValue(s3_cells[[1,1]], "List of SNPs that make up the differences between sample and background for all SNPs.")
  setCellValue(s3_cells[[2,1]], paste("Sample is ", sample_to_keep_CV, " from MR38.", sep =''))
  SNP_DF = data.frame(SNP = temp[[2]])
  addDataFrame(SNP_DF, sheet = s3, row.names = TRUE, startRow = 3)
  
  s4_rows = createRow(s4, rowIndex = 1:2)
  s4_cells = createCell(s4_rows, colIndex = 1)
  setCellValue(s4_cells[[1,1]], "Count Data from 01-Find_1st_2nd_Background_search_clusters_v1.R")
  setCellValue(s4_cells[[2,1]], paste("Sample is ", sample_to_keep_CV, " from MR38. Differences for differing SNPs from 1st background and all backgrounds.", sep =''))
  addDataFrame(temp[[8]], sheet = s4, row.names = TRUE, startRow = 3)
  
  addDataFrame(temp[[5]], sheet = s5, row.names = TRUE, startRow = 1)
  
  s6_rows = createRow(s6, rowIndex = 1:3)
  s6_cells = createCell(s6_rows, colIndex = 1)
  setCellValue(s6_cells[[1,1]], "# Description: Output from 01-Find_1st_2nd_Background_search_clusters_v1")
  setCellValue(s6_cells[[2,1]], paste("# Details: Finding 1st and 2nd background for ", sample_to_keep_CV, " from MR38", sep =''))
  setCellValue(s6_cells[[3,1]], "# Type: Data")
  
  
  smpl_ident_STR = chartr(old = "/", new = " ", sample_to_keep_CV)
  raw_data_file_STR = paste(data_dir, smpl_ident_STR, "_pooled.xlsx", sep ='' )
  saveWorkbook(wb, raw_data_file_STR)
}

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

# Create dataframe to store analysis summary
analysis_DF = data.frame(Name = character(),
                         First_BG = character(),
                         First_BG_diff = character(),
                         Second_BG = character(),
                         Second_BG_diff = character())


# Loop through Samples and find 1st and 2nd BG for every one. Store Raw data in Excel files and append summary data to Results-MiniMUGA

sample_to_keep_LST = list(c("2968 DNA"),
                          c("C57BL/6N (Chr)","C57BL/6N (Chr)A"),
                          c("C57BL/6N (FEM)","C57BL/6N (FEM)A"),
                          c("C57BL/6N (Janvier)","C57BL/6N (Janvier)A"),
                          c("C57BL6/J","C57BL6/JA"),
                          c("C57BL6/N","C57BL6/NA"),
                          c("Casp 1","Casp 1A","Casp 1B"),
                          c("Casp 11","Casp 11A","Casp 11B"),
                          c("GSDMD","GSDMDA"),
                          c("IL1a","IL1aA","IL1aB"),
                          c("IL1b","IL1bA","IL1bB"),
                          c("Mut439","Mut439A","Mut439B"),
                          c("P2rx7","P2rx7A"),
                          c("Slc26a6")
                          )

# for (nme in unique(combined_DF$Sample.ID)){
for (nme in sample_to_keep_LST){  
  sample_to_keep_CV = c(nme)
  temp = Find_Backgrounds_FUN(combined_DF, track_SNPs_DF, sample_to_keep_CV)
  create_report_FUN(MR38_dir_STR, sample_to_keep_CV, temp)
  
  First_BG_STR = names(temp[[3]])
  First_BG_STR = paste("1.", First_BG_STR[2], " 2.", First_BG_STR[3], " 3.", First_BG_STR[4], " 4.", First_BG_STR[5], " 5.", First_BG_STR[6], sep ='')
  
  last_row = dim(temp[[3]])[1]
  First_BG_diff_STR = ''
  for (i in seq(2,6,1)){
    First_BG_diff_STR = paste(First_BG_diff_STR, i-1,":",temp[[3]][last_row-2,i],"/",temp[[3]][last_row,i]," ", sep='')
  }
  
  Second_BG_STR = names(temp[[8]])
  Second_BG_STR = paste("1.", Second_BG_STR[2], " 2.", Second_BG_STR[3], " 3.", Second_BG_STR[4], " 4.", Second_BG_STR[5], " 5.", Second_BG_STR[6], sep ='')
  
  last_row = dim(temp[[8]])[1]
  Second_BG_diff_STR = ''
  for (i in seq(2,6,1)){
    Second_BG_diff_STR = paste(Second_BG_diff_STR, i-1,":",temp[[8]][last_row-2,i],"/",temp[[8]][last_row,i]," ", sep='')
  }
  
  temp_DF = data.frame(Name = sample_to_keep_CV,
                       First_BG = First_BG_STR,
                       First_BG_diff = First_BG_diff_STR,
                       Second_BG = Second_BG_STR,
                       Second_BG_diff = Second_BG_diff_STR
  )
  analysis_DF = rbind(analysis_DF, temp_DF)
  print (paste(sample_to_keep_CV, "finished!"))
}


Results_MiniMUGA_DF = read.csv2(Results_MiniMUGA_file_STR, skip=3)
Results_MiniMUGA_DF = left_join(Results_MiniMUGA_DF, analysis_DF, by = c('ID'='Name'))

write.xlsx(Results_MiniMUGA_DF, report_file_STR)









