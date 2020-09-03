library("dplyr")
library("tidyr")
library("ggplot2")
# Compare the background of various BL6 strains from MR38 with each other

# Set working directory ####
#setwd("S:/AG/AG-Knauf/Daten/AG Knauf/Gabriel Kirchgraber/GK-01.00-MUGA")
setwd("C:/Users/sweet/Documents/MUGA")
# Set File Paths ####
final_report_STR = 'S:/AG/AG-Knauf/Daten/AG Knauf/Gabriel Kirchgraber/GK-01.00-MUGA/.Data/Neogene_Files/Charite_Berlin_MURCOMV02_20200608_FinalReport.txt'
consensus_calls_STR = 'S:/AG/AG-Knauf/Daten/AG Knauf/Gabriel Kirchgraber/GK-01.00-MUGA/.Data/Neogene_Files/MiniMUGAv2 Consensus Calls.txt'
# final_report_STR = 'C:/Users/sweet/Documents/MUGA/DATA/Original/Charite_Berlin_MURCOMV02_20200608_FinalReport.txt'
# consensus_calls_STR = 'C:/Users/sweet/Documents/MUGA/DATA/Original/MiniMUGAv2 Consensus Calls.txt'

# Load data ####
## read data into tibble
final_report_DF = read.delim(final_report_STR, skip=9)
consensus_DF = read.delim(consensus_calls_STR)
final_report_TB = as_tibble(final_report_DF)
consensus_TB = as_tibble(consensus_DF)
## rename columns
final_report_TB = rename(final_report_TB, A1 = Allele1...Forward, A2 = Allele2...Forward)
consensus_TB = rename(consensus_TB, Position = Position..b38.)
## drop unused columns
final_report_TB = select(final_report_TB, -c('X', 'Y', 'Theta', 'X.Raw', 'Y.Raw', 'R'))
consensus_TB = select(consensus_TB, c('Marker', 'Chromosome', 'Position'))
## Merge consensus calls with final_report
combined_TB = left_join(final_report_TB, consensus_TB, by = c('SNP.Name' = 'Marker'))
## Clear memory
rm (consensus_DF, consensus_TB, final_report_DF, final_report_TB, consensus_calls_STR, final_report_STR)

# Clean Data ####

## Create tibble to track dropped SNPs
track_SNPs_TB = tibble(Action = 'Initialize', Dropped_SNPs = 0, Remaining_SNPs = length(unique(combined_TB$SNP.Name)))

## Replace '-' with N
levels(combined_TB$A1) = c(levels(combined_TB$A1), 'N')
levels(combined_TB$A2) = c(levels(combined_TB$A2), 'N')
combined_TB = combined_TB %>% mutate(A1 = replace(A1, A1 == '-', 'N'))
combined_TB = combined_TB %>% mutate(A2 = replace(A2, A2 == '-', 'N'))

## Make all nucleotides uppercase in consensus calls
#consensus_TB$C57BL.6NCrl = toupper(consensus_TB$C57BL.6NCrl)

## Drop duplicated markers
SNPs_beforeINT = length(unique(combined_TB$SNP.Name))
### Protect markers with NA for Chromosome and Position from distinct. 
### All SNPS where Chromosome is NA get a unique negative Number for position
combined_TB[is.na(combined_TB$Chromosome),]$Position = seq(-1, -dim(combined_TB[is.na(combined_TB$Chromosome),])[1], -1)
### Remove duplicated 
combined_TB = distinct(combined_TB, Sample.ID, Position, Chromosome, .keep_all = TRUE)
### Change neg. Position values back to NA
combined_TB[combined_TB$Position<0,]$Position = NA
SNPsafterINT = length(unique(combined_TB$SNP.Name))
track_SNPs_TB = add_row(track_SNPs_TB, Action = 'Drop duplicates', Dropped_SNPs = SNPs_beforeINT - SNPsafterINT, Remaining_SNPs = SNPsafterINT)

## Drop SNPs where every Allele is 'N' (GC.Score == 0)
SNPs_beforeINT = length(unique(combined_TB$SNP.Name))
GC_TB = combined_TB %>% group_by(SNP.Name) %>% summarise(gcscore = sum(GC.Score))
GC_TB = filter(GC_TB, gcscore > 0)
combined_TB = combined_TB[combined_TB$SNP.Name %in% GC_TB$SNP.Name,]
SNPsafterINT = length(unique(combined_TB$SNP.Name))
track_SNPs_TB = add_row(track_SNPs_TB, Action = 'Drop SNPs with all N', Dropped_SNPs = SNPs_beforeINT - SNPsafterINT, Remaining_SNPs = SNPsafterINT)
rm(GC_TB)

## Drop SNPs where at least one SNP has N
SNPs_beforeINT = length(unique(combined_TB$SNP.Name))
### 1. Create copy of combined_DF to hold list auf SNPs with at least one N call
###    Drop unused columns. A2 can be dropped because if A1 = N, A2 is always N
N_DF = select(combined_TB, SNP.Name, Sample.ID, A1)
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
combined_TB = combined_TB[!(combined_TB$SNP.Name %in% N_DF),]
SNPsafterINT = length(unique(combined_TB$SNP.Name))
track_SNPs_TB = add_row(track_SNPs_TB, Action = 'Drop SNPs with at least one N', Dropped_SNPs = SNPs_beforeINT - SNPsafterINT, Remaining_SNPs = SNPsafterINT)
rm(N_DF)

# Add information about Genotype/ Grouping to Dataframe ####
GT_DF = data.frame(Genotype = c('ko', 'ko', 'ko', 'ko', 'ko', 'wt', 'wt', 'wt' ,'wt', 'wt'),
                   Sample.ID = c(189, 190, 191, 192, 193, 194, 195, 196, 197, 198))
combined_TB = left_join(combined_TB, GT_DF)

# Add information about group Allele. All animals in Group with same nucleotid, get G/A/T/C, when heterozygous get h #### 
# only works if no more Ns in dataset
combined_TB$Allele = NA
combined_TB$Allele[combined_TB$A1 == combined_TB$A2] = filter(combined_TB, A1==A2)$A1
combined_TB$Allele[is.na(combined_TB$Allele)] = 'h'

# Add information about Allele combination of group ####
## Walk through dataframe. Use vectors below to get ko and wt Alleles and save in ko and wt columns as string
ko_select_VB = GT_DF$Genotype == 'ko'
wt_select_VB = GT_DF$Genotype == 'wt'
step_size_INT = length(ko_select_VB)
combined_TB = arrange(combined_TB, SNP.Name)
combined_TB$wt = NA
combined_TB$ko = NA
combined_TB$wt1 = NA
combined_TB$ko1 = NA
combined_TB$different = NA

for (i in seq(0, SNPsafterINT*step_size_INT-step_size_INT, step_size_INT)){
  combined_TB[(i+1):(i+step_size_INT),]$ko = paste(combined_TB[(i+1):(i+step_size_INT),]$Allele[ko_select_VB], collapse = '')
  combined_TB[(i+1):(i+step_size_INT),]$wt = paste(combined_TB[(i+1):(i+step_size_INT),]$Allele[wt_select_VB], collapse = '')
}

#save(combined_TB, file = 'S:/AG/AG-Knauf/Daten/AG Knauf/Gabriel Kirchgraber/GK-01.00-MUGA/.Analysis/Output/Compare_SNPs_WT_KO.rdata')
save(combined_TB, file = 'C:/Users/sweet/Documents/MUGA/DATA/Output/Compare_SNPs_WT_KO.rdata')

df = combined_TB

df = df %>% group_by(SNP.Name, ko, wt, Chromosome, GC.SCore, Position, different) %>% summarise()

for (i in 1:nrow(df)){
  df$wt1[i] = strsplit(df$wt[[i]], split = '')
  df$ko1[i] = strsplit(df$ko[[i]], split = '')
  df$different[i] = sum(df$wt1[[i]] %in% df$ko1[[i]])
}

save(df, file = 'C:/Users/sweet/Documents/MUGA/DATA/Output/Compare_SNPs_WT_KO.rdata')
load('C:/Users/sweet/Documents/MUGA/DATA/Output/Compare_SNPs_WT_KO.rdata')
df1 = df1 %>% group_by(SNP.Name, ko, wt, Chromosome, Position, different) %>% summarise()

df1 = df[c(1:4,7:9)]
write.csv2(df1, 'C:/Users/sweet/Documents/MUGA/DATA/Output/Compare_SNPs_WT_KO.csv')
write.csv2(GT_DF, 'C:/Users/sweet/Documents/MUGA/DATA/Output/Group mapping.csv')
write.csv(track_SNPs_TB, 'C:/Users/sweet/Documents/MUGA/DATA/Output/Dropped SNPs.csv')


combined_TB = read.csv2('S:/AG/AG-Knauf/Daten/AG Knauf/Gabriel Kirchgraber/GK-01.00-MUGA/.Analysis/SNPs_comparison_KO_WT.csv')
SNP_annotation_DF = read.csv('S:/AG/AG-Knauf/Daten/AG Knauf/Gabriel Kirchgraber/GK-01.00-MUGA/.Data/mart_export.txt')
combined_TB = left_join(combined_TB, SNP_annotation_DF, by = c("Position" = "Chromosome.scaffold.position.start..bp."))
combined_TB = combined_TB[c(1:7, 9, 11:13)]


combined_TB = distinct(combined_TB, Position, .keep_all = TRUE)

write.csv2(combined_TB, 'S:/AG/AG-Knauf/Daten/AG Knauf/Gabriel Kirchgraber/GK-01.00-MUGA/.Analysis/Output/SNPs_comparison_KO_WT.csv')







