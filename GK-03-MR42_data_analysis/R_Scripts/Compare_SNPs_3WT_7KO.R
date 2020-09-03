library("dplyr")
library("tidyr")
library("biomaRt")

# Compare the background of various BL6 strains from MR38 with each other

# Set working directory ####
setwd("S:/AG/AG-Knauf/Daten/AG Knauf/Gabriel Kirchgraber/GK-01.00-MUGA/GK-01.30-MR42_data_analysis/")
#setwd("C:/Users/sweet/Documents/MUGA")
# Set File Paths ####
final_report_STR = 'S:/AG/AG-Knauf/Daten/AG Knauf/Gabriel Kirchgraber/GK-01.00-MUGA/.Data/Neogene_Files/Charite_Berlin_MURCOMV02_20200608_FinalReport.txt'
consensus_calls_STR = 'S:/AG/AG-Knauf/Daten/AG Knauf/Gabriel Kirchgraber/GK-01.00-MUGA/.Data/Neogene_Files/MiniMUGAv2 Consensus Calls.txt'
# final_report_STR = 'C:/Users/sweet/Documents/MUGA/DATA/Original/Charite_Berlin_MURCOMV02_20200608_FinalReport.txt'
# consensus_calls_STR = 'C:/Users/sweet/Documents/MUGA/DATA/Original/MiniMUGAv2 Consensus Calls.txt'

# Load data ####
## read data into dataframe
final_report_DF = read.delim(final_report_STR, skip=9)
consensus_DF = read.delim(consensus_calls_STR)
## rename columns
final_report_DF = rename(final_report_DF, A1 = Allele1...Forward, A2 = Allele2...Forward)
consensus_DF = rename(consensus_DF, Position = Position..b38.)
## drop unused columns
final_report_DF = subset(final_report_DF, select = -c(X, Y, Theta, X.Raw, Y.Raw, R))
consensus_DF = subset(consensus_DF, select = c(Marker, Chromosome, Position))
## Merge consensus calls with final_report
combined_DF = left_join(final_report_DF, consensus_DF, by = c('SNP.Name' = 'Marker'))
## Clear memory
rm (consensus_DF, final_report_DF, consensus_calls_STR, final_report_STR)

# Clean Data ####

## Create tibble to track dropped SNPs
track_SNPs_TB = tibble(Action = 'Initialize', Dropped_SNPs = 0, Remaining_SNPs = length(unique(combined_DF$SNP.Name)))

## Replace '-' with N
levels(combined_DF$A1) = c(levels(combined_DF$A1), 'N')
levels(combined_DF$A2) = c(levels(combined_DF$A2), 'N')
combined_DF = combined_DF %>% mutate(A1 = replace(A1, A1 == '-', 'N'))
combined_DF = combined_DF %>% mutate(A2 = replace(A2, A2 == '-', 'N'))

## Make all nucleotides uppercase in consensus calls
#consensus_TB$C57BL.6NCrl = toupper(consensus_TB$C57BL.6NCrl)

## Drop duplicated markers
SNPs_beforeINT = length(unique(combined_DF$SNP.Name))
### Protect markers with NA for Chromosome and Position from distinct. 
### All SNPS where Chromosome is NA get a unique negative Number for position
combined_DF[is.na(combined_DF$Chromosome),]$Position = seq(-1, -dim(combined_DF[is.na(combined_DF$Chromosome),])[1], -1)
### Remove duplicated 
combined_DF = distinct(combined_DF, Sample.ID, Position, Chromosome, .keep_all = TRUE)
### Change neg. Position values back to NA
combined_DF[combined_DF$Position<0,]$Position = NA
SNPsafterINT = length(unique(combined_DF$SNP.Name))
track_SNPs_TB = add_row(track_SNPs_TB, Action = 'Drop duplicates', Dropped_SNPs = SNPs_beforeINT - SNPsafterINT, Remaining_SNPs = SNPsafterINT)

## Drop SNPs where every Allele is 'N' (GC.Score == 0)
SNPs_beforeINT = length(unique(combined_DF$SNP.Name))
GC_TB = combined_DF %>% group_by(SNP.Name) %>% summarise(gcscore = sum(GC.Score))
GC_TB = filter(GC_TB, gcscore > 0)
combined_DF = combined_DF[combined_DF$SNP.Name %in% GC_TB$SNP.Name,]
SNPsafterINT = length(unique(combined_DF$SNP.Name))
track_SNPs_TB = add_row(track_SNPs_TB, Action = 'Drop SNPs with all N', Dropped_SNPs = SNPs_beforeINT - SNPsafterINT, Remaining_SNPs = SNPsafterINT)
rm(GC_TB)

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
rm(N_DF)

# Add information about Genotype/ Grouping to Dataframe ####
GT_DF = data.frame(Genotype = c('ko', 'ko', 'ko', 'ko', 'ko', 'wt', 'wt', 'wt' ,'ko', 'ko'),
                   Sample.ID = c(189, 190, 191, 192, 193, 194, 195, 196, 197, 198))
combined_DF = left_join(combined_DF, GT_DF)

# Add information about group Allele. All animals in Group with same nucleotid, get G/A/T/C, when heterozygous get h #### 
# only works if no more Ns in dataset
combined_DF$Allele = NA
combined_DF$Allele[combined_DF$A1 == combined_DF$A2] = filter(combined_DF, A1==A2)$A1
combined_DF$Allele[is.na(combined_DF$Allele)] = 'h'

# Add information about Allele combination of group und difference score ####
## Walk through dataframe. Use vectors below to get ko and wt Alleles and save in ko and wt columns as string
ko_select_VB = GT_DF$Genotype == 'ko'
wt_select_VB = GT_DF$Genotype == 'wt'
step_size_INT = length(ko_select_VB)
combined_DF = arrange(combined_DF, SNP.Name)
combined_DF$wt = 'NA'
combined_DF$ko = 'NA'
combined_DF$wt1 = c('NA', 'NA')
combined_DF$ko1 = c('NA', 'NA')
combined_DF$difference_score = -1
## Add Allele information for grouping as string
for (i in seq(0, SNPsafterINT*step_size_INT-step_size_INT, step_size_INT)){
  combined_DF[(i+1):(i+step_size_INT),]$ko = paste(combined_DF[(i+1):(i+step_size_INT),]$Allele[ko_select_VB], collapse = '')
  combined_DF[(i+1):(i+step_size_INT),]$wt = paste(combined_DF[(i+1):(i+step_size_INT),]$Allele[wt_select_VB], collapse = '')
}
## Group dataframe so that loop is shorter
df = combined_DF
df = df %>% group_by(SNP.Name, ko, wt, Chromosome, wt1, ko1, Position, difference_score) %>% summarise(GC.Score = mean(GC.Score))
## Add a difference score 
for (i in 1:nrow(df)){
  df$wt1[i] = strsplit(df$wt[[i]], split = '')
  df$ko1[i] = strsplit(df$ko[[i]], split = '')
  df$difference_score[i] = sum(df$wt1[[i]] %in% df$ko1[[i]])
}

## removed ko1 and wt1 column
df = select(df, -c('wt1', 'ko1'))

# Export Report ####
## Download annotation files from Ensmbl
ensembl = useEnsembl(biomart = 'snps', dataset = 'mmusculus_snp')

## Set filter type 
bio_filters = c('chromosomal_region')

## Set filter values for SNPs to download. All SNPs where not ALL SNPs were different [< max value]
bio_values_VC = filter(df, difference_score< max(df$difference_score))
bio_values_VC$list = paste(bio_values_VC$Chromosome, ':', bio_values_VC$Position, ':', bio_values_VC$Position, sep = "")
bio_values_VC = bio_values_VC$list

## Set output information
bio_Attributes = c('refsnp_id', 'chr_name', 'chrom_start', 'chrom_end', 'allele', 'ensembl_gene_stable_id', 'reg_feature_stable_id')
bio_values = bio_values_VC

## Start query 
annotation_DF = getBM(attributes = bio_Attributes, filters = bio_filters, values = bio_values, mart = ensembl, uniqueRows = TRUE)

## Combine annotations with df
df = left_join(df, annotation_DF, by = c('Position' = 'chrom_start'))

## Drop columns
df = subset(df, select = -c(wt1, ko1, chr_name, chrom_end))

## Write files
output_directory_STR = 'S:/AG/AG-Knauf/Daten/AG Knauf/Gabriel Kirchgraber/GK-01.00-MUGA/GK-01.30-MR42_data_analysis/R_Scripts/Output/'

write.csv2(df, 'S:/AG/AG-Knauf/Daten/AG Knauf/Gabriel Kirchgraber/GK-01.00-MUGA/GK-01.30-MR42_data_analysis/R_Scripts/Output/Compare_SNPs_3WT_7KO.csv', na = '')
write.csv2(GT_DF, 'S:/AG/AG-Knauf/Daten/AG Knauf/Gabriel Kirchgraber/GK-01.00-MUGA/GK-01.30-MR42_data_analysis/R_Scripts/Output/Compare_SNPs_3WT_7KO_Group_mapping.csv', na = '')
write.csv2(track_SNPs_TB, 'S:/AG/AG-Knauf/Daten/AG Knauf/Gabriel Kirchgraber/GK-01.00-MUGA/GK-01.30-MR42_data_analysis/R_Scripts/Output/Compare_SNPs_3WT_7KO_Dropped_SNPs.csv', na = '')





