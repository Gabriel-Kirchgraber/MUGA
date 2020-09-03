library("dplyr")
library("tidyr")
library("biomaRt")
library("ggplot2")
options(scipen=999)
# Compare the background of various BL6 strains from MR38 with each other

# Set working directory ####
setwd("S:/AG/AG-Knauf/Daten/AG Knauf/Gabriel Kirchgraber/GK-01.00-MUGA/GK-01.30-MR42_data_analysis/")
# setwd("C:/Users/sweet/Documents/MUGA")

# Set File Paths ####
final_report_STR = 'S:/AG/AG-Knauf/Daten/AG Knauf/Gabriel Kirchgraber/GK-01.00-MUGA/.Data/Neogene_Files/Charite_Berlin_MURCOMV02_20191231_FinalReport.txt'
consensus_calls_STR = 'S:/AG/AG-Knauf/Daten/AG Knauf/Gabriel Kirchgraber/GK-01.00-MUGA/.Data/Neogene_Files/MiniMUGAv2 Consensus Calls.txt'
# final_report_STR = 'C:/Users/sweet/Documents/MUGA/DATA/Original/Charite_Berlin_MURCOMV02_20191231_FinalReport.txt'
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
consensus_DF = subset(consensus_DF, select = c(Marker, Chromosome, Position, C57BL.6NCrl))

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
combined_DF$C57BL.6NCrl = toupper(combined_DF$C57BL.6NCrl)

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

# Remove samples ####
samples_to_keep_VC = c("C57BL/6N (Chr)","C57BL/6N (Chr)A","C57BL/6N (FEM)","C57BL/6N (FEM)A","C57BL/6N (Janvier)","C57BL/6N (Janvier)A")
combined_DF = filter(combined_DF, Sample.ID %in% samples_to_keep_VC)

# Add information about group Allele. All animals in Group with same nucleotid, get G/A/T/C, when heterozygous get h #### 
# only works if no more Ns in dataset
combined_DF$Allele = NA
combined_DF$Allele[combined_DF$A1 == combined_DF$A2] = filter(combined_DF, A1==A2)$A1
combined_DF$Allele[is.na(combined_DF$Allele)] = 'H'

# Compare sample Allele to background allele ####
df = combined_DF
df$equal_to_background = NA
df$equal_to_background = df$Allele == df$C57BL.6NCrl

# Plot ideogram of SNPs that are different
SNPs_of_interest_VC = unique(filter(df, equal_to_background == FALSE)$SNP.Name)
plot_DF = filter(df, SNP.Name %in% SNPs_of_interest_VC)
order_plot_DF = data.frame(Sample.ID = unique(plot_DF$Sample.ID), Order = seq(1,length(unique(plot_DF$Sample.ID)),1))
plot_DF = left_join(plot_DF, order_plot_DF)
plot_DF$Chromosome = factor(plot_DF$Chromosome, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "X", "Y"))

gg1 = ggplot()
gg1 = gg1 + geom_segment(data = plot_DF,
                         mapping = aes(x=Order-0.25, y=Position, colour=equal_to_background, xend = Order+0.25, yend=Position))
gg1 = gg1 + scale_x_continuous(breaks = plot_DF$Order, labels = plot_DF$Sample.ID)
gg1 = gg1 + scale_y_continuous(breaks = seq(0,150000000,50000000))
gg1 = gg1 + facet_wrap(~Chromosome)
gg1 = gg1 + theme_classic()
gg1 = gg1 + theme(axis.text.x = element_text(angle = 90))
gg1
ggsave(plot = gg1, filename = 'C:/Users/sweet/Dropbox/Arbeit/MUGA/DATA/Processed/Background_comparison_ideogram.jpg', scale = 3)
rm(plot_DF, gg1, order_plot_DF)

# Export Report ####
df_wide = subset(df, select = -c(A1, A2, GC.Score, equal_to_background))
df_wide = filter(df_wide, SNP.Name %in% SNPs_of_interest_VC)
df_wide = pivot_wider(df_wide, names_from = Sample.ID, values_from = Allele)
names(df_wide) = make.names(names(df_wide))
df_wide = df_wide %>% group_by(SNP.Name, Chromosome, Position, C57BL.6NCrl, C57BL.6N..Chr.A, C57BL.6N..Chr., C57BL.6N..Janvier.A, C57BL.6N..Janvier., C57BL.6N..FEM.A, C57BL.6N..FEM.)

df_filtered = filter(df, SNP.Name %in% SNPs_of_interest_VC)


# output_directory_STR = 'S:/AG/AG-Knauf/Daten/AG Knauf/Gabriel Kirchgraber/GK-01.00-MUGA/GK-01.30-MR42_data_analysis/R_Scripts/Output/'

# write.csv2(df, 'S:/AG/AG-Knauf/Daten/AG Knauf/Gabriel Kirchgraber/GK-01.00-MUGA/GK-01.30-MR42_data_analysis/R_Scripts/Output/Compare_SNPs_3WT_7KO.csv', na = '')
# write.csv2(GT_DF, 'S:/AG/AG-Knauf/Daten/AG Knauf/Gabriel Kirchgraber/GK-01.00-MUGA/GK-01.30-MR42_data_analysis/R_Scripts/Output/Compare_SNPs_3WT_7KO_Group_mapping.csv', na = '')
# write.csv2(track_SNPs_TB, 'S:/AG/AG-Knauf/Daten/AG Knauf/Gabriel Kirchgraber/GK-01.00-MUGA/GK-01.30-MR42_data_analysis/R_Scripts/Output/Compare_SNPs_3WT_7KO_Dropped_SNPs.csv', na = '')
write.csv(df, 'C:/Users/sweet/Dropbox/Arbeit/MUGA/DATA/Processed/Background_comparison_all.csv', na = '')
write.csv(df_wide, 'C:/Users/sweet/Dropbox/Arbeit/MUGA/DATA/Processed/Background_comparison_wide.csv', na = '')
write.csv(df_filtered, 'C:/Users/sweet/Dropbox/Arbeit/MUGA/DATA/Processed/Background_comparison_filtered.csv', na = '')





