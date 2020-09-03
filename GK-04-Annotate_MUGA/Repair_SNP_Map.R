# Description: R Script to repair wrong information in SNP_Map.
# Details: 88 MT SNPs and 3 PAR SNPs are wrongly labelled to be on Chromosome 0 in the original SNP_Map.txt file. This is corrected and output is a new SNP_Map file named SNP_Map_v2.csv
# Type: Script
library(tidyr)
library(dplyr)
library(janitor)

# Set working directory ####
# setwd("S:/AG/AG-Knauf/Daten/AG Knauf/Gabriel Kirchgraber/GK-01-MUGA/")
setwd("C:/Users/sweet/Dropbox/Arbeit/GK-01-MUGA/")
# Set file Paths ####
data_Neo_dir_STR = "./_Data/Neogene_Files/"

SNP_map_STR = "SNP_Map.txt"
output_STR = "SNP_Map_v2.csv"
consensus_calls_STR = "MiniMUGAv2 Consensus Calls.txt"
annotations_STR = "miniMUGA-Marker-Annotations November 2018.csv"
final_report_STR = "Charite_Berlin_MURCOMV02_20191231_FinalReport.txt"

# Load data ####
SNP_map_file_STR = paste(data_Neo_dir_STR, SNP_map_STR, sep = '')
SNP_map_DF = read.delim(SNP_map_file_STR)

consensus_calls_file_STR = paste(data_Neo_dir_STR, consensus_calls_STR, sep = '')
consensus_DF = read.delim(consensus_calls_file_STR)

annotations_file_STR = paste(data_Neo_dir_STR, annotations_STR, sep = '')
annotations_DF = read.csv(annotations_file_STR)

final_report_file_STR = paste(data_Neo_dir_STR, final_report_STR, sep ='')
final_report_DF = read.delim(final_report_file_STR, skip=9)

# Clean data ####

## The two files provided by Neogene reagarding SNP information include information for different SNPs.
## Combine information from those two files and then map them to the SNPs included in our miniMUGA.
## Check if there is conflicting information in the files.
## Save the final annotated file for later use.

### Drop Index
SNP_map_DF = dplyr::select(SNP_map_DF, !1)

### Only need three columns: Marker, Chromosome, Position
consensus_DF = dplyr::select(consensus_DF, 1:3)
annotations_DF = dplyr::select(annotations_DF, 1:3)

### Only need first column with SNP.Names
final_report_DF = dplyr::select(final_report_DF, 1)
final_report_DF= distinct(final_report_DF)

### Some sanity checks
annotations_DF[!(annotations_DF$Marker %in% final_report_DF$SNP.Name)]
consensus_DF[!(consensus_DF$Marker %in% final_report_DF$SNP.Name)]
sum(final_report_DF$SNP.Name %in% SNP_map_DF$Name)
unique(consensus_DF$Chromosome)
unique(annotations_DF$Chromosome)
unique(SNP_map_DF$Chromosome)

SNP_map_chr_count_DF = SNP_map_DF %>% group_by(Chromosome) %>% summarise(n=n())
annotations_chr_count_DF = annotations_DF %>% group_by(Chromosome) %>% summarise(n=n())
consensus_chr_count_DF = consensus_DF %>% group_by(Chromosome) %>% summarise(n=n())

SNP_map_chr_count_DF = rename(SNP_map_chr_count_DF, n_SNP_map=n)
annotations_chr_count_DF = rename(annotations_chr_count_DF, n_annotations=n)
consensus_chr_count_DF = rename(consensus_chr_count_DF, n_consensus=n)

chr_count_DF = full_join(SNP_map_chr_count_DF, annotations_chr_count_DF, by = 'Chromosome')
chr_count_DF = full_join(chr_count_DF, consensus_chr_count_DF, by = 'Chromosome')
chr_count_DF$Chromosome = factor(chr_count_DF$Chromosome, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "MT", "PAR", "X", "Y", "0"))
chr_count_DF = arrange(chr_count_DF, Chromosome)
chr_count_DF = adorn_totals(chr_count_DF)
### SNP_Map has the same 11125 markers as final_report.
### If you look at chr_count_Df it is clear that some of the info in SNP_Map is not correct. The consensus_calls file has all SNPs that are not constructs.
### This means that 88 MT SNPs and 3 PAR SNPs have misslabelled the Chromsome as 0 in SNP_Map.

# Repair data ####

## There are 88 SNPs on MT in SNP_Map that are falsley labbeled as beeing on Chromosome 0.
SNP_map_v2_DF = SNP_map_DF
SNP_map_v2_DF[SNP_map_v2_DF$Name %in% filter(consensus_DF, Chromosome == 'MT')$Marker,]$Chromosome = 'MT'

## There are 3 SNPs on PAR in SNP_Map that are falsley labbeled as beeing on Chromosome 0. 
SNP_map_v2_DF[SNP_map_v2_DF$Name %in% filter(consensus_DF, Chromosome == 'PAR')$Marker,]$Chromosome = 'PAR'

## Count up new Sample_map and add to previous df
SNP_map_v2_chr_count_DF = SNP_map_v2_DF %>% group_by(Chromosome) %>% summarise(n=n())
SNP_map_v2_chr_count_DF = rename(SNP_map_v2_chr_count_DF, n_SNP_map_v2=n)

chr_count_DF = full_join(chr_count_DF, SNP_map_v2_chr_count_DF, by = 'Chromosome')
chr_count_DF$n_SNP_map_v2[25] = sum(chr_count_DF[1:24,]$n_SNP_map_v2)

# Save data
output_file_STR = paste(data_Neo_dir_STR, output_STR, sep ='')
write.csv(SNP_map_v2_DF, output_file_STR, row.names = FALSE)
