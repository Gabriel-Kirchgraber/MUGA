# Description: R Script to download information for all SNPs from Biomart.
# Details:
# Type: Script
library(tidyr)
library(dplyr)
library(biomaRt)

# Set working directory ####
# setwd("S:/AG/AG-Knauf/Daten/AG Knauf/Gabriel Kirchgraber/GK-01-MUGA/")
setwd("C:/Users/sweet/Dropbox/Arbeit/GK-01-MUGA/")
# Set file Paths ####
output_dir_STR = "./_Data/Mouse_Genome/"
data_Neo_dir_STR = "./_Data/Neogene_Files/"
SNP_map_STR = "SNP_Map.txt"
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
## Check if there are conflicting information in the files.
## Save the final annotated file for later use.

### Only need three columns: Marker, Chromosome, Position
consensus_DF = select(consensus_DF, 1:3)
annotations_DF = select(annotations_DF, 1:3)
SNP_map_DF = select(SNP_map_DF,2:4)

### Only need first column with SNP.Names
final_report_DF = select(final_report_DF, 1)
final_report_DF= distinct(final_report_DF)

### rename columns
# consensus_DF = rename(consensus_DF, PositionConsensus=Position..b38.)
# consensus_DF = rename(consensus_DF, ChromosomeConsensus=Chromosome)
# consensus_DF = rename(consensus_DF, MarkerConsensus=Marker)
# 
# annotations_DF = rename(annotations_DF, PositionAnnotation = Position..b38.)
# annotations_DF = rename(annotations_DF, ChromosomeAnnotation = Chromosome)
# annotations_DF = rename(annotations_DF, MarkerAnnotation = Marker)

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


# Download annotation information ####
## Connect to ensmbl database
ensembl = useEnsembl(biomart = 'snps', dataset = 'mmusculus_snp')

## Function to download annotations
get_annotation_FUN = function(bio_filters, bio_Attributes, bio_values, ensembl){
  df = getBM(attributes = bio_Attributes, filters = bio_filters, values = bio_values, mart = ensembl, uniqueRows = TRUE)
  return(df)
}

## Set up download parameters
bio_filters = c('chromosomal_region')
bio_Attributes = c('refsnp_id', 'chr_name', 'chrom_start', 'chrom_end', 'allele', 'ensembl_gene_stable_id', 'reg_feature_stable_id')

## Because of biomart server problems, only load small chunks of the data at once
chunk_size_INT = 50

## column to document failed chunks
consensus_DF$annotated = NA

## create empty dataframe to hold annotations
annotation_df = data.frame(refsnp_id=character(0),
                           chr_name=character(0),
                           chrom_start=numeric(0),
                           chrom_end=numeric(0),
                           allele=character(0),
                           ensembl_gene_stable_id=character(0),
                           reg_feature_stable_id=character(0))

## Loop to download data
for (i in seq(1,dim(consensus_DF)[1]-10719,chunk_size_INT)){
  bio_values = c()
  for (j in seq(i, i+chunk_size_INT-1,1)){
    chr_reg_STR = paste(consensus_DF[j,2], ":", consensus_DF[j,3], ":", consensus_DF[j,3], ":1", sep ="")
    bio_values = append(bio_values, chr_reg_STR)
  }
  df = get_annotation_FUN(bio_filters, bio_Attributes, bio_values, ensembl)
  print (dim(df)[1])
  annotation_df = rbind(annotation_df, df)
}

bio_values = c('1:6645188:6645188:1', '1:13203373:13203373:1')
a = get_annotation_FUN(bio_filters, bio_Attributes, bio_values, ensembl)




