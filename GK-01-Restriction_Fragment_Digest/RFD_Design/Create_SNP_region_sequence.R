# Description: R Script to get the context of SNPs and format them to help with primer design
# Details: Out put is chr8_region_SNP_sequences.csv
# Type: Script

library("DECIPHER")
library("seqinr")
library("tidyverse")
Sys.setenv(LANG = "en")

# Create sequences of region surrounding our SNPs of interest with other SNPs as 'N'
# Create sequences for A1 and A2

file_dna         = 'S:/AG/AG-Knauf/Daten/AG Knauf/Gabriel Kirchgraber/Data analysis/Common_datasets/ncbi_dataset/data/GCF_000001635.26/chr8.fna'
file_SNPs        = 'S:/AG/AG-Knauf/Daten/AG Knauf/Gabriel Kirchgraber/Data analysis/MUGA/DATA/processed/chr8_region_SNPs.csv'

df_SNP = read.csv2(file_SNPs)
dna = read.fasta(file_dna)

#Change SNP locations to N
dna_N = dna
dna_N$NC_000074.6[df_SNP$Start] = 'n'

#Fill sequence fields
region_size = 200
df_SNP$Sequence = NA
df_SNP$In_Region = NA
df_SNP$Sequence_A1 = NA
df_SNP$Sequence_A2 = NA

for (row in 1:nrow(df_SNP)){
  if (df_SNP$SNP_of_Interest[row]){
    SNP_position = df_SNP$Start[row]
    df_SNP$Sequence[row] = toupper(paste(dna_N$NC_000074.6[(SNP_position-region_size):(SNP_position+region_size)], collapse = ''))
    
    str = df_SNP$Sequence[row]
    A1 = substr(df_SNP$Alleles[row],1,1)
    A2 = substr(df_SNP$Alleles[row],3,3)
    substr(str, region_size+1, region_size+1) = A1
    df_SNP$Sequence_A1[row] = str
    substr(str, region_size+1, region_size+1) = A2
    df_SNP$Sequence_A2[row] = str
    }
}

write.csv2(df_SNP, 'S:/AG/AG-Knauf/Daten/AG Knauf/Gabriel Kirchgraber/Data analysis/MUGA/DATA/processed/chr8_region_SNP_sequences.csv')






























