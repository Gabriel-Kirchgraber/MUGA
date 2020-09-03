# Description: R Script to download information for all SNPs from Biomart.
# Details: Add info to SNP_map_v2.csv and save as SNP_map_v3.csv
# Type: Script
library(tidyr)
library(dplyr)
library(biomaRt)

# Set working directory ####
# setwd("S:/AG/AG-Knauf/Daten/AG Knauf/Gabriel Kirchgraber/GK-01-MUGA/")
setwd("C:/Users/sweet/Dropbox/Arbeit/GK-01-MUGA/")
# Set file Paths ####
data_Neo_dir_STR = "./_Data/Neogene_Files/"
output_dir_STR = "./GK-04-Annotate_MUGA/Output/"
SNP_map_STR = "SNP_Map_v2.csv" # Use repaired version
annotations_STR = "annotations.Rdata"

# Load data ####
SNP_map_file_STR = paste(data_Neo_dir_STR, SNP_map_STR, sep = '')
SNP_map_DF = read.csv(SNP_map_file_STR)

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
chunk_size_INT = 25

## create empty dataframe to hold annotations
annotation_DF = data.frame(refsnp_id=character(0),
                           chr_name=character(0),
                           chrom_start=numeric(0),
                           chrom_end=numeric(0),
                           allele=character(0),
                           ensembl_gene_stable_id=character(0),
                           reg_feature_stable_id=character(0))

## Loop to download data
for (i in seq(1,dim(SNP_map_DF)[1],chunk_size_INT)){
  bio_values = c()
  for (j in seq(i, i+chunk_size_INT-1,1)){
    chr_reg_STR = paste(SNP_map_DF[j,2], ":", SNP_map_DF[j,3], ":", SNP_map_DF[j,3], ":1", sep ="")
    bio_values = append(bio_values, chr_reg_STR)
  }
  try({
    df = get_annotation_FUN(bio_filters, bio_Attributes, bio_values, ensembl)
    annotation_DF = rbind(annotation_DF, df)
    print (paste(i,":"))
    print (df)
  })
}

## runtime at home was about 6h. Save dataframe.
# save(annotation_DF, file= paste(output_dir_STR, annotations_STR, sep =''))
load(paste(output_dir_STR, annotations_STR, sep =''))
annotation_DF = annotation_df
rm(annotation_df)

## run again for SNPs with no annotations
### remove duplicates
annotation_DF = annotation_DF %>% distinct(chrom_start, .keep_all = TRUE)

### 8283 SNPs are single nucleotide variations that only concern one base
sum(annotation_DF$chrom_start == annotation_DF$chrom_end)
sum(annotation_DF$chrom_start != annotation_DF$chrom_end)

### Join all those 8283 to SNP_map and create SNP_map_v3
SNP_map_v3_DF = SNP_map_DF
SNP_map_v3_DF = left_join(SNP_map_v3_DF, filter(annotation_DF, chrom_start == chrom_end), by = c('Position'='chrom_start'), keep=TRUE)

### For the 75 SNPs where the SNP concerns more than one base
### Compare them only to the SNPs in SNP_Map that don't have information yet
### If there is a hit add Position to multi_SNP_DF 
no_annot_DF = SNP_map_v3_DF %>% filter(is.na(chrom_end))
multi_SNP_DF = annotation_DF %>% filter(chrom_start != chrom_end)
multi_SNP_DF$chrom_start < multi_SNP_DF$chrom_end
multi_SNP_DF$Position = NA
for (i in seq(1,dim(multi_SNP_DF)[1],1)){
  c_s = multi_SNP_DF[i, 'chrom_start']
  c_e = multi_SNP_DF[i, 'chrom_end']
  c  = multi_SNP_DF[i, 'chr_name']
  hit = no_annot_DF[(no_annot_DF$Chromosome==c)&(between(no_annot_DF$Position, c_s, c_e)),]
  if (dim(hit)[1]==1){
    multi_SNP_DF[i,'Position'] = hit$Position
  }
}
### Add info for multi_SNPs to new SNP_map

setDT(SNP_map_v3_DF)
setDT(multi_SNP_DF)
SNP_map_v3_DF[multi_SNP_DF, on = c("Position"="Position"), ':=' (refsnp_id = i.refsnp_id, 
                                                                     chr_name = i.chr_name, 
                                                                     chrom_start =i.chrom_start,
                                                                     chrom_end = i.chrom_end,
                                                                     allele = i.allele,
                                                                     ensembl_gene_stable_id = i.ensembl_gene_stable_id,
                                                                     reg_feature_stable_id = i.reg_feature_stable_id)]

### Run SNP download again for remaining SNPs without info without constructs (chr == 0)
rm(annotation_DF, check_multi_DF, df, hit, multi_SNP_DF, no_annot_DF, SNP_map_DF)
rm(c, c_e, c_s, chr_reg_STR, chunk_size_INT, failed_chunk_CV,i,j,pos)

## create empty dataframe to hold annotations
annotation_DF = data.frame(refsnp_id=character(0),
                           chr_name=character(0),
                           chrom_start=numeric(0),
                           chrom_end=numeric(0),
                           allele=character(0),
                           ensembl_gene_stable_id=character(0),
                           reg_feature_stable_id=character(0))

SNPs_no_info_DF = filter(SNP_map_v3_DF, (Chromosome != 0)&(is.na(chrom_start))) 
chunk_size_INT = 6
for (i in seq(1,dim(SNPs_no_info_DF)[1],chunk_size_INT)){
  bio_values = c()
  for (j in seq(i, i+chunk_size_INT-1,1)){
    chr_reg_STR = paste(SNPs_no_info_DF[j,2], ":", SNPs_no_info_DF[j,3], ":", SNPs_no_info_DF[j,3], ":1", sep ="")
    bio_values = append(bio_values, chr_reg_STR)
  }
    df = get_annotation_FUN(bio_filters, bio_Attributes, bio_values, ensembl)
    annotation_DF = rbind(annotation_DF, df)
    print (paste(i,":"))
    print (df)
}

save(annotation_DF, file= paste(output_dir_STR, 'annotations_2.Rdata', sep =''))

### remove duplicates
annotation_DF = annotation_DF %>% distinct(chrom_start, .keep_all = TRUE)
setDT(SNP_map_v3_DF)
setDT(annotation_DF)
SNP_map_v3_DF[annotation_DF, on = c("Position"="chrom_start"), ':=' (refsnp_id = i.refsnp_id, 
                                                                 chr_name = i.chr_name, 
                                                                 chrom_start =i.chrom_start,
                                                                 chrom_end = i.chrom_end,
                                                                 allele = i.allele,
                                                                 ensembl_gene_stable_id = i.ensembl_gene_stable_id,
                                                                 reg_feature_stable_id = i.reg_feature_stable_id)]

write.csv(SNP_map_v3_DF, paste(data_Neo_dir_STR, "SNP_Map_v3.csv", sep = ''))


