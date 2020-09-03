library("dplyr")
library("tidyr")
library("biomaRt")
library("ggplot2")
options(scipen=999)
# Compare the background of BL6N to BL6J

# Set working directory ####
setwd("S:/AG/AG-Knauf/Daten/AG Knauf/Gabriel Kirchgraber/GK-01.00-MUGA/GK-01.30-MR42_data_analysis/")
# setwd("C:/Users/sweet/Documents/MUGA")

# Set File Paths ####
consensus_calls_STR = 'S:/AG/AG-Knauf/Daten/AG Knauf/Gabriel Kirchgraber/GK-01.00-MUGA/.Data/Neogene_Files/MiniMUGAv2 Consensus Calls.txt'
# consensus_calls_STR = 'C:/Users/sweet/Documents/MUGA/DATA/Original/MiniMUGAv2 Consensus Calls.txt'

# Load data ####
## read data into dataframe
consensus_DF = read.delim(consensus_calls_STR)

## rename columns
consensus_DF = rename(consensus_DF, Position = Position..b38.)

## drop unused columns
consensus_DF = subset(consensus_DF, select = c(Marker, Chromosome, Position, C57BL.6NCrl, C57BL.6J))

## Drop SNPs. Use only those from NLPR3 analysis.
load('S:/AG/AG-Knauf/Daten/AG Knauf/Gabriel Kirchgraber/GK-01.00-MUGA/GK-01.20-MR38_data_analysis/R_Scripts/Output/SNP_list_for_NLPR2.rdata')
consensus_DF = filter(consensus_DF, consensus_DF$Marker %in% keep_VC)

## Make all nucleotides uppercase in consensus calls
consensus_DF[, c(4:dim(consensus_DF)[2])] = mutate_all(consensus_DF[, c(4:dim(consensus_DF)[2])], .funs=toupper)

## Clear memory
rm (consensus_calls_STR)

# Clean Data ####
combined_DF = consensus_DF
rm(consensus_DF)

## Replace '-' with N
levels(combined_DF$C57BL.6NCrl) = c(levels(combined_DF$A1), 'N')
levels(combined_DF$C57BL.6J) = c(levels(combined_DF$A2), 'N')
combined_DF = combined_DF %>% mutate(C57BL.6NCrl = replace(C57BL.6NCrl, C57BL.6NCrl == '-', 'N'))
combined_DF = combined_DF %>% mutate(C57BL.6J = replace(C57BL.6J, C57BL.6J == '-', 'N'))

## Drop duplicated markers
combined_DF = distinct(combined_DF, Position, .keep_all = TRUE)

## Drop SNPs where any Allele is 'N'
combined_DF = filter(combined_DF, (C57BL.6J != 'N')&(C57BL.6NCrl != 'N'))

# Compare background alleles ####
combined_DF$Different = combined_DF$C57BL.6J != combined_DF$C57BL.6NCrl
combined_DF$Chromosome = factor(combined_DF$Chromosome, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "X", "Y"))
df = combined_DF %>% group_by(Chromosome) %>% summarise(sum(Different))

# Plot ideogram of SNPs that are different
plot_DF = filter(combined_DF, Different == TRUE)
plot_DF = pivot_longer(plot_DF, values_to = c('Allele'), cols = c('C57BL.6NCrl', 'C57BL.6J'))
order_plot_DF = data.frame(name = unique(plot_DF$name), Order = seq(1,length(unique(plot_DF$name)),1))
plot_DF = left_join(plot_DF, order_plot_DF)

gg1 = ggplot()
gg1 = gg1 + geom_segment(data = plot_DF,
                         mapping = aes(x=Order-0.25, y=Position, colour=Allele, xend = Order+0.25, yend=Position))
gg1 = gg1 + scale_x_continuous(breaks = plot_DF$Order, labels = plot_DF$name)
gg1 = gg1 + facet_wrap(~Chromosome)
gg1 = gg1 + theme_classic()
gg1 = gg1 + theme(axis.text.x = element_text(angle = 90))
gg1

#ggsave(plot = gg1, filename = 'C:/Users/sweet/Dropbox/Arbeit/MUGA/DATA/Processed/Background_comparison_ideogram.jpg', scale = 3)
ggsave(plot=gg1, file = 'S:/AG/AG-Knauf/Daten/AG Knauf/Gabriel Kirchgraber/GK-01.00-MUGA/GK-01.20-MR38_data_analysis/R_Scripts/Plots/Background_comparison_6N_J_ideogram.jpg', scale = 3)

# Export Report ####
write.csv2(combined_DF, 'S:/AG/AG-Knauf/Daten/AG Knauf/Gabriel Kirchgraber/GK-01.00-MUGA/GK-01.20-MR38_data_analysis/R_Scripts/Output/Background_comparison_6N_J_all.csv', na = '')
#write.csv2(combined_DF, 'S:/AG/AG-Knauf/Daten/AG Knauf/Gabriel Kirchgraber/GK-01.00-MUGA/GK-01.20-MR38_data_analysis/R_Scripts/Output/Background_comparison_6N_J_all.csv', na = '')
# write.csv(df, 'C:/Users/sweet/Dropbox/Arbeit/MUGA/DATA/Processed/Background_comparison_all.csv', na = '')
# write.csv(df_wide, 'C:/Users/sweet/Dropbox/Arbeit/MUGA/DATA/Processed/Background_comparison_wide.csv', na = '')
# write.csv(df_filtered, 'C:/Users/sweet/Dropbox/Arbeit/MUGA/DATA/Processed/Background_comparison_filtered.csv', na = '')





