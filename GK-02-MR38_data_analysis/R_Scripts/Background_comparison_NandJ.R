# Description: Compare various groups of mice to a background.
# Details: NandJ
# Type: Script

options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx8192m"))
library("dplyr")
library("tidyr")
library("ggplot2")
library("janitor")
library("xlsx")
options(scipen = 999)

# Set working directory ####
setwd("S:/AG/AG-Knauf/Daten/AG Knauf/Gabriel Kirchgraber/GK-01-MUGA/")
# setwd("C:/Users/sweet/Dropbox/Arbeit/GK-01-MUGA/")

# Set File Paths ####
output_dir_STR = "./GK-02-MR38_data_analysis/R_Scripts/Output/"
plots_dir_STR = "./GK-02-MR38_data_analysis/R_Scripts/Plots/"
data_Neo_dir_STR = "./_Data/Neogene_Files/"

final_report_STR = "Charite_Berlin_MURCOMV02_20191231_FinalReport.txt"
SNP_Map_STR = "SNP_Map_v3.csv"
consensus_STR = "MiniMUGAv2 Consensus Calls.txt"
fun_STR = "MUGA_functions.R"
ideogram_complex_STR = "Background_comparison_NandJ_ideogram_complex.jpg"
ideogram_simple_STR = "Background_comparison_NandJ_ideogram_simple.jpg"
report_STR = "Background_comparison_NandJ.xlsx"

final_report_file_STR = paste(data_Neo_dir_STR, final_report_STR, sep = '')
SNP_map_file_STR = paste(data_Neo_dir_STR, SNP_Map_STR, sep = '')
consensus_calls_file_STR = paste(data_Neo_dir_STR, consensus_STR, sep = '')
fun_file = paste('./', fun_STR, sep = '')
ideogram_complex_file_STR = paste(plots_dir_STR, ideogram_complex_STR, sep = '')
ideogram_simple_file_STR = paste(plots_dir_STR, ideogram_simple_STR, sep = '')
report_file_STR   = paste(output_dir_STR, report_STR, sep = '')

# Load MUGA functions and set parameters ####
source(fun_file)
sample_to_keep_CV = c(
  "C57BL6/J",
  "C57BL6/JA",
  "C57BL6/N",
  "C57BL6/NA"
)
BG_to_keep_CV = c('C57BL.6NCrl')
SNP_to_keep_CV = c()

# Load data ####
temp = load_clean_data_FUN(
  final_report_file_STR = final_report_file_STR,
  consensus_calls_file_STR = consensus_calls_file_STR,
  sample_to_keep_CV = sample_to_keep_CV,
  BG_to_keep_CV = BG_to_keep_CV,
  SNP_map_file_STR = SNP_map_file_STR
)
combined_DF = temp[[1]]
track_SNPs_DF = temp[[2]]
rm(temp)

# Compare sample Allele to background allele ####

combined_DF$equal_to_background = NA
combined_DF$equal_to_background = combined_DF[, 'Allele'] == combined_DF[, BG_to_keep_CV]

# Plot ideogram of SNPs that are different from background in at least one sample
SNPs_of_interest_VC = unique(filter(combined_DF, equal_to_background == FALSE)$SNP.Name)
plot_DF = filter(combined_DF, SNP.Name %in% SNPs_of_interest_VC)
order_plot_DF = data.frame(Sample.ID = unique(plot_DF$Sample.ID),
                           Order = seq(1, length(unique(plot_DF$Sample.ID)), 1))
plot_DF = left_join(plot_DF, order_plot_DF)

# Plotting ####
## Create Sample.ID for background so it can be plotted along with the others
plot_smpid_DF = plot_DF %>% group_by(
  SNP.Name,
  Chromosome,
  Position,
  C57BL.6NCrl,
  refsnp_id,
  ensembl_gene_stable_id,
  reg_feature_stable_id
) %>% summarise()
plot_smpid_DF$Order = length(unique(plot_DF$Sample.ID)) + 1 
plot_smpid_DF$Sample.ID = "Background(C57BL.6NCrl)"
plot_smpid_DF = rename(plot_smpid_DF, 'Allele' = C57BL.6NCrl)

plot_DF = bind_rows(plot_DF, plot_smpid_DF)

# Create complex plot ####
gg1 = ggplot()
gg1 = gg1 + geom_segment(data = plot_DF,mapping = aes(x=Order-0.25, y=Position, color=Allele, xend = Order+0.25, yend=Position))
gg1 = gg1 + geom_text(data = plot_DF, mapping = aes(x=5.8, y=Position, label = refsnp_id), size = 1.5)
gg1 = gg1 + scale_x_continuous(breaks = plot_DF$Order, labels = plot_DF$Sample.ID, limits = c(0,6))
gg1 = gg1 + scale_y_continuous(breaks = seq(0,150000000,50000000))
gg1 = gg1 + facet_wrap(~Chromosome, ncol =3)
gg1 = gg1 + theme_classic()
gg1 = gg1 + theme(axis.text.x = element_text(angle = 90))
gg1 = gg1 + scale_color_manual(values=c("yellow", "blue", "green", "black", "red", "orange"))
ggsave(plot = gg1, filename = ideogram_complex_file_STR, width = 16, height = 40)

# Create simple plot ####
plot_DF_1 = filter(plot_DF, !is.na(equal_to_background))
plot_DF_1$equal_str = as.character("")
plot_DF_1[plot_DF_1$equal_to_background == TRUE,]$equal_str = "TRUE"
plot_DF_1[plot_DF_1$equal_to_background == FALSE,]$equal_str = "FALSE"
plot_DF_1[plot_DF_1$Allele == 'N',]$equal_str = "N"

# Create complex report
gg2 = ggplot()
gg2 = gg2 + geom_segment(data = plot_DF_1, mapping = aes(x=Order-0.25, y=Position, colour=equal_str, xend = Order+0.25, yend=Position))
gg2 = gg2 + scale_x_continuous(breaks = plot_DF$Order, labels = plot_DF$Sample.ID)
gg2 = gg2 + scale_y_continuous(breaks = seq(0,150000000,50000000))
gg2 = gg2 + facet_wrap(~Chromosome, ncol = 3)
gg2 = gg2 + theme_classic()
gg2 = gg2 + theme(axis.text.x = element_text(angle = 90))
ggsave(plot = gg2, filename = ideogram_simple_file_STR,  width = 16, height = 40)

# Export Report ####

df_wide = subset(plot_DF, select = -c(A1, A2, GC.Score, C57BL.6NCrl, equal_to_background, Order))
df_wide = filter(df_wide, SNP.Name %in% SNPs_of_interest_VC)
df_wide = pivot_wider(df_wide, names_from = Sample.ID, values_from = Allele)

wb = createWorkbook()
s1 = createSheet(wb, 'Raw Data')
s2 = createSheet(wb, 'Filtered_Data')
s3 = createSheet(wb, 'Filtered_Data_wide')
s4 = createSheet(wb, 'Dropped_SNPs')
s5 = createSheet(wb, 'Ideogram_simple')
s6 = createSheet(wb, 'Ideogram_complex')
s7 = createSheet(wb, 'Description')

addDataFrame(combined_DF, sheet = s1, row.names = FALSE)
addDataFrame(plot_DF, sheet = s2, row.names = FALSE)
addDataFrame(df_wide, sheet = s3, row.names = TRUE)
addDataFrame(track_SNPs_DF, sheet = s4, row.names = FALSE)

saveWorkbook(wb, report_file_STR)





