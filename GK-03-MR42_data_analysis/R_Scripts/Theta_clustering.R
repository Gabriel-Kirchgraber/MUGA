library ("argyle")
library("dplyr")
library("tidyr")
library("ggplot2")

setwd("S:/AG/AG-Knauf/Daten/AG Knauf/Gabriel Kirchgraber/Data analysis/MUGA")

load("./DATA/SNP_Map.Rda")
final_report = read.delim("./DATA/Charite_Berlin_MURCOMV02_20191231_FinalReport.txt", skip=9)
final_report = read.delim("./DATA/Charite_Berlin_MURCOMV02_20200608_FinalReport.txt", skip=9)


geno <- read.beadstudio(prefix = 'Charite_Berlin_MURCOMV02_20191231', in.path = "./DATA", snps = SNP_Map, keep.intensity = TRUE)
geno_recode = argyle::recode(geno, "01")

head(geno)
summary(geno)
head(samples(geno))
head(markers(geno))


final_report$het= !(final_report$Allele1...Forward == final_report$Allele2...Forward)
summary(final_report$het)

# df1 = final_report[final_report$Sample.ID=="Casp 1A",]
# df1 = final_report
# df1 = final_report[final_report$SNP.Name == "DIV070020526",]
# df1 = filter(final_report, (het==TRUE)&(Sample.ID == "Casp 1A")&(GC.Score>0.15))
# df2 = filter(final_report, (het==FALSE)&(Sample.ID == "Casp 1A")&(GC.Score>0.15))
# df1 = filter(final_report, (het==TRUE)&(Sample.ID == "Casp 1A"))
# df2 = filter(final_report, (het==FALSE)&(Sample.ID == "Casp 1A"))
# 
df1 = filter(final_report, (het==TRUE)&(GC.Score>0.15))
df2 = filter(final_report, (het==FALSE)&(GC.Score>0.15))
# df1 = filter(final_report, (het==TRUE))
# df2 = filter(final_report, (het==FALSE))

# gg1 = ggplot(df1, aes(x = Theta, y = R, colour = het))
# gg1 = gg1 + geom_point(size = 0.5, alpha = 2)
# #gg1 = gg1 + geom_text(aes(label=Allele2...Forward), vjust = 0, hjust = -2)
# gg1 = gg1 + xlim(0,1) + ylim(0,8)
# gg1


gg1 = ggplot() +
      geom_point(data = df1, mapping = aes(x = Theta, y = R), size = 1, alpha = 1, color = "red") +
      geom_point(data = df2, mapping = aes(x = Theta, y = R), size = 0.5, alpha = 1/5, color = "blue") +
      ylim(0,8) +
      facet_wrap(vars(Sample.ID))

ggsave(filename = "20200608ThetaRClustering_mit_Filter.jpg", plot=gg1, scale = 2, dpi = 600 )









