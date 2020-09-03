library("dplyr")
library("tidyr")
library("ggplot2")

# Set working directory
setwd("S:/AG/AG-Knauf/Daten/AG Knauf/Gabriel Kirchgraber/GK-01.00-MUGA")

# Load data
consensus    = read.delim("./.Data/Neogene_Files/MiniMUGAv2 Consensus Calls.txt")
final_report = read.delim("./.Data/Neogene_Files/Charite_Berlin_MURCOMV02_20200608_FinalReport.txt", skip=9)

# Create TRUE/FALSE for heterozigosity column
final_report$het= !(final_report$Allele1...Forward == final_report$Allele2...Forward)

# Replace - with N
final_report$Allele1...Forward = as.character(final_report$Allele1...Forward)
final_report$Allele2...Forward = as.character(final_report$Allele2...Forward)

final_report$Allele1...Forward[final_report$Allele1...Forward=='-']='N'
final_report$Allele2...Forward[final_report$Allele2...Forward=='-']='N'

# Create dataframe with Consensus calls and SNPs
cmpSNP = merge(final_report, consensus, by.x = "SNP.Name", by.y = "Marker" )

# Add Allele column [not longer needed]

# fltr_cmpSNP = filter(cmpSNP, (Sample.ID==192)&(Chromosome==8), (GC.Score>0.15))
# 
# cmpSNP$Allele = NA
# for (row in 1:nrow(cmpSNP)){
#   if (cmpSNP[row,"Allele1...Forward"]==cmpSNP[row,"Allele2...Forward"]){
#     cmpSNP[row, "Allele"]=cmpSNP[row,"Allele1...Forward"]
#   }
#   if (cmpSNP[row,"Allele1...Forward"]!=cmpSNP[row,"Allele2...Forward"]){
#     cmpSNP[row, "Allele"]="H"
#   }
# }

# Create list of SNPs on chr8, note if they are het/hom for Bl6/129)

## Drop unused columns
fltr_cmpSNP = cmpSNP[c(1:4,7,13:14,24, 105)]
## Drop GC.Score==0
fltr_cmpSNP = filter(fltr_cmpSNP, (GC.Score>0))
## Drop where X129 == N or H
fltr_cmpSNP = filter(fltr_cmpSNP, !((X129X1.SvJ=='H')|(X129X1.SvJ=='N')))
## Drop where C57Bl6 == N or H
fltr_cmpSNP = filter(fltr_cmpSNP, !((C57BL.6J=='H')|(C57BL.6J=='N')))
## Only look on chr 8
fltr_cmpSNP = filter(fltr_cmpSNP, (Chromosome==8))
## Drop duplicated markers
fltr_cmpSNP = distinct(fltr_cmpSNP, Sample.ID, Position..b38., .keep_all = TRUE)

## For each row classify parentage according to the following scheme:
## A1  A2	129	BL6	hom/Bl6	hom/129	Bl6or129 het	unexplained
## G   G	A	  G	  x
## G	 G  G	  A		        x
## G	 G	G	  G	                  x
## G	 G	A	  A				                           x
## A	 G	A	  G			                        x
## A	 G	G	  A			                        x
## A	 G	A	  A			            	               x
## A	 G	G	  G				                           x
## A	 G	A	  A				                           x

fltr_cmpSNP$Allele = NA

for (row in 1:nrow(fltr_cmpSNP)){
  A1      = fltr_cmpSNP[row, 'Allele1...Forward']
  A2      = fltr_cmpSNP[row, 'Allele2...Forward']
  X129    = fltr_cmpSNP[row, 'X129X1.SvJ']
  BL6     = fltr_cmpSNP[row, 'C57BL.6J']
  
  if (A1==A2){
    if ((A1==X129)&(A1!=BL6)){
      fltr_cmpSNP[row, 'Allele'] = 'hom.X129'
    }
    else if ((A1==BL6)&(A1!=X129)){
      fltr_cmpSNP[row, 'Allele'] = 'hom.BL6'
    }
    else if ((A1==X129)&(A1==BL6)){
      fltr_cmpSNP[row, 'Allele'] = 'hom.BL6.or.X129'
    }
    else if ((A1!=X129)&(A1!=BL6)){
      fltr_cmpSNP[row, 'Allele'] = 'unexplained'
    }
  }
  else if (A1!=A2){
    if (X129!=BL6){
      fltr_cmpSNP[row, 'Allele'] = 'het'
    }
    else if (BL6!=129){
      fltr_cmpSNP[row, 'Allele'] = 'unexplained'
    }
  }
}

# Plot 

## Drop unwanted columns
plot_cmpSNP = fltr_cmpSNP[c(1:2, 7, 10)]
## df for chromosome outline
rect_cmpSNP = data.frame(X = unique(plot_cmpSNP$Sample.ID))
rect_cmpSNP$Y = 65000000

gg1 = ggplot()
gg1 = gg1 + geom_segment(data = filter(plot_cmpSNP, (Allele != 'hom.BL6')&(Allele != 'hom.BL6.or.X129')),
                         mapping = aes(x=Sample.ID-0.25, y=Position..b38., colour=Allele, xend = Sample.ID+0.25, yend=Position..b38.))
# gg1 = gg1 + geom_segment(data = filter(plot_cmpSNP, Allele == 'hom.X129'), 
#                          mapping = aes(x=Sample.ID-0.25, y=Position..b38., xend = Sample.ID+0.25, yend=Position..b38.), 
#                          alpha = 1, colour = 'red')
# gg1 = gg1 + geom_segment(data = filter(plot_cmpSNP, Allele == 'het'), 
#                          mapping = aes(x=Sample.ID-0.25, y=Position..b38., xend = Sample.ID+0.25, yend=Position..b38.), 
#                          alpha = 1, colour = 'blue')
# gg1 = gg1 + geom_segment(data = filter(plot_cmpSNP, Allele == 'unexplained'), 
#                          mapping = aes(x=Sample.ID-0.25, y=Position..b38., xend = Sample.ID+0.25, yend=Position..b38.), 
#                          alpha = 1, colour = 'green')
gg1 = gg1 + scale_y_continuous(labels = scales::comma, breaks = c(0,25000000, 50000000, 75000000, 100000000, 125000000))
# gg1 = gg1 + scale_y_continuous(labels = scales::comma, breaks = c(80000000, 90000000, 100000000, 110000000, 120000000, 130000000))
gg1 = gg1 + scale_x_continuous(breaks =c(189,190,191,192,193,194,195,196,197,198))
gg1 = gg1 + geom_tile(data  = rect_cmpSNP,mapping = aes(x=X, y = Y), 
                      fill = NA, colour = 'black', size = 1, width = 0.5, height = 130000000)
gg1 = gg1 + labs(x = "Sample ID", y = "Position", title = "Chromosome 8: Primary C57BL6/J, Secondary 129X1/SvJ")
gg1 = gg1 + theme(
  plot.title = element_text(hjust = 0.5)
)
# gg1 = gg1 + ylim(80000000, 111000000)

gg1

ggsave(plot = gg1, filename = './Plots/Background_all.jpg', scale = 2)

# Save SNP classification data
pivot_cmpSNP = fltr_cmpSNP[c(1:2,7,10)]
pivot_cmpSNP = pivot_wider(pivot_cmpSNP, names_from = 'Sample.ID', values_from = 'Allele')
pivot_cmpSNP = pivot_cmpSNP[c( "SNP.Name", "Position..b38.", "189", "190", "191", "192", "193", "194", "195", "196", "197", "198")]         
  
write.csv2(pivot_cmpSNP, './.Analysis/Output/chr8_classifiedSNP.csv')


# Create list of SNPs on chr2, note if they are het/hom for Bl6/129)

## Drop unused columns
fltr_cmpSNP = cmpSNP[c(1:4,7,13:14,24, 105)]
## Drop GC.Score==0
fltr_cmpSNP = filter(fltr_cmpSNP, (GC.Score>0))
## Drop where X129 == N or H
fltr_cmpSNP = filter(fltr_cmpSNP, !((X129X1.SvJ=='H')|(X129X1.SvJ=='N')))
## Drop where C57Bl6 == N or H
fltr_cmpSNP = filter(fltr_cmpSNP, !((C57BL.6J=='H')|(C57BL.6J=='N')))
## Only look on chr 2
fltr_cmpSNP = filter(fltr_cmpSNP, (Chromosome==2))
## Drop duplicated markers
fltr_cmpSNP = distinct(fltr_cmpSNP, Sample.ID, Position..b38., .keep_all = TRUE)

## For each row classify parentage according to the following scheme:
## A1  A2	129	BL6	hom/Bl6	hom/129	Bl6or129 het	unexplained
## G   G	A	  G	  x
## G	 G  G	  A		        x
## G	 G	G	  G	                  x
## G	 G	A	  A				                           x
## A	 G	A	  G			                        x
## A	 G	G	  A			                        x
## A	 G	A	  A			            	               x
## A	 G	G	  G				                           x
## A	 G	A	  A				                           x

fltr_cmpSNP$Allele = NA

for (row in 1:nrow(fltr_cmpSNP)){
  A1      = fltr_cmpSNP[row, 'Allele1...Forward']
  A2      = fltr_cmpSNP[row, 'Allele2...Forward']
  X129    = fltr_cmpSNP[row, 'X129X1.SvJ']
  BL6     = fltr_cmpSNP[row, 'C57BL.6J']
  
  if (A1==A2){
    if ((A1==X129)&(A1!=BL6)){
      fltr_cmpSNP[row, 'Allele'] = 'hom.X129'
    }
    else if ((A1==BL6)&(A1!=X129)){
      fltr_cmpSNP[row, 'Allele'] = 'hom.BL6'
    }
    else if ((A1==X129)&(A1==BL6)){
      fltr_cmpSNP[row, 'Allele'] = 'hom.BL6.or.X129'
    }
    else if ((A1!=X129)&(A1!=BL6)){
      fltr_cmpSNP[row, 'Allele'] = 'unexplained'
    }
  }
  else if (A1!=A2){
    if (X129!=BL6){
      fltr_cmpSNP[row, 'Allele'] = 'het'
    }
    else if (BL6!=129){
      fltr_cmpSNP[row, 'Allele'] = 'unexplained'
    }
  }
}

# Plot 

## Drop unwanted columns
plot_cmpSNP = fltr_cmpSNP[c(1:2, 7, 10)]
## df for chromosome outline
rect_cmpSNP = data.frame(X = unique(plot_cmpSNP$Sample.ID))
rect_cmpSNP$Y = 182113224/2

gg1 = ggplot()
gg1 = gg1 + geom_segment(data = filter(plot_cmpSNP, (Allele == 'hom.BL6')|(Allele == 'hom.BL6.or.X129')),
                         mapping = aes(x=Sample.ID-0.25, y=Position..b38., colour=Allele, xend = Sample.ID+0.25, yend=Position..b38.))
# gg1 = gg1 + geom_segment(data = filter(plot_cmpSNP, (Allele != 'hom.BL6')&(Allele != 'hom.BL6.or.X129')),
#                          mapping = aes(x=Sample.ID-0.25, y=Position..b38., colour=Allele, xend = Sample.ID+0.25, yend=Position..b38.))
# gg1 = gg1 + geom_segment(data = plot_cmpSNP,
#                          mapping = aes(x=Sample.ID-0.25, y=Position..b38., colour=Allele, xend = Sample.ID+0.25, yend=Position..b38.))
# gg1 = gg1 + scale_y_continuous(labels = scales::comma, breaks = c(0,25000000, 50000000, 75000000, 100000000, 125000000, 150000000, 175000000))
gg1 = gg1 + scale_y_continuous(labels = scales::comma, breaks = seq(110000000,140000000, 5000000),
                               limits=(c(109000000, 138500000) ))
gg1 = gg1 + scale_x_continuous(breaks =c(189,190,191,192,193,194,195,196,197,198))
gg1 = gg1 + geom_tile(data  = rect_cmpSNP,mapping = aes(x=X, y = Y),
                      fill = NA, colour = 'black', size = 1, width = 0.5, height = 182113224)
gg1 = gg1 + labs(x = "Sample ID", y = "Position", title = "Chromosome 2: Primary C57BL6/J, Secondary 129X1/SvJ")
gg1 = gg1 + theme(
  plot.title = element_text(hjust = 0.5)
)

gg1

# ggsave(plot = gg1, filename = './Plots/chr2_Background_all.jpg', scale = 2)
# ggsave(plot = gg1, filename = './Plots/chr2_Background_without_hom.jpg', scale = 2)
# ggsave(plot = gg1, filename = './Plots/chr2_Background_only_hom.jpg', scale = 2)

# ggsave(plot = gg1, filename = './Plots/chr2_Background_all_zoom.jpg', scale = 2)
# ggsave(plot = gg1, filename = './Plots/chr2_Background_without_homBl6_zoom.jpg', scale = 2)
#ggsave(plot = gg1, filename = './Plots/chr2_Background_only_hom_zoom.jpg', scale = 2)

# Save SNP classification data
pivot_cmpSNP = fltr_cmpSNP[c(1:2,6,7,10)]
pivot_cmpSNP = pivot_wider(pivot_cmpSNP, names_from = 'Sample.ID', values_from = 'Allele')
pivot_cmpSNP = pivot_cmpSNP[c( "SNP.Name", "Position..b38.", "Chromosome", "189", "190", "191", "192", "193", "194", "195", "196", "197", "198")]         

write.csv2(pivot_cmpSNP, './DATA/Processed/chr2_classifiedSNP.csv')






















