# Description: Simulate the gel picture of a restriction digest. 
# Details: Needs as input a .csv with band size, instensity, alpha values etc. Output is a .jpg with the plot.
# Type: Script

library("ggplot2")

# Create a plot of the expected fragemnts sizes as seen on the gel. Data is exported from a list in Restriction_digest_design.xlsx.

# Set working directory
setwd("S:/AG/AG-Knauf/Daten/AG Knauf/Gabriel Kirchgraber/GK-01.00-MUGA")

# Load data
RD_fragment_DF = read.csv2('S:/AG/AG-Knauf/Daten/AG Knauf/Gabriel Kirchgraber/GK-01.00-MUGA/GK-01.10-Restriction_Fragment_Digest/R_Scripts/Data/RD_fragment_size.csv')


# Create Plot
gg1 = ggplot()
gg1 = gg1 + geom_segment(data = RD_fragment_DF,
                         mapping = aes(x=Pos-0.25, y=Band.size, xend=Pos+0.25, yend=Band.size, alpha = Intensity), size = 1.1)
gg1 = gg1 + geom_segment(mapping = aes(x = 0.75, y = 500, xend = 1.25, yend=500), size = 2)
gg1 = gg1 + scale_y_continuous(breaks = seq(0,1000,100))
gg1 = gg1 + scale_x_continuous(breaks = RD_fragment_DF$Pos, labels =RD_fragment_DF$Name)
gg1 = gg1 + theme_minimal()
gg1 = gg1 + theme(axis.text.x = element_text(angle = 90),
                  legend.position = 'none')
gg1 = gg1 + xlab("")
gg1 = gg1 +ylim(c(0,600))
gg1


ggsave(plot = gg1, filename = 'S:/AG/AG-Knauf/Daten/AG Knauf/Gabriel Kirchgraber/GK-01.00-MUGA/GK-01.10-Restriction_Fragment_Digest/RD_expected_bands.jpg')
 


