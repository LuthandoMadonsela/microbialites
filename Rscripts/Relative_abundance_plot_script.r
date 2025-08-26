#GENERATIING BARPLOT OF ALL OTUs FROM ALL SAMPLES: 2019-2022 (FIG. 5A)

#Load libraries
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(vegan)
library(grid)

#Read in data
data <- read.csv(file = "All_OTUs.txt", header=T, sep="\t")

#Convert data frame from a "wide" format to a "long" format
melted_data = melt(data, id = c("Sample"))

#Keep samples ordered the way they're ordered in the data
melted_data$Sample <- factor(melted_data$Sample,levels=unique(melted_data$Sample))

#Generate a palette with distinctive colours
library(randomcoloR)
library(viridis)
n <- 800
palette <- distinctColorPalette(n)


#Generating plot
barplot = ggplot(melted_data, aes(x = Sample, fill = forcats::fct_rev(variable), y = value)) + 
  geom_bar(stat = "identity") + 
  theme(panel.background = element_blank(), axis.text.x = element_text(angle = 90, size = 8, colour = "black", vjust = 0.5, hjust = 1), 
        axis.title.y = element_text(size = 10), legend.title = element_text(size = 16), 
        axis.text.y = element_text(colour = "black", size = 10, face = "bold")) + 
  scale_y_continuous(expand = c(0,0)) + 
  theme(legend.position = "none") +
  scale_fill_manual(values = rep(palette,20))+
  labs(x = "", y = "Relative Abundance (%)", fill = "OTU")

#Show plot
barplot