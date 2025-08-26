#GENERATING BARPLOT FOR CONSERVED OTUs: 2019 (FIG.3A)

#Load libraries
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(vegan)
library(grid)

#Read in data
conserved_data <- read.csv(file = "Conserved_OTUs_2019.txt", header=T, sep="\t")

#Convert data frame from a "wide" format to a "long" format
conserved_melted_data = melt(conserved_data, id = c("Sample"))

#Keep samples ordered the way they're ordered in the data
conserved_melted_data$Sample <- factor(conserved_melted_data$Sample,levels=unique(conserved_melted_data$Sample))

#Generating a palette with distinctive colours
library(randomcoloR)
n <- 21
palette <- distinctColorPalette(n)

#Generating plot
conserved_barplot = ggplot(conserved_melted_data, aes(x = Sample, fill = forcats::fct_rev(variable), y = value)) + 
  geom_bar(stat = "identity") + 
  theme(panel.background = element_blank(), axis.text.x = element_text(angle = 90, size = 7, colour = "black", vjust = 0.5, hjust = 1), 
        axis.title.y = element_text(size = 10), legend.title = element_text(size = 10), 
        axis.text.y = element_text(colour = "black", size = 10, face = "bold")) + 
  scale_y_continuous(expand = c(0,0)) + 
  scale_fill_manual(values = rep(palette,20))+
  labs(x = "", y = "Relative Abundance (%)", fill = "OTU")

conserved_barplot

#GENERATING BARPLOT FOR CONSERVED OTUs: 2019-2022 (Fig.6)

#Load libraries
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(vegan)
library(grid)

#Read in data
conserved_data <- read.csv(file = "Conserved_OTUs_All.txt", header=T, sep="\t")

#Convert data frame from a "wide" format to a "long" format
conserved_melted_data = melt(conserved_data, id = c("Sample"))

#Keep samples ordered the way they're ordered in the data
conserved_melted_data$Sample <- factor(conserved_melted_data$Sample,levels=unique(conserved_melted_data$Sample))

#Generating a palette with distinctive colours
library(randomcoloR)
n <- 7
palette <- distinctColorPalette(n)

#Generating plot
conserved_barplot = ggplot(conserved_melted_data, aes(x = Sample, fill = forcats::fct_rev(variable), y = value)) + 
  geom_bar(stat = "identity") + 
  theme(panel.background = element_blank(), axis.text.x = element_text(angle = 90, size = 7, colour = "black", vjust = 0.5, hjust = 1), 
        axis.title.y = element_text(size = 10), legend.title = element_text(size = 10), 
        axis.text.y = element_text(colour = "black", size = 10, face = "bold")) + 
  scale_y_continuous(expand = c(0,0)) + 
  scale_fill_manual(values = rep(palette,20))+
  labs(x = "", y = "Relative Abundance (%)", fill = "OTU")

conserved_barplot