Chapter 2 R scripts:


#### GENERATING BARPLOT OF ALL PHYLA FROM ALL SAMPLES 

#Read in data
data <- read.csv(file = "OTUs_PHYLA.txt", header=T, sep="\t")

#convert data frame from a "wide" format to a "long" format
melted_data = melt(data, id = c("Sample"))

#Keep samples ordered the way they're ordered in the data
melted_data$Sample <- factor(melted_data$Sample,levels=unique(melted_data$Sample))

#Making a palette with distinctive colours
library(randomcoloR)
n <- 22
palette <- distinctColorPalette(n)

#Generating plot
barplot = ggplot(conserved_melted_data, aes(x = Sample, fill = forcats::fct_rev(variable), y = value)) + 
  geom_bar(stat = "identity") + 
  theme(panel.background = element_blank(), axis.text.x = element_text(angle = 90, size = 7, colour = "black", vjust = 0.5, hjust = 1), 
        axis.title.y = element_text(size = 10), legend.title = element_text(size = 10), 
        axis.text.y = element_text(colour = "black", size = 10, face = "bold")) + 
  scale_y_continuous(expand = c(0,0)) + 
  scale_fill_manual(values = rep(palette,20))+
  labs(x = "", y = "Relative Abundance (%)", fill = "OTU")

#Display plot
barplot
---------------------------------------------------------------------------------------------------------------------------------------
#GENERATING BARPLOT FOR DOMINANT ORDERS and Conserved OTUs

library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(vegan)
library(grid)

#Read in data
data <- read.csv(file = "Orders.txt", header=T, sep="\t")   #### (or Conserved_OTUs)

#convert data frame from a "wide" format to a "long" format
melted_data = melt(data, id = c("Sample"))

#Keep samples ordered the way they're ordered in the data
melted_data$Sample <- factor(melted_data$Sample,levels=unique(melted_data$Sample))

#Making a palette with distinctive colours
n <- 21
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

#Generating plot

barplot = ggplot(melted_data, aes(x = Sample, fill = forcats::fct_rev(variable), y = value)) + 
  geom_bar(stat = "identity", colour = "black") + 
  theme(panel.background = element_blank(),axis.text.x = element_text(angle = 90, size = 14, colour = "black", vjust = 0.5, hjust = 1), 
        axis.title.y = element_text(size = 16), legend.title = element_text(size = 16), 
        axis.text.y = element_text(colour = "black", size = 12, face = "bold")) + 
  scale_y_continuous(expand = c(0,0)) + 
  scale_fill_manual(values = col_vector) +
  labs(x = "", y = "Relative Abundance (%)", fill = "OTU")

barplot
---------------------------------------------------------------------------------------------------------------------------------------
#### GENERATING NMDS PLOT (2D)

#Load libraries
library(vegan)
library(ggplot2)
library(ggrepel)   # optional, for better labels

#Import data
abund = read.csv("ALL_OTUs.txt", header=TRUE, sep="\t")

#Community matrix - abundance columns start at 3
nums = abund[,3:ncol(abund)]
nums_matrix = as.matrix(nums)

#Run NMDS (2D)
set.seed(123)
nmds = metaMDS(nums_matrix, k=2, distance="bray")

#Extract NMDS scores
data.scores = as.data.frame(scores(nmds, "sites"))

#Add metadata
data.scores$Sample = abund$Sample
data.scores$Location = abund$Location

#2D NMDS plot with ggplot2 (no background lines)
plot_2D <- ggplot(data.scores, aes(x=NMDS1, y=NMDS2, color=Location)) +
  geom_point(size=4) +
  geom_text_repel(aes(label=Sample), size=3) +   # sample names, spread out
  theme_classic() +                              # removes grid lines
  labs(title="2D NMDS Plot", x="NMDS1", y="NMDS2")

#Show plot
plot_2D
---------------------------------------------------------------------------------------------------------------------------------------
#### CALCULATE PAIRWISE ANOSIM (GROUPED BY LOCATION) - use ALL_OTUs table

#Find unique groups in data
groups = unique(abund[c("Location")])
group_list = as.list(groups$Location)

#Calculate pairwise ANOSIM statistics
for (i in group_list) {
  for (j in group_list) {
    if (i==j) next
    sub_df1 = subset(abund, subset=(Location==i))
    sub_df2 = subset(abund, subset=(Location==j))
    sub_df = rbind(sub_df1, sub_df2)
    nums = sub_df[,3:ncol(sub_df)]
    nums_matrix = as.matrix(nums)
    set.seed(123)
    ano = anosim(nums_matrix, sub_df$Location, distance = "bray", permutations = 9999)
    Rvalue = ano$statistic
    pvalue = ano$signif
    cat(i,"\t",j ,"\tR-value:\t", Rvalue, "\tP-value:\t", pvalue, '\n')
    
  }
}
---------------------------------------------------------------------------------------------------------------------------------------
#### CALCULATE ANOSIM (MICROBIALITES V WATER)

library(vegan)

# Read in data
abund <- read.csv("ALL_OTUs.txt", header = TRUE, sep = "\t")

# Collapse microbialite locations into one group
abund$Group <- ifelse(abund$Location %in% c("West","South","East"),
                      "Microbialite", "Water")

# Extract only numeric OTU columns
nums <- abund[, 3:ncol(abund)]

# Drop any non-numeric columns automatically
nums <- nums[, sapply(nums, is.numeric)]

# Replace NAs with zeros (common in OTU tables)
nums[is.na(nums)] <- 0

# Convert to matrix
nums_matrix <- as.matrix(nums)

# Run ANOSIM
set.seed(123)
ano <- anosim(nums_matrix, abund$Group, distance = "bray", permutations = 9999)

cat("Microbialites vs Water", "\tR-value:\t", ano$statistic,
    "\tP-value:\t", ano$signif, "\n")
---------------------------------------------------------------------------------------------------------------------------------------
#### GENERATE SHANNON PLOT (DIVERSITY INDEX)

# Load libraries
library(vegan)
library(ggplot2)

# Import abundance table
abund <- read.csv("All_OTUs.txt", header=TRUE, sep="\t")

# Extract abundance columns (metadata in first 2 cols, abund starts at col 3)
nums <- abund[,3:ncol(abund)]
nums_matrix <- as.matrix(nums)

# Calculate Shannon diversity
shan <- diversity(nums_matrix, index="shannon")

# Build dataframe with metadata + Shannon scores
scores_df <- abund[,c(1,2)]   # assumes Sample + Location metadata
scores_df$Shannon <- shan

# Plot Shannon diversity boxplots (grouped by Location)
shan_plot <- ggplot(scores_df, aes(x=Location, y=Shannon, fill=Location)) +
  geom_boxplot(colour="black") +
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(colour="black")) +
  labs(title="Shannon Diversity by Location", y="Shannon Index")

shan_plot
---------------------------------------------------------------------------------------------------------------------------------------
#### CALCULATE SHANNON INDEX (PER SAMPLE)

# Load libraries
library(vegan)

# Import abundance table
abund <- read.csv("All_OTUs.txt", header=TRUE, sep="\t")

# Extract abundance columns (metadata in first 2 cols, abund starts at col 3)
nums <- abund[,3:ncol(abund)]
nums_matrix <- as.matrix(nums)

# Calculate Shannon diversity per sample
shan <- diversity(nums_matrix, index="shannon")

# Build dataframe with metadata + Shannon scores
scores_df <- abund[,c(1,2)]   # assumes Sample + Location metadata
scores_df$Shannon <- shan

# View Shannon index per sample
print(scores_df)
---------------------------------------------------------------------------------------------------------------------------------------
# Spearman's Rank Correlation
# Relationship between nubbin mass and Shannon diversity

# Load data


data <- read.delim("Mass_Shannon.txt",
                   header = TRUE,
                   stringsAsFactors = FALSE)


# Check the structure of the data
str(data)

# Ensure variables are numeric
data$Mass_g <- as.numeric(data$Mass_g)
data$Shannon <- as.numeric(data$Shannon)

# Test for normality (Shapiro-Wilk)

# Mass
mass_normality <- shapiro.test(data$Mass_g)

# Shannon diversity
shannon_normality <- shapiro.test(data$Shannon)

# Display results
cat("\n")
cat("=====================================\n")
cat("Shapiro-Wilk Normality Tests\n")
cat("=====================================\n")

cat("\nMass_g\n")
print(mass_normality)

cat("\nShannon\n")
print(shannon_normality)


# Spearman's Rank Correlation

spearman_result <- cor.test(data$Mass_g,
                            data$Shannon,
                            method = "spearman",
                            exact = FALSE)

cat("\n")
cat("=====================================\n")
cat("Spearman's Rank Correlation\n")
cat("=====================================\n")

print(spearman_result)

# Display only the statistics needed for reporting
cat("\n")
cat("Spearman's rho =", round(spearman_result$estimate, 3), "\n")
cat("P-value =", signif(spearman_result$p.value, 3), "\n")
---------------------------------------------------------------------------------------------------------------------------------------
