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
#GENERATING CCA PLOT

library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(vegan)
library(grid)

bv.step <- function(fix.mat, var.mat, 
                    fix.dist.method="bray", var.dist.method="euclidean", correlation.method="spearman",
                    scale.fix=FALSE, scale.var=TRUE,
                    max.rho=0.95,
                    min.delta.rho=0.001,
                    random.selection=TRUE,
                    prop.selected.var=0.2,
                    num.restarts=10,
                    var.always.include=NULL,
                    var.exclude=NULL,
                    output.best=10
){
  
  if(dim(fix.mat)[1] != dim(var.mat)[1]){stop("fixed and variable matrices must have the same number of rows")}
  if(sum(var.always.include %in% var.exclude) > 0){stop("var.always.include and var.exclude share a variable")}
  require(vegan)
  
  if(scale.fix){fix.mat<-scale(fix.mat)}else{fix.mat<-fix.mat}
  if(scale.var){var.mat<-scale(var.mat)}else{var.mat<-var.mat}
  
  fix.dist <- vegdist(as.matrix(fix.mat), method=fix.dist.method)
  
  #an initial removal phase
  var.dist.full <- vegdist(as.matrix(var.mat), method=var.dist.method)
  full.cor <- suppressWarnings(cor.test(fix.dist, var.dist.full, method=correlation.method))$estimate
  var.comb <- combn(1:ncol(var.mat), ncol(var.mat)-1)
  RES <- data.frame(var.excl=rep(NA,ncol(var.comb)), n.var=ncol(var.mat)-1, rho=NA)
  for(i in 1:dim(var.comb)[2]){
    var.dist <- vegdist(as.matrix(var.mat[,var.comb[,i]]), method=var.dist.method)
    temp <- suppressWarnings(cor.test(fix.dist, var.dist, method=correlation.method))
    RES$var.excl[i] <- c(1:ncol(var.mat))[-var.comb[,i]]
    RES$rho[i] <- temp$estimate
  }
  delta.rho <- RES$rho - full.cor
  exclude <- sort(unique(c(RES$var.excl[which(abs(delta.rho) < min.delta.rho)], var.exclude)))
  
  if(random.selection){
    num.restarts=num.restarts
    prop.selected.var=prop.selected.var
    prob<-rep(1,ncol(var.mat))
    if(prop.selected.var< 1){
      prob[exclude]<-0
    }
    n.selected.var <- min(sum(prob),prop.selected.var*dim(var.mat)[2])
  } else {
    num.restarts=1
    prop.selected.var=1  
    prob<-rep(1,ncol(var.mat))
    n.selected.var <- min(sum(prob),prop.selected.var*dim(var.mat)[2])
  }
  
  RES_TOT <- c()
  for(i in 1:num.restarts){
    step=1
    RES <- data.frame(step=step, step.dir="F", var.incl=NA, n.var=0, rho=0)
    attr(RES$step.dir, "levels") <- c("F","B")
    best.comb <- which.max(RES$rho)
    best.rho <- RES$rho[best.comb]
    delta.rho <- Inf
    selected.var <- sort(unique(c(sample(1:dim(var.mat)[2], n.selected.var, prob=prob), var.always.include)))
    while(best.rho < max.rho & delta.rho > min.delta.rho & RES$n.var[best.comb] < length(selected.var)){
      #forward step
      step.dir="F"
      step=step+1
      var.comb <- combn(selected.var, RES$n.var[best.comb]+1, simplify=FALSE)
      if(RES$n.var[best.comb] == 0){
        var.comb.incl<-1:length(var.comb)
      } else {
        var.keep <- as.numeric(unlist(strsplit(RES$var.incl[best.comb], ",")))
        temp <- NA*1:length(var.comb)
        for(j in 1:length(temp)){
          temp[j] <- all(var.keep %in% var.comb[[j]]) 
        }
        var.comb.incl <- which(temp==1)
      }
      
      RES.f <- data.frame(step=rep(step, length(var.comb.incl)), step.dir=step.dir, var.incl=NA, n.var=RES$n.var[best.comb]+1, rho=NA)
      for(f in 1:length(var.comb.incl)){
        var.incl <- var.comb[[var.comb.incl[f]]]
        var.incl <- var.incl[order(var.incl)]
        var.dist <- vegdist(as.matrix(var.mat[,var.incl]), method=var.dist.method)
        temp <- suppressWarnings(cor.test(fix.dist, var.dist, method=correlation.method))
        RES.f$var.incl[f] <- paste(var.incl, collapse=",")
        RES.f$rho[f] <- temp$estimate
      }
      
      last.F <- max(which(RES$step.dir=="F"))
      RES <- rbind(RES, RES.f[which.max(RES.f$rho),])
      best.comb <- which.max(RES$rho)
      delta.rho <- RES$rho[best.comb] - best.rho 
      best.rho <- RES$rho[best.comb]
      
      if(best.comb == step){
        while(best.comb == step & RES$n.var[best.comb] > 1){
          #backward step
          step.dir="B"
          step <- step+1
          var.keep <- as.numeric(unlist(strsplit(RES$var.incl[best.comb], ",")))
          var.comb <- combn(var.keep, RES$n.var[best.comb]-1, simplify=FALSE)
          RES.b <- data.frame(step=rep(step, length(var.comb)), step.dir=step.dir, var.incl=NA, n.var=RES$n.var[best.comb]-1, rho=NA)
          for(b in 1:length(var.comb)){
            var.incl <- var.comb[[b]]
            var.incl <- var.incl[order(var.incl)]
            var.dist <- vegdist(as.matrix(var.mat[,var.incl]), method=var.dist.method)
            temp <- suppressWarnings(cor.test(fix.dist, var.dist, method=correlation.method))
            RES.b$var.incl[b] <- paste(var.incl, collapse=",")
            RES.b$rho[b] <- temp$estimate
          }
          RES <- rbind(RES, RES.b[which.max(RES.b$rho),])
          best.comb <- which.max(RES$rho)
          best.rho<- RES$rho[best.comb]
        }
      } else {
        break()
      }
      
    }
    
    RES_TOT <- rbind(RES_TOT, RES[2:dim(RES)[1],])
    print(paste(round((i/num.restarts)*100,3), "% finished"))
  }
  
  RES_TOT <- unique(RES_TOT[,3:5])
  
  
  if(dim(RES_TOT)[1] > output.best){
    order.by.best <- RES_TOT[order(RES_TOT$rho, decreasing=TRUE)[1:output.best],]
  } else {
    order.by.best <-  RES_TOT[order(RES_TOT$rho, decreasing=TRUE), ]
  }
  rownames(order.by.best)<-NULL
  
  order.by.i.comb <- c()
  for(i in 1:length(selected.var)){
    f1 <- which(RES_TOT$n.var==i)
    f2 <- which.max(RES_TOT$rho[f1])
    order.by.i.comb <- rbind(order.by.i.comb, RES_TOT[f1[f2],])
  }
  rownames(order.by.i.comb)<-NULL
  
  if(length(exclude)<1){var.exclude=NULL} else {var.exclude=exclude}
  out <- list(
    order.by.best=order.by.best,
    order.by.i.comb=order.by.i.comb,
    best.model.vars=paste(colnames(var.mat)[as.numeric(unlist(strsplit(order.by.best$var.incl[1], ",")))], collapse=","),
    best.model.rho=order.by.best$rho[1],
    var.always.include=var.always.include,
    var.exclude=var.exclude
  )
  out
  
}

abund_table = read.csv("2019_CCA_spe.txt", row.names=1, header= TRUE, sep="\t")
#Transpose the data to have sample names on rows
abund_table<-t(abund_table)
meta_table = read.csv("2019_CCA_env.txt", row.names=1, header= TRUE, sep="\t")
#Get grouping information
grouping_info = read.csv("2019_grouping_info.txt", row.names=1, header= TRUE, sep="\t")
head(grouping_info)

#Parameters
cmethod<-"pearson" #Correlation method to use: pearson, pearman, kendall
fmethod<-"bray" #Fixed distance method: euclidean, manhattan, gower, altGower, canberra, bray, kulczynski, morisita,horn, binomial, and cao
vmethod<-"bray" #Variable distance method: euclidean, manhattan, gower, altGower, canberra, bray, kulczynski, morisita,horn, binomial, and cao
nmethod<-"bray" #NMDS distance method:  euclidean, manhattan, gower, altGower, canberra, bray, kulczynski, morisita,horn, binomial, and cao

res.bv.step.biobio <- bv.step(wisconsin(abund_table), wisconsin(abund_table), 
                              fix.dist.method=fmethod, var.dist.method=vmethod,correlation.method=cmethod,
                              scale.fix=FALSE, scale.var=FALSE, 
                              max.rho=0.95, min.delta.rho=0.001,
                              random.selection=TRUE,
                              prop.selected.var=0.3,
                              num.restarts=10,
                              output.best=10,
                              var.always.include=NULL) 

#Get the 10 best subset of taxa
taxaNames<-colnames(abund_table)
bestTaxaFit<-""
for(i in (1:length(res.bv.step.biobio$order.by.best$var.incl)))
{
  bestTaxaFit[i]<-paste(paste(taxaNames[as.numeric(unlist(strsplit(res.bv.step.biobio$order.by.best$var.incl[i], split=",")))],collapse=' + '), " = ",res.bv.step.biobio$order.by.best$rho[i],sep="")
}
bestTaxaFit<-data.frame(bestTaxaFit)
colnames(bestTaxaFit)<-"Best combination of taxa with similarity score"

bestTaxaFit

#Generate NMDS plot
MDS_res=metaMDS(abund_table, distance = nmethod, k = 2, trymax = 50)

bio.keep <- as.numeric(
  unlist(
    strsplit(
      res.bv.step.biobio$order.by.best$var.incl[2],
      ","
    )
  )
)

# Check selected taxa
print(colnames(abund_table)[bio.keep])

# Fit selected taxa onto NMDS
bio.fit <- envfit(
  MDS_res,
  abund_table[, bio.keep, drop = FALSE],
  perm = 999
)

# View envfit statistics
bio.fit

# Extract site scores
df <- scores(MDS_res, display = "sites")

# Add grouping information
df <- data.frame(
  df,
  Type = grouping_info[rownames(df), 1]
)

# Extract vector scores
df_biofit <- scores(
  bio.fit,
  display = "vectors"
)

# Scale arrows
df_biofit <- df_biofit * vegan:::ordiArrowMul(df_biofit)

df_biofit <- as.data.frame(df_biofit)

# Plot NMDS with taxa vectors
fitted_NMDS <- ggplot() +
  
  geom_point(
    data = df,
    aes(NMDS1, NMDS2, colour = Type),
    size = 3
  ) +
  
  geom_segment(
    data = df_biofit,
    aes(
      x = 0,
      y = 0,
      xend = NMDS1,
      yend = NMDS2
    ),
    arrow = arrow(length = unit(0.2, "cm")),
    colour = "#808080"
  ) +
  
  geom_text(
    data = as.data.frame(df_biofit * 1.1),
    aes(
      NMDS1,
      NMDS2,
      label = rownames(df_biofit)
    ),
    colour = "black",
    size = 4
  ) +
  
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  )

# Display plot
fitted_NMDS
