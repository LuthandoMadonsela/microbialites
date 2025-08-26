#GENERATE 3D NMDS PLOT FROM ALL OTUs: 2019 (FIG. 3B)

#Read in data
abund = read.csv("2019_data.txt", header= TRUE, sep="\t")

#Make community matrix - extract columns with abundance information
nums = abund[,3:ncol(abund)]

#Turn abundance data frame into a matrix
nums_matrix = as.matrix(nums)

#Create 3D plot
set.seed(123)
nmds = metaMDS(nums_matrix, k=3, distance = "bray")
scores(nmds, "sites")

#Extract NMDS scores (x,y and z coordinates)
data.scores = as.data.frame(scores(nmds, "sites"))

#Add columns to data frame 
data.scores$Sample = abund$Sample
data.scores$Location = abund$Location

head(data.scores)

#Load libraries
library(plotly)

#Generate plot
plot_3D <- plot_ly(data.scores, x = ~NMDS1, y = ~NMDS2, z = ~NMDS3, color = ~Location, marker = list(size = 5), text = ~paste('<br>Sample:', Sample, '<br>Location:', Location)) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'NMDS1'),
                      yaxis = list(title = 'NMDS2'),
                      zaxis = list(title = 'NMDS3')))
#Show plot
plot_3D

#ANOSIM ANALYSIS: 2019 (TABLE_S3)

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


#GENERATE 3D NMDS PLOT FROM ALL OTUs: 2019-2022 (FIG. 5C)

#Read in data
abund = read.csv("NMDS_data.txt", header= TRUE, sep="\t")

#Make community matrix - extract columns with abundance information
nums = abund[,4:ncol(abund)]

#Turn abundance data frame into a matrix
nums_matrix = as.matrix(nums)

#Create 3D plot
set.seed(123)
nmds = metaMDS(nums_matrix, k=3, distance = "bray")
scores(nmds, "sites")

#Extract NMDS scores (x,y and z coordinates)
data.scores = as.data.frame(scores(nmds, "sites"))

#Add columns to data frame 
data.scores$Sample = abund$Sample
data.scores$Time = abund$Time
data.scores$Location = abund$Location

head(data.scores)

#Load libraries
library(plotly)

#Generate plot - points coloured by time and shaped by location
plot_3D2 <- plot_ly(data.scores, x = ~NMDS1, y = ~NMDS2, z = ~NMDS3, color = ~Time, symbol = ~Location, symbols = c('circle','diamond','square','cross'), marker = list(size = 5), text = ~paste('<br>Sample:', Sample, '<br>Location:', Location)) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'NMDS1'),
                      yaxis = list(title = 'NMDS2'),
                      zaxis = list(title = 'NMDS3')))
#Show plot
plot_3D2

#ANOSIM ANALYSIS: 2019-2022 (TABLE_S5)
#Calculate pairwise ANOSIM scores, grouped by location

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
    nums = sub_df[,4:ncol(sub_df)]
    nums_matrix = as.matrix(nums)
    set.seed(123)
    ano = anosim(nums_matrix, sub_df$Location, distance = "bray", permutations = 9999)
    Rvalue = ano$statistic
    pvalue = ano$signif
    cat(i,"\t",j ,"\tR-value:\t", Rvalue, "\tP-value:\t", pvalue, '\n')
    
  }
}

#Calculate pairwise ANOSIM scores, grouped by time

#Find unique groups in data
groups = unique(abund[c("Time")])
group_list = as.list(groups$Time)

#Calculate pairwise ANOSIM statistics
for (i in group_list) {
  for (j in group_list) {
    if (i==j) next
    sub_df1 = subset(abund, subset=(Time==i))
    sub_df2 = subset(abund, subset=(Time==j))
    sub_df = rbind(sub_df1, sub_df2)
    nums = sub_df[,4:ncol(sub_df)]
    nums_matrix = as.matrix(nums)
    set.seed(123)
    ano = anosim(nums_matrix, sub_df$Time, distance = "bray", permutations = 9999)
    Rvalue = ano$statistic
    pvalue = ano$signif
    cat(i,"\t",j ,"\tR-value:\t", Rvalue, "\tP-value:\t", pvalue, '\n')
    
  }
}