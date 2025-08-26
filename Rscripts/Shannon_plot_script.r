#GENERATING DIVERSITY BOXPLOTS (Fig. 5C)

#Read in data
abund = read.csv("NMDS_data.txt", header= TRUE, sep="\t")

#Make community matrix - extract columns with abundance information
nums = abund[,4:ncol(abund)]

#Turn abundance data frame into a matrix
nums_matrix = as.matrix(nums)

#Calculate shannon diversity
shan = diversity(nums_matrix, index="shannon")
shan.scores = as.data.frame(scores(shan))
colnames(shan.scores) = c('Shannon')

#Create main dataframe with metadata
scores_df <- abund[c(1,2,3)]

#Add in shannon and simpson diversity scores
scores_df$Shannon = shan.scores$Shannon

box_data = melt(scores_df, id = "Location")
sub_box = subset(box_data,variable !='Sample')
sub_box$value = as.numeric(as.character(sub_box$value))

#Plot Shannon diversity box plots
shan_box_plot = ggplot(scores_df, aes(x = Location, y = Shannon, fill = Time)) + 
  geom_boxplot(colour = "black", position = position_dodge(1)) +
  facet_wrap(~Location,ncol = 4, scales = "free_x") +
  theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                    panel.background = element_blank(), axis.line = element_line(colour = "black"))

#Show plot
shan_box_plot