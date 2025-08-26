#NUTRIENT CONCENTRATIONS ANALYSIS (FIG.S3)

library(ggplot2)

#Read in data
data <- read.csv(file = "Nutrients.txt", header = TRUE, sep = "\t")

#Create linegraph
ggplot(data, aes(x = Collection_date, y = Concentration, color = Nutrient, group = Nutrient)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = Concentration - SD, ymax = Concentration + SD),width = 0.2) +
  facet_wrap(~Nutrient, nrow = 1, scales = "free_y") + 
  labs(y = "Nutrient Concentration (uM/L)", color = "Nutrient") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.line = element_line(colour = "black"))