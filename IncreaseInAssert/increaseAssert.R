setwd("/Users/mkamdar/Desktop/PhD/lod_query/increaseInAssert")
library(ggplot2)
library(RColorBrewer)

onto_stats <- read.csv("cbound_agg.tsv", sep="\t", header=TRUE) #increase.tsv .. 
df <- onto_stats
is.finite.data.frame <- function(obj){
  sapply(obj,FUN = function(x) all(is.finite(x)))
}

colnames(df) <- c("index", "V1", "V2", "V3")

p10 <- ggplot(df, aes(x = V2, y = V1, fill = V3)) +
  geom_boxplot(alpha=0.7, outlier.size = 5, outlier.shape = 20) +
  scale_y_log10(name = "Total Assertion Count (Log Scale) ---->") +
  scale_x_discrete(name = "Linked Data Source added to PhLeGrA Querying Scheme ---->") +
  #ggtitle("Information Increase (Overload?) on Querying More Data Sources") +
  theme_bw() + theme(legend.justification = c(1.1, -0.1), legend.position = c(1, 0), legend.background = element_rect(colour="black", fill="white")) +
  theme(plot.title = element_text(size = 24, family = "Tahoma", face = "bold"),
        text = element_text(size = 24, family = "Tahoma"),
        axis.title = element_text(size=32,  family = "Verdana"),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 40, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 50)),
        axis.text.x = element_text(size = 24, face="bold"),
        legend.text = element_text(size = 24),
        legend.title = element_text(size = 24, face="bold"),
        axis.text.y = element_text(size = 24, face="bold")) +
  scale_fill_brewer(palette = "Accent") +
  #geom_label(data = outliers, aes(label=index), position = position_jitter(w = 0.3, h=0.05), size=5, hjust = 0) +
  labs(fill = "Entity Type", face="bold") 
p10

