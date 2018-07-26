setwd("/Users/mkamdar/Desktop/PhD/lod_query/ReuseAcrossLD")
library(ggplot2)
library(RColorBrewer)
library(grid)
library(gtable)

pmin_pmax_clip <- function(x, a, b) 0.7
onto_stats <- read.csv("compSet_df_a.tsv", sep="\t", header=TRUE)
obo_stats <- onto_stats[onto_stats["type"] == "establishedOnto",]
obon_stats <- onto_stats[onto_stats["type"] == "establishedOnto-N",]
umls_stats <- onto_stats[onto_stats["type"] == "establishedVocab",]
other_stats <- onto_stats[onto_stats["type"] == "unpublishedVocab",]
sel_set1 <- umls_stats

p1 <- ggplot()+
  geom_point(data=umls_stats, aes(x=dsets, y=clcount), size=umls_stats$clSize*6, na.rm=TRUE, fill = "blue", colour = "black", stroke = 1, pch=21, alpha = 0.5) +
  geom_point(data=obo_stats, aes(x=dsets, y=clcount), size=obo_stats$clSize*8, na.rm=TRUE, fill = "yellow", colour = "black", stroke = 1, pch=22, alpha = 0.5) +
  geom_point(data=obon_stats, aes(x=dsets, y=clcount), size=obon_stats$clSize*8, na.rm=TRUE, fill = "orange", colour = "black", stroke = 1, pch=22, alpha = 0.5) +
  geom_point(data=other_stats, aes(x=dsets, y=clcount), position = position_jitter(w = 0.1, h = 0.1), size=other_stats$clSize*4, na.rm=TRUE, fill = "green", colour = "black", stroke = 1, pch=23, alpha = 0.5) +
  #geom_text(data=onto_stats, aes(x=onto_stats$ui/onto_stats$size, y=onto_stats$api/onto_stats$size, label = onto_stats$index), size=4, hjust = 0, nudge_x = 0.01, check_overlap = TRUE) +
  geom_label(data = obo_stats, aes(x=obo_stats$dsets, y=obo_stats$clcount, label=obo_stats$index, fill="red"), size=4, colour = "white", fontface = "bold", hjust = 0, nudge_x = 0.01, alpha = 0.8) +
  geom_label(data = obon_stats, aes(x=obon_stats$dsets, y=obon_stats$clcount, label=obon_stats$index, fill="red"), size=4, colour = "white", fontface = "bold", hjust = 0, nudge_x = 0.01, alpha = 0.8) +
  geom_label(data = sel_set1, aes(x=sel_set1$dsets, y=sel_set1$clcount, label=sel_set1$index), size=4, hjust = 0, nudge_x = 0.01, alpha = 0.5, fontface = "bold") +
  xlab("Number of Linked Data Sources using the Ontology/Vocabulary (Log Scale)") + ylab("Schema Elements Reused (Log Scale)") + 
  theme_bw() + scale_y_log10() + scale_x_log10() +
  ggtitle("Ontologies and Vocabularies reused across Biomedical Linked Data sources") +
  #geom_label(data = outliers, aes(label=index), position = position_jitter(w = 0.3, h=0.05), size=5, hjust = 0) +
  theme(plot.title = element_text(size = 20, family = "Tahoma", face = "bold"),
        axis.title.y = element_text(size = rel(1.0), angle = 90, face="bold"), 
        axis.title.x = element_text(size = rel(1.0), angle = 00, face="bold"), text = element_text(size=20),
        panel.background = element_rect(colour = "black"), legend.position = "none")

p1

df <- read.csv("plot_dsethist_a.tsv", sep="\t", header=TRUE)
#rep("c", dim(df)[1]))
p2 <- ggplot(df, aes(x=dsetnumber, y=sharedschema, fill="Reds"))+
  geom_bar(stat="identity", color="black", fill="sienna1")  + scale_y_log10() +
  xlab("Sharing Linked Data Sources") + ylab("Schema Elements (Log Scale)") + 
  theme_minimal() + theme(legend.position="none") +
  #scale_fill_brewer(palette="Reds", direction=-1)+
  #ggtitle("Number of Linked Data Sources that Share a Set of Schema Elements") +
  theme(plot.title = element_text(size = 20, family = "Tahoma", face = "bold"),
        axis.text.x = element_text(size = rel(1.0), face="bold"),
        axis.text.y = element_text(size = rel(1.0), face="bold"),
        axis.title.y = element_text(size = rel(1.0), angle = 90, face="bold", margin = margin(t = 0, r = 20, b = 0, l = 50)), 
        axis.title.x = element_text(size = rel(1.0), angle = 00, face="bold", margin = margin(t = 20, r = 0, b = 50, l = 0)), text = element_text(size=20),
        panel.background = element_rect(colour = "black"), legend.position = "none")

p2
#size 24

gt1 <- ggplot_gtable(ggplot_build(p1))
gt2 <- ggplot_gtable(ggplot_build(p2))

maxWidth <- unit.pmax(gt1$widths[2:3], gt2$widths[2:3])

# Set the maximums in the gtables for gt1, gt2 and gt3
gt1$widths[2:3] <- as.list(maxWidth)
gt2$widths[2:3] <- as.list(maxWidth)

# Combine the scatterplot with the two marginal boxplots
# Create a new gtable
gt <- gtable(widths = unit(c(7, 2), "null"), height = unit(c(3, 7), "null"))

# Instert gt1, gt2 and gt3 into the new gtable
gt <- gtable_add_grob(gt, gt1, 2, 1)
gt <- gtable_add_grob(gt, gt2, 1, 1)

# And render the plot
grid.newpage()
grid.draw(gt)