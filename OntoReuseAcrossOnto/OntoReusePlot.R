setwd("/Users/mkamdar/Desktop/PhD/lod_query/OntoReuseAcrossOnto")
library(ggplot2)

df2 <- read.csv("reuseontoc.tsv", sep="\t", header=TRUE)
df2$reusepercent <- as.numeric(as.vector(df2$reusepercent))
ggplot(data=df2, aes(x=reusepercent, y=reusecount, fill=reusetype)) +
  geom_bar(stat="identity", position=position_dodge(), color="black")+ 
  #coord_cartesian(ylim = c(0, 100)) +
  scale_fill_brewer(name="Reuse Type", palette="Reds", direction=-1)+ 
  scale_x_continuous(breaks=c(0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100))+
  scale_y_continuous(breaks=c(0,10,25,50,100,200,300))+
  #scale_y_log10() +
  xlab("% of Terms explicitly reused from other Ontologies") + ylab("Number of Ontologies") + 
  theme_minimal() + theme(legend.justification = c(1.2, 1.2), legend.position = c(1, 1), legend.background = element_rect(colour="black", fill="white")) +
  theme(plot.title = element_text(size = 20, family = "Tahoma", face = "bold"),
        axis.text.x = element_text(size = 24, face="bold"),
        axis.text.y = element_text(size = 24, face="bold"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20, face="bold"),
        axis.title.y = element_text(size = rel(2.0), angle = 90, face="bold", margin = margin(t = 0, r = 20, b = 0, l = 50)), 
        axis.title.x = element_text(size = rel(2.0), angle = 00, face="bold", margin = margin(t = 20, r = 0, b = 50, l = 0)), text = element_text(size=20),
        panel.background = element_rect(colour = "black")) 
