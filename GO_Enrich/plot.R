library(ggplot2);
library(reshape2);
library(grid);
mf <- read.table("data/mf.ggplot.tab",header=T,sep="\t");

#bp <- read.table("bp.ggplot.tab",header=T,sep="\t");
bp1 <- read.table("data/bp1.ggplot.tab",header=T,sep="\t");
bp2 <- read.table("data/bp2.ggplot.tab",header=T,sep="\t");
cc <- read.table("cc.ggplot.tab",header=T,sep="\t");

pdf("GO_plots_percent.pdf")

#MF
mf$Samples <- factor(mf$Samples, levels = c("Rozella", "Ecun","Cneo",
                                   "Ncra", "Anid", "Spom", "Scer","Calb"))

cc$Samples <- factor(cc$Samples, levels = c("Rozella", "Ecun","Cneo",
                                   "Ncra", "Anid", "Spom", "Scer","Calb"))

bp1$Samples <- factor(bp1$Samples, levels = c("Rozella", "Ecun","Cneo",
                                   "Ncra", "Anid", "Spom", "Scer","Calb"))

bp2$Samples <- factor(bp2$Samples, levels = c("Rozella", "Ecun","Cneo",
                                   "Ncra", "Anid", "Spom", "Scer","Calb"))

p <- ggplot(mf, aes(GO_Term, Percentage, fill = Samples)) +
  geom_bar(position="dodge",stat="identity") + coord_flip() +  
  xlab(NULL) + #ylab(eNULL)
  labs(title="Molecular Function") + #scale_fill_brewer(palette="Dark2") +
  theme_bw() +
  theme(axis.text.y=element_text(angle=0,size=10,hjust=0),
        plot.title=element_text(lineheight=.8, face="bold"),
        legend.position = c(-1,-0.05),# "bottom"
        legend.key.size = unit(0.25,"cm"),
        panel.grid.major = element_line(colour = NA), 
        panel.grid.minor = element_line(colour = NA), 
        legend.direction="horizontal",legend.title=element_blank(),
        plot.title=element_text(lineheight=.8, face="bold"))

print(p)

#legend.position = c(-1,-0.02),# "bottom"
#        legend.key.size = unit(0.25,"cm"),
#        legend.direction="horizontal",legend.title=element_blank(),
#BP1
p <- ggplot(bp1, aes(GO_Term, Percentage, fill = Samples)) +
  geom_bar(position="dodge",stat="identity") + coord_flip() +
  xlab(NULL) + #ylab(NULL)
  theme_bw() +
  labs(title="Biological Process") +
  theme(axis.text.y=element_text(angle=0,size=10,hjust=0),
        legend.position = c(-0.3,-0.06),# "bottom"
        legend.key.size = unit(0.25,"cm"),
        panel.grid.major = element_line(colour = NA), 
        panel.grid.minor = element_line(colour = NA), 
        legend.direction="horizontal",legend.title=element_blank(),
        plot.title=element_text(lineheight=.8, face="bold"))

print(p)


#BP2 
p <- ggplot(bp2, aes(GO_Term, Percentage, fill = Samples)) +
  geom_bar(position="dodge",stat="identity") + coord_flip() +
  xlab(NULL) + #ylab(NULL)
  theme_bw() +
  theme(legend.position = c(-0.6,-0.06),# "bottom"
        legend.key.size = unit(0.25,"cm"),
        legend.direction="horizontal",legend.title=element_blank(),
        axis.text.y=element_text(angle=0,size=10,hjust=0),
        panel.grid.major = element_line(colour = NA), 
        panel.grid.minor = element_line(colour = NA), 
        plot.title=element_text(lineheight=.8, face="bold")) +
  labs(title="Biological Process") 

print(p)

#CC
p <- ggplot(cc, aes(GO_Term,Percentage, fill = Samples)) +
  geom_bar(position="dodge",stat="identity") + coord_flip() +
  xlab(NULL) + #ylab(NULL)
  theme_bw() +
  theme(#legend.position =  "bottom",
        legend.key.size = unit(0.25,"cm"),
        legend.position = c(-0.15,-0.06),# "bottom"
        legend.direction="horizontal",legend.title=element_blank(),
        axis.text.y=element_text(angle=0,size=10,hjust=0),
        panel.grid.major = element_line(colour = NA), 
        panel.grid.minor = element_line(colour = NA), 
        plot.title=element_text(lineheight=.8, face="bold")  ) +
  labs(title="Cellular Compartment")
print(p)

pdf("GO_plots_count.pdf")

#MF stacked
#p <- ggplot(mf, aes(GO_Term, Percentage, fill = Samples)) +
#  geom_bar(stat="identity") + coord_flip() +
#  xlab(NULL) + #ylab(NULL)
#  theme_bw() +
#  theme(legend.position = c(-1,-0.02),# "bottom"
#        legend.key.size = unit(0.25,"cm"),
#        legend.direction="horizontal",legend.title=element_blank(),
#        axis.text.y=element_text(angle=0,size=10,hjust=0),
#        plot.title=element_text(lineheight=.8, face="bold")) +
#  labs(title="Molecular Function")
#print(p)

#BP1 stacked
## p <- ggplot(bp1, aes(GO_Term, Count, fill = Samples)) +
##   geom_bar(stat="identity") + coord_flip() +
##   xlab(NULL) + #ylab(NULL)
##   theme(axis.text.y=element_text(angle=0,size=10,hjust=0)) +
##   labs(title="Biological Process") +
##   theme(legend.position = c(-0.3,-0.05),# "bottom"
##         legend.key.size = unit(0.25,"cm"),
##         legend.direction="horizontal",legend.title=element_blank()) +
##   theme(plot.title=element_text(lineheight=.8, face="bold"))

## print(p)

## #BP2 stacked
## p <- ggplot(bp2, aes(GO_Term, Count, fill = Samples)) +
##   geom_bar(stat="identity") + coord_flip() +
##   xlab(NULL) + #ylab(NULL)
##   theme(legend.position = c(-0.6,-0.05),# "bottom"
##         legend.key.size = unit(0.25,"cm"),
##         legend.direction="horizontal",legend.title=element_blank()) +
##   theme(axis.text.y=element_text(angle=0,size=10,hjust=0)) +
##   labs(title="Biological Process") +
##   theme(plot.title=element_text(lineheight=.8, face="bold"))

## print(p)



## #CC stacked
## p <- ggplot(cc, aes(GO_Term, Count, fill = Samples)) +
##   geom_bar(stat="identity") + coord_flip() +
##   xlab(NULL) + #ylab(NULL)
##   theme(legend.position = "bottom",
##         legend.key.size = unit(0.25,"cm"),
##         legend.direction="horizontal",legend.title=element_blank()) +
##   theme(axis.text.y=element_text(angle=0,size=10,hjust=0)) +
##   ggtitle("Cellular Compartment") +
##   theme(plot.title=element_text(lineheight=.8, face="bold"))

## print(p)
