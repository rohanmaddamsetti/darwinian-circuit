##cellDNAanalysis.R by Rohan Maddamsetti.
##This script analyzes the raw flow cytometry data in .fcs format
##that represents DNA content per cell in REL606 and in REL8593A, which is a
##clone isolated from Ara-1 at 20,000 generations.
##I want to find the fold-increase in DNA content per cell after 20,000
##generations of evolution.

##This script uses packages from Bioconductor that are amenable for flow
##cytometry analysis.

## Load packages.

library(flowCore)
library(flowStats)
library(flowQ)
library(flowViz) # for flow data visualization.
library(ggplot2)

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


##########################
### PICOGREEN ANALYSIS
##########################

## I gated the data by size and fluorescent signal using BDFACSDiva software.
## I collected 30000 events passing these criteria.

pg.606.data1 <- read.FCS("./FACS020515/RM020515/REL606_P5.fcs")
pg.20K.data1 <- read.FCS("./FACS020515/RM020515/REL8593a_P5.fcs")

pg.606.data2 <- read.FCS("./FACS021615/RM021615/REL606_P5.fcs")
pg.20K.data2 <- read.FCS("./FACS021615/RM021615/REL8593a_P5.fcs")

pg.606.data3 <- read.FCS("./FACS021715/RM021715/REL606_P5.fcs")
pg.20K.data3 <- read.FCS("./FACS021715/RM021715/REL8593a_P5.fcs")

pg.606.data4 <- read.FCS("./FACS021815/RM021815/REL606.fcs")
pg.20K.data4 <- read.FCS("./FACS021815/RM021815/REL8593a.fcs")

pg.606.data5 <- read.FCS("./FACS021915/RM021915/REL606_P5.fcs")
pg.20K.data5 <- read.FCS("./FACS021915/RM021915/REL8593a_P5.fcs")

pg.606.data6 <- read.FCS("./FACS022015/RM022015/REL606_P5.fcs")
pg.20K.data6 <- read.FCS("./FACS022015/RM022015/REL8593a_P5.fcs")

##Plot data.

anc1 <- as.vector(exprs(pg.606.data1$"FITC-A")[,1])
evo1 <- as.vector(exprs(pg.20K.data1$"FITC-A")[,1])

anc2 <- as.vector(exprs(pg.606.data2$"FITC-A")[,1])
evo2 <- as.vector(exprs(pg.20K.data2$"FITC-A")[,1])

anc3 <- as.vector(exprs(pg.606.data3$"FITC-A")[,1])
evo3 <- as.vector(exprs(pg.20K.data3$"FITC-A")[,1])

anc4 <- as.vector(exprs(pg.606.data4$"FITC-A")[,1])
evo4 <- as.vector(exprs(pg.20K.data4$"FITC-A")[,1])

anc5 <- as.vector(exprs(pg.606.data5$"FITC-A")[,1])
evo5 <- as.vector(exprs(pg.20K.data5$"FITC-A")[,1])

anc6 <- as.vector(exprs(pg.606.data6$"FITC-A")[,1])
evo6 <- as.vector(exprs(pg.20K.data6$"FITC-A")[,1])


plot.data <- data.frame("Fluorescence"=c(anc1,anc2,anc3,anc4,anc5,anc6,evo1,evo2,evo3,evo4,evo5,evo6),"Replicate"=c(rep("606 1", 30000),rep("606 2", 30000),rep("606 3", 30000), rep("606 4", 30000), rep("606 5", 30000), rep("606 6", 30000), rep("20K 1", 30000), rep("20K 2", 30000),rep("20K 3", 30000),rep("20K 4", 30000),rep("20K 5", 30000),rep("20K 6", 30000)))

p <- ggplot(plot.data, aes(x=Fluorescence, colour=Replicate)) + geom_density()

p1 <- ggplot(subset(plot.data,Replicate=="606 1"|Replicate== "20K 1"), aes(x=Fluorescence, fill=Replicate)) + geom_histogram(binwidth=1000) + theme_bw() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()) + ylab("Cell count") + scale_y_continuous(limits=c(0,2500)) + xlab("Fluorescence (arbitrary units)")

p2 <- ggplot(subset(plot.data,Replicate=="606 2"|Replicate== "20K 2"), aes(x=Fluorescence, fill=Replicate)) + geom_histogram(binwidth=1000) + theme_bw() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()) + ylab("Cell count") + scale_y_continuous(limits=c(0,2500)) + xlab("Fluorescence (arbitrary units)")


p3 <- ggplot(subset(plot.data,Replicate=="606 3"|Replicate== "20K 3"), aes(x=Fluorescence, fill=Replicate)) + geom_histogram(binwidth=1000) + theme_bw() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()) + ylab("Cell count") + scale_y_continuous(limits=c(0,2500)) + xlab("Fluorescence (arbitrary units)")


p4 <- ggplot(subset(plot.data,Replicate=="606 4"|Replicate== "20K 4"), aes(x=Fluorescence, fill=Replicate)) + geom_histogram(binwidth=1000) + theme_bw() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()) + ylab("Cell count") + scale_y_continuous(limits=c(0,2500)) + xlab("Fluorescence (arbitrary units)")


p5 <- ggplot(subset(plot.data,Replicate=="606 5"|Replicate== "20K 5"), aes(x=Fluorescence, fill=Replicate)) + geom_histogram(binwidth=1000) + theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()) + ylab("Cell count") + scale_y_continuous(limits=c(0,2500)) + xlab("Fluorescence (arbitrary units)")


p6 <- ggplot(subset(plot.data,Replicate=="606 6"|Replicate== "20K 6"), aes(x=Fluorescence, fill=Replicate)) + geom_histogram(binwidth=1000) + theme_bw() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()) + ylab("Cell count") + scale_y_continuous(limits=c(0,2500)) + xlab("Fluorescence (arbitrary units)")


multiplot(p1, p2, p3, p4, p5, p6,layout=matrix(c(1,2,3,4,5,6), nrow=3, byrow=TRUE))

## figure out fold increase.

mean.anc.fluor1 <- each_col(pg.606.data1$"FITC-A", mean)
sd.anc.fluor1 <- each_col(pg.606.data1$"FITC-A", sd)
std.err.anc.fluor1 <- sd.anc.fluor1/sqrt(nrow(pg.606.data1))

mean.evo.fluor1 <- each_col(pg.20K.data1$"FITC-A", mean)
sd.evo.fluor1 <- each_col(pg.20K.data1$"FITC-A", sd)
std.err.evo.fluor1 <- sd.evo.fluor1/sqrt(nrow(pg.20K.data1))

mean.anc.fluor2 <- each_col(pg.606.data2$"FITC-A", mean)
sd.anc.fluor2 <- each_col(pg.606.data2$"FITC-A", sd)
std.err.anc.fluor2 <- sd.anc.fluor2/sqrt(nrow(pg.606.data2))

mean.evo.fluor2 <- each_col(pg.20K.data2$"FITC-A", mean)
sd.evo.fluor2 <- each_col(pg.20K.data2$"FITC-A", sd)
std.err.evo.fluor2 <- sd.evo.fluor2/sqrt(nrow(pg.20K.data2))

mean.anc.fluor3 <- each_col(pg.606.data3$"FITC-A", mean)
sd.anc.fluor3 <- each_col(pg.606.data3$"FITC-A", sd)
std.err.anc.fluor3 <- sd.anc.fluor3/sqrt(nrow(pg.606.data3))

mean.evo.fluor3 <- each_col(pg.20K.data3$"FITC-A", mean)
sd.evo.fluor3 <- each_col(pg.20K.data3$"FITC-A", sd)
std.err.evo.fluor3 <- sd.evo.fluor3/sqrt(nrow(pg.20K.data3))

mean.anc.fluor4 <- each_col(pg.606.data4$"FITC-A", mean)
sd.anc.fluor4 <- each_col(pg.606.data4$"FITC-A", sd)
std.err.anc.fluor4 <- sd.anc.fluor4/sqrt(nrow(pg.606.data4))

mean.evo.fluor4 <- each_col(pg.20K.data4$"FITC-A", mean)
sd.evo.fluor4 <- each_col(pg.20K.data4$"FITC-A", sd)
std.err.evo.fluor4 <- sd.evo.fluor4/sqrt(nrow(pg.20K.data4))

mean.anc.fluor5 <- each_col(pg.606.data5$"FITC-A", mean)
sd.anc.fluor5 <- each_col(pg.606.data5$"FITC-A", sd)
std.err.anc.fluor5 <- sd.anc.fluor5/sqrt(nrow(pg.606.data5))

mean.evo.fluor5 <- each_col(pg.20K.data5$"FITC-A", mean)
sd.evo.fluor5 <- each_col(pg.20K.data5$"FITC-A", sd)
std.err.evo.fluor5 <- sd.evo.fluor5/sqrt(nrow(pg.20K.data5))

mean.anc.fluor6 <- each_col(pg.606.data6$"FITC-A", mean)
sd.anc.fluor6 <- each_col(pg.606.data6$"FITC-A", sd)
std.err.anc.fluor6 <- sd.anc.fluor6/sqrt(nrow(pg.606.data6))

mean.evo.fluor6 <- each_col(pg.20K.data6$"FITC-A", mean)
sd.evo.fluor6 <- each_col(pg.20K.data6$"FITC-A", sd)
std.err.evo.fluor6 <- sd.evo.fluor6/sqrt(nrow(pg.20K.data6))

##divide 20K by 606, and calculate the standard error.
fold.increase1 <- mean.evo.fluor1/mean.anc.fluor1
fold.increase2 <- mean.evo.fluor2/mean.anc.fluor2
fold.increase3 <- mean.evo.fluor3/mean.anc.fluor3
fold.increase4 <- mean.evo.fluor4/mean.anc.fluor4
fold.increase5 <- mean.evo.fluor5/mean.anc.fluor5
fold.increase6 <- mean.evo.fluor6/mean.anc.fluor6

folds = c(fold.increase1, fold.increase2, fold.increase3, fold.increase4, fold.increase5, fold.increase6)

## Calculate mean and confidence interval.
fold.estimate <- mean(folds) #1.67 fold increase.
fold.sd <- sd(folds)
n <- length(folds)
error <- qt(0.975,df=n-1)*fold.sd/sqrt(n)
left.conf <- fold.estimate - error  ## left: 1.45
right.conf <- fold.estimate + error ## right: 1.90
