#! /usr/bin/Rscript

args = commandArgs(trailingOnly=TRUE);

library(ggplot2)
library(readr)

datafile<-read_csv("tempViolin.csv") # should have a column called value and a column called id

datafile$value<-log10(datafile$value+0.2)

datafile$id<-as.factor(datafile$id)

plot<-ggplot(datafile, aes(x=id, y=value, fill=id)) + geom_violin() + labs(y = "Log10 transformed expression levels", x="")

plot<-plot+scale_fill_manual(values=c("royalblue1", "orange1")) # blue = control, orange = case

plot <-plot + geom_boxplot(width=0.05) + stat_summary(fun.y=mean, geom="point", size=2, color="white") + theme_bw()+theme(legend.position = "none")

# plot+coord_cartesian(ylim = c(-2,5))

ggsave(args[1])

#########################################################################################################################################################

barFile1<-read_csv("tempBar1.csv") #should have column 1 as H2B names and column 2 as H2B values with header called values
row.names(barFile1)=barFile1$UID
barFile1$UID=NULL

barColors1<-read_csv("tempColors1.csv") # should have column 1 as H2B names and column 2 as colors with header called colors
row.names(barColors1)=barColors1$UID
barColors1$UID=NULL

barFile1$values<-log10(barFile1$values+1.02)

q1 = ggplot(barFile1) + scale_fill_manual(values = barColors1$colors) + geom_col(aes(x = row.names(barFile1), y = barFile1$values, fill=row.names(barFile1))) + labs(y = "Log10 transformed expression levels", x="")

q1  + theme(axis.text.x = element_text(angle = 90, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none") + coord_cartesian(ylim = c(0,4))


ggsave(args[2])

#########################################################################################################################################################

barFile2<-read_csv("tempBar2.csv") #should have column 1 as H2B names and column 2 as H2B values with header called values
row.names(barFile2)=barFile2$UID
barFile2$UID=NULL

barColors2<-read_csv("tempColors2.csv") # should have column 1 as H2B names and column 2 as colors with header called colors
row.names(barColors2)=barColors2$UID
barColors2$UID=NULL

barFile2$values<-log10(barFile2$values+1.02)

q2 = ggplot(barFile2) + scale_fill_manual(values = barColors2$colors) + geom_col(aes(x = row.names(barFile2), y = barFile2$values, fill=row.names(barFile2))) + labs(y = "Log10 transformed expression levels", x="")

q2  + theme(axis.text.x = element_text(angle = 90, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none") + coord_cartesian(ylim = c(0,4))

ggsave(args[3])