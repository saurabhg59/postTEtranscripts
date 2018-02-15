#! /usr/bin/Rscript

args = commandArgs(trailingOnly=TRUE);

library(ggplot2)
library(readr)
library(magrittr)
library(ggpubr)


datafile<-read_csv("tempViolin.csv") # should have a column called value and a column called id

datafile$value<-log10(datafile$value+0.2)

datafile$id<-factor(datafile$id,labels = c("Wild-Type L1","Mutant L1","Wild-Type Alu","Mutant Alu"))
#########################################################################################################################################################
# VIOLIN PLOT PART
#########################################################################################################################################################
# plot<-ggplot(datafile, aes(x=id, y=value, fill=id)) + geom_violin() + labs(y = "Log10 transformed expression levels", x="")

# plot<-plot+scale_fill_manual(values=c("royalblue1", "orange1")) # blue = control, orange = case

# plot <-plot + geom_boxplot(width=0.05, outlier.shape = NA) + stat_summary(fun.y=mean, geom="point", size=2, color="white") + theme_bw()+theme(legend.position = "none")

# plot+coord_cartesian(ylim = c(-2,5))
#########################################################################################################################################################
# BOXPLOT PART
#########################################################################################################################################################

df=data.frame(datafile[1:145,1],datafile[146:290,1])

a = ks.test(df$value,df$value.1)

plot1<-ggboxplot(datafile,x="id",y="value",color="black",fill="id",palette=c("royalblue1","orange1","royalblue4","orange4"),shape="id", xlab=FALSE, ylab="Log10 transformed expression levels", show.legend=FALSE)

# plot1+geom_label(aes(x=1.5 , y = 4, label = paste("p-value = ",a$p.value,sep = "")), size = 5)

plot1<-plot1+geom_text(aes(x=1.5 , y = 4, label = paste("P = ",formatC(a$p.value, format = "e", digits = 2),sep = "")), size = 5)+geom_segment(mapping=aes(x=1,y=3.8,xend=2,yend=3.8))+geom_segment(mapping=aes(x=1,y=3.8,xend=1,yend=3.5))+geom_segment(mapping=aes(x=2,y=3.8,xend=2,yend=3.5))

### add part to calculate ks-test values for Alu's and then update the variable in the below line.

plot1<-plot1+geom_text(aes(x=3.5 , y = 0, label = paste("P = ",a$p.value,sep = "")), size = 5)+geom_segment(mapping=aes(x=3,y=0.2,xend=4,yend=0.2))+geom_segment(mapping=aes(x=3,y=0.2,xend=3,yend=0.5))+geom_segment(mapping=aes(x=4,y=0.2,xend=4,yend=0.5))

ggsave(args[1])

#########################################################################################################################################################

barFile<-read_csv("tempBar.csv") #should have column 1 as H2B names and column 2 as H2B values and column 3 as colors

barFile$values<-log10(barFile$values+1.02)

plot2<-ggbarplot(barFile, x = "UID", y = "values",fill = "UID",color = "black",palette = barFile$colors,sort.by.groups = FALSE,x.text.angle = 90,xlab=FALSE,ylab="Log10 transformed expression values",show.legend=FALSE)

ggsave(args[2])

# barFile1<-read_csv("tempBar1.csv") #should have column 1 as H2B names and column 2 as H2B values with header called values
# row.names(barFile1)=barFile1$UID
# barFile1$UID=NULL

# barColors1<-read_csv("tempColors1.csv") # should have column 1 as H2B names and column 2 as colors with header called colors
# row.names(barColors1)=barColors1$UID
# barColors1$UID=NULL

# barFile1$values<-log10(barFile1$values+1.02)

# q1 = ggplot(barFile1) + scale_fill_manual(values = barColors1$colors) + geom_col(aes(x = row.names(barFile1), y = barFile1$values, fill=row.names(barFile1))) + labs(y = "Log10 transformed expression levels", x="")

# q1  + theme(axis.text.x = element_text(angle = 90, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none") + coord_cartesian(ylim = c(0,4))


# ggsave(args[2])

# #########################################################################################################################################################

# barFile2<-read_csv("tempBar2.csv") #should have column 1 as H2B names and column 2 as H2B values with header called values
# row.names(barFile2)=barFile2$UID
# barFile2$UID=NULL

# barColors2<-read_csv("tempColors2.csv") # should have column 1 as H2B names and column 2 as colors with header called colors
# row.names(barColors2)=barColors2$UID
# barColors2$UID=NULL

# barFile2$values<-log10(barFile2$values+1.02)

# q2 = ggplot(barFile2) + scale_fill_manual(values = barColors2$colors) + geom_col(aes(x = row.names(barFile2), y = barFile2$values, fill=row.names(barFile2))) + labs(y = "Log10 transformed expression levels", x="")

# q2  + theme(axis.text.x = element_text(angle = 90, hjust = 1),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position = "none") + coord_cartesian(ylim = c(0,4))

# ggsave(args[3])

#########################################################################################################################################################

##### add KS test commands here, and write the KS test values to a file

# variable<-ks.test(run1_without0_normal_l1_only$`4395.bam.T`,run1_without0_normal_l1_only$`4403.bam.C`,alternative = "less")
# write.table(variable$p.value,file="test.jpeg",quote = F,col.names = F,row.names = F)