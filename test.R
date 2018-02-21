#! /usr/bin/Rscript

args = commandArgs(trailingOnly=TRUE);

options(warn=-1)

library(ggplot2)
library(readr)
library(magrittr)
library(ggpubr)


suppressMessages(datafile<-read_csv("tempViolin.csv")) # should have a column called value and a column called id

df=data.frame(datafile[1:145,1],datafile[146:290,1])

aluDF=data.frame(datafile[291:331,1],datafile[332:372,1])

a = ks.test(df$value,df$value.1,alternative="l")

b = ks.test(aluDF$value,aluDF$value.1,alternative="l")

# l1WilcoxTwosidedBeforeLog = wilcox.test(df$value-df$value.1,alternative="two.sided")

# aluWilcoxTwosidedBeforeLog = wilcox.test(aluDF$value-aluDF$value.1,alternative="two.sided")

# l1WilcoxLessBeforeLog = wilcox.test(df$value-df$value.1,alternative="less")

# aluWilcoxLessBeforeLog = wilcox.test(aluDF$value-aluDF$value.1,alternative="less")

# l1WilcoxGreaterBeforeLog = wilcox.test(df$value-df$value.1,alternative="greater")

# aluWilcoxGreaterBeforeLog = wilcox.test(aluDF$value-aluDF$value.1,alternative="greater")

# l1KsLessBeforeLog = ks.test(df$value,df$value.1,alternative="less")

# aluKsLessBeforeLog = ks.test(aluDF$value,aluDF$value.1,alternative="less")

# l1KsGreaterBeforeLog = ks.test(df$value,df$value.1,alternative="greater")

# aluKsGreaterBeforeLog = ks.test(aluDF$value,aluDF$value.1,alternative="greater")

datafile$value<-log10(datafile$value+0.2)

# l1KsTwosidedAfterLog = ks.test(df$value,df$value.1,alternative="two.sided")

# aluKsTwosidedAfterLog = ks.test(aluDF$value,aluDF$value.1,alternative="two.sided")

# l1WilcoxTwosidedAfterLog = wilcox.test(df$value-df$value.1,alternative="two.sided")

# aluWilcoxTwosidedAfterLog = wilcox.test(aluDF$value-aluDF$value.1,alternative="two.sided")

# l1WilcoxLessAfterLog = wilcox.test(df$value-df$value.1,alternative="less")

# aluWilcoxLessAfterLog = wilcox.test(aluDF$value-aluDF$value.1,alternative="less")

# l1WilcoxGreaterAfterLog = wilcox.test(df$value-df$value.1,alternative="greater")

# aluWilcoxGreaterAfterLog = wilcox.test(aluDF$value-aluDF$value.1,alternative="greater")

# l1KsLessAfterLog = ks.test(df$value,df$value.1,alternative="less")

# aluKsLessAfterLog = ks.test(aluDF$value,aluDF$value.1,alternative="less")

# l1KsGreaterAfterLog = ks.test(df$value,df$value.1,alternative="greater")

# aluKsGreaterAfterLog = ks.test(aluDF$value,aluDF$value.1,alternative="greater")

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

# plot1<-ggboxplot(datafile,x="id",y="value",color="black",fill="id",palette=c("royalblue1","orange1","royalblue4","orange4"),shape="id", xlab=FALSE, ylab="Log10 transformed expression levels", show.legend=FALSE)
plot1<-ggboxplot(datafile,x="id",y="value",color="black",fill="id",palette=c("royalblue1","orange1","royalblue4","orange4"),shape="id", xlab=FALSE, ylab=FALSE, show.legend=FALSE)
if(a$p.value>=0.0001 && a$p.value<=1){
	plot1<-plot1+geom_text(aes(x=1.5 , y = 4, label = paste("P = ",round(a$p.value,digits=4),sep = ""),fontface="italic"), size = 2)+geom_segment(mapping=aes(x=1,y=3.8,xend=2,yend=3.8))+geom_segment(mapping=aes(x=1,y=3.8,xend=1,yend=3.5))+geom_segment(mapping=aes(x=2,y=3.8,xend=2,yend=3.5))
}else if(a$p.value>0 && a$p.value<0.0001){ 
	plot1<-plot1+geom_text(aes(x=1.5 , y = 4, label = paste("P = ",formatC(a$p.value, format = "e", digits = 2),sep = ""),fontface="italic"), size = 2)+geom_segment(mapping=aes(x=1,y=3.8,xend=2,yend=3.8))+geom_segment(mapping=aes(x=1,y=3.8,xend=1,yend=3.5))+geom_segment(mapping=aes(x=2,y=3.8,xend=2,yend=3.5))
}else{
	plot1<-plot1+geom_text(aes(x=1.5 , y = 4, label = paste("P < ","2.22e-16",sep = ""),fontface="italic"), size = 2)+geom_segment(mapping=aes(x=1,y=3.8,xend=2,yend=3.8))+geom_segment(mapping=aes(x=1,y=3.8,xend=1,yend=3.5))+geom_segment(mapping=aes(x=2,y=3.8,xend=2,yend=3.5))
}

# plot1+geom_label(aes(x=1.5 , y = 4, label = paste("p-value = ",a$p.value,sep = "")), size = 5)

### add part to calculate ks-test values for Alu's and then update the variable in the below line.

if(b$p.value>=0.0001 && b$p.value<=1){
	plot1<-plot1+geom_text(aes(x=3.5 , y = 0, label = paste("P = ",round(b$p.value,digits=4),sep = ""),fontface="italic"), size = 2)+geom_segment(mapping=aes(x=3,y=0.2,xend=4,yend=0.2))+geom_segment(mapping=aes(x=3,y=0.2,xend=3,yend=0.5))+geom_segment(mapping=aes(x=4,y=0.2,xend=4,yend=0.5))
}else if(b$p.value>0 && b$p.value<0.0001){
	plot1<-plot1+geom_text(aes(x=3.5 , y = 0, label = paste("P = ",formatC(b$p.value, format = "e", digits = 2),sep = ""),fontface="italic"), size = 2)+geom_segment(mapping=aes(x=3,y=0.2,xend=4,yend=0.2))+geom_segment(mapping=aes(x=3,y=0.2,xend=3,yend=0.5))+geom_segment(mapping=aes(x=4,y=0.2,xend=4,yend=0.5))
}else{
	plot1<-plot1+geom_text(aes(x=3.5 , y = 0, label = paste("P < ","2.22e-16",sep = ""),fontface="italic"), size = 2)+geom_segment(mapping=aes(x=3,y=0.2,xend=4,yend=0.2))+geom_segment(mapping=aes(x=3,y=0.2,xend=3,yend=0.5))+geom_segment(mapping=aes(x=4,y=0.2,xend=4,yend=0.5))
}

#plot1<-plot1+geom_text(aes(x=3.5 , y = 0, label = paste("P = ",a$p.value,sep = "")), size = 5)+geom_segment(mapping=aes(x=3,y=0.2,xend=4,yend=0.2))+geom_segment(mapping=aes(x=3,y=0.2,xend=3,yend=0.5))+geom_segment(mapping=aes(x=4,y=0.2,xend=4,yend=0.5))

suppressMessages(ggsave(args[1],width=36,height=50,units="mm"))

# suppressMessages(ggsave(args[1]))

#########################################################################################################################################################

suppressMessages(barFile<-read_csv("tempBar.csv")) #should have column 1 as H2B names and column 2 as H2B values and column 3 as colors

barFile$values<-log10(barFile$values+1.02)

# plot2<-ggbarplot(barFile, x = "UID", y = "values",fill = "UID",color = "black",palette = barFile$colors,sort.by.groups = FALSE,x.text.angle = 90,xlab=FALSE,ylab="Log10 transformed expression values",show.legend=FALSE) + coord_cartesian(ylim = c(0, 4))

plot2<-ggbarplot(barFile, x = "UID", y = "values",fill = "UID",color = "black",palette = barFile$colors,sort.by.groups = FALSE,label=FALSE,xlab=FALSE,ylab=FALSE,show.legend=FALSE) + coord_cartesian(ylim = c(0, 4))

suppressMessages(ggsave(args[2],width=18,height=50,units="mm"))

# suppressMessages(ggsave(args[2]))

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

# table=c(a$p.value,
# b$p.value,
# l1WilcoxTwosidedBeforeLog$p.value,
# aluWilcoxTwosidedBeforeLog$p.value,
# l1WilcoxLessBeforeLog$p.value,
# aluWilcoxLessBeforeLog$p.value,
# l1WilcoxGreaterBeforeLog$p.value,
# aluWilcoxGreaterBeforeLog$p.value,
# l1KsLessBeforeLog$p.value,
# aluKsLessBeforeLog$p.value,
# l1KsGreaterBeforeLog$p.value,
# aluKsGreaterBeforeLog$p.value,
# l1KsTwosidedAfterLog$p.value,
# aluKsTwosidedAfterLog$p.value,
# l1WilcoxTwosidedAfterLog$p.value,
# aluWilcoxTwosidedAfterLog$p.value,
# l1WilcoxLessAfterLog$p.value,
# aluWilcoxLessAfterLog$p.value,
# l1WilcoxGreaterAfterLog$p.value,
# aluWilcoxGreaterAfterLog$p.value,
# l1KsLessAfterLog$p.value,
# aluKsLessAfterLog$p.value,
# l1KsGreaterAfterLog$p.value,
# aluKsGreaterAfterLog$p.value)

# table=as.data.frame(table)

# colnames(table)=args[3]
# rownames(table)=c('l1KsTwosidedBeforeLog',
# 'aluKsTwosidedBeforeLog',
# 'l1WilcoxTwosidedBeforeLog',
# 'aluWilcoxTwosidedBeforeLog',
# 'l1WilcoxLessBeforeLog',
# 'aluWilcoxLessBeforeLog',
# 'l1WilcoxGreaterBeforeLog',
# 'aluWilcoxGreaterBeforeLog',
# 'l1KsLessBeforeLog',
# 'aluKsLessBeforeLog',
# 'l1KsGreaterBeforeLog',
# 'aluKsGreaterBeforeLog',
# 'l1KsTwosidedAfterLog',
# 'aluKsTwosidedAfterLog',
# 'l1WilcoxTwosidedAfterLog',
# 'aluWilcoxTwosidedAfterLog',
# 'l1WilcoxLessAfterLog',
# 'aluWilcoxLessAfterLog',
# 'l1WilcoxGreaterAfterLog',
# 'aluWilcoxGreaterAfterLog',
# 'l1KsLessAfterLog',
# 'aluKsLessAfterLog',
# 'l1KsGreaterAfterLog',
# 'aluKsGreaterAfterLog')

# write.table(table,file=paste(args[3],"_pvalues.csv",sep=""),quote = F,row.names = F,sep=",",eol="\n")