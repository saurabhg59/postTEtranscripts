#! /usr/bin/Rscript

args = commandArgs(trailingOnly=TRUE);

options(warn=-1)

library(magrittr)
library(ggplot2)
library(readr)
library(ggpubr)


suppressMessages(datafile<-read_csv(args[2])) # should have a column called value and a column called id

datafile$value<-log10(datafile$value+1.1)

l1baseWT<-subset(datafile,datafile$id==1)

l1baseM<-subset(datafile,datafile$id==2)

youngWT<-subset(datafile,datafile$id==3)

youngM<-subset(datafile,datafile$id==4)

middleWT<-subset(datafile,datafile$id==5)

middleM<-subset(datafile,datafile$id==6)

oldWT<-subset(datafile,datafile$id==7)

oldM<-subset(datafile,datafile$id==8)

a = ks.test(l1baseWT$value,l1baseM$value,alternative="g")

b = ks.test(youngWT$value,youngM$value,alternative="g")

c = ks.test(middleWT$value,middleM$value,alternative="g")

d = ks.test(oldWT$value,oldM$value,alternative="g")

l1basePos=max(max(l1baseWT$value),max(l1baseM$value))+1
youngPos=max(max(youngWT$value),max(youngM$value))+1
middlePos=max(max(middleWT$value),max(middleM$value))+1
oldPos=max(max(oldWT$value),max(oldM$value))+1

datafile$id<-factor(datafile$id,labels = c("L1Base WT","L1Base Mutant","Young L1 WT","Young L1 Mutant","Middle L1 WT","Middle L1 Mutant","Old L1 WT","Old L1 Mutant"))

plot<-ggboxplot(datafile,x="id",y="value",color="black",fill="id",palette=c("royalblue1","orange1","royalblue1","orange1","royalblue1","orange1","royalblue1","orange1"),shape="id", xlab=FALSE, ylab="Log10 transformed expression levels", show.legend=FALSE,outlier.shape=NA)

plot<-plot+rotate_x_text(angle=90)

if(a$p.value<=0.05){
	if(a$p.value>=0.0001 && a$p.value<=1){
		plot<-plot+geom_text(aes(x=1.5 , y = l1basePos, label = paste("P = ",round(a$p.value,digits=4),sep = ""),fontface="italic",family = "Helvetica"), size = 2)+geom_segment(mapping=aes(x=1,y=l1basePos-0.2,xend=2,yend=l1basePos-0.2))+geom_segment(mapping=aes(x=1,y=l1basePos-0.2,xend=1,yend=l1basePos-0.5))+geom_segment(mapping=aes(x=2,y=l1basePos-0.2,xend=2,yend=l1basePos-0.5))
	}else if(a$p.value>0 && a$p.value<0.0001){ 
		plot<-plot+geom_text(aes(x=1.5 , y = l1basePos, label = paste("P = ",formatC(a$p.value, format = "e", digits = 2),sep = ""),fontface="italic",family = "Helvetica"), size = 2)+geom_segment(mapping=aes(x=1,y=l1basePos-0.2,xend=2,yend=l1basePos-0.2))+geom_segment(mapping=aes(x=1,y=l1basePos-0.2,xend=1,yend=l1basePos-0.5))+geom_segment(mapping=aes(x=2,y=l1basePos-0.2,xend=2,yend=l1basePos-0.5))
	}else{
		plot<-plot+geom_text(aes(x=1.5 , y = l1basePos, label = paste("P < ","2.22e-16",sep = ""),fontface="italic",family = "Helvetica"), size = 2)+geom_segment(mapping=aes(x=1,y=l1basePos-0.2,xend=2,yend=l1basePos-0.2))+geom_segment(mapping=aes(x=1,y=l1basePos-0.2,xend=1,yend=l1basePos-0.5))+geom_segment(mapping=aes(x=2,y=l1basePos-0.2,xend=2,yend=l1basePos-0.5))
	}	
}

if(b$p.value<=0.05){
	if(b$p.value>=0.0001 && b$p.value<=1){
		plot<-plot+geom_text(aes(x=3.5 , y = youngPos, label = paste("P = ",round(b$p.value,digits=4),sep = ""),fontface="italic",family = "Helvetica"), size = 2)+geom_segment(mapping=aes(x=3,y=youngPos-0.2,xend=4,yend=youngPos-0.2))+geom_segment(mapping=aes(x=3,y=youngPos-0.2,xend=3,yend=youngPos-0.5))+geom_segment(mapping=aes(x=4,y=youngPos-0.2,xend=4,yend=youngPos-0.5))
	}else if(b$p.value>0 && b$p.value<0.0001){ 
		plot<-plot+geom_text(aes(x=3.5 , y = youngPos, label = paste("P = ",formatC(b$p.value, format = "e", digits = 2),sep = ""),fontface="italic",family = "Helvetica"), size = 2)+geom_segment(mapping=aes(x=3,y=youngPos-0.2,xend=4,yend=youngPos-0.2))+geom_segment(mapping=aes(x=3,y=youngPos-0.2,xend=3,yend=youngPos-0.5))+geom_segment(mapping=aes(x=4,y=youngPos-0.2,xend=4,yend=youngPos-0.5))
	}else{
		plot<-plot+geom_text(aes(x=3.5 , y = youngPos, label = paste("P < ","2.22e-16",sep = ""),fontface="italic",family = "Helvetica"), size = 2)+geom_segment(mapping=aes(x=3,y=youngPos-0.2,xend=4,yend=youngPos-0.2))+geom_segment(mapping=aes(x=3,y=youngPos-0.2,xend=3,yend=youngPos-0.5))+geom_segment(mapping=aes(x=4,y=youngPos-0.2,xend=4,yend=youngPos-0.5))
	}
}

if(c$p.value<=0.05){
	print(round(c$p.value,digits=4))
	plot<-plot+geom_text(aes(x=5.5 , y = youngPos, label = "P = ",fontface="italic",family = "Helvetica"), size = 2)+geom_segment(mapping=aes(x=5,y=youngPos-0.2,xend=6,yend=youngPos-0.2))+geom_segment(mapping=aes(x=5,y=youngPos-0.2,xend=5,yend=youngPos-0.5))+geom_segment(mapping=aes(x=6,y=youngPos-0.2,xend=6,yend=youngPos-0.5))
}

if(d$p.value<=0.05){
	if(d$p.value>=0.0001 && d$p.value<=1){
		plot<-plot+geom_text(aes(x=7.5 , y = oldPos, label = paste("P = ",round(d$p.value,digits=4),sep = ""),fontface="italic",family = "Helvetica"), size = 2)+geom_segment(mapping=aes(x=7,y=oldPos-0.2,xend=8,yend=oldPos-0.2))+geom_segment(mapping=aes(x=7,y=oldPos-0.2,xend=7,yend=oldPos-0.5))+geom_segment(mapping=aes(x=8,y=oldPos-0.2,xend=8,yend=oldPos-0.5))
	}else if(d$p.value>0 && d$p.value<0.0001){ 
		plot<-plot+geom_text(aes(x=7.5 , y = oldPos, label = paste("P = ",formatC(d$p.value, format = "e", digits = 2),sep = ""),fontface="italic",family = "Helvetica"), size = 2)+geom_segment(mapping=aes(x=7,y=oldPos-0.2,xend=8,yend=oldPos-0.2))+geom_segment(mapping=aes(x=7,y=oldPos-0.2,xend=7,yend=oldPos-0.5))+geom_segment(mapping=aes(x=8,y=oldPos-0.2,xend=8,yend=oldPos-0.5))
	}else{
		plot<-plot+geom_text(aes(x=7.5 , y = oldPos, label = paste("P < ","2.22e-16",sep = ""),fontface="italic",family = "Helvetica"), size = 2)+geom_segment(mapping=aes(x=7,y=oldPos-0.2,xend=8,yend=oldPos-0.2))+geom_segment(mapping=aes(x=7,y=oldPos-0.2,xend=7,yend=oldPos-0.5))+geom_segment(mapping=aes(x=8,y=oldPos-0.2,xend=8,yend=oldPos-0.5))
	}	
}

ggsave(args[1])