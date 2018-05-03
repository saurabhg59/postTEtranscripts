#! /usr/bin/Rscript

args = commandArgs(trailingOnly=TRUE);

library(ggplot2)
library(readr)
library(ggpubr)
library(magrittr)

datafile<-read_csv("temp.csv") # should have a column called value and a column called id

datafile$value<-log10(datafile$value+0.2)

datafile$id<-factor(datafile$id,labels = c("L1Base WT","L1Base Mutant","Young L1 WT","Young L1 Mutant","Middle L1 WT","Middle L1 Mutant","Old L1 WT","Old L1 Mutant"))

plot<-ggboxplot(datafile,x="id",y="value",color="black",fill="id",palette=c("royalblue1","orange1","royalblue4","orange4","royalblue4","orange4","royalblue4","orange4"),shape="id", xlab=FALSE, ylab="Log10 transformed expression levels", show.legend=FALSE)

ggsave(args[1])