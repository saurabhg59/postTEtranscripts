args = commandArgs(trailingOnly=TRUE);

library(ggplot2)
library(readr)

datafile<-read_csv("temp.csv")

datafile$value<-log10(datafile$value+0.1)

datafile$id<-as.factor(datafile$id)

plot<-ggplot(datafile, aes(x=id, y=value, fill=id)) + geom_violin() + labs(y = "Log10 transformed expression levels", x="")

plot<-plot+scale_fill_manual(values=c("royalblue4","royalblue1", "orange1"))

plot <-plot + geom_boxplot(width=0.05) + stat_summary(fun.y=mean, geom="point", size=2, color="white") + theme_bw()+theme(legend.position = "none")

plot+coord_cartesian(ylim = c(-2,5))

ggsave("args[1]")