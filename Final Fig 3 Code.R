### The following packages only need to be installed once!
install.packages("ggplot2")
install.packages("plyr")
install.packages("grid")

### Code begins here!
read.csv("OverviewSelCoeff_BachelerFilter.csv") -> CSV
library(ggplot2)
library(plyr)
library(grid)
library(seqinr)

#all green points for the left synonomous site graph
a <- frequenciesOfSynAmutsNONCP <- CSV[which(((CSV[, 4] == "syn" )|(CSV[, 4] == "overlap")) & (CSV[, 7] == "a")),3]
b <- frequenciesOfSynTmutsNONCP <- CSV[which(((CSV[, 4] == "syn" )|(CSV[, 4] == "overlap")) & (CSV[, 7] == "t")),3]
c <- frequenciesOfSynGmutsNONCP <- CSV[which(((CSV[, 4] == "syn" )|(CSV[, 4] == "overlap")) & (CSV[, 7] == "g")),3]
d <- frequenciesOfSynCmutsNONCP <- CSV[which(((CSV[, 4] == "syn" )|(CSV[, 4] == "overlap")) & (CSV[, 7] == "c")),3]

#blue dots of left graph that show a CpG-forming mutation
e <- frequenciesOfSynAmutsCP <- CSV[which(((CSV[, 4] == "syn" )|(CSV[, 4] == "overlap")) & (CSV[, 7] == "a")& (CSV[, 15] == "1")),3]
f <- frequenciesOfSynTmutsCP <- CSV[which(((CSV[, 4] == "syn" )|(CSV[, 4] == "overlap")) & (CSV[, 7] == "t")& (CSV[, 15] == "1")),3]

#all green points for the right NON-synonomous site graph
g <- frequenciesOfNONSynAmutsNONCP <- CSV[which((CSV[, 4] == "nonsyn" ) & (CSV[, 7] == "a")& (CSV[, 14] == "0")),3]
h <- frequenciesOfNONSynTmutsNONCP <- CSV[which((CSV[, 4] == "nonsyn" ) & (CSV[, 7] == "t")& (CSV[, 14] == "0")),3]
i <- frequenciesOfNONSynGmutsNONCP <- CSV[which((CSV[, 4] == "nonsyn" ) & (CSV[, 7] == "g")& (CSV[, 14] == "0")),3]
j <- frequenciesOfNONSynCmutsNONCP <- CSV[which((CSV[, 4] == "nonsyn" ) & (CSV[, 7] == "c")& (CSV[, 14] == "0")),3]

#all yellow points for the right NON-SYN site graph
k <- frequenciesOfNONSynAmutsDRASTIC <- CSV[which((CSV[, 4] == "nonsyn" ) & (CSV[, 7] == "a")& (CSV[, 14] == "1")),3]
l <- frequenciesOfNONSynTmutsDRASTIC <- CSV[which((CSV[, 4] == "nonsyn" ) & (CSV[, 7] == "t")& (CSV[, 14] == "1")),3]
m <- frequenciesOfNONSynGmutsDRASTIC <- CSV[which((CSV[, 4] == "nonsyn" ) & (CSV[, 7] == "g")& (CSV[, 14] == "1")),3]
n <- frequenciesOfNONSynCmutsDRASTIC <- CSV[which((CSV[, 4] == "nonsyn" ) & (CSV[, 7] == "c")& (CSV[, 14] == "1")),3]

#red points on the right
o <- frequenciesOfNONSynAmutsCPDRASTIC <- CSV[which(((CSV[, 4] == "nonsyn" ) & (CSV[, 7] == "a")& (CSV[, 14] == "1") &(CSV[, 15] == "1"))),3]
p <- frequenciesOfNONSynTmutsCPDRASTIC <- CSV[which(((CSV[, 4] == "nonsyn" ) & (CSV[, 7] == "t")& (CSV[, 14] == "1") &(CSV[, 15] == "1"))),3]

mylist <- list (a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p)

maxlngth = length (mylist[[2]])
for (i in 1:16){
  if (length(mylist[[i]])>maxlngth){
    maxlngth = length(mylist[[i]])
  }
}

means<-c()
ucls<-c()
lcls<-c()
dfval<-data.frame(matrix(nrow=maxlngth))
for (i in 1:16){
  mn<-mean(mylist[[i]],na.rm = TRUE)
  means<-c(means,mn)
  sdev<-sd(mylist[[i]],na.rm = TRUE)
  upl<-mn + sdev
  ucls<-c(ucls,upl)
  lwl<-mn-sdev
  lcls<-c(lcls,lwl)
  length(mylist[[i]]) = maxlngth
  dfval<-cbind(dfval,data.frame(mylist[[i]]))
}
mxstats<-rbind(means,lcls)
mxstats<-rbind(mxstats,ucls)
dfstats<-data.frame(mxstats)

dfstats["means","X1"]

names (dfval)[1] <- "xxx"
for (i in 2:17){
  names (dfval) [i] = paste ("val" , i-1, sep = "_")
}

dfval<-rename(dfval,c("val_1"="val_a","val_5"="val_b","val_2"="val_c","val_6"="val_d","val_4"="val_e","val_3"="val_f","val_7"="val_g","val_11"="val_h","val_15"="val_i","val_8"="val_j","val_16"="val_k","val_12"="val_l","val_10"="val_m","val_14"="val_n","val_9"="val_o","val_13"="val_p"))
dfstats<-rename(dfstats,c("X1"="mut_a","X5"="mut_b","X2"="mut_c","X6"="mut_d","X4"="mut_e","X3"="mut_f","X7"="mut_g","X11"="mut_h","X15"="mut_i","X8"="mut_j","X16"="mut_k","X12"="mut_l","X10"="mut_m","X14"="mut_n","X9"="mut_o","X13"="mut_p"))

dfstats

DFmylist <- data.frame(mylist)
class(DFmylist)
View(DFmylist)


Synplt <- ggplot() +
  
  geom_point( data = dfval, mapping = aes(x="a", y= val_a), colour = "red", size = 0.1) +xlab("Mutation Type") + ylab("Samples' Average Frequencies of Mutations") +
  geom_point(data = dfstats, mapping = aes (x = "a", y = dfstats["means","mut_a"]), colour = "red",size = 5.0) +
  geom_errorbar(data = dfstats,aes(x = "a", ymin= dfstats["lcls","mut_a"], ymax= dfstats["ucls","mut_a"], color = "red",width=.5)) +
  
  annotate("text", x = 1.5, y = 0.00001, label = "A -> G") +
  
  geom_point( data = dfval, mapping = aes(x="b", y= val_b), colour = "red", size = 0.1) +
  geom_point(data = dfstats, mapping = aes (x = "b", y = dfstats["means","mut_b"]), colour = "red",size = 5.0) +
  geom_errorbar(data = dfstats,aes(x = "b", ymin= dfstats["lcls","mut_b"], ymax= dfstats["ucls","mut_b"], color = "red",width=.5)) +
  
  geom_vline(aes(linetype=1, colour="black"),xintercept =c(2.5)) +
  
  geom_point( data = dfval, mapping = aes(x="c", y= val_c), colour = "red", size = 0.1) +
  geom_point(data = dfstats, mapping = aes (x = "c", y = dfstats["means","mut_c"]), colour = "red",size = 5.0) +
  geom_errorbar(data = dfstats,aes(x = "c", ymin= dfstats["lcls","mut_c"], ymax= dfstats["ucls","mut_c"], color = "red",width=.5)) +
  
  annotate("text", x = 3.5, y = 0.00001, label = "T -> C") +
  
  geom_point( data = dfval, mapping = aes(x="d", y= val_d), colour = "red", size = 0.1) +
  geom_point(data = dfstats, mapping = aes (x = "d", y = dfstats["means","mut_d"]), colour = "red",size = 5.0) +
  geom_errorbar(data = dfstats,aes(x = "d", ymin= dfstats["lcls","mut_d"], ymax= dfstats["ucls","mut_d"], color = "red",width=.5)) +
  
  geom_vline(aes(linetype=1, colour="black"),xintercept =c(4.5)) +
  
  geom_point( data = dfval, mapping = aes(x="e", y= val_e), colour = "red", size = 0.1) +
  geom_point(data = dfstats, mapping = aes (x = "e", y = dfstats["means","mut_e"]), colour = "red",size = 5.0) +
  geom_errorbar(data = dfstats,aes(x = "e", ymin= dfstats["lcls","mut_e"], ymax= dfstats["ucls","mut_e"], color = "red",width=.5)) +
  
  geom_point(data = dfstats, mapping = aes (x = "e1", y = 0.0), colour = "red",size = 0.0) +
  
  annotate("text", x = 5.5, y = 0.00001, label = "C -> T") +
  
  geom_vline(aes(linetype=1, colour="black"),xintercept =c(6.5)) +
  
  geom_point( data = dfval, mapping = aes(x="f", y= val_f), colour = "red", size = 0.1) +
  geom_point(data = dfstats, mapping = aes (x = "f", y = dfstats["means","mut_f"]), colour = "red",size = 5.0) +
  geom_errorbar(data = dfstats,aes(x = "f", ymin= dfstats["lcls","mut_f"], ymax= dfstats["ucls","mut_f"], color = "red",width=.5)) +
  
  geom_point(data = dfstats, mapping = aes (x = "f1", y = 0.0), colour = "red",size = 0.0) +
  
  annotate("text", x = 7.5, y = 0.00001, label = "G -> A") +
  
  scale_x_discrete(labels=c("a" = "", "b" = "", "c" = "", "d" = "","e" = "", "e1" = "","f" ="", "f1" = "")) +
  
  #ggtitle("HIV Genetic Mutations Study - Frequency vs Type") +
  scale_y_log10() +
  theme(legend.position="none") +
  expand_limits(y = c(0.00001, 0.1)) +
  theme(plot.margin = unit(c(1,1,2.5,1), "cm")) +
  theme (panel.border = element_rect(colour = "black", fill=NA, size=3),plot.title = element_text(hjust = 0.5))

NonSynplt <- ggplot() +
  
  geom_point( data = dfval, mapping = aes(x="g", y= val_g), colour = "red", size = 0.1) +xlab("Mutation Type") + ylab("Samples' Average Frequencies of Mutations") +
  geom_point(data = dfstats, mapping = aes (x = "g", y = dfstats["means","mut_g"]), colour = "red",size = 5.0) +
  geom_errorbar(data = dfstats,aes(x = "g", ymin= dfstats["lcls","mut_g"], ymax= dfstats["ucls","mut_g"], color = "red",width=.5)) +
  
  annotate("text", x = 2.0, y = 0.00001, label = "A -> G") +
  
  geom_point( data = dfval, mapping = aes(x="h", y= val_h), colour = "red", size = 0.1) +
  geom_point(data = dfstats, mapping = aes (x = "h", y = dfstats["means","mut_h"]), colour = "red",size = 5.0) +
  geom_errorbar(data = dfstats,aes(x = "h", ymin= dfstats["lcls","mut_h"], ymax= dfstats["ucls","mut_h"], color = "red",width=.5)) +
  
  geom_point( data = dfval, mapping = aes(x="i", y= val_i), colour = "red", size = 0.1) +
  geom_point(data = dfstats, mapping = aes (x = "i", y = dfstats["means","mut_i"]), colour = "red",size = 5.0) +
  geom_errorbar(data = dfstats,aes(x = "i", ymin= dfstats["lcls","mut_i"], ymax= dfstats["ucls","mut_i"], color = "red",width=.5)) +
  
  geom_vline(aes(linetype=1, colour="black"),xintercept =c(3.5)) +
  
  geom_point( data = dfval, mapping = aes(x="j", y= val_j), colour = "red", size = 0.1) +
  geom_point(data = dfstats, mapping = aes (x = "j", y = dfstats["means","mut_j"]), colour = "red",size = 5.0) +
  geom_errorbar(data = dfstats,aes(x = "j", ymin= dfstats["lcls","mut_j"], ymax= dfstats["ucls","mut_j"], color = "red",width=.5)) +
  
  annotate("text", x = 5, y = 0.00001, label = "T -> C") +
  
  geom_point( data = dfval, mapping = aes(x="k", y= val_k), colour = "red", size = 0.1) +
  geom_point(data = dfstats, mapping = aes (x = "k", y = dfstats["means","mut_k"]), colour = "red",size = 5.0) +
  geom_errorbar(data = dfstats,aes(x = "k", ymin= dfstats["lcls","mut_k"], ymax= dfstats["ucls","mut_k"], color = "red",width=.5)) +
  
  geom_point( data = dfval, mapping = aes(x="l", y= val_l), colour = "red", size = 0.1) +
  geom_point(data = dfstats, mapping = aes (x = "l", y = dfstats["means","mut_l"]), colour = "red",size = 5.0) +
  geom_errorbar(data = dfstats,aes(x = "l", ymin= dfstats["lcls","mut_l"], ymax= dfstats["ucls","mut_l"], color = "red",width=.5)) +
  
  geom_vline(aes(linetype=1, colour="black"),xintercept =c(6.5)) +
  
  geom_point( data = dfval, mapping = aes(x="m", y= val_m), colour = "red", size = 0.1) +
  geom_point(data = dfstats, mapping = aes (x = "m", y = dfstats["means","mut_m"]), colour = "red",size = 5.0) +
  geom_errorbar(data = dfstats,aes(x = "m", ymin= dfstats["lcls","mut_m"], ymax= dfstats["ucls","mut_m"], color = "red",width=.5)) +
  
  annotate("text", x = 7.5, y = 0.00001, label = "C -> T") +
  
  geom_point( data = dfval, mapping = aes(x="n", y= val_n), colour = "red", size = 0.1) +
  geom_point(data = dfstats, mapping = aes (x = "n", y = dfstats["means","mut_n"]), colour = "red",size = 5.0) +
  geom_errorbar(data = dfstats,aes(x = "n", ymin= dfstats["lcls","mut_n"], ymax= dfstats["ucls","mut_n"], color = "red",width=.5)) +
  
  geom_vline(aes(linetype=1, colour="black"),xintercept =c(8.5)) +
  
  geom_point( data = dfval, mapping = aes(x="o", y= val_o), colour = "red", size = 0.1) +
  geom_point(data = dfstats, mapping = aes (x = "o", y = dfstats["means","mut_o"]), colour = "red",size = 5.0) +
  geom_errorbar(data = dfstats,aes(x = "o", ymin= dfstats["lcls","mut_o"], ymax= dfstats["ucls","mut_o"], color = "red",width=.5)) +
  
  annotate("text", x = 9.5, y = 0.00001, label = "G -> A") +
  
  geom_point( data = dfval, mapping = aes(x="p", y= val_p), colour = "red", size = 0.1) +
  geom_point(data = dfstats, mapping = aes (x = "p", y = dfstats["means","mut_p"]), colour = "red",size = 5.0) +
  geom_errorbar(data = dfstats,aes(x = "p", ymin= dfstats["lcls","mut_p"], ymax= dfstats["ucls","mut_p"], color = "red",width=.5)) +
  
  scale_x_discrete(labels=c("g" = "", "h" = "", "i" = "", "j" = "","k" = "","l" ="","m" ="", "n"="","o" = "","p" = "")) +
  
  #ggtitle("HIV Genetic Mutations Study - Frequency vs Type") +
  scale_y_log10() +
  theme(legend.position="none") +
  expand_limits(y = c(0.00001, 0.1)) +
  theme(plot.margin = unit(c(1,1,2.5,1), "cm")) +
  theme (panel.border = element_rect(colour = "black", fill=NA, size=3),plot.title = element_text(hjust = 0.5))

#, vp=viewport(width=1.0, height=0.97)
library("gridExtra", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4")
require(grid)
require(gridExtra)
title1=textGrob("
                Fig. 3: type of site vs frequency, highlighting mean and 1 sd range", gp=gpar(fontface="bold", fontsize = 16, cex = 1))
grid.arrange( top =title1,Synplt + ggtitle('Synonymous Plots'), NonSynplt + ggtitle('Non-synonymous Plots'), nrow=1)

grid.text("
          Text for margin note annotation goes here. Looks like you can just type a paragraph
          and it will show where specified.",
          x = unit(3, "cm"), y = unit(1.00,"cm"), just = "left", vjust = unit(0.0,"cm"))