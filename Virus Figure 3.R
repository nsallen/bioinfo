#The following packages only need to be installed once:
install.packages("seqinr")
install.packages("ggplot2")

#Actual code begins here:
library(seqinr)
library(ggplot2)

read.csv("OverviewSelCoeff_BachelerFilter.csv") -> CSV

#all green points for the left synonomous site graph
a <- frequenciesOfSynAmutsNONCP <- CSV[which(((CSV[, 4] == "syn"  )|(CSV[, 4] ==  "overlap")) & (CSV[, 7] == "a")),3]
b <- frequenciesOfSynTmutsNONCP <- CSV[which(((CSV[, 4] == "syn"  )|(CSV[, 4] ==  "overlap")) & (CSV[, 7] == "t")),3]
c <- frequenciesOfSynGmutsNONCP <- CSV[which(((CSV[, 4] == "syn"  )|(CSV[, 4] ==  "overlap")) & (CSV[, 7] == "g")),3]
d <- frequenciesOfSynCmutsNONCP <- CSV[which(((CSV[, 4] == "syn"  )|(CSV[, 4] ==  "overlap")) & (CSV[, 7] == "c")),3]

#blue dots of left graph that show a CpG-forming mutation
e <- frequenciesOfSynAmutsCP <- CSV[which(((CSV[, 4] == "syn"  )|(CSV[, 4] ==  "overlap")) & (CSV[, 7] == "a")& (CSV[, 15] == "1")),3]
f <- frequenciesOfSynTmutsCP <- CSV[which(((CSV[, 4] == "syn"  )|(CSV[, 4] ==  "overlap")) & (CSV[, 7] == "t")& (CSV[, 15] == "1")),3]

#all green points for the right NON-synonomous site graph
g <- frequenciesOfNONSynAmutsNONCP <- CSV[which((CSV[, 4] == "nonsyn"  ) & (CSV[, 7] == "a")& (CSV[, 14] == "0")),3]
h <- frequenciesOfNONSynTmutsNONCP <- CSV[which((CSV[, 4] == "nonsyn"  ) & (CSV[, 7] == "t")& (CSV[, 14] == "0")),3]
i <- frequenciesOfNONSynGmutsNONCP <- CSV[which((CSV[, 4] == "nonsyn"  ) & (CSV[, 7] == "g")& (CSV[, 14] == "0")),3]
j <- frequenciesOfNONSynCmutsNONCP <- CSV[which((CSV[, 4] == "nonsyn"  ) & (CSV[, 7] == "c")& (CSV[, 14] == "0")),3]

#all yellow points for the right NON-SYN site graph
k <- frequenciesOfNONSynAmutsDRASTIC <- CSV[which((CSV[, 4] == "nonsyn"  ) & (CSV[, 7] == "a")& (CSV[, 14] == "1")),3]
l <- frequenciesOfNONSynTmutsDRASTIC <- CSV[which((CSV[, 4] == "nonsyn"  ) & (CSV[, 7] == "t")& (CSV[, 14] == "1")),3]
m <- frequenciesOfNONSynGmutsDRASTIC <- CSV[which((CSV[, 4] == "nonsyn"  ) & (CSV[, 7] == "g")& (CSV[, 14] == "1")),3]
n <- frequenciesOfNONSynCmutsDRASTIC <- CSV[which((CSV[, 4] == "nonsyn"  ) & (CSV[, 7] == "c")& (CSV[, 14] == "1")),3]

#red points on the right
o <- frequenciesOfNONSynAmutsCPDRASTIC <- CSV[which(((CSV[, 4] == "nonsyn"  ) & (CSV[, 7] == "a")& (CSV[, 14] == "1") &(CSV[, 15] == "1"))),3]
p <- frequenciesOfNONSynTmutsCPDRASTIC <- CSV[which(((CSV[, 4] == "nonsyn"  ) & (CSV[, 7] == "t")& (CSV[, 14] == "1") &(CSV[, 15] == "1"))),3]

a
class(a)
pugs <- data.frame(a)
pugs
class(pugs)
tst = c(1,5,10,20)

plot(tst)
plot(pugs)


NUMS<-a
ddf = data.frame(NUMS, GRP = "a",replace=T)
library(lattice)

stripplot(NUMS,data=ddf, jitter.data=T)

class(NUMS)

ggplot(ddf, aes(x="A to T", y=NUMS)) +
  
#avoid plotting outliers twice
geom_jitter(position=position_jitter(width=.0025, height=0))

#p<-ggplot(mapping = aes(x = "a", y = "values")) + geom_point(data = data.frame(tst))
#p + geom_point(data = data.frame(tst))

amn <- mean(a)

ggplot() +
  #Graph=ggplot(data=df,mapping=aes(x=x,y=y,color=C),size=0.1)+geom_point()
  geom_point( data = data.frame(amn), mapping = aes(x="T to C", y=c, color= "red"), size = 2.0) +
  geom_point( data = data.frame(c), mapping = aes(x="T to C", y=c, color = "red"), size = 0.1, position=position_jitter(width=.0025, height=0)) +
  geom_point( data = data.frame(mean(d)), mapping = aes(x="A to T", y=mean.d., color = "blue"), size = 2.0) +
  geom_point( data = data.frame(d), aes(x=labels[2], y=d, color ="blue"), size = 0.7, position=position_jitter(width=.0025, height=0))+
  geom_point( data = data.frame(mean(b)), mapping = aes(x="A to gg", y=mean.b., color = "green"), size = 2.0) +
  geom_point( data = data.frame(b), aes(x="A to gg", y=b, color ="green"), size = 0.7, position=position_jitter(width=.0025, height=0))
b
mylist <- list(a,b)
mylist
a
b

mean(b)
mean(b)

n <- max(length(a), length(b), length(c), length(d), length(e), length(f), length(g), length(h), length(i), length(j))
length(a) <- n
length(c) <- n
length(b) <- n
length(d) <- n
length(e) <- n
length(f) <- n
length(g) <- n
length(h) <- n
length(i) <- n
length(j) <- n

vbmn <- rep (bmn,n)
vbmn
vamn <- rep (amn,n)
vamn
df1 <- data.frame(a,vamn,b,vbmn,c,d,e,f,g,h,i,j)
df1
zz <- 0.03

ggplot() +
  #Graph=ggplot(data=df,mapping=aes(x=x,y=y,color=C),size=0.1)+geom_point()
  geom_point( data = df1, mapping = aes(x="T to C", y=a, color = "red"), size = 2.0) +
  geom_point(data = df1, mapping = aes (x = "T to C", y = vamn, color ="blue"), size = 5.0) +
  geom_point( data = df1, mapping = aes(x="A to T", y=b, color ="green"), size = 1.) +
  geom_point(data = df1, mapping = aes (x = "A to T", y = vbmn, color ="blue"), size = 4.0) +
  geom_point( data = df1, mapping = aes(x="P         ", y=c, color ="red"), size = 2.0) +
  geom_point( data = df1, mapping = aes(x="Sir Pu   gs", y=d,color ="yellow"), size = 1.)+
  geom_point( data = df1, mapping = aes(x="    ", y=e,color ="red"), size = 1.)+
  geom_point( data = df1, mapping = aes(x="A  to T", y=f, color ="green"), size = 1.) +
  geom_point( data = df1, mapping = aes(x="P           ", y=f, color ="red"), size = 2.0) +
  geom_point( data = df1, mapping = aes(x="S ir Pugs", y=g,color ="yellow"), size = 1.)+
  geom_point( data = df1, mapping = aes(x=" ", y=i,color ="red"), size = 1.)+
  geom_point( data = df1, mapping = aes(x="    ", y=j,color ="red"), size = 1.)

labels <- c("T to C","A to T","C to T","G to A", " T to C ")



#rename(df1,c("a"="a","b"="e","c"="b","d"="f","e"="d","f"="c","g"="g","h"="k","i"="o","j"="h","k"="p","l"="l","m"="j","n"="n","o"="i","p"="m"))
col1 <- c(0.0)
amn <- mean(a)
bmn <- mean(b, na.rm = TRUE)
bmn
b

n <- max(length(a), length(b), length(c), length(d), length(e), length(f), length(g), length(h), length(i), length(j))
length(a) <- n
length(c) <- n
length(b) <- n
length(d) <- n
length(e) <- n
length(f) <- n
length(g) <- n
length(h) <- n
length(i) <- n
length(j) <- n

vbmn <- rep (bmn,n)
vbmn
vamn <- rep (amn,n)
vamn

c <- rep(col1,n)
df1<-data.frame(a,vamn,b,vbmn,c,d,e,f,g,h,i,j)
df1
c
se = 0.002

ggplot() +
  #Graph=ggplot(data=df,mapping=aes(x=x,y=y,color=C),size=0.1)+geom_point()
  geom_point( data = df1, mapping = aes(x="a", y=a), colour = "red", size = 0.1) +xlab("Mutation Type") + ylab("Samples' Average Frequencies of Mutations") +
  geom_point(data = df1, mapping = aes (x = "a", y = vamn), colour = "red",size = 5.0) +
  geom_errorbar(aes(x = "a",ymin=amn-se, ymax=amn + se), color = "red",width=.5) +
  geom_point( data = df1, mapping = aes(x="b", y=b), color ="green", size = 1.) +
  geom_point(data = df1, mapping = aes (x = "b", y = vbmn), color ="blue", size = 4.0) +
  geom_errorbar(aes(x = "b",ymin=amn-se, ymax=bmn + se), color = "red",width=.5) +
  geom_point( data = df1, mapping = aes(x="c", y=c), size = 0.0) +
  scale_x_discrete(labels=c("a" = "T -> G","b" = "T -> G", "c" = "")) +
  geom_point( data = df1, mapping = aes(x="Sir Pu   gs", y=d),color ="yellow", size = 1.)+
  geom_point( data = df1, mapping = aes(x=" V   ", y=e),color ="red", size = 1.)+
  geom_point( data = df1, mapping = aes(x="wA  to T", y=f), color ="green", size = 1.) +
  geom_point( data = df1, mapping = aes(x="P           ", y=f), color ="red", size = 2.0) +
  geom_point( data = df1, mapping = aes(x="S ir Pugs", y=g),color ="yellow", size = 1.)+
  geom_point( data = df1, mapping = aes(x="T ", y=i),color ="red", size = 1.)+
  geom_point( data = df1, mapping = aes(x="U    ", y=j),color ="red", size = 1.)+
  
  ggtitle("Genetic Divergence vs Elapsed Generations") +
  scale_y_log10() +
  theme(legend.position="none") +
  
  #scale_x_continuous()
  geom_vline(aes(xintercept=as.numeric(x[c(3, 5)])),linetype=1, colour="black") +
  geom_line(data = c(0,0), c(1,1))
            expand_limits(y = c(0.00001, 0.1)) +
              annotate("text", x  = 1.75, y = 0.00001, label = "T -> G") +
              theme(plot.margin = unit(c(1,1,2.5,1), "cm")) +
              annotate("text", x  = 1.75, y = 0.00001, label = "T -> G") +
              theme (panel.border = element_rect(colour = "black", fill=NA, size=3))
            require(grid)
            grid.text("
                      Text for margin note annotation goes here. Looks like you can just type a paragraph
                      and it will show where specified.",
                      x = unit(3, "cm"), y = unit(1.00,"cm"), just = "left", vjust = unit(0.0,"cm"))
            grid.lines(x = unit(c(1, 5), "cm"),
                       y = unit(c(1, 5), "cm"),
                       default.units = "cm",
                       arrow = NULL, name = NULL,
                       gp=gpar(), draw = TRUE, vp = NULL)
            -as.numeric(x[3])
            d <- data.frame(c(rep("A",5), rep("B",5)),
                            sample(c(1:10), 10, replace=TRUE),    
                            sample(c(1:10), 10, replace=TRUE),    
                            sample(c(1:10), 10, replace=TRUE),  
                            sample(c(1:10), 10, replace=TRUE),    
                            sample(c(1:10), 10, replace=TRUE))
            colnames(d) <- c("col1","col2","col3","col4","col5","col6" )
            cols_to_plot <- c("col2","col4","col6")
            
            for (i in seq_along(cols_to_plot)) {
              print(ggplot(data=d, aes_string(x = cols_to_plot[i], fill= "col1")) +
                      geom_density(alpha = 0.5))
            }
            
            d[,1]
            d
            
            
#left plot (Synonymous Sites)
            
par(mfrow=c(1,4))
plot(frequenciesOfSynAmutsNONCP,type="p",col="dark green",
                 xlab="A->G", ylab="Mutation Frequency",
                 xaxt="n")
            points(frequenciesOfSynAmutsCP, col="darkblue")
            plot(frequenciesOfSynTmutsNONCP,type="p",col="dark green",
                 xlab="T->C", ylab="Mutation Frequency",
                 xaxt="n")
            points(frequenciesOfSynTmutsCP, col="darkblue")
            plot(frequenciesOfSynCmutsNONCP, type="p", col="dark green",
                 xlab= "C->T", ylab="Mutation Frequency",
                 xaxt="n")
            plot(frequenciesOfSynGmutsNONCP, type="p", col="dark green",
                 xlab= "G->A", ylab="Mutation Frequency",
                 xaxt="n")
            
#right plot (Non-Synonymous Sites) -- would need to add blue points
par(mfrow=c(1,4))
plot(frequenciesOfNONSynAmutsNONCP,type="p",col="dark green",
                 xlab="A->G", ylab="Mutation Frequency",
                 xaxt="n")
            points(frequenciesOfNONSynAmutsDRASTIC, col="darkyellow")
            points(frequenciesOfNONSynAmutsCPDRASTIC, col="red")
plot(frequenciesOfSynTmutsNONCP,type="p",col="dark green",
                 xlab="T->C", ylab="Mutation Frequency",
                 xaxt="n")
            points(frequenciesOfNONSynTmutsDRASTIC, col="darkyellow")
            points(frequenciesOfNONSynTmutsCPDRASTIC, col="red")
plot(frequenciesOfSynCmutsNONCP, type="p", col="dark green",
                 xlab= "C->T", ylab="Mutation Frequency",
                 xaxt="n")
            points(frequenciesOfNONSynCmutsDRASTIC, col="darkyellow")
plot(frequenciesOfSynGmutsNONCP, type="p", col="dark green",
                 xlab= "G->A", ylab="Mutation Frequency",
                 xaxt="n")
            points(frequenciesOfNONSynGmutsDRASTIC, col="darkyellow")