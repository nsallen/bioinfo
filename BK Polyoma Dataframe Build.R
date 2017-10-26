library(seqinr)
library(ape)

# Reads BK Polyoma alignment and stores as variable
PolyDF <- read.alignment("BKpolyomavirus_VP1.fasta.mu.trim05", format = "fasta")
class(PolyDF)
PolyDF


# Prepares BK Polyoma sequence as a matrix
BKseq <- read.dna("BKpolyomavirus_VP1.fasta.mu.trim05", format = "fasta", as.character=TRUE)
class(BKseq)
BKseq

detach("package:ape", unload=TRUE)
# Determines the consensus of BK Polyoma alignment
PolyCon <- consensus(PolyDF)
length(PolyCon)

# Counts the total number of sequences in BK Polyoma data
numCons <- length(which(BKseq[,1]==PolyCon[1]))
numCons


transition <- function(nt){
if(nt=="a") {return("g")}
if(nt=="g") {return("a")}
if(nt=="c") {return("t")}
if(nt=="t") {return("c")}}

numTrans<-c()
for (i in 1:ncol(BKseq)){
  numTrans<-c(numTrans,(length(which(BKseq[,i]==transition(PolyCon[i])))/ncol(BKseq)))
}

locNum<-c(1:length(PolyCon))
df<-data.frame(locNum)
df<-cbind(df,PolyCon)
df<-cbind(df,numTrans)
dffin<-df[-which(df$numTrans ==0.0),]
dffin
View(dffin)
