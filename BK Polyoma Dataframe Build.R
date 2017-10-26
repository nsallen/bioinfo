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

# Detaches "ape" package to avoid conflict with "seqinr" package functions
detach("package:ape", unload=TRUE)

# Determines the consensus of BK Polyoma alignment
WTnt <- consensus(PolyDF)
length(WTnt)

# Counts the total number of sequences in BK Polyoma data
numCons <- length(which(BKseq[,1]==WTnt[1]))
numCons

# Writes new function to create transition mutations in nucleotides
transition <- function(nt){
  if(nt=="a") {return("g")}
  if(nt=="g") {return("a")}
  if(nt=="c") {return("t")}
  if(nt=="t") {return("c")}}

# For-loop to calculate mean frequency of transition mutations for each nucleotide:
MeanFreq<-c()
for (i in 1:ncol(BKseq)){
  MeanFreq<-c(MeanFreq,(length(which(BKseq[,i]==transition(WTnt[i])))/ncol(BKseq)))
}

# Constructing final data frame:
num<-c(1:length(WTnt))
df<-data.frame(num)
df<-cbind(df,WTnt)
df<-cbind(df,MeanFreq)
dffin<-df[-which(df$MeanFreq ==0.0),]
dffin
View(dffin)
