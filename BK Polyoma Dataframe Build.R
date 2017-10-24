library(seqinr)
library(ape)

# Reads BK Polyoma alignment and stores as variable
PolyDF <- read.alignment("BKpolyomavirus_VP1.fasta.mu.trim05", format = "fasta")
class(PolyDF)

# Determines the consensus of BK Polyoma alignment
PolyCon <- consensus(PolyDF)
class(PolyCon)

# Prepares BK Polyoma sequence as a matrix
BKseq <- read.dna("BKpolyomavirus_VP1.fasta.mu.trim05", format = "fasta", as.character=TRUE)
class(BKseq)
View(BKseq)

# Counts the total number of sequences in BK Polyoma data
numCons <- length(which(BKseq[,1]==PolyCon[1]))
numCons

# Determines transition mutations from the consensus
numTrans<-length(which(BKseq[,1]=="g"))
numTrans

transition <- function(nt){
  if(nt=="a") {return("g")}
  if(nt=="g") {return("a")}
  if(nt=="c") {return("t")}
  if(nt=="t") {return("c")}
  }

BKtrans <- transition(PolyCon)

numTrans<-length(which(BKseq[,12]==transition(PolyCon[12])))
numTrans
