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

# Building a transition function
transition <- function(nt){
  if(nt=="a") {return("g")}
  if(nt=="g") {return("a")}
  if(nt=="c") {return("t")}
  if(nt=="t") {return("c")}
  }

# Creates a new variable to make an empty data frame with exact number of rows
pos<-c(1:length(PolyCon))

# Creates the data frame with 1089 rows 
BKdata<-data.frame(pos)

# Creates 2 new columns 
BKdata$nt=""
BKdata$tr=""
View(BKdata)

# Places BK polyoma consensus nucleotides into position
PolyCon -> BKdata$nt

# Determines transition nucleotides from consensus matrix
BKtrans<-c()
for(i in 1:length(PolyCon)){
  trans_letter<-transition(PolyCon[i])
  BKtrans<-c(BKtrans,trans_letter)
}

BKtrans

### WHAT'S LEFT TO DO: Calculate frequency of transition mutations from each sequence for each nucleotide position (% out of 116 possible) compared to consensus sequence
