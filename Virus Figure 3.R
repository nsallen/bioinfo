#The following packages only need to be installed once:
install.packages("seqinr")
install.packages("ggplot2")

#Actual code begins here:
library(seqinr)
library(ggplot2)


read.csv("OverviewSelCoeff_BachelerFilter.csv") -> CSV

frequenciesOfSynAmuts <- CSV[which(((CSV[, 4] == "syn"  )|(CSV[, 4] ==  "overlap")) & (CSV[, 7] == "a")),3]
frequenciesOfSynTmuts <- CSV[which(((CSV[, 4] == "syn"  )|(CSV[, 4] ==  "overlap")) & (CSV[, 7] == "t")),3]
frequenciesOfSynGmuts <- CSV[which(((CSV[, 4] == "syn"  )|(CSV[, 4] ==  "overlap")) & (CSV[, 7] == "g")),3]
frequenciesOfSynCmuts <- CSV[which(((CSV[, 4] == "syn"  )|(CSV[, 4] ==  "overlap")) & (CSV[, 7] == "c")),3]


frequenciesOfNONSynAmuts <- CSV[which(CSV[, 4] == "nonsyn"  ) & (CSV[, 7] == "a",3]
frequenciesOfNONSynTmuts <- CSV[which(CSV[, 4] == "nonsyn"  ) & (CSV[, 7] == "t"),3]
frequenciesOfNONSynGmuts <- CSV[which(CSV[, 4] == "nonsyn"  ) & (CSV[, 7] == "g"),3]
frequenciesOfNONSynCmuts <- CSV[which(CSV[, 4] == "nonsyn"  ) & (CSV[, 7] == "c"),3]

frequenciesOfNONSynAmuts <- CSV[which(CSV[, 4] == "nonsyn"  ) & (CSV[, 7] == "a",3]
frequenciesOfNONSynTmuts <- CSV[which(CSV[, 4] == "nonsyn"  ) & (CSV[, 7] == "t"),3]
frequenciesOfNONSynGmuts <- CSV[which(CSV[, 4] == "nonsyn"  ) & (CSV[, 7] == "g"),3]
frequenciesOfNONSynCmuts <- CSV[which(CSV[, 4] == "nonsyn"  ) & (CSV[, 7] == "c"),3]


mean(frequenciesOfNONSynAmuts)
mean(frequenciesOfNONSynTmuts)
mean(frequenciesOfNONSynGmuts)
mean(frequenciesOfNONSynCmuts)
mean(frequenciesOfSynAmuts)
mean(frequenciesOfSynTmuts)
mean(frequenciesOfSynGmuts)
mean(frequenciesOfSynCmuts)







