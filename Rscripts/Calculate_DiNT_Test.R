# calculate dinucleotide frequencies
library(msa)
library(ape)
library(car)
library(seqinr)
#library(bios2mds)
#library(pegas)


# Compare two ways to calcualte
# 1. read alignment together
align1<-read.dna(paste0("FASTA.Aligned.Sep2019/",Fastafiles[i]),format = "fasta")
# calculate dinucleotide frequency
count(align1, 2)
#aa    ac    ag    at    ca    cc    cg    ct    ga    gc    gg    gt    ta    tc    tg    tt 
#79678    97   197    90   108 47300    38   139   184    80 56979    67    85   108    93 56289 

# 2. read as alignment and calcualte dinucleotide freq for each seq 
i=1
fas1<-read.alignment(paste0("FASTA.Aligned.Sep2019/",Fastafiles[i]),format = "fasta")
freq<-list()
for (j in 1: length(fas1[[3]])){
    seq<-fas1[[3]][[j]]
    seq<-substring(seq, 1:nchar(seq),1:nchar(seq))
    
    #record dinucleotide freq
    freq[[j]]<-as.data.frame(seqinr::count(seq,2))
    names(freq)[j]<-fas1[[2]][j]
    #freq$rho<-rho(seq,2)
    #freq$Zscore<-zscore(seq, modele="base")
}

for (j in 1:length(freq)) {colnames(freq[[j]])<-c("diNT",paste0(names(freq[j])))}

DiNtFreq<-freq%>% purrr::reduce(full_join, by='diNT')
DiNtFreq<-as.numeric(DiNt)
rowSums(DiNtFreq[,2:223])
#[1] 30003 13645 21266 14933 16557 16303   450 14276 17529 10526 16733 12536 15751  7115 18867 14837

# Different-> use #2 method