library(ape)
library(seqinr)

al<-read.fasta("BKpolyomavirus_VP1.fasta",as.string=TRUE)
al<-read.fasta("humanparainfluenzavirus3.fasta", as.string=T)
pos<-c()
for (i in 1: length(al)){
        seq<-paste(al[i])
        #split into 3 bases from the start (if frame1)
        seq<-substring(seq,seq(1,nchar(seq),3), seq(3, nchar(seq),3) )
        p<-which(seq =="taa"|seq=="tga"|seq=="tag")
        pos<-c(pos,p)
        
}

poss<-unique(pos)
if (length(poss)==0) pos.remove<-0
if (length(poss)!=0) pos.remove<-c(poss*3-2, poss*3-1,poss*3)
keep<-length(as.list(algn1))

algn1<-read.dna("BKpolyomavirus_VP1.fasta", format = "fasta")
algn1<-read.dna("humanparainfluenzavirus3.fasta", format = "fasta")
keep<-c(1:length(as.list(algn1)[[1]]))-pos.remove
algn.trimmed<-algn1[,keep]
a<-as.list(algn1)
length(a[[1]])
write.FASTA(algn.trimmed, file="Trimmed.fasta" )

strsplit(al.removed[[1]][1])

######################################

pos<-c() 
stop.seq.no<-c()
for (i in 1: length(al)){
        seq<-paste(al[i])
        #split into 3 bases from the start (if frame1)
        seq<-substring(seq,seq(1,nchar(seq),3), seq(3, nchar(seq),3) )
        p<-which(seq =="taa"|seq=="tga"|seq=="tag")
        
        for (j in 1:length(p)){
                if (p[j]==length(seq)) pos<-c(pos,p[j])
                else  stop.seq.no<-c(stop.seq.no,i)
        }
}

al.removed<-al[-stop.seq.no]
write.fasta(al.removed,names(al.removed),"BKpolyomavirus.trimmed.fasta" )
names(al.removed,names)

