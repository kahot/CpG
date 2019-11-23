library(ape)
library(seqinr)

al<-read.fasta("HCV_Alignment_CDS.fasta",as.string=TRUE)
al2<-al
pos<-c()
for (i in 1: length(al)){
        seq<-paste(al[i])
        #split into 3 bases from the start (if frame1)
        seq<-substring(seq,seq(1,nchar(seq),3), seq(3, nchar(seq),3) )
        p<-which(seq =="taa"|seq=="tga"|seq=="tag")
        seq[p]<-"nnn"
        newseq<-capture.output(cat(seq, sep="")) 
        al2[i]<-newseq
}

write.fasta(al2, names=names(al2), file.out="hcv-no.stop.codons.fasta")

#double checking
for (i in 1: length(al2)){
        seq<-paste(al2[i])
        #split into 3 bases from the start (if frame1)
        seq<-substring(seq,seq(1,nchar(seq),3), seq(3, nchar(seq),3) )
        p<-which(seq =="taa"|seq=="tga"|seq=="tag")
        pos<-c(pos,p)
}


