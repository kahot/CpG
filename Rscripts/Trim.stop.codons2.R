library(ape)
library(seqinr)

fastas<-list.files("fasta_files/")

i=2

al<-read.fasta(paste0("fasta_files/",fastas[i]),as.string=T)

#al<-read.fasta("Data/HPIV1_CDS1.fasta", as.string=T)
pos1<-c() #record the position numbers where Stop codon exists
pos2<-c()  # should record if the stop codon exists at the very end of sequences 
stop.seq.no<-c() # The sequence number which has stop codon(s) in the middle of sequence
stop.pos<-list()
k=1
for (i in 1: length(al)){
        seq<-paste(al[i])
        #split into 3 bases from the start (if frame1)
        seq<-substring(seq,seq(1,nchar(seq),3), seq(3, nchar(seq),3) )
        p<-which(seq =="taa"|seq=="tga"|seq=="tag")
        pos1<-c(pos1,p)
        if (length(p)!=0){
                stop.pos[[k]]<-p
                names(stop.pos)[k]<-i
                k=k+1
                for (j in 1:length(p)){
                        if (p[j]==length(seq)) pos2<-c(pos2,p[j])
                        else  stop.seq.no<-c(stop.seq.no,i)
                } 
        }
}

#which seqeunce has stop codon(s) in where?
stop.pos

#Where do stop codons appear?
poss<-unique(pos1)

#codon position(s)
poss

#nucleotide position(s)
c(poss*3-2, poss*3-1, poss*3)

#which sequence(s) have stop codons?
stop.seq.no



pos<-c()
stop.seq.no<-c()
for (i in 1: length(al)){
        seq<-paste(al[i])
        #split into 3 bases from the start (if frame1)
        seq<-substring(seq,seq(1,nchar(seq),3), seq(3, nchar(seq),3) )
        p<-which(seq =="taa"|seq=="tga"|seq=="tag")
        if (length(p)!=0){
                for (j in 1:length(p)){
                        if (p[j]==length(seq)) pos<-c(pos,p[j])
                        else  stop.seq.no<-c(stop.seq.no,i)
                } 
        }
        
}
poss<-unique(pos)
poss
unique(stop.seq.no)
# if stop.seq.no is empty, al[-stop.seq.no] will erase the whole file. So let it go through the if loops. At the end, 
# you will have "xxxxx.trimmed.fasta" regardless whether you trimmed/removed seqeunce(s)

if (length(stop.seq.no)!=0){
        al<-al[-stop.seq.no]
        write.fasta(al,names(al),"humanparainfluenzavirus3.filtered.fasta" )

        if (length(poss)!=0){
                al.filtered<-read.dna("humanparainfluenzavirus3.filtered.fasta", format = "fasta")
                al.trimmed<-al.filtered[,-c(poss*3-2, poss*3-1,poss*3)]
                write.FASTA(al.trimmed, file="humanparainfluenzavirus3.trimmed.fasta" )}
        else  write.FASTA(al, file="humanparainfluenzavirus3.trimmed.fasta" )
}

if (length(stop.seq.no)==0){
        if (length(poss)!=0){
                al.filtered<-read.dna("humanparainfluenzavirus3.fasta", format = "fasta")
                al.trimmed<-al.filtered[,-c(poss*3-2, poss*3-1,poss*3)]
                write.FASTA(al.trimmed, file="humanparainfluenzavirus3.trimmed.fasta" )}
        else file.copy("humanparainfluenzavirus3.fasta", "humanparainfluenzavirus3.trimmed.fasta" )
}

