#install.packages("seqinr")
library(seqinr)
library(rentrez)

### fasta handling

genbankQueryFasta <- function(accCodeChunk){
  rawFasta <- entrez_fetch(db = "nucleotide", id = accCodeChunk, rettype = "fasta")
  semiParsedFasta <- read.fasta(file = textConnection(rawFasta), as.string = T)
  
  return(semiParsedFasta)
}

### sticking together fastas
genbankFasta <- function(accCodes) {
  
  accCodes <- as.character(accCodes)
  
  fastas <- c()
  
  if(length(accCodes) > 500) {     #chopUp
    
    numChunks <- ceiling(length(accCodes)/500)
    maxCodes <- length(accCodes)
    
    for(k in 1:numChunks) {  
      if(k*500 > maxCodes) {
        #handle last case specially
        accCodeChunk <- accCodes[(k*500-499):maxCodes]
      } else {
        #send 500
        accCodeChunk <- accCodes[(k*500-499):(k*500)]
      }
      #collect result
      fastas <- c(fastas, genbankQueryFasta(accCodeChunk))
    }  
  } else {
    fastas <- genbankQueryFasta(accCodes)  
  } 
  
  return(fastas)
}

separateSeqs <- function(originalNames, fastaFile) {
  seqsOnly <- list()
  for(k in 1:length(originalNames)) {
    fastaElem <- fastaFile[names(fastaFile) == originalNames[k]]
    fastaElem[1]
  }
  
  outputSeqs <- list()
  for(i in 1:numberOfoutputSeqs){
    # Extract the fasta headers 
    inputSeqsIndexNames[i] <- attr(fastaFile[[originalNames[k]]],"name")
    # Extract the sequence
    outputSeqs[[i]] <- inputSeqs[[inputSeqsIndex[i]]][initialPos:((initialPos+lengthOutputSeqs)-1)]
  }
}