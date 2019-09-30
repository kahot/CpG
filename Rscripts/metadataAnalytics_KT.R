viruses <- read.csv("Data/virusInstructionFile.csv")

#transparency function
transp <- function(col, alpha=.5){
  res <- apply(col2rgb(col),2, function(c) rgb(c[1]/255, c[2]/255, c[3]/255, alpha))
  return(res)
}

#plot tha makes a histogram and cumulative percent of sequence length
distributionOfLength <- function(i) {
  exampleVirus <- read.csv( paste0("metadata/metadata-", gsub(" ", "", viruses[i,1]), ".csv"))
  test <- hist(exampleVirus$seqLength, plot=F, breaks=25)
  test.h <- cumsum(test$counts/(sum(test$counts))) *100
  plot(test, col="dodgerblue3", main=paste(viruses[i,1], "\ndistribution of seqence length"),
       las=2, ylab="frequency", xlab="", border="dodgerblue3")
  par(new = T)
  plot(x=test$mids, y=test.h, type="l", axes=F, col=transp("firebrick",0.8), lwd=3,  xlab="length (bp)", ylab="", ylim=c(0,100))
  box()
  axis(4, seq(0, 100, 20), las=1)
  mtext("cumulative % of sequences", 4, line=2.5, las=3, cex=0.8)

}

containsMetadata <- function() {
  
  seqsWithDate <- c()
  seqsWithLoc <- c()
  seqsWithDateAndLoc <- c()
  
  for(i in 2:19) {
    exampleVirus <- read.csv( paste0("metadata/metadata-", gsub(" ", "", viruses[i,1]), ".csv"))
    #collate fraction of sequences with date
    seqsWithDate[i-1] <- sum(is.na(exampleVirus$collectionDate)==FALSE)
    #collate fraction of sequences with location
    seqsWithLoc[i-1] <- sum(is.na(exampleVirus$location)==FALSE)
    #fraction of sequences with date and location
    seqsWithDateAndLoc[i-1] <- nrow(exampleVirus[is.na(exampleVirus$collectionDate)==F & is.na(exampleVirus$location)==F, ]) 
  }

  #using subset for now
  viruses2 <- viruses[2:19,]
  virusNames=viruses2[,1]
  virusesWithoutDate <- data.frame(seqsWithDateAndLoc,  seqsWithLoc, seqsWithDate, totSeqs=viruses2[,3]  )
   
  virusesWithoutDate <- data.frame(seqsWithDate, totSeqs=viruses2[,3]-seqsWithDate )
  virusesWithoutLoc <- data.frame(seqsWithLoc, totSeqs=viruses2[,3]-seqsWithLoc )
  virusesWithoutDateAndLoc <- data.frame(seqsWithDateAndLoc, totSeqs=viruses2[,3]-seqsWithDateAndLoc )
  
  par(mfrow=c(1,4), las=1)
  par(mar=c(4,13,3,0))
  barplot(t(virusesWithoutDate), names.arg=virusNames, col="white", horiz=T, border=NA, axes=F)
  par(mar=c(4,0,3,1))
  #make cumulative barplot
  barplot(t(virusesWithoutDate),  horiz=T, main="Seqs with Date", col=c("dodgerblue4", "dodgerblue2" ), border=NA) #names.arg=virusNames,
  barplot(t(virusesWithoutLoc), horiz=T, main="Seqs with Location", col=c("dodgerblue4", "dodgerblue2" ), border=NA)
  barplot(t(virusesWithoutDateAndLoc),  horiz=T, main="Seqs with Date and Location", col=c("dodgerblue4", "dodgerblue2" ), border=NA)
  legend("topright", legend=c("has info", "missing info"),
         col=c("dodgerblue4", "dodgerblue2" ), border=NA,  pch=15, pt.cex=2, bty="n")
}

removeShortSeqs <- function(metaDataTable, lengthLimit) {
  metaDataTable <- as.data.frame(metaDataTable)
  metaDataTable <- metaDataTable[as.numeric(as.character(metaDataTable$seqLength)) > lengthLimit,] 
  return(metaDataTable)
}

removeNoLocation <- function(metaDataTable) {
  metaDataTable <- as.data.frame(metaDataTable)
  metaDataTable <- metaDataTable[is.na(metaDataTable$location)==F, ]
  return(metaDataTable)
}

removeNoDate <- function(metaDataTable) {
  metaDataTable <- as.data.frame(metaDataTable)
  metaDataTable <- metaDataTable[is.na(metaDataTable$collectionDate)==F, ]
  return(metaDataTable)
}

continentLookUp <- function(metaDataTable) {
  metaDataTable <- as.data.frame(metaDataTable)
  countries <- read.csv("continentLookup.csv", header=T, colClasses="character", na.strings="none")
  metaDataTable <- data.frame(metaDataTable, continent=countries[match(metaDataTable$location, countries$country), 2])
  return(metaDataTable)
}

require(lubridate)
#uses lubridate function to format the miscellaneous date formats, and just return the year. 
#NB, shouldnt be used for anything else with date, because where Month or day are unknown, it just gives 1st of the 1st
#this is fine for extracting date though
addYearOnly <- function(metaDataTable) {
  testDates <- as.character(metaDataTable$collectionDate)
  years <- year(parse_date_time(testDates, c("dby", "dbY", "by", "y", "Ymd", "Ym", "mdY")) )
  metaDataTable <- data.frame(metaDataTable, year=years)
}
