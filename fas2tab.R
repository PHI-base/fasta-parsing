#MU 2015-05-10     Read fasta sequences
##ubuntu 25 machine connected to Owncloud_MU
#file: 2_makePHIfasta_file.R
#purpose: takes output of uniprot_retrieve.pl called "seqs" and coverts into 2-column dataframe
#input: reads mutiple uniprot fasta file
#output: saves   "phi3.fas" file, seqs
#================================================================================================
#--Environment----
#install.packages("dplyr")
library(magrittr)
library(plyr)
library(dplyr)
library(reshape)
setwd("D:/OwnCloud_MU/VIMdb/R/phi4pt1")
#--var declaration---
fastaProtfile="phi41_add.txt"
#http://stackoverflow.com/questions/26843995/r-read-fasta-files-into-data-frame-using-base-r-not-biostrings-and-the-like
#-- Function
ReadFasta<-function(file) {
   # Read the file line by line
   fasta<-readLines(file)
   # Identify header lines
   ind<-grep(">", fasta)
   # Identify the sequence lines
   s<-data.frame(ind=ind, from=ind+1, to=c((ind-1)[-1], length(fasta)))
   # Process sequence lines
   seqs<-rep(NA, length(ind))
   for(i in 1:length(ind)) {
      seqs[i]<-paste(fasta[s$from[i]:s$to[i]], collapse="")
   }
   # Create a data frame 
   DF<-data.frame(name=gsub(">", "", fasta[ind]), sequence=seqs)
   # Return the data frame as a result object from the function
   return(DF)
}
#--MAIN---------------------------------------------------------
seqs<-ReadFasta(fastaProtfile)
colnames(seqs)=c("Accession","seq")
nrow(seqs)

seqs$Accession=sapply(seqs$Accession,gsub,pattern="[a-z][a-z]\\|(.*)\\|.*$", replacement= "\\1",perl=TRUE)
tail(seqs$Accession,5)
seqs=tbl_df(seqs)
i <- sapply(seqs, is.factor)   #put all factor variables into character variables
seqs[i] <- lapply(seqs[i], as.character)
write.csv(seqs, file="seqs.csv", row.names=FALSE, na="")
