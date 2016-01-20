#   CLEAN UP SCRIPTS for the master spread sheet:  phi4==CURATION==CURRENT
#   FILE NAME: merger_phenotypes.R
#   purpose: transform data for fasta extraction from Uniprot and generate FASTA.file
#
#   input files: phi4-master spreadsheet in format csv
#   output files:  phi4.fasta file named "phi4.fas"

# ==install packages==
# install.packages("XLConnect")
# install.packages("dplyr")
# install.packages("magrittr")
# install.packages("tidyr")
# install.packages("stringr")
# install.packages("stringi)

# ==load libraries==
library("dplyr")
library("plyr")
library("dplyr2")
library("magrittr")
library("tidyr")
library("stringr")
library("stringi")
#==set working directory==
setwd("D:/Owncloud_MU/VIMdb/R/phi4pt1")

#==assign path to  phi4 curation spread sheet and read phi4 sheet==
excel.file<-file.path("D:/OwnCloud_MU/VIMdb/R/phi4pt1/phi4.1.csv")

#convert factor variables into character variables to be able to substitute na with ""
phi4=read.csv(excel.file, stringsAsFactors=FALSE, strip.white=TRUE, header=TRUE)
#only read in NEW rows to keep memory requirements low
phi4=phi4[6479:7433,] #for phi4pt1 new accessions: PHI:5023, NO GENE TESTED to PHI:4257, ATR1

#==clean phi_tmp==
#==== strip \r and \n to remove carriage returns and new lines within strings using #gsub("[\r\n]", " ", x)t   
phi_tmp=as.data.frame(apply(phi4,2,function(x)gsub("[\r\n]", " ",x))) #2= columns , 1= rows; 
phi_tmp["Phenotype_of_mutant"] = as.data.frame(sapply(phi_tmp["Phenotype_of_mutant"], tolower))


#---MERGE MIXED PHENOTYPES---
phi3m=aggregate(Phenotype_of_mutant ~ PHI_MolConn_ID, phi_tmp, FUN=function(x){paste(unique(x), collapse="|")})
phi3m$Phenotype_of_mutant=as.factor(phi3m$Phenotype_of_mutant)

                                 
#-- replace aggregated phenotypes "reduced virulence | reduced virulence"  with "mixed outcome" in one column
#-- below needs to be rewritten to allow multiple substitutes:

phi3m$Phenotype_of_mutant=sapply(phi3m$Phenotype_of_mutant, gsub,pattern=".*\\|.*", replacement= "mixed outcome")

#---replace aggregated Phenotype_of_mutant column with the same column in df 'phi3'
#---for this first remove duplicate 'PHI_base_accession' rows in phi3
phi3u=phi_tmp %>% group_by(PHI_MolConn_ID) %>% distinct(PHI_MolConn_ID)
phi3u$Phenotype_of_mutant= phi3m$Phenotype_of_mutant
#delete last empty row, don't know where it comes from
phi3u=head(phi3u, -1)
#---drop columns--i

phi3u[8]=NULL  #drop pesky column "phi3u$AA.sequence..no.EMBL."
#drops <- c("phi3u$AA.sequence..no.EMBL.")
#phi3u[,!(names(phi3u) %in% drops)]
seqs=read.csv("seqs.csv")
#---now merge uniprot sequences with phi3u using vlookup with the plyr package 
phi3u<- join(phi3u, seqs, by ="ProteinID")
write.csv(phi3u, file="D:/OwnCloud_MU/VIMdb/R/phi4pt1/phi4.1_add.csv", row.names=FALSE, na="")

#---all done, need only to drop and reorder columns ---
#do with: subset command. You don't have to provide anything for the 

#subset argument (i.e., you'll keep all the rows) but you can re-order the 
#columns by listing them in the order you want, within the select argument, 
#like this 
#=================MARTIN TO CONTINUE CHANGE HEADERS TO PHI4 headers==============================

include=c("Record ID", "PHI_MolConn_ID", "ProteinID","Gene_name",  "Pathogen_NCBI_species_Taxonomy.ID","Pathogen_species", "Phenotype_of_mutant", "seq")
phi3out<-subset(phi3u,select=names(phi3u) %in% include)
phi3out$id=NULL
library("reshape2")
#within(phi3out,id=paste(PHI.base.accession.no, PHI.base.accession, Gene.name, Accession, Pathogen.NCBI.Taxonomy.ID, Phenotype.of.mutant), sep=" | ", na.rm=TRUE)
phi3out=as.data.frame(na.omit(phi3out))
phi3out$id=paste(phi3out$PHI_ID, phi3out$PHI_MolConn_ID, phi3out$ProteinID, phi3out$Gene_name,
                 phi3out$"Pathogen_NCBI_species_Taxonomy.ID", phi3out$Pathogen_species,phi3out$Phenotype_of_mutant, sep="|")

#========================Save remove superfluous columns and save as fasta file
phi3outtab=as.data.frame(phi3out[,-c(1:5)])  #drop columns 1-6
names(phi3outtab)=c("row","seq","id")
phi3outtab=as.data.frame(phi3outtab[c("id","seq")])    #swap columns
#---remove white space from Taxonomy.id colum
phi3outtab$id=sapply(phi3outtab$id,gsub,pattern="\\|\\s+(.*)\\|", replacement= "\\|\\1\\|",perl=TRUE)


#---sprintf is a lot more powerful than plain concatenation (less "" etc)---
phi3outtab$id=sapply(phi3outtab$id, sub,pattern="\\|", replacement= "")
phi3outtab$id=sapply(phi3outtab$id, gsub,pattern="^", replacement= ">\\1")
#---next convert phi3outtab with two columns to fasta format---
#  awk and python scripts not necessary but possible alterantives
#awk '{print ">"$2"\n"$1}' tabseq.tsv > seqs.fa]
#python tab2fasta.py tabseq.tsv 5 1 2 4  > tabseq.fa
## remove PHI-base accessions without valid uniprot sequences
phi3outtab=na.omit(phi3outtab)
phi3outtab=matrix(rbind(phi3outtab$id,as.character(phi3outtab$seq)))
write.table(phi3outtab, file="phi_add.fas", quote=FALSE, eol="\n", sep="",row.names=FALSE, col.names=FALSE, na="")

