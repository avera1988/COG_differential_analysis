#!/usr/bin/env RScript
##########################################################
#	Script for ploting COG categories into a heatmap
#	Author Aruto Vera
#	May 2020
#	avera@ccg.unam.mx
######################################################################

if (!require("tidyverse")){
  install.packages("tidyverse")
  library(tidyverse)
}
if (!require('pheatmap')){
  install.packages('pheatmap',
                   repos = "http://cran.us.r-project.org")
  library(pheatmap)
}
if (!require("viridis")){
  install.packages("viridis",
                  repos = "http://cran.us.r-project.org")
  library(viridis)
} 

#Reading tab file with COG annot for each genome

file.names <- list.files(pattern = ".*catego.tab$")
names <- gsub("\\..*","",file.names)
my_read_delim <- function(path){
  readr::read_delim(path, "\t", escape_double = FALSE, 
                    trim_ws = TRUE)
}

#Generate a list of these files
list.tables <- purrr::map(file.names,my_read_delim)
names(list.tables) <- names

#Reading COG annotations

COGTSV.file <- list.files(pattern = "*.tsv$")

#Reading all COG annotaitons
COG <- read.delim(COGTSV.file,header=F,col.names = c("COG","Annot"))

#Reducing all tables into a single Dataframe
reducing <- function(x){
  rtab <- left_join(COG,x,by="COG") %>%
    select(matches("Gene"))
}
TotalCOG <- purrr::reduce(map(list.tables,reducing),cbind) %>% 
  cbind(COG) %>%
  select(matches("COG|Annot"),everything())
colnames(TotalCOG) <- gsub(".fasta..*","",gsub("Number.Of.Genes.","",colnames(TotalCOG)))

#Replace Na for 0

TotalCOG <- TotalCOG %>% replace(is.na(.),0)

#Sorting and arranging 


TotalCOG <- TotalCOG %>%
  arrange(Annot)

#Heatmap

#Obtaining the matrix
TotCOGmat <- TotalCOG %>%
  unite(COG_Annot,c("COG","Annot"),sep="_") %>%
  column_to_rownames("COG_Annot")



COGPh <- pheatmap(TotCOGmat)
ggsave(COGPh, file=paste0(COGCat,".heatmap.pdf"))