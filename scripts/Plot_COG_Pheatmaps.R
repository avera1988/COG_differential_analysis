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
if (!require("pheatmap")){
  install.packages("pheatmap")
  library(pheatmap)
}
if (!require("viridis")){
  install.packages("viridis")
  library(viridis)
} 

#Reading tab file with COG annot for each genome

file.names <- list.files(pattern = ".*tab$")
names <- gsub(".fasta..*","",file.names)
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
TotCOGforMAt <- TotalCOG %>%
  unite(COG_Annot,c("COG","Annot"),sep="_")

rownames(TotCOGforMAt) <- TotCOGforMAt$COG_Annot
TotCOGmat <- TotCOGforMAt[,2:18]

pheatmap(TotCOGmat)

#Generating a DataFrame for annotation with the source

ph_annot <- data.frame(Genome=colnames(TotalCOG)[3:19])
ph_annot <- ph_annot %>%
  mutate(Source = if_else(grepl("Dorea",Genome),"Public","P.americana"))
ph_annotMod <- data.frame(Source=ph_annot$Source,row.names =ph_annot$Genome )

#Adding Colors

anot_colors <- list(Source=c(Public="#2166AC",
                             P.americana="#92C5DE"))
Colors <-RColorBrewer::brewer.pal(n=9,"Greys")
COGPh <- pheatmap(TotCOGmat, 
                  color= Colors,
                  cluster_rows = F,
                  cluster_cols = F,
                  treeheight_row = 0,
                  annotation_col  = ph_annotMod,
                  annotation_colors = anot_colors,
                  breaks = c(0,0.25,0.5,1,2,4,6,8,10))
COGCat <- gsub(".tsv","",COGTSV.file)

ggsave(COGPh, file=paste0(COGCat,".heatmap.pdf"),width = 15,height = 15)
pdf("scale.pdf")
plot(NULL, xlim=c(0,length(Colors)), ylim=c(0,1), 
     xlab="", ylab="", xaxt="n", yaxt="n")
rect(0:(length(Colors)-1), 0, 1:length(Colors), 1, col=Colors)
dev.off()
