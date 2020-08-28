# COG barplot: compare Dorea vs PAL genome COG datasets
# Using tab delimited "func_stats.txt" files produced from COG output for each genome
# all "PublicDoreaGenomes" obtained from RefSeq database
# see COG_BarPlot.rmd for source markdown document

# load libraries:

library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(stringr)

#make function to read tables

my_read_delim <- function(path,name){
  readr::read_delim(path, "\t", escape_double = FALSE, 
                    trim_ws = TRUE,col_names = c("COG","Func",name))
}

# list all .txt files and shorten to "Names" w/o extensions (after 1st ".")

list.filenames <- list.files(pattern = "tab$")
Names <- gsub(".out","",gsub(".func_stats.tab.*","",list.filenames))

# apply new my_read_delim function to obtain a list with all data frames

list.Data <- purrr::map2(list.filenames,Names,my_read_delim)

# name all data frames in the list

names(list.Data) <- Names
head(list.Data$PAL_113)

# remove unecessary "COG" column 

list.data2 <- lapply(list.Data,function(x){x[,c(2,3)]})
head(list.data2$PAL_113)

# merge all dataframes into 1 big df by the "Func" column

CogValues <- purrr::reduce(list.data2,full_join,by="Func")
head(CogValues)

# compare no. of genes (normalize by rel proportion of genes)

sapply(names(CogValues)[-1], function(x) {
  CogValues[paste0(x,".pct")] <<- (CogValues[x]*100/sum(CogValues[x]))
})

# select proportion of genes ("pct" columns) 
# new df contains rel freq percent values instead of raw counts

CogValuesMod <- dplyr::select(CogValues, matches("Func|pct"))
head(CogValuesMod)

# filter out rows for "extracellular" and "nuclear" cog functions

CogValuesMod <- dplyr::filter(CogValuesMod,str_detect(Func,"Nuclear", negate = TRUE)) %>%
  filter(str_detect(Func,"Extracellular", negate = TRUE))

# copy CogValuesMod df into duplicate df for manipulation

CV <- CogValuesMod

# separate objects into 2 new dfs based on genome names to compare
# we are comparing "Dorea" and "PAL" data

CV_Dorea <- dplyr::select(CV,matches("(Dorea|Func)"))
CV_PAL <- dplyr::select(CV,matches("(Func|PAL)")) 

# shift tables from wide to long format; add column categories & genome source
# for PAL:

names_PAL <- colnames(CV_PAL)[-1]
CV_PAL_mod <- tidyr::gather(CV_PAL,key="Genome",value="Numberofgenes",names_PAL)
CV_PAL_mod$Source <- "P.americana"
head(CV_PAL_mod)

# for other Dorea:

names_Dorea <- colnames(CV_Dorea)[-1]
CV_Dorea_mod <- tidyr::gather(CV_Dorea,key="Genome",value="Numberofgenes",names_Dorea)
CV_Dorea_mod$Source <- "PublicDoreaGenomes"
head(CV_Dorea_mod)

# bind both PAL and Dorea dfs into 1 joint df

CVtotal <- dplyr::bind_rows(CV_PAL_mod,CV_Dorea_mod)

# group by COG function categories and genome source
# get avg and stdev of genes in each category for both source groups (PAL and Dorea)

CVtotalStats <- CVtotal %>% 
  group_by(Func,Source) %>%
  dplyr::summarise(Average =mean(Numberofgenes),
                   stde=sd(Numberofgenes)/sqrt(length(Numberofgenes))) %>%
  dplyr::rename( "Numberofgenes"=Average) #rename not working?
head(CVtotalStats)

# now to make the 2 color barplot with ggplot2

BarPlot <- ggplot(data=CVtotalStats,aes(x=reorder(Func, +Numberofgenes),
                                         y=Numberofgenes,fill=Source)) +
  geom_bar(stat="identity",position=position_dodge())+
  labs(x="COG",y="Proportion of Genes") + 
  geom_errorbar(aes(ymin=Numberofgenes-stde, ymax=Numberofgenes+stde), width=.5,
                position=position_dodge(.6))+
  scale_fill_brewer(palette = "Paired")+
  coord_flip()+
  theme_bw()

BarPlot

ggsave(BarPlot,file="COG.Dorea.pdf")

colSums(CogValues[,2:17])

Model <- lm(Numberofgenes ~ Func+Source +Func:Source, data = CVtotal)
res.aov <- aov(Model)
tuck <- TukeyHSD(res.aov)
Sum.AOV <- summary(res.aov)

totoTotal <- CVtotal
levfunci<-unique(totoTotal$Func)
TTest <- function(funcion){
  D <- filter(totoTotal,grepl(funcion,Func))
  Pt <- pairwise.t.test(D$Numberofgenes,D$Source, p.adjust.method = "none")
  return(Pt)
} 
TtestTable <- lapply(levfunci,function(x){
  TTest(x)
})
names(TtestTable) <- levfunci
cat(capture.output(print(TtestTable),file="Ttestwandwiout.final.txt"))
