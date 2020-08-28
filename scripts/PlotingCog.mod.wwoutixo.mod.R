#! /usr/bin/env Rscript
#############################################################################
# From Spiroplasma ixodetis of Dactylopius project
#This script produce a barplot to compare the differences in COG categoires 
#  Using tab delimited "func_stats.txt" files produced from COG output for each genome
######################################################################################


library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(stringr)
library(DescTools)

list.filenames <- list.files(pattern = "tab$")
list.data<-list()
for (i in 1:length(list.filenames)){
  list.data[[i]]<-read.table(list.filenames[i],sep="\t",header=F,quote="")
}
Names <- gsub("\\..*","",list.filenames)
names(list.data) <- Names
for (i in 1:length(list.filenames)){
  colnames(list.data[[i]]) <- c("COG","Func",Names[i])
}

list.data3 <- lapply(list.data,function(x){x[,c(2,3)]})

CogValues2 <- Reduce(merge,list.data3)

sapply(names(CogValues2)[-1], function(x) {
  CogValues2[paste0(x,"pct")] <<- (CogValues2[x]/sum(CogValues2[x]))
})
CogValuesMod2 <- dplyr::select(CogValues2, matches("Func|pct"))
CogValuesMod2 <- dplyr::filter(CogValuesMod2,str_detect(Func,"Nuclear", negate = TRUE)) %>%
  filter(str_detect(Func,"Extracellular", negate = TRUE))
Names2 <- colnames(CogValuesMod2)[-1]

toto <- CogValuesMod2
totoDacty <- dplyr::select(toto,matches("(DC|DO|Func)"))
totoPublic <- dplyr::select(toto,matches("(Func|Spi)")) %>%
  select(-matches("WSS|Danaus"))

totoIxo <- dplyr::select(toto,matches("(DC|DO|WSS|Danaus|Func)"))



#write.table(CogValues2,"COG.spiroplasma.tab",sep="\t",row.names=F,quote=F)
write.csv(CogValues2,"COG.Spiroplasma.csv")

##### New idea adding SD and error bar####
NamesDacty <- colnames(totoDacty)[-1]

totoDactymod <- gather(totoDacty,key="Genome",value="Numberofgenes",NamesDacty)
totoDactymod$Source <- "S.ixodetisDactylopius"#"Dacty"

NamesNoDacty <- colnames(totoPublic)[-1]

totoNoDacty <- gather(totoPublic,key="Genome",value="Numberofgenes",NamesNoDacty)
totoNoDacty$Source <- "OtherSpiroplasma" #"NoDacty"

NamesIxo <- colnames(totoIxo)[-1]

totoIxomod <- gather(totoIxo,key="Genome",value="Numberofgenes",NamesIxo)
totoIxomod$Source <- "SpiroplasmaIxodetis"

totoTotal <- dplyr::bind_rows(totoDactymod,totoNoDacty,totoIxomod)

totoTolaStats <- totoTotal %>% group_by(Func) %>% 
  group_by(Func,Source) %>%
  dplyr::summarise(Average =mean(Numberofgenes),
                   stde=sd(Numberofgenes)/sqrt(length(Numberofgenes))) %>%
  dplyr::rename("Numberofgenes"=Average)


COgSortedPlot <- ggplot(data=totoTolaStats,aes(x=reorder(Func, +Numberofgenes),
                                           y=Numberofgenes,fill=Source)) +
  geom_bar(stat="identity",position=position_dodge())+
  labs(x="COG",y="Proportion of Genes") + 
  geom_errorbar(aes(ymin=Numberofgenes-stde, ymax=Numberofgenes+stde), width=.5,
                position=position_dodge(.6))+
  scale_fill_brewer(palette = "Set2")+
  coord_flip()+
  theme_bw()
ggsave(COgSortedPlot,filename = "COG.Spiroplasma.w.errorbars.color.wandwiout.pdf")
write.table(totoTolaStats,"COG.Statswandwiout.tab",sep="\t",quote=F,row.names = F)
write.csv(totoTolaStats,"COG.Statswandwiout.csv",quote=F,row.names = F)
Model <- lm(Numberofgenes ~ Func+Source +Func:Source, data = totoTotal)
res.aov <- aov(Model)
tuck <- TukeyHSD(res.aov)
Sum.AOV <- summary(res.aov)
write.table(tuck$`Func:Source`,"COG.Tuckey.Func.Sourcewandwiout.tsv",sep="\t",quote = F)
#Applying paired t-test

totoTotal.mod <- totoTotal
levfunci<-unique(levels(totoTotal$Func))
TTest <- function(funcion){
  D <- filter(totoTotal,grepl(funcion,Func))
  Pt <- pairwise.t.test(D$Numberofgenes,D$Source, p.adjust.method = "none")
  return(Pt)
} 
TtestTable <- lapply(levfunci,function(x){
  TTest(x)
})
names(TtestTable) <- levfunci
TtestTable$`Amino acid transport and metabolism`
cat(capture.output(print(TtestTable),file="Ttestwandwiout.final.txt"))
