## Graphics 
# NOVEMBER 2020
rm(list = ls(all = TRUE))
gc()
## plot 

baseDir <- "~/Documents/ROSARIO/fichiers/RCA_virus"
SourceDir <- "~/Documents/ROSARIO/fichiers/RCA_virus/src"
OutputDir <- "~/Documents/ROSARIO/fichiers/RCA_virus/plot"
DataDir <- "~/Documents/ROSARIO/fichiers/RCA_virus/data"

# baseDir <- "/home/lohmanne/Documents/RCA_virus"
# SourceDir <- "/home/lohmanne/Documents/RCA_virus/src"
# OutputDir <- "/home/lohmanne/Documents/RCA_virus/plot"
# DataDir <- "/home/lohmanne/Documents/RCA_virus/data"


setwd(baseDir)
library(optparse)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(reshape2)
library(scales)
library(tidyverse)
library(rlist)

option_list = list(
  make_option(c("-p", "--pattern"), type="character", default="csv", help="pattern for file names [default= %default]", metavar="character")
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

csvFile  <- grep(opt$pattern,list.files(DataDir),value=TRUE) # find names of every data files 

Data<- list()
for (el in csvFile){
  el_name <- sapply(strsplit(as.character(el), "\\."), "[[", 1)
  Data[[el_name]] <- read.delim(paste(DataDir,el, sep = "/"), row.names=1)
}
ncbi_virus_lineages <- read.csv(paste(baseDir,"lineage","ncbi_lineages-virus.csv", sep = "/"))

# j'enleve les colonnes qui ne contiennes que des NA 
### Colonne contenant que des NA
lineage <- ncbi_virus_lineages[, apply(ncbi_virus_lineages, 2, function(x) !all(is.na(x)))]
# write.table(lineage, file=paste(baseDir,"lineage","ncbi_lineages-virus_short.csv", sep = "/"), col.names = TRUE , row.names=FALSE, sep=",")

x <- grep(paste(rownames(Data[["Abundance_reduced_taxid"]]),  collapse="|"),lineage$tax_id) # lineage position en récupérant par TAXID
table_analyse<-lineage[x,][which(lineage[x,]$tax_id%in%rownames(Data[["Abundance_reduced_taxid"]])),] # lineage en ne gardant que les ligne d'intéret 



df <- Data[["Abundance_reduced_taxid"]]
nom_colonnes <- colnames(table_analyse)
## filtre de la table d'analyse lineage 
#  [1] "tax_id"       "superkingdom" "phylum"       "class"        "order"       
#  [6] "family"       "genus"        "species"      "kingdom"      "no.rank"     
# [11] "no.rank1"     "no.rank2"     "no.rank22"    "no.rank3"     "subclass"    
# [16] "subfamily"    "subkingdom"   "suborder"     "subphylum"    "subspecies"  
# [21] "subtribe"     "superfamily" 

for (i in colnames(table_analyse)){
  table_analyse[i][which(table_analyse[i]==""),] <- NA
}


table_analyse$taxonomy<-NA
for (el in 1:nrow(table_analyse)){
  name="Viruses"
  for (i in c("family","subkingdom","genus","species")){

  name=paste(name,table_analyse[el,i],sep=";")
  
  }
  table_analyse$taxonomy[el]<-name
}


df$tax_id <- rownames(df)
df<- transform(df, tax_id = as.numeric(tax_id))
df <-df[order(df$tax_id),]
table_analyse <- table_analyse[order(table_analyse$tax_id),]

dd<- bind_cols(df, table_analyse)


vect<-which(colnames(dd)%in%c("Row.names" ,"tax_id","phylum","superkingdom", "class", "kingdom", "no.rank2", "no.rank22", "no.rank3", "subclass", "subfamily", "suborder" , "subphylum", "subspecies", "subtribe" ,"superfamily" ))
#"Row.names" ,"tax_id", "phylum", "class", "kingdom", "no.rank2", "no.rank22", "no.rank3", "subclass", "subfamily", "suborder" , "subphylum", subspecies", "subtribe" ,"superfamily" 
dd<-dd[,-vect]
colnames(dd)


tmpdd <-dd[order(-dd$RCA.Skin.5,-dd$RCA.Skin.7,-dd$RCA.Skin.10,-dd$RCA.Skin.3,-dd$RCA.Skin.8,-dd$RCA.Skin.12,-dd$RCA.Skin.1,-dd$RCA.Skin.4,-dd$RCA.Skin.6,-dd$RCA.Skin.11,-dd$RCA.Skin.2,-dd$RCA.Skin.9),]
tmpNamedtable <-Data[["Abundance_atleast2seq"]][order(-Data[["Abundance_atleast2seq"]]$RCA.Skin.5 ,-Data[["Abundance_atleast2seq"]]$RCA.Skin.7 ,-Data[["Abundance_atleast2seq"]]$RCA.Skin.10 ,-Data[["Abundance_atleast2seq"]]$RCA.Skin.3 ,-Data[["Abundance_atleast2seq"]]$RCA.Skin.8 ,-Data[["Abundance_atleast2seq"]]$RCA.Skin.12 ,-Data[["Abundance_atleast2seq"]]$RCA.Skin.1 ,-Data[["Abundance_atleast2seq"]]$RCA.Skin.4,-Data[["Abundance_atleast2seq"]]$RCA.Skin.6 ,-Data[["Abundance_atleast2seq"]]$RCA.Skin.11 ,-Data[["Abundance_atleast2seq"]]$RCA.Skin.2 ,-Data[["Abundance_atleast2seq"]]$RCA.Skin.9 ),]

dd1 <- cbind(tmpNamedtable,tmpdd)

# colnames(dd1[c(0:24,52,56:58,67,73)])
#  [1] "RCA.Skin.5"   "RCA.Vulva.5"  "RCA.Vulva.8"  "RCA.Skin.7"   "RCA.Vulva.3" 
#  [6] "RCA.Vulva.7"  "RCA.Skin.10"  "RCA.Skin.3"   "RCA.Vulva.4"  "RCA.Skin.8"  
# [11] "RCA.Vulva.9"  "RCA.Skin.12"  "RCA.Skin.1"   "RCA.Vulva.1"  "RCA.Vulva.6" 
# [16] "RCA.Skin.4"   "RCA.Skin.6"   "RCA.Vulva.10" "RCA.Vulva.11" "RCA.Vulva.12"
# [21] "RCA.Skin.11"  "RCA.Skin.2"   "RCA.Vulva.2"  "RCA.Skin.9"   "superkingdom"
# [26] "family"       "genus"        "species"      "subkingdom"   "taxonomy"    


dd2 <- dd1[c(0:24,52:55,58,59)]

dd2$ESPECES <- rownames(dd2)


write.csv(dd2,paste(DataDir,"table_compte_more_than2.csv",sep="/"))
#|||||||||||||||||||||||||||||||||||||||||||||||||-


var_Data <- colnames( dd2[grepl("^RCA",colnames(dd2))])


restructuredData = melt( dd2, # table with all info
  id = c("ESPECES","family","subkingdom","genus","species"),
    measure.vars = var_Data, 
    variable_name = "number" )



levels(restructuredData$variable) <- c( "RCA.Skin.1","RCA.Skin.2", "RCA.Skin.3","RCA.Skin.4","RCA.Skin.5", "RCA.Skin.6" ,"RCA.Skin.7", "RCA.Skin.8","RCA.Skin.9",
"RCA.Skin.10", "RCA.Skin.11","RCA.Skin.12"  , "RCA.Vulva.1" ,"RCA.Vulva.2" ,"RCA.Vulva.3","RCA.Vulva.4","RCA.Vulva.5" , "RCA.Vulva.6" ,"RCA.Vulva.7" , "RCA.Vulva.8" , "RCA.Vulva.9",       
    "RCA.Vulva.10", "RCA.Vulva.11" ,"RCA.Vulva.12"  )

restructuredData$from_ <- NA
g <- grepl("Skin", restructuredData$variable)
restructuredData$from_[which(g)] <- "Skin"
restructuredData$from_[-which(g)] <- "Vulva"


## table without all the 0 values 
reducedrestructuredData <- restructuredData[!restructuredData$value==0,]
rownames(reducedrestructuredData) <- 1:nrow(reducedrestructuredData)


##################
Data[["skin"]] <- dd2[grep("Skin|tax_id|order|family|genus|species|ESPECES",colnames(dd2))]
Data[["vulva"]] <- dd2[grep("Vulva|tax_id|order|family|genus|species|ESPECES",colnames(dd2))]

  
Data[["skinreduced"]] <- Data[["skin"]][rowSums(Data[["skin"]]==0, na.rm=TRUE)<ncol(Data[["skin"]]), ] # table with no virus == 0 
Data[["vulvareduced"]] <- Data[["vulva"]][rowSums(Data[["vulva"]]==0, na.rm=TRUE)<ncol(Data[["vulva"]]), ]

Data[["skinreducedsup10"]] <- Data[["skinreduced"]][which(rowSums(Data[["skinreduced"]][,0:12])>10),] # compte > 10
Data[["vulvareducedsup10"]] <- Data[["vulvareduced"]][which(rowSums(Data[["vulvareduced"]][,0:12])>10),]

write.table(Data[["skinreducedsup10"]],paste(DataDir,"skin_treshold10.csv",sep="/"),row.names = F,col.names = T)
write.table(Data[["vulvareducedsup10"]],paste(DataDir,"vulva_treshold10.csv",sep="/"),row.names = F,col.names = T)
write.table(reducedrestructuredData,paste(DataDir,"OUTPUT_to_ExPanD.csv",sep="/"),row.names = F,col.names = T)

#######
# reproduction rmd
df <- dd2[c(0:24,30)]

# Tableau de comptage final
taxonomy_split <- list()

for (line in df$taxonomy){
  element_combi <- ""
  elements <- unlist(str_split(line, ";"))
  elements <- elements[ elements != ""]
  for (element in elements){
    if (element_combi == ""){
      element_combi <- element
      taxonomy_split <- list.append(taxonomy_split, element_combi)
    }
    else{
      element_combi <- paste(element_combi, element, sep=";")
      taxonomy_split <- list.append(taxonomy_split, element_combi)
    }
  }
}
taxonomy_split <- unique(unlist(taxonomy_split))

Table_count <- data.frame(matrix(0L,ncol=25, nrow=length(taxonomy_split)))
colnames(Table_count) <- c("taxonomy",levels(as.factor(restructuredData$variable)))
Table_count$taxonomy <- taxonomy_split


for (line in Table_count$taxonomy){
  str_to_search <- str_replace_all(noquote(line), "[(]","[(]")
  str_to_search <- str_replace_all(str_to_search, "[)]","[)]")
  df$get <- str_detect(df$taxonomy, str_to_search)
  line[2:25] <- c(colSums(df[df$get == "TRUE",1:24]))
  Table_count[Table_count$taxonomy == line[1],] <- line
  #grace au vecteur found on sait quelle ligne sont a recuperer dans df et a additionner dans Table_count
  
}

write.csv(Table_count, paste(DataDir,"matrice_comptage_final.txt",sep="/"))

## TEST plots
# > names(reducedrestructuredData)
# [1] "ESPECES"    "family"     "subkingdom" "genus"      "species"   
# [6] "variable"   "value"      "from_"   

p<-ggplot(data=reducedrestructuredData, aes(x=from_, y=value, fill=variable)) +
  geom_bar(stat="identity")
ggplot(data=reducedrestructuredData, aes(x=from_, y=value, fill=genus)) +
  geom_bar(stat="identity")
ggplot(data=reducedrestructuredData, aes(x=from_, y=value, fill=family)) +
  geom_bar(stat="identity")
ggplot(data=reducedrestructuredData, aes(x=from_, y=value, fill=ESPECES)) +
  geom_bar(stat="identity")


p<-ggplot(data=reducedrestructuredData, aes(x=variable, y=value, fill=genus)) +
  geom_bar(stat="identity") + scale_y_continuous(labels = percent)



reducedrestructuredData %>%
     mutate(text = text %>% forcats::fct_reorder(variable)) %>%
     ggplot(aes(x = text, y = count, fill = count)) + 
     geom_col() +
     facet_wrap(~ group,  scales = "free_y") +
     coord_flip() 



reducedrestructuredData %>% arrange( desc(from_))%>%
     ggplot(aes(x=variable, y=value, fill=family)) +
  geom_bar(position = position_fill(), stat = "identity") +
   theme(axis.text.x = element_text(angle = 45, hjust=1))+
       scale_y_continuous(labels = scales::percent_format())
## plots
# > names(reducedrestructuredData)
# [1] "ESPECES"  "family"   "genus"    "variable" "value"    "from_"  


reducedrestructuredData %>% arrange( desc(from_))%>%
     ggplot(aes(x=genus, y=value, fill=from_)) +
  geom_bar(position = position_fill(), stat = "identity") +
   theme(axis.text.x = element_text(angle = 45, hjust=1))+
       scale_y_continuous(labels = scales::percent_format())

reducedrestructuredData %>% arrange( desc(from_))%>%
     ggplot(aes(x=genus, y=value, fill=from_)) +
  geom_bar( stat = "identity") +
   theme(axis.text.x = element_text(angle = 45, hjust=1))


