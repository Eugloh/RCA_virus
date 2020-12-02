## Graphics 
# NOVEMBER 2020
rm(list = ls(all = TRUE))
gc()
## plot 

# baseDir <- "~/Documents/ROSARIO/fichiers/RCA_virus"
# SourceDir <- "~/Documents/ROSARIO/fichiers/RCA_virus/src"
# OutputDir <- "~/Documents/ROSARIO/fichiers/RCA_virus/plot"
# DataDir <- "~/Documents/ROSARIO/fichiers/RCA_virus/data"

baseDir <- "/home/lohmanne/Documents/RCA_virus"
SourceDir <- "/home/lohmanne/Documents/RCA_virus/src"
OutputDir <- "/home/lohmanne/Documents/RCA_virus/plot"
DataDir <- "/home/lohmanne/Documents/RCA_virus/data"


setwd(baseDir)
library(optparse)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(reshape2)

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

# write.table(table_analyse, file=paste(baseDir,"lineage","ncbi_lineages_only_analysed.csv", sep = "/"), col.names = TRUE , row.names=FALSE, sep=",")

#colnames(Data[["Abundance_reduced_taxid"]])<- str_replace(colnames(Data[["Abundance_reduced_taxid"]]), "RCA.", "")
#colnames(Data[["Abundance_atleast2seq"]])<- str_replace(colnames(Data[["Abundance_atleast2seq"]]), "RCA.", "")
Data[["Abundance_reduced_taxid"]]$SUM <- rowSums(Data[["Abundance_reduced_taxid"]])
Data[["Abundance_reduced_taxid"]]$taxid <- rownames(Data[["Abundance_reduced_taxid"]])



df <-Data[["Abundance_reduced_taxid"]][order(-Data[["Abundance_reduced_taxid"]]$SUM,-Data[["Abundance_reduced_taxid"]]$RCA.Skin.5,-Data[["Abundance_reduced_taxid"]]$RCA.Skin.7,-Data[["Abundance_reduced_taxid"]]$RCA.Skin.10,-Data[["Abundance_reduced_taxid"]]$RCA.Skin.3,-Data[["Abundance_reduced_taxid"]]$RCA.Skin.8,-Data[["Abundance_reduced_taxid"]]$RCA.Skin.12,-Data[["Abundance_reduced_taxid"]]$RCA.Skin.1,-Data[["Abundance_reduced_taxid"]]$RCA.Skin.4,-Data[["Abundance_reduced_taxid"]]$RCA.Skin.6,-Data[["Abundance_reduced_taxid"]]$RCA.Skin.11,-Data[["Abundance_reduced_taxid"]]$RCA.Skin.2,-Data[["Abundance_reduced_taxid"]]$RCA.Skin.9),]
df2 <-Data[["Abundance_atleast2seq"]][order(-Data[["Abundance_atleast2seq"]]$SUM,-Data[["Abundance_atleast2seq"]]$RCA.Skin.5 ,-Data[["Abundance_atleast2seq"]]$RCA.Skin.7 ,-Data[["Abundance_atleast2seq"]]$RCA.Skin.10 ,-Data[["Abundance_atleast2seq"]]$RCA.Skin.3 ,-Data[["Abundance_atleast2seq"]]$RCA.Skin.8 ,-Data[["Abundance_atleast2seq"]]$RCA.Skin.12 ,-Data[["Abundance_atleast2seq"]]$RCA.Skin.1 ,-Data[["Abundance_atleast2seq"]]$RCA.Skin.4,-Data[["Abundance_atleast2seq"]]$RCA.Skin.6 ,-Data[["Abundance_atleast2seq"]]$RCA.Skin.11 ,-Data[["Abundance_atleast2seq"]]$RCA.Skin.2 ,-Data[["Abundance_atleast2seq"]]$RCA.Skin.9 ),]

df3 <- cbind(df2,df)
df4 <- df3[c(0:25,51)]

df4$ESPECES <- rownames(df4)

rownames(table_analyse) <- table_analyse$tax_id
rownames(df4)<- df4$taxid

dd <- merge(df4,table_analyse,by=0,all=TRUE)

dd$tax_id<-as.character(dd$tax_id)

summary(dd)

names(dd)
#  [1] "Row.names"    "RCA.Skin.5"   "RCA.Vulva.5"  "RCA.Vulva.8"  "RCA.Skin.7"  
#  [6] "RCA.Vulva.3"  "RCA.Vulva.7"  "RCA.Skin.10"  "RCA.Skin.3"   "RCA.Vulva.4" 
# [11] "RCA.Skin.8"   "RCA.Vulva.9"  "RCA.Skin.12"  "RCA.Skin.1"   "RCA.Vulva.1" 
# [16] "RCA.Vulva.6"  "RCA.Skin.4"   "RCA.Skin.6"   "RCA.Vulva.10" "RCA.Vulva.11"
# [21] "RCA.Vulva.12" "RCA.Skin.11"  "RCA.Skin.2"   "RCA.Vulva.2"  "RCA.Skin.9"  
# [26] "SUM"          "taxid"        "ESPECES"      "tax_id"       "superkingdom"
# [31] "phylum"       "class"        "order"        "family"       "genus"       
# [36] "species"      "kingdom"      "no.rank"      "no.rank1"     "no.rank2"    
# [41] "no.rank22"    "no.rank3"     "subclass"     "subfamily"    "subkingdom"  
# [46] "suborder"     "subphylum"    "subspecies"   "subtribe"     "superfamily" 

vect<-which(colnames(dd)%in%c("Row.names" ,"tax_id","phylum","superkingdom", "class", "kingdom", "no.rank2", "no.rank22", "no.rank3", "subclass", "subfamily", "suborder" , "subphylum", "subspecies", "subtribe" ,"superfamily" ))
#"Row.names" ,"tax_id", "phylum", "class", "kingdom", "no.rank2", "no.rank22", "no.rank3", "subclass", "subfamily", "suborder" , "subphylum", subspecies", "subtribe" ,"superfamily" 
Data[["Table"]]<-dd[,-vect]
colnames(Data[["Table"]])
#  [1] "RCA.Skin.5"   "RCA.Vulva.5"  "RCA.Vulva.8"  "RCA.Skin.7"   "RCA.Vulva.3" 
#  [6] "RCA.Vulva.7"  "RCA.Skin.10"  "RCA.Skin.3"   "RCA.Vulva.4"  "RCA.Skin.8"  
# [11] "RCA.Vulva.9"  "RCA.Skin.12"  "RCA.Skin.1"   "RCA.Vulva.1"  "RCA.Vulva.6" 
# [16] "RCA.Skin.4"   "RCA.Skin.6"   "RCA.Vulva.10" "RCA.Vulva.11" "RCA.Vulva.12"
# [21] "RCA.Skin.11"  "RCA.Skin.2"   "RCA.Vulva.2"  "RCA.Skin.9"   "SUM"         
# [26] "taxid"        "ESPECES"      "order"        "family"       "genus"       
# [31] "species"      "no.rank"      "no.rank1"     "subkingdom" 

rownames(Data[["Table"]])<- Data[["Table"]]$ESPECES


for (i in colnames(Data[["Table"]])){
  Data[['Table']][i][which(Data[['Table']][i]==""),] <- NA
}

var_Data <- colnames(
  Data[['Table']][grepl("^RCA",
                   colnames(Data[['Table']]))])


restructuredData = melt(Data[['Table']], # table with all info
  id = c("ESPECES","family","genus"), 
    measure.vars = var_Data, 
    variable_name = "number" )

restructuredData$from_ <- NA
g <- grepl("Skin", restructuredData1$variable)
restructuredData$from_[which(g)] <- "Skin"
restructuredData$from_[-which(g)] <- "Vulva"


## table without all the 0 values 
reducedrestructuredData <- restructuredData[!restructuredData$value==0,]



##################
Data[["skin"]] <- Data[["Table"]][grep("Skin|tax_id|order|family|genus|species|ESPECES",colnames(Data[["Table"]]))]
Data[["vulva"]] <- Data[["Table"]][grep("Vulva|tax_id|order|family|genus|species|ESPECES",colnames(Data[["Table"]]))]

  
Data[["skinreduced"]] <- Data[["skin"]][rowSums(Data[["skin"]]==0, na.rm=TRUE)<ncol(Data[["skin"]]), ] # table with no virus == 0 
Data[["vulvareduced"]] <- Data[["vulva"]][rowSums(Data[["vulva"]]==0, na.rm=TRUE)<ncol(Data[["vulva"]]), ]

Data[["skinreducedsup10"]] <- Data[["skinreduced"]][which(rowSums(Data[["skinreduced"]][,0:12])>10),] # compte > 10
Data[["vulvareducedsup10"]] <- Data[["vulvareduced"]][which(rowSums(Data[["vulvareduced"]][,0:12])>10),]

setwd(DataDir)
write.table(Data[["skinreducedsup10"]],"skin_treshold10.csv",row.names = F,col.names = T)
write.table(Data[["vulvareducedsup10"]],"vulva_treshold10.csv",row.names = F,col.names = T)
# 

