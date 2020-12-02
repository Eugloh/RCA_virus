## Graphics 
# NOVEMBER 2020
rm(list = ls(all = TRUE))
gc()
## plot 

baseDir <- "~/Documents/ROSARIO/fichiers/RCA_virus"
SourceDir <- "~/Documents/ROSARIO/fichiers/RCA_virus/src"
OutputDir <- "~/Documents/ROSARIO/fichiers/RCA_virus/plot"
DataDir <- "~/Documents/ROSARIO/fichiers/RCA_virus/data"

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
# [1] "Row.names"    "Skin.5"       "Vulva.5"      "Vulva.8"      "Skin.7"       "Vulva.3"      "Vulva.7"      "Skin.10"      "Skin.3"       "Vulva.4"     
# [11] "Skin.8"       "Vulva.9"      "Skin.12"      "Skin.1"       "Vulva.1"      "Vulva.6"      "Skin.4"       "Skin.6"       "Vulva.10"     "Vulva.11"    
# [21] "Vulva.12"     "Skin.11"      "Skin.2"       "Vulva.2"      "Skin.9"       "tax_id"       "superkingdom" "phylum"       "class"        "order"       
# [31] "family"       "genus"        "species"      "kingdom"      "no.rank"      "no.rank1"     "no.rank2"     "no.rank22"    "no.rank3"     "subclass"    
# [41] "subfamily"    "subkingdom"   "suborder"     "subphylum"    "subspecies"   "subtribe"     "superfamily" 

vect<-which(colnames(dd)%in%c("phylum","superkingdom", "class", "kingdom", "no.rank2", "no.rank22", "no.rank3", "subclass", "subfamily", "suborder" , "subphylum", "subspecies", "subtribe" ,"superfamily" ))
# "phylum", "class", "kingdom", "no.rank2", "no.rank22", "no.rank3", "subclass", "subfamily", "suborder" , "subphylum", subspecies", "subtribe" ,"superfamily" 
Data[["Table"]]<-dd[,-vect]
colnames(Data[["Table"]])
# [1] "Row.names"    "RCA.Skin.5"   "RCA.Vulva.5"  "RCA.Vulva.8"  "RCA.Skin.7"   "RCA.Vulva.3"  "RCA.Vulva.7"  "RCA.Skin.10" 
# [9] "RCA.Skin.3"   "RCA.Vulva.4"  "RCA.Skin.8"   "RCA.Vulva.9"  "RCA.Skin.12"  "RCA.Skin.1"   "RCA.Vulva.1"  "RCA.Vulva.6" 
# [17] "RCA.Skin.4"   "RCA.Skin.6"   "RCA.Vulva.10" "RCA.Vulva.11" "RCA.Vulva.12" "RCA.Skin.11"  "RCA.Skin.2"   "RCA.Vulva.2" 
# [25] "RCA.Skin.9"   "SUM"          "taxid"        "ESPECES"      "tax_id"       "order"        "family"       "genus"       
# [33] "species"      "no.rank"      "no.rank1"     "subkingdom"  
########################################################

rownames(Data[["Table"]])<- Data[["Table"]]$ESPECES


for (i in colnames(Data[["Table"]])){
  Data[['Table']][i][which(Data[['Table']][i]==""),] <- NA
}

var_Data <- colnames(
  Data[['Table']][grepl("^RCA",
                   colnames(Data[['Table']]))])
Data_vertical <- melt( Data[['Table']],
                            id.vars = "ESPECES",
                            measure.vars = var_Data)

Data_vertical$Sample <- NA

g <- grepl("Skin", Data_vertical$variable)
Data_vertical$Sample[which(g)] <- "Skin"
Data_vertical$Sample[-which(g)] <- "Vulva"

Data_vertical$order <- NA
Data_vertical$family <- NA
Data_vertical$genus <- NA
Data_vertical$species <- NA
Data_vertical$subkingdom <- NA

for (el in c("order","family","genus","species","subkingdom")){
  for (name in unique(Data_vertical$ESPECES)){
    for (item in unique(Data_vertical$variable)){
      Data_vertical
    }
    
  }
}
Data[['Table']][28:36]

Data_vertical$sample<- 

##################
Data[["skin"]] <- Data[["Table"]][grep("Skin|tax_id|order|family|genus|species|ESPECES",colnames(Data[["Table"]]))]
Data[["vulva"]] <- Data[["Table"]][grep("Vulva|tax_id|order|family|genus|species|ESPECES",colnames(Data[["Table"]]))]

  
Data[["skinreduced"]] <- Data[["skin"]][rowSums(Data[["skin"]]==0, na.rm=TRUE)<ncol(Data[["skin"]]), ] # table with no virus == 0 
Data[["vulvareduced"]] <- Data[["vulva"]][rowSums(Data[["vulva"]]==0, na.rm=TRUE)<ncol(Data[["vulva"]]), ]

Data[["skinreducedsup10"]] <- Data[["skinreduced"]][which(rowSums(Data[["skinreduced"]][,0:12])>10),] # compte > 10
Data[["vulvareducedsup10"]] <- Data[["vulvareduced"]][which(rowSums(Data[["vulvareduced"]][,0:12])>10),]

tibble::as_tibble(Data[["skinreducedsup10"]])

write.table(Data[["skinreducedsup10"]],"skin_treshold10.csv",row.names = F,col.names = T)
write.table(Data[["vulvareducedsup10"]],"vulva_treshold10.csv",row.names = F,col.names = T)
# 




# Barplots empilés avec plusieurs groupes
plot(Data[["vulvareducedsup10"]]$ESPECES)
barplot(vertical$variable~vertical$value)






mensuelles <- data.frame(produit=LETTERS[1:4],
                         ventes_2019_01=c(100,120,80,95),
                         ventes_2019_02=c(95,135,82,95),
                         ventes_2019_03=c(101,140,85,90),
                         ventes_2019_04=c(100,100,70,85),
                         ventes_2019_05=c(122,81,82,100),
                         ventes_2019_06=c(150,50,90,100)
)


# étape 1 : on met la table à la verticale avec melt

# on récupère les noms des colonnes "ventesXXX" à transposer
var_ventes <- colnames(
  mensuelles[grepl("^ventes",
                   colnames(mensuelles))])
mensuelles_vertical <- melt(mensuelles,
                            id.vars = "produit",
                            measure.vars = var_ventes)


Data[["vulvareducedsup10"]][1:5]
var <- colnames(Data[["vulvareducedsup10"]][1:12])
vertical <- melt(Data[["vulvareducedsup10"]][1:13],
                 id.vars = "ESPECES",
                 measure.vars = var)

barplot(vertical)
