## Graphics 
# NOVEMBER 2020
rm(list = ls(all = TRUE))
gc()
## plot 

#baseDir <- "~/Documents/ROSARIO/fichiers/RCA_virus"
#SourceDir <- "~/Documents/ROSARIO/fichiers/RCA_virus/src"
#OutputDir <- "~/Documents/ROSARIO/fichiers/RCA_virus/plot"
#DataDir <- "~/Documents/ROSARIO/fichiers/RCA_virus/data"

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
library(scales)
library(tidyverse)
library(rlist)
library(plotly)

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
# je mets des NA partout ici 

for (i in colnames(ncbi_virus_lineages)){
  ncbi_virus_lineages[i][which(ncbi_virus_lineages[i]==""),] <- NA
}


# j'enleve les colonnes qui ne contiennes que des NA 
### Colonne contenant que des NA
lineage <- ncbi_virus_lineages[, apply(ncbi_virus_lineages, 2, function(x) !all(is.na(x)))]
# write.table(lineage, file=paste(baseDir,"lineage","ncbi_lineages-virus_short.csv", sep = "/"), col.names = TRUE , row.names=FALSE, sep=",")
x_<-which(colnames(lineage)%in%c("subkingdom","phylum","no.rank","no.rank1","class", "kingdom", "no.rank2", "no.rank22", "no.rank3", "subclass", "subfamily", "suborder" , "subphylum", "subspecies", "subtribe" ,"superfamily" ))
lineage <- lineage[-x_]
x <- grep(paste(rownames(Data[["Abundance_reduced_taxid"]]),  collapse="|"),lineage$tax_id) # lineage position en récupérant par TAXID

######### TABLE DE TAXONOMIE 
table_analyse<-lineage[x,][which(lineage[x,]$tax_id%in%rownames(Data[["Abundance_reduced_taxid"]])),] # lineage en ne gardant que les ligne d'intéret 

 
table_analyse$taxonomy<-NA
for (el in 1:nrow(table_analyse)){
  name=""

  for (i in c("superkingdom","order","family","genus","species")){
    if (name==""){
      name=table_analyse[el,i]
    }
    else {
       name=paste(name,table_analyse[el,i],sep=";")
    }
  }
  table_analyse$taxonomy[el]<-name
}


############ table contenant les comptes 
df <- Data[["Abundance_reduced_taxid"]] # version avec les TAXiD

df$tax_id <- rownames(df)
table_analyse<- transform(table_analyse, tax_id = as.character(tax_id))
table_analyse <-table_analyse[order(table_analyse$tax_id),]
df <-df[order(df$tax_id),]

#for(i in 1:nrow(df)){ print(df$tax_id[i]==table_analyse$tax_id[i]) }
dd<- merge(df, table_analyse)


dd$tax_id <- as.integer(dd$tax_id)


tmpdd <-dd[order(-dd$RCA.Skin.5,-dd$RCA.Skin.7,-dd$RCA.Skin.10,-dd$RCA.Skin.3,-dd$RCA.Skin.8,-dd$RCA.Skin.12,-dd$RCA.Skin.1,-dd$RCA.Skin.4,-dd$RCA.Skin.6,-dd$RCA.Skin.11,-dd$RCA.Skin.2,-dd$RCA.Skin.9),]
tmpNamedtable <-Data[["Abundance_atleast2seq"]][order(-Data[["Abundance_atleast2seq"]]$RCA.Skin.5 ,-Data[["Abundance_atleast2seq"]]$RCA.Skin.7 ,-Data[["Abundance_atleast2seq"]]$RCA.Skin.10 ,-Data[["Abundance_atleast2seq"]]$RCA.Skin.3 ,-Data[["Abundance_atleast2seq"]]$RCA.Skin.8 ,-Data[["Abundance_atleast2seq"]]$RCA.Skin.12 ,-Data[["Abundance_atleast2seq"]]$RCA.Skin.1 ,-Data[["Abundance_atleast2seq"]]$RCA.Skin.4,-Data[["Abundance_atleast2seq"]]$RCA.Skin.6 ,-Data[["Abundance_atleast2seq"]]$RCA.Skin.11 ,-Data[["Abundance_atleast2seq"]]$RCA.Skin.2 ,-Data[["Abundance_atleast2seq"]]$RCA.Skin.9 ),]
tmpNamedtable$Description <- rownames(tmpNamedtable)

dd1 <- merge(tmpNamedtable,tmpdd)

# write.csv(dd1,paste(DataDir,"table_compte_more_than2.csv",sep="/"))


var_Data <- colnames( dd1[grepl("^RCA",colnames(dd1))])

restructuredData = melt( dd1,
  id = c("Description","order","family","genus","species"),
    measure.vars = var_Data, 
    variable_name = "number" )



levels(restructuredData$variable) <- c( "RCA.Skin.1","RCA.Skin.2", "RCA.Skin.3","RCA.Skin.4","RCA.Skin.5",
                                        "RCA.Skin.6" ,"RCA.Skin.7", "RCA.Skin.8","RCA.Skin.9","RCA.Skin.10", 
                                        "RCA.Skin.11","RCA.Skin.12" ,"RCA.Vulva.1","RCA.Vulva.2","RCA.Vulva.3",
                                        "RCA.Vulva.4","RCA.Vulva.5" , "RCA.Vulva.6" ,"RCA.Vulva.7","RCA.Vulva.8",
                                        "RCA.Vulva.9","RCA.Vulva.10", "RCA.Vulva.11" ,"RCA.Vulva.12" )

restructuredData$from_ <- "Skin"
restructuredData$from_[which(grepl("Vulva", restructuredData$variable))] <- "Vulva"



## table without all the 0 values 
reducedrestructuredData <- restructuredData[!restructuredData$value==0,]



#######
# reproduction rmd
df <- dd1[,c(0:24,length(dd1))]

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

Table_count <- data.frame(matrix(0L,ncol=25, nrow=length(taxonomy_split))) # creation de la table de compte 
colnames(Table_count) <- c("taxonomy",levels(as.factor(restructuredData$variable)))
Table_count$taxonomy <- taxonomy_split

# je range df
df <- df[c( "RCA.Skin.1","RCA.Skin.2", "RCA.Skin.3","RCA.Skin.4","RCA.Skin.5",
                                        "RCA.Skin.6" ,"RCA.Skin.7", "RCA.Skin.8","RCA.Skin.9","RCA.Skin.10", 
                                        "RCA.Skin.11","RCA.Skin.12" ,"RCA.Vulva.1","RCA.Vulva.2","RCA.Vulva.3",
                                        "RCA.Vulva.4","RCA.Vulva.5" , "RCA.Vulva.6" ,"RCA.Vulva.7","RCA.Vulva.8",
                                        "RCA.Vulva.9","RCA.Vulva.10", "RCA.Vulva.11" ,"RCA.Vulva.12","taxonomy" )]



for (line in Table_count$taxonomy){
  str_to_search <- str_replace_all(noquote(line), "[(]","[(]")
  str_to_search <- str_replace_all(str_to_search, "[)]","[)]")
  df$get <- str_detect(df$taxonomy, str_to_search)
  line[2:25] <- c(colSums(df[df$get == "TRUE",1:24]))
  Table_count[Table_count$taxonomy == line[1],][-1] <- as.integer(line[-1])
  #grace au vecteur found on sait quelle ligne sont a recuperer dans df et a additionner dans Table_count
}

##############################################
# Graphique barres empilees, repartition globale
##############################################

shinyApp(

  ui = fluidPage(
  #Widgets 
    sidebarPanel(
      h2(strong("Options :")),
      selectInput("taxon_radio", width = "100%",
                  h4("Niveau taxonomique:"),
                  choices=list("Description"=0,"order"=1,"family"=2,"genus"=3,"species"=4),selected = 0),
        checkboxGroupInput("hist_checkbox",
                         h4("Echantillons desires:"), width = "100%",
                         names(Table_count[,2:25]), selected = colnames(Table_count[,2:25]), inline=F),
      uiOutput("slider_num"),
      width=3),
  #Graphe barres empilees
    mainPanel(plotlyOutput("histogramme_c"),width=9)
  ),

  server = function(input, output) {
  #Widget change seuil max = nombre individus du re-echantillonnage
    output$slider_num <- renderUI({sliderInput("hist_num", width = "100%",
                  h4("Nombre d'individus minimum (seuil):"),
                  min=1, max=input$sample_num, step=1, value=20)})
  #Graphique creation
    output$histogramme_c <- renderPlotly({
      #recuperer echantillons
      graphe <- data.frame(subset.data.frame(Table_count, select = input$hist_checkbox))
      graphe$taxonomy <- Table_count$taxonomy
      
      #recuperer Taxon
      regular <- paste("^([^\\;]*\\;){",input$taxon_radio,"}[^\\;]*$",sep="")
      graphe_trans <- graphe[str_detect(graphe$taxonomy, regular),] 
      graphe_t <- graphe_trans[str_detect(graphe_trans$taxonomy, "(unclassified)$",negate=T),]
      graphe_t$taxonomy <- str_extract(graphe_t$taxonomy,"[^\\;]*$")
      
  
      #Graphique barres empilees
      data_graphe <- data.frame(species=graphe_t$taxonomy, melt(graphe_t, id.vars=c("taxonomy"))) %>% select(-taxonomy)
  g <- data_graphe %>% filter(value>=input$hist_num) %>% ggplot(aes(x=variable,y=value,fill=as.factor(species))) + 
       geom_bar(position="fill",stat="identity") +
       theme_classic()+theme(axis.text.x = element_text(angle=90, hjust=1), plot.title = element_text(color="Black", size=14, 
                                                                                                      face="bold"),
             axis.title = element_text(face="bold",hjust = 0.5) , legend.direction="horizontal", legend.position = "bottom",
             legend.title= element_text(size=13,face="bold") ) +
       guides(fill=guide_legend(title="Taxons")) +
       xlab("Echantillons") + ylab("Abondance") + ggtitle("Composition des differents echantillons") 
      ggplotly(g)%>% layout(height = 850)})
      
  },
  options = list(height = 950, width = "100%")
)
