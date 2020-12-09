
#recuperer echantillons
graphe <- data.frame(subset.data.frame(Table_count, select = colnames(Table_count[c(2:25)])))
graphe$taxonomy <- Table_count$taxonomy
checkboxNa <- TRUE
checkboxStacked <- FALSE
Log <- TRUE
hist_num <- 1
#recuperer Taxon
regular <- paste("^([^\\;]*\\;){",0,"}[^\\;]*$",sep="")
# "Description"=0,"family"=1,"genus"=2,"species"=3
graphe_trans <- graphe[str_detect(graphe$taxonomy, regular),] 

if (checkboxNa == FALSE){
  graphe_t <- graphe_trans[str_detect(graphe_trans$taxonomy, "NA",negate=T),]
  graphe_t$taxonomy <- str_extract(graphe_t$taxonomy,"[^\\;]*$")
} else { graphe_t <- graphe_trans
graphe_t$taxonomy <- str_extract(graphe_t$taxonomy,"[^\\;]*$") }

#Graphique barres empilees
data_graphe <- data.frame(species=graphe_t$taxonomy, melt(graphe_t, id.vars=c("taxonomy"))) %>% select(-taxonomy)

g <- data_graphe %>% filter(value>=hist_num) %>% ggplot(aes(x=variable,y=value,fill=as.factor(species))) + 
  
  {if(checkboxStacked)geom_bar(position="fill",stat="identity")}+
  {if(checkboxStacked==FALSE)geom_bar(stat="identity")}+
  
  {if(Log==TRUE) scale_y_continuous(trans = log10_trans(),
                                          breaks = trans_breaks("log10", function(x) 10^x),
                                          labels = trans_format("log10", math_format(10^.x))) }+
  
  theme_classic()+theme(axis.text.x = element_text(angle=90, hjust=1), plot.title = element_text(color="Black", size=14, 
                                                                                                 face="bold"),
                        axis.title = element_text(face="bold",hjust = 0.5) , legend.direction="horizontal", legend.position = "bottom",
                        legend.title= element_text(size=13,face="bold") ) +
  guides(fill=guide_legend(title="Type")) +
  xlab("Samples") + { if (Log==TRUE) ylab("Log Abundance") } +  { if (Log==FALSE) ylab("Abundance") }  + ggtitle(paste0("Composition of the different samples minimum threshold at :", hist_num, " sequences")) 
g



Papillo <- Papillo[Papillo$taxonomy %like% "apilloma", ]
graphe <- data.frame(subset.data.frame(Papillo, select = colnames(Papillo[c(2:25)])))
graphe$taxonomy <- Papillo$taxonomy
checkboxNa <- TRUE
checkboxStacked <- FALSE
Log <- TRUE
hist_num <- 2
#recuperer Taxon
regular <- paste("^([^\\;]*\\;){",2,"}[^\\;]*$",sep="")
# "Description"=0,"family"=1,"genus"=2,"species"=3
graphe_trans <- graphe[str_detect(graphe$taxonomy, regular),] 

if (checkboxNa == FALSE){
  graphe_t <- graphe_trans[str_detect(graphe_trans$taxonomy, "NA",negate=T),]
  graphe_t$taxonomy <- str_extract(graphe_t$taxonomy,"[^\\;]*$")
} else { graphe_t <- graphe_trans
graphe_t$taxonomy <- str_extract(graphe_t$taxonomy,"[^\\;]*$") }

#Graphique barres empilees
data_graphe <- data.frame(species=graphe_t$taxonomy, melt(graphe_t, id.vars=c("taxonomy"))) %>% select(-taxonomy)

g <- data_graphe %>% filter(value>=hist_num) %>% ggplot(aes(x=variable,y=value,fill=as.factor(species))) + 
  
  {if(checkboxStacked)geom_bar(position="fill",stat="identity")}+
  {if(checkboxStacked==FALSE)geom_bar(stat="identity")}+
  
  {if(Log==TRUE) scale_y_continuous(trans = log10_trans(),
                                    breaks = trans_breaks("log10", function(x) 10^x),
                                    labels = trans_format("log10", math_format(10^.x))) }+
  
  theme_classic()+theme(axis.text.x = element_text(angle=90, hjust=1), plot.title = element_text(color="Black", size=14, 
                                                                                                 face="bold"),
                        axis.title = element_text(face="bold",hjust = 0.5) , legend.direction="horizontal", legend.position = "bottom",
                        legend.title= element_text(size=13,face="bold") ) +
  guides(fill=guide_legend(title="Genus")) +
  xlab("Samples") + { if (Log==TRUE) ylab("Log Abundance") } +  { if (Log==FALSE) ylab("Abundance") }  + ggtitle(paste0("Composition of the different samples minimum threshold at :", hist_num, " sequences \n without \"Human Papillomavirus\"")) 
g
