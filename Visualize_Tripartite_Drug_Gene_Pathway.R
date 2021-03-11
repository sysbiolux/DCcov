#library(ggrepel)
#library(igraph)
library(networkD3)
library(dplyr)
library(ggbipart)
library(magrittr)
library(htmlwidgets)
library(htmltools)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
data <- read.csv('SKO_ALL_Studies.csv',row.names = 1)
#colnames(data) <- c('Genes','Genes_Symbols','Phenotype Essentiality', 'Phenotype Safety')
png(filename="Figs/SKO_Scatter_Plot.png", units="in", width=4, height=6, res=300)
ggplot(data, aes(Phenotype_Safety, Phenotype_Essentiality, label = rownames(data))) +
  geom_text_repel(max.overlaps = 23, min.segment.length = 1,force_pull = 0.5,force = 2 ) +
  geom_point(color = 'red') +
  #geom_label_repel()+
  theme_classic(base_size = 12)+ xlab("Safety Score") +
  ylab("Essentiality Score") 
dev.off()

ggsave("Figs/SKO_Scatter_Plot.png",dpi = 300)

#plot(Phenotype_Essentiality ~Phenotype_Safety, col="lightblue", pch=19, cex=2,data=data)
#text(Phenotype_Essentiality ~Phenotype_Safety, labels=rownames(data),data=data, cex=0.9, font=1)
#devtools::install_github("pedroj/bipartite_plots")

# Read tripartite table for DKO drugs-genes-pathways
dko_drugs <- read.csv('KO_data/DKO_Drugs_OneDrug_Pathways_tripartite.csv',row.names = 1)

# A connection data frame is a list of flows with intensity for each flow
m1<- as.matrix(dko_drugs[,c(1,2)])
m2 <-as.matrix(dko_drugs[,c(2,6)])
m <- rbind(m1,m2)

links <- data.frame(
  source=m[,1], 
  target=m[,2], 
  value=rep(0.2,length(m[,2])),
  edge_color=c(m[,2])
)

# From these flows we need to create a node data frame: it lists every entities involved in the flow
names_d =c(as.character(dko_drugs$Drug)) %>% unique()
names_g =c(as.character(dko_drugs$Gene_Pair)) %>% unique()
names_p =c(as.character(dko_drugs$Pathways)) %>% unique()
group <-c(rep('Drug',length(names_d)),rep('Gene',length(names_g)),rep('Pathway',length(names_p)))
names <- c(names_d,names_g, names_p)
nodes <- data.frame(
  name=names,
  color=group,
  stringsAsFactors= FALSE
)


# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1
#my_color <- 'd3.scaleOrdinal() .domain(["edge_color"]) .range(["edge_color"])'


# Make the Network
p <- sankeyNetwork(Links = links, Nodes = nodes,fontSize = 12,fontFamily = 'color',#nodeWidth = #colourScale =my_color_scale,# JS("d3.scaleOrdinal(links$color);"),
                   LinkGroup="edge_color", # colourScale=my_color,
                   NodeGroup = "color",
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", 
                   sinksRight=TRUE)
p #%>% 
  #htmlwidgets::prependContent(htmltools::tags$h1("Title")) %>%
  #saveNetwork(file = 'm.html')
  
htmlwidgets::saveWidget(p, file="m.html")
webshot::webshot("m.html", "Figs/DKO_OneDrug_tripartite_graph_3colors.png",zoom=2)


# Read tripartite table for DKO drugs-genes-pathways
dko_drugs <- read.csv('KO_data/DKO_Drug_Repurposing_Reduced_Pathways_tripartite.csv',row.names = 1)

# A connection data frame is a list of flows with intensity for each flow
m1<- as.matrix(dko_drugs[,c(1,3)])
m2 <-as.matrix(dko_drugs[,c(2,3)])
m3 <-as.matrix(dko_drugs[,c(3,7)])
#m4 <-as.matrix(dko_drugs[,c(1,2)])
m <- rbind(m1,m2,m3)
len= length(m[,2])/3
color <- cbind(t(dko_drugs$Drugs1),t(dko_drugs$Drugs2),t(rep('n',length(m[,2])/3)))
#color <- cbind(t(rep('1',len)),t(rep('2',len)),t(rep('3',len)))

x <-colors()
y <-as.factor(color)


links <- data.frame(
  source=m[,1], 
  target=m[,2], 
  value=rep(0.7,length(m[,2])),
  #color=as.numeric(color),#m[,1],
  edge_color=c(m[,2]),
  stringsAsFactors= FALSE
)

# From these flows we need to create a node data frame: it lists every entities involved in the flow
names_d =c(as.character(dko_drugs$Drugs1),  as.character(dko_drugs$Drugs2)) %>% unique()
names_g =c(as.character(dko_drugs$Gene_Pair)) %>% unique()
names_p =c(as.character(dko_drugs$Pathways)) %>% unique()
group <-c(rep('Drug',length(names_d)),rep('Gene',length(names_g)),rep('Pathway',length(names_p)))
names <- c(names_d,names_g, names_p)
nodes <- data.frame(
  name=names,
  color=group,
  stringsAsFactors= FALSE
  #name=c(as.character(links$source), 
  #       as.character(links$target)) %>% unique()
  #, color=as.character(color) %>% unique()#m[,1]
)

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1
#links$IDcolor <- match(links$color, nodes$name)-1

# Make the Network
p <- sankeyNetwork(Links = links,fontFamily = 'color',
                    #nodeWidth =2,#nodeWidth = #colourScale =my_color_scale,# JS("d3.scaleOrdinal(links$color);"),
                   #LinkGroup = 'IDcolor',
                   NodeGroup = "color",
                   LinkGroup="edge_color",
                   Nodes = nodes,fontSize = 11,#colourScale = 'color',
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", 
                   sinksRight=FALSE)
p
htmlwidgets::saveWidget(p, file="m.html")
webshot::webshot("m.html", "Figs/DKO_Reduced_list_tripartite_graph_3colors.png",zoom=2)

my_color_scale <- 'd3.scaleOrdinal() .domain(uniq_drugs) .range(my_color)'
uniq_drugs <- unique(color[1,])
my_color <- list()
for(i in 1:length(uniq_drugs)){
  drug_color = colors()[1+i]
  my_color <- append(my_color, drug_color)
}
my_color=as.character(my_color)

# Read tripartite table for SKO drugs-genes-pathways
sko_drugs <- read.csv('KO_data/SKO_db_Drugs_Unique_Pathways_tripartite.csv',row.names = 1)

# A connection data frame is a list of flows with intensity for each flow
m1<- as.matrix(sko_drugs[,c(1,2)])
m2 <-as.matrix(sko_drugs[,c(2,5)])
m <- rbind(m1,m2)
value <-cbind(t(rep(1,length(m[,2])/2)),t(sko_drugs$Phenotype_Essentiality))

links <- data.frame(
  source=m[,1], 
  target=m[,2], 
  value=rep(0.1,length(m[,1])),#t(value)
  edge_color=c(m[,2])
)

names_d =c(as.character(sko_drugs$Drugs)) %>% unique()
names_g =c(as.character(sko_drugs$Gene_Symbol)) %>% unique()
names_p =c(as.character(sko_drugs$Pathways)) %>% unique()
group <-c(rep('Drug',length(names_d)),rep('Gene',length(names_g)),rep('Pathway',length(names_p)))
names <- c(names_d,names_g, names_p)

nodes <- data.frame(
  name=names,
  color=group,
  stringsAsFactors= FALSE
  #name=c(as.character(links$source), 
  #       as.character(links$target)) %>% unique()
  #, color=as.character(color) %>% unique()#m[,1]
)


# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1

library(networkD3)

colors <- paste(links$nodes$colors, collapse = '", "')
colorJS <- paste('d3.scaleOrdinal(["', colors, '"])')
# Make the Network
p <- sankeyNetwork(#Links = links$links, Nodes = links$nodes, 
                     Links = links, Nodes = nodes,#nodePadding = 2,
                    fontSize = 12,fontFamily = 'color',#nodeWidth = #colourScale =my_color_scale,# JS("d3.scaleOrdinal(links$color);"),
                   #LinkGroup = 'IDcolor',
                   NodeGroup = "color",
                   Source = "IDsource", Target = "IDtarget",
                    #Source = 'source', Target = 'target',
                   LinkGroup="edge_color",
                   Value = "value", NodeID = "name", 
                   sinksRight=TRUE)
p

htmlwidgets::saveWidget(p, file="m.html")
webshot::webshot("m.html", "Figs/SKO_Drugs_tripartite_graph_3colors.png",zoom=2.4)

#convert +append Figs/SKO_Scatter_Plot.png Figs/SKO_Drugs_tripartite_graph_3colors.png Figs/Fig3.png
#convert +append Figs/DKO_Reduced_list_tripartite_graph_3colors.png Figs/DKO_OneDrug_tripartite_graph_3colors.png  Figs/DKO_image_combined.png
