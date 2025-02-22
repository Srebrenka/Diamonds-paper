
################
# Diamonds 1   # 
################
# script 1

# Read in the data - all networks

#########################
#                       #
#     Load libraries    #
#                       #
#########################

library(dplyr)
library(openxlsx)
library(readxl)
library(igraph)
library(purrr)



# Reading in the data - undirected networks


# 1. Monastery 
net1<-read_graph("Diamonds/Monastery interactions.graphml", format = "graphml")


gdf <-get.data.frame(net1)
E(net1)$sign <- gdf$weight
V(net1)$name <- as.character(seq(1, gorder(net1)))


net1_neg <- net1 - edge(which(E(net1)$sign == 1))
net1_pos <- net1 - edge(which(E(net1)$sign == -1))


# 2. Eastern College
net2<-read_graph("Diamonds/Eastern College.graphml", format = "graphml")

gdf <-get.data.frame(net2)
E(net2)$sign <- gdf$weight
V(net2)$name <- as.character(seq(1, gorder(net2)))

net2_neg <- net2 - edge(which(E(net2)$sign == 1))
net2_pos <- net2 - edge(which(E(net2)$sign == -1))

# 3. Fraternity
net3<-read_graph("Diamonds/Fraternity.graphml", format = "graphml")

gdf <-get.data.frame(net3)
E(net3)$sign <- gdf$weight

V(net3)$name <- as.character(seq(1, gorder(net3)))

net3_neg <- net3 - edge(which(E(net3)$sign == 1))
net3_pos <- net3 - edge(which(E(net3)$sign == -1))

# 4. Highland tribes

net4<-read_graph("Diamonds/Highland tribes.graphml", format = "graphml")


gdf <-get.data.frame(net4)
E(net4)$sign <- gdf$weight

V(net4)$name <- as.character(seq(1, gorder(net4)))

net4_neg <- net4 - edge(which(E(net4)$sign == 1))
net4_pos <- net4 - edge(which(E(net4)$sign == -1))

# 5. International networks
net5 <- read_graph("Diamonds/N46.net", format = "pajek")

gdf <-get.data.frame(net5)
E(net5)$sign <- gdf$weight
V(net5)$name <-as.character(seq(1, gorder(net5)))

net5_neg <- net5 - edge(which(E(net5)$sign == 1))
net5_pos <- net5 - edge(which(E(net5)$sign == -1))


# put it in the list

nets<- list(net1, net2, net3, net4, net5)
nets_pos<- list(net1_pos, net2_pos, net3_pos, net4_pos, net5_pos)
nets_neg<- list(net1_neg, net2_neg, net3_neg, net4_neg, net5_neg)

# Organizations
load("Diamonds/recens_networks_modified.Rdata")
#ls()
rec_net <- recens_ud_strong


# add names to nodes
# inspect N of nets
Ns_rec <- map(rec_net, gorder) # they are in the different order
Ns_rec_den <- map(rec_net, edge_density)
rec_net <- rec_net[- 8]  

rec_pos <- list()
rec_neg <- list()
for(i in 1:length(rec_net)){
  rec_neg[[i]] <- rec_net[[i]] - edge(which(E(rec_net[[i]])$sign == 1))
  rec_pos[[i]] <- rec_net[[i]] - edge(which(E(rec_net[[i]])$sign == -1))
}


# villages
load("Diamonds/villages.Rdata")
villages <- village.net
# add names to nodes
for(i in 1:length(villages)){
  V(villages[[i]])$name <- as.character(seq(1, gorder(villages[[i]])))
}

vil_pos <- list()
vil_neg <- list()
for(i in 1:length(villages)){
  vil_neg[[i]] <- villages[[i]] - edge(which(E(villages[[i]])$sign == 1))
  vil_pos[[i]] <- villages[[i]] - edge(which(E(villages[[i]])$sign == -1))
}
# inspect N of nets
Ns_vil <- map(villages, gorder)
Ns_vil <- unlist(Ns_vil)
Ns_vil 

