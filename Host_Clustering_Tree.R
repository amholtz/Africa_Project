####################################################
##                           Host Species Cluster in Tree 
## ####   Looking a Hosts per node to determine host by host clustering
## ####   Based off of Previous R Script from Andrew Holtz 
## author: Andrew Holtz
## creation date: 2024/July/10
###################################################

rm(list=ls())
library(dplyr)
library(cepiigeodist)
library(countrycode)
library(seqinr)
library(readr)
library(dplyr)
library(data.table)
library(optparse)
library(ape)
library(reshape)
library(data.table)
library(DT)





## Read Nexus Tree

tree_full <- read.nexus('/Users/aholtz/Dropbox/rabies/Africa_Project/data/Africa_CI_OutRem.date_pastml/named.tree_Africa_CI_OutRem.date.nexus')
africa_metadata_hostDog <- read.delim("~/Dropbox/rabies/Africa_Project/data/metadata_final_wisolate_DOG.tab")

#Create simple host_accession table
host_table <- africa_metadata_hostDog %>% select(Accession, dog)


####
select.tip.or.node <- function(element, tree) {
  ifelse(element < Ntip(tree)+1,
         tree$tip.label[element],
         tree$node.label[element-Ntip(tree)])
}


edge_table_full <- data.frame(
  "child" = tree_full$edge[,2],
  "node" = sapply(tree_full$edge[,2],
                  select.tip.or.node,
                  tree = tree_full))

edge_table_full$tips<- sapply(edge_table_full$child, function(x){
  phangorn::Descendants(tree_full, x, type = 'tips')
})

node_full<-edge_table_full
for(i in 1:nrow(edge_table_full)) {
  node_full$key_lists[i] <- list(edge_table_full$node[unlist(edge_table_full$tips[i])])
}

names(node_full)[1] <- 'node_key'

#Join with metadata for host
host_table <- africa_metadata_hostDog %>% select(Accession, dog)
node_full <- left_join(node_full, host_table, by = c('node'='Accession'))

## Now group
setDT(node_full)[, l:=sapply(key_lists,length)]

# Function to get top N hosts and their counts
get_top_hosts_with_counts <- function(node_ids, meta_host, N = 5) {
  # Filter meta_host for relevant node_ids
  hosts <- meta_host %>%
    filter(Accession %in% node_ids) %>%
    group_by(dog) %>%
    summarise(count = n()) %>%
    arrange(desc(count))
  
  # Get the top N hosts and their counts
  top_hosts <- hosts %>%
    head(N)
  
  # Extract hosts and counts
  top_hosts_list <- top_hosts$dog
  top_counts_list <- top_hosts$count
  
  # Pad with NA if fewer than N hosts
  length(top_hosts_list) <- N
  length(top_counts_list) <- N
  
  return(list(top_hosts_list, top_counts_list))
}

# Apply the function to each row in node_full
top_hosts_with_counts_list <- lapply(node_full$key_lists, get_top_hosts_with_counts, meta_host = host_table, N = 5)

# Separate the hosts and counts into different lists
top_hosts_list <- lapply(top_hosts_with_counts_list, `[[`, 1)
top_counts_list <- lapply(top_hosts_with_counts_list, `[[`, 2)

# Convert the lists to data frames
top_hosts_df <- do.call(rbind, lapply(top_hosts_list, function(x) as.data.frame(t(x))))
top_counts_df <- do.call(rbind, lapply(top_counts_list, function(x) as.data.frame(t(x))))

# Name the columns appropriately
colnames(top_hosts_df) <- paste0('host_major', 1:5)
colnames(top_counts_df) <- paste0('N_', 1:5)

# Combine the host and count data frames with the original node_full dataframe
host_cluster_result <- cbind(node_full, top_hosts_df, top_counts_df) %>% 
  select(node,l, host_major1, N_1, host_major2, N_2, host_major3,
         N_3,host_major4, N_4,host_major5, N_5) 
colnames(host_cluster_result)[2] <- 'Ntips'
