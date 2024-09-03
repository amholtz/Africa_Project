# This code was written by Samuel Hong.
# Samuel Hong will be part of the co-authors of this article.

library(ape)
library(glue)
library(dplyr)
library(readr)

# Helper functions

# Gets loc from name of strain.
get_loc <- function(string,delimiter,index,reverse=FALSE){
    if(reverse==TRUE){
        splitted = strsplit(string,split=delimiter,fixed=T)
        reved = lapply(splitted,rev)
        return(lapply(reved,"[[",index))
    }
    else{
    return(lapply(strsplit(string,split=delimiter,fixed=T),"[[",index))
    }
}


# This function returns the location of sequence n.
get_loc_from_meta <- function(n,metadf){
    loc<- (metadf %>% filter(Accession==n))$loc
    return(as.character(loc))
}

# This function gets the location trait for each sequence in name_array.
# The metadata file must have a 'loc' and a 'name' column.
get_loc_array_meta <- function(name_array,meta){
    locs<-sapply(name_array,get_loc_from_meta,metadf=meta)
    return(as.character(locs))
}

# Returns the two direct descendant nodes.
getDescendants<-function(tree,node){
    if(node<=length(tree$tip.label)){
        return(-1)
    } else{
        return( tree$edge[which(tree$edge[,1]==node),2] )
    }
}

getAllDescendants<-function(tree,node,curr=NULL){
    if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
    if(is.null(curr)) curr<-vector()
    daughters<-tree$edge[which(tree$edge[,1]==node),2]
    curr<-c(curr,daughters)
    if(length(curr)==0&&node<=Ntip(tree)) curr<-node
    w<-which(daughters>Ntip(tree))
    if(length(w)>0) for(i in 1:length(w)) 
        curr<-getAllDescendants(tree,daughters[w[i]],curr)
    return(curr)
}

getTipsFromNode<-function(tree,node){
    desc <- getAllDescendants(tree,node,curr=NULL)
    return(desc[desc<length(tree$tip.label)])
}


getNodeAnnots <- function(tree,tip_annot){
    tree<-reorder(tree,"postorder") # For each internal node, we will first iterate descendants.
    # prev: "pruningwise". According to reorder documentation, "postorder" is more efficient.
    num_nodes <- (length(tree$tip.label)+tree$Nnode)
    node_annot_list <- rep(0, length = num_nodes)
    rootnode <- length(tree$tip.label)+1
    for(node in c(as.array(tree$edge[,2]),rootnode)){ #Iterate over all nodes
        
        if(node<=length(tree$tip.label)){ # If node is tip, annot is loc.
            node_annot_list[node] <- tip_annot[node]
        }
        
        else{
            children_nodes <- getDescendants(tree,node)
            locs_in_children <- c()
            for (c in children_nodes){
                locs_in_children <- c(locs_in_children, node_annot_list[c])
            }

            locs_in_children <- unique(unlist(locs_in_children))
            
            # If all descendant tips have the same location, label node with this location, MIXED otherwise.
            if( ("MIXED" %in% locs_in_children) | (length(locs_in_children)>1) ){
                node_annot_list[node] = "MIXED"
            }
            else{
                node_annot_list[node] = node_annot_list[c]
            }
        }
    
    }
    return(node_annot_list)
}


# Returns the biggest non-mixed clades.
getClusters <- function(tree,node_annot){
    non_mixed <- which(node_annot!="MIXED" & node_annot!=loc_to_keep)
    children_nodes <- c() # This vector will contain each non-mixed node, which is child of another non-mixed node.
    for( internal_node in non_mixed[non_mixed>length(tree$tip.label)] ){
        children_nodes <- c(children_nodes, getDescendants(tree,internal_node))
    }
    return(setdiff(non_mixed[non_mixed>length(tree$tip.label)], children_nodes)) 
    # Only return non-mixed node which are not children --> Returns the big non-mixed clade and not the corresponding descending subclades.
}

# This function creates a df: [taxa, cluster, remove]
# taxa: tip's name
# cluster: Parent node of the biggest homogeneous cluster the tip belongs to.
# remove [bool]: Do we remove this tip? Based on a random selection.
pruneClusters<-function(tree,cluster_nodes){
    toremove_df <- as.data.frame(tree$tip.label)
    colnames(toremove_df) <- c("taxa")
    toremove_df$cluster <- 0
    toremove_df$remove <- 0
    for (c in cluster_nodes){
        subtree_tips <- tree$tip.label[getTipsFromNode(tree, c)] # Here again, we already go through all tips when going backward --> There's probably a way to save them in advance.
        filtered <- sample(subtree_tips, length(subtree_tips)-1) # Randomly select all tips but one.
        toremove_df$cluster[which(toremove_df$taxa %in% subtree_tips )] = c
        toremove_df$remove[which(toremove_df$taxa %in% filtered )] = 1
    }
    return(toremove_df)
}


# I didn't adapt this function.
downsampleTree <- function(tree,split_char,split_index,reverse){
    num_taxa<-length(tree$tip.label)
    # get tip annotations by delimiter and index, use reverse to count from the back
    tip_annot <- get_loc(tre$tip.label, split_char, split_index,reverse)
    node_annots<-getNodeAnnots(tree,tip_annot)
    clusters<-getClusters(tree,node_annots)
    toremove_df <- pruneClusters(tree,clusters)
    num_removed <- length(toremove_df[which(toremove_df$remove==1),'taxa'])
    print(glue("Tree pruned from {num_taxa} taxa to {num_taxa-num_removed} after collapsing {length(clusters)} clusters"))
    return(list(toremove_df,node_annots))
}

# This function outputs:
#   - df with tips that should be removed, based on the heuristic.
#   - A list of labels for each node: whether MIXED or with the trait if the cluster is homogeneous.
# from metadata file with location annotated as loc
downsampleTree_meta <- function(tree,meta){
    num_taxa<-length(tree$tip.label)
    # get tip annotations by delimiter and index, use reverse to count from the back
    tip_annot <- get_loc_array_meta(tree$tip.label,meta)
    node_annots<-getNodeAnnots(tree,tip_annot)
    clusters<-getClusters(tree,node_annots)
    toremove_df <- pruneClusters(tree,clusters)
    num_removed <- length(toremove_df[which(toremove_df$remove==1),'taxa'])
    print(glue("Tree pruned from {num_taxa} taxa to {num_taxa-num_removed} after collapsing {length(clusters)} clusters"))
    return(list(toremove_df,node_annots))
}

# I didn't adapt this function.
pruneClusters_N<-function(tree,cluster_nodes,N){
    toremove_df <- as.data.frame(tree$tip.label)
    colnames(toremove_df) <- c("taxa")
    toremove_df$cluster <- 0
    toremove_df$remove <- 0
    for (c in cluster_nodes){
        subtree_tips <- tree$tip.label[getTipsFromNode(tree, c)]
        toremove_df$cluster[which(toremove_df$taxa %in% subtree_tips )] = c
        subtree_size<-length(subtree_tips)
        if(subtree_size>N){
            filtered <- sample(subtree_tips, subtree_size-N)
            toremove_df$remove[which(toremove_df$taxa %in% filtered )] = 1
        }
        else{
            filtered <- NA
        }
    }
    return(toremove_df)
}

# I didn't adapt this function.
# keep only N per cluster
downsampleTree_N <- function(tree,N,split_char,split_index,reverse){
    num_taxa <- length(tree$tip.label)
    # get tip annotations by delimiter and index, use reverse to count from the back
    tip_annot <- get_loc(tree$tip.label, split_char, split_index,reverse)
    node_annots <- getNodeAnnots(tree,tip_annot)
    clusters <- getClusters(tree,node_annots)
    toremove_df <- pruneClusters_N(tree,clusters,N)
    num_removed <- length(toremove_df[which(toremove_df$remove==1),'taxa'])
    print(glue("Tree pruned from {num_taxa} taxa to {num_taxa-num_removed} after collapsing {length(clusters)} clusters.\n{num_removed} taxa removed."))
    return(list(toremove_df,node_annots))
}

#analysis
# 100seq example
loc_to_keep <- "ct" # No sequence with this corresponding 'loc' value within metadata file will be excluded from the tree.
tre <- read.tree("/Users/aholtz/Dropbox/rabies/Africa_Project/data/Africa_RABV_July2024.nwk_collapsed_0.5.nwk")
metadata <- read.delim2("/Users/aholtz/Dropbox/rabies/Africa_Project/data/meta_edited.tab") 
metadata <- metadata %>% filter(Accession %in% tre$tip.label)
colnames(metadata)[5] <- 'loc'
sampleClusters <- downsampleTree_meta(tre, meta = metadata)
annotated_nodes = sampleClusters[[2]]
clustered_df = sampleClusters[[1]]
# note, cluster = 0 means the taxa does not fall in a monophyletic clade w sequences from the same location

#write.table(clustered_df, "cluster_sampling.csv", sep = ',', row.names=FALSE)

keep<-(clustered_df %>% filter(remove==0))$taxa
drop<-(clustered_df %>% filter(remove==1))$taxa
ktre<-drop.tip(tre,drop)
write.tree(ktre,'subsampled_tree.tre')

tips_orig <- tre$tip.label
meta_country <- metadata %>% left_join(clustered_df, by = c('Accession'='taxa'))
meta_country <- meta_country %>% filter(Accession %in% tips_orig)
write.table(meta_country, "meta_sampling_clustering.csv", sep = ',', row.names=FALSE)

meta_removed <- metadata %>% filter(Accession %in% keep)
country_diversity_removed <- meta_removed %>% group_by(loc) %>% tally()
country_diversity <- metadata %>% group_by(loc) %>% tally()

colnames(country_diversity)[2] <- 'n_orig'
colnames(country_diversity_removed)[2] <- 'n_post'

country_div_comp <- left_join(country_diversity, country_diversity_removed, by = 'loc')
write.table(country_div_comp, "removed_seq_country.csv", sep = ',', row.names=FALSE)

sub_meta <- meta_country %>% filter(remove == 0) %>% select(-remove)
write_delim(sub_meta, 'sub_meta.tab', delim = '\t', quote = 'none')
