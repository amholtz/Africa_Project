####################################################
##                           CLEAN ALL RABV METADATA
## ####   Removing Sequences with no country, date data, and small fragments
## ####   Defining Simple Clade information
## ####   Defining host definitions
## ####   Defining the genetic region of each sequence
## ####   Based off of Previous R Script from Andrew Holtz 
## author: Andrew Holtz
## creation date: 2024/January/30
###################################################

rm(list=ls())
library(dplyr)
library(cepiigeodist)
library(countrycode)
library(ggplot2)
library(seqinr)
library(readr)
library(dplyr)
library(stringr)
library(data.table)
library(taxize)
library(optparse)
library(ape)
library(reshape)
library(data.table)
library(DT)


meta <- read.csv("~/Dropbox/rabies/Africa_Project/data/ncbi_meta.csv")
glue <- read.delim("~/Dropbox/rabies/Africa_Project/data/RABVGlue_data.txt")

meta$Geo_Location <- sub(":.*", "", meta$Geo_Location)

meta <- left_join(meta, glue, by = c('Accession'='sequence.sequenceID'))
meta <- meta %>% select(Accession, Length, alignment.name,Isolate,sequence.m49_country.display_name, sequence.collection_year, sequence.host, Geo_Location, Host, Collection_Date)

names(meta)[3] <- 'Clade'
names(meta)[4] <- 'Isolate'
names(meta)[5] <- 'Country'
names(meta)[6] <- 'Collection_Date'
names(meta)[7] <- 'Host'
names(meta)[8] <- 'Location_GenBank'
names(meta)[9] <- 'Host_GenBank'
names(meta)[10] <- 'Collection_Date_Genbank'

meta$Country <- ifelse(is.na(meta$Country), meta$Location_GenBank, meta$Country)
meta$Country <- gsub("USA", "United States", meta$Country)
meta$Host <- ifelse(is.na(meta$Host), meta$Host_GenBank, meta$Host)
meta$Collection_Date <- ifelse(is.na(meta$Collection_Date), substr(meta$Collection_Date_Genbank, 1, 4), meta$Collection_Date)

meta <- meta %>% select(-Location_GenBank, -Host_GenBank, -Collection_Date_Genbank)

#write.table(meta, file = 'data/metadata.tab', quote = FALSE, sep = '\t', row.names =FALSE)

## Creating a simpler column for clade definitions
meta$clade_simple <- ifelse(str_detect(meta$Clade, "Bats"), "Bat_Clade", 
                            ifelse(str_detect(meta$Clade, "Africa_2"), "Africa2_Clade", 
                                   ifelse(str_detect(meta$Clade, "Africa_3"), "Africa3_Clade", 
                                          ifelse(str_detect(meta$Clade, "Arctic"), "Arctic_Clade", 
                                                 ifelse(str_detect(meta$Clade, "Cosmo"), "Cosmopolitan_Clade", 
                                                        ifelse(str_detect(meta$Clade, "Asian"), "Asian_Clade", meta$Clade))))))



## NEED TO ADD COUNTRY INFORMATION FROM KISSY1995 and TALBI 2010

#import both csv with data from isolate number and country (note: you need isolate column *see above for edit)
isolate_country <- read.csv("~/Dropbox/rabies/Africa_Project/data/isolate_countrycode.csv")

# Extract the three-letter country code from the 'Isolate' column
meta$country_code <- sub(".*([A-Z]{3}).*", "\\1", meta$Isolate)


# Map the extracted country code to the corresponding country name
meta <- meta %>%
  left_join(isolate_country, by = c("country_code" = "COUNTRY_CODE")) %>%
  select(-country_code)

#combine country from isolate now with original country
meta$Country <- ifelse(meta$Country == '-', meta$COUNTRY_KISSY, meta$Country)

meta <- meta %>% select(-COUNTRY_KISSY)

## Removing sequences with data quality issues
meta <- meta %>% filter(!is.na(Country)) %>% filter(!is.na(Collection_Date)) %>% 
  filter(Collection_Date > 1971) %>% filter(Length > 99) %>% filter(Country != 'not collected')

meta$region23<- countryname(meta$Country, destination = "region23")

## Good from here for filtering 


## Creating text files for sequence accession numbers that contain nucleotides in a subgenomic region

aln = read.fasta(opt$aln)

# Convert alignment
aln = as.alignment(
  nb = length(aln), 
  nam = sapply(aln, function(x) attributes(x)$name),
  seq = lapply(aln, function(x) {
    attributes(x) = NULL
    x
  })
)

# Trim alignment to positions each gene position (or WGS)
aln$seqWG = lapply(aln$seq, function(x) x[72:min(length(x), 11000)])
aln$seqN = lapply(aln$seq, function(x) x[72:min(length(x), 1424)])
aln$seqP = lapply(aln$seq, function(x) x[1512:min(length(x), 2405)])
aln$seqM = lapply(aln$seq, function(x) x[2494:min(length(x), 3102)])
aln$seqG = lapply(aln$seq, function(x) x[3315:min(length(x), 4889)])
aln$seqL = lapply(aln$seq, function(x) x[5408:min(length(x), 11000)])




# Verify that all sequences contain the N gene
Ngenes = unlist(sapply(aln$nam, function(x) { 
  s = aln$seqN[[x]]
  if (sum(s %in% "-") < 1300) {
    gsub("_.*$", "", x)
  } else {
    NULL
  }
}
))

# Verify that all sequences contain the P gene
Pgenes = unlist(sapply(aln$nam, function(x) { 
  s = aln$seqP[[x]]
  if (sum(s %in% "-") < 890) {
    gsub("_.*$", "", x)
  } else {
    NULL
  }
}
))

# Verify that all sequences contain the M gene
Mgenes = unlist(sapply(aln$nam, function(x) { 
  s = aln$seqM[[x]]
  if (sum(s %in% "-") < 605) {
    gsub("_.*$", "", x)
  } else {
    NULL
  }
}
))

# Verify that all sequences contain the G gene
Ggenes = unlist(sapply(aln$nam, function(x) { 
  s = aln$seqG[[x]]
  if (sum(s %in% "-") < 1570) {
    gsub("_.*$", "", x)
  } else {
    NULL
  }
}
))

# Verify that all sequences contain the L gene
Lgenes = unlist(sapply(aln$nam, function(x) { 
  s = aln$seqL[[x]]
  if (sum(s %in% "-") < 5590) {
    gsub("_.*$", "", x)
  } else {
    NULL
  }
}
))



#######
meta$fragment <- ifelse(meta$Length > 10000, 'WGS', 
                        ifelse(meta$Accession %in% Ngenes, 'N', 
                               ifelse(meta$Accession %in% Pgenes, 'P',
                                      ifelse(meta$Accession %in% Mgenes, 'M',
                                             ifelse(meta$Accession %in% Ggenes, 'G',
                                                    ifelse(meta$Accession %in% Lgenes, 'L', 'removed'))))))

write_delim(meta, 'data/meta_edited_July2024.tab', delim = '\t', quote = 'none')
###END UPDATE 24 JULY 2024


meta_n <- meta %>% filter(meta$fragment == 'N') %>% select(Accession) %>% 
  write.csv(opt$out_n_text,
            row.names = FALSE, quote = FALSE, col.names = NA)
meta_p <- meta %>% filter(meta$fragment == 'P') %>% select(Accession) %>% 
  write.csv(opt$out_p_text,
            row.names = FALSE, quote = FALSE, col.names = FALSE)
meta_m <- meta %>% filter(meta$fragment == 'M') %>% select(Accession) %>% 
  write.csv(opt$out_m_text,
            row.names = FALSE, quote = FALSE, col.names = FALSE)
meta_g <- meta %>% filter(meta$fragment == 'G') %>% select(Accession) %>% 
  write.csv(opt$out_g_text,
            row.names = FALSE, quote = FALSE, col.names = FALSE)
meta_l <- meta %>% filter(meta$fragment == 'L') %>% select(Accession) %>% 
  write.csv(opt$out_l_text,
            row.names = FALSE, quote = FALSE, col.names = FALSE)
meta_wgs <- meta %>% filter(meta$fragment == 'WGS') %>% select(Accession) %>% 
  write.csv(opt$out_wgs_text,
            row.names = FALSE, quote = FALSE, col.names = FALSE)



## Remove sequences from G cut alignment for just the G and WGS genes (repeat)
## Hard coded
aln_n = read.fasta('/Users/aholtz/Dropbox/rabies/Africa_Project/data/sequence_alignments/gene_specific_analysis/cutalign_N.fa')
N_wgs <- read.table("~/Dropbox/rabies/Africa_Project/data/sequence_alignments/gene_specific_analysis/N_wgs.txt", quote="\"", comment.char="")
n_wgs_alignment <- aln_n[which(names(aln_n) %in% N_wgs$V1)]

aln_m = read.fasta('/Users/aholtz/Dropbox/rabies/Africa_Project/data/sequence_alignments/gene_specific_analysis/cutalign_M.fa')
M_wgs <- read.table("~/Dropbox/rabies/Africa_Project/data/sequence_alignments/gene_specific_analysis/M_wgs.txt", quote="\"", comment.char="")
m_wgs_alignment <- aln_m[which(names(aln_m) %in% M_wgs$V1)]

aln_p = read.fasta('/Users/aholtz/Dropbox/rabies/Africa_Project/data/sequence_alignments/gene_specific_analysis/cutalign_P.fa')
P_wgs <- read.table("~/Dropbox/rabies/Africa_Project/data/sequence_alignments/gene_specific_analysis/P_wgs.txt", quote="\"", comment.char="")
p_wgs_alignment <- aln_p[which(names(aln_p) %in% P_wgs$V1)]


aln_g = read.fasta('/Users/aholtz/Dropbox/rabies/Africa_Project/data/sequence_alignments/gene_specific_analysis/cutalign_G.fa')
G_wgs <- read.table("~/Dropbox/rabies/Africa_Project/data/sequence_alignments/gene_specific_analysis/G_wgs.txt", quote="\"", comment.char="")
g_wgs_alignment <- aln_g[which(names(aln_g) %in% G_wgs$V1)]

aln_l = read.fasta('/Users/aholtz/Dropbox/rabies/Africa_Project/data/sequence_alignments/gene_specific_analysis/cutalign_L.fa')
L_wgs <- read.table("~/Dropbox/rabies/Africa_Project/data/sequence_alignments/gene_specific_analysis/L_wgs.txt", quote="\"", comment.char="")
l_wgs_alignment <- aln_l[which(names(aln_l) %in% L_wgs$V1)]


write.fasta(sequences = n_wgs_alignment, names = names(n_wgs_alignment), file.out = "cutalign_N_only.fa")
write.fasta(sequences = m_wgs_alignment, names = names(m_wgs_alignment), file.out = "cutalign_M_only.fa")
write.fasta(sequences = p_wgs_alignment, names = names(p_wgs_alignment), file.out = "cutalign_P_only.fa")
write.fasta(sequences = g_wgs_alignment, names = names(g_wgs_alignment), file.out = "cutalign_G_only.fa")
write.fasta(sequences = l_wgs_alignment, names = names(l_wgs_alignment), file.out = "cutalign_L_only.fa")

########
########
########
########


## It was now deceided to remove sequences besides for all African origin (country) and 
## part of the African clade. Plus, 5 sequences from each additional clade to give context
## to help with rooting


###Updated below July 24, 2024
meta_edited <- meta

africa_clades <- meta_edited %>% 
  filter(str_detect(Clade, 'Africa|AF1') | str_detect(region23, 'Africa'))
random_other_clades <- meta_edited %>% filter(fragment == 'WGS') %>% 
  group_by(Clade) %>%
  sample_n(size = 4, replace = TRUE) %>%
  ungroup()

africa_total <- rbind(africa_clades, random_other_clades)

#ALIGNMENT FROM MAFFT
aln = read.fasta("/Users/aholtz/Dropbox/rabies/Africa_Project/data/concat_seq_genes.fasta")

africa_alignment <- aln[which(names(aln) %in% africa_total$Accession)]

write.fasta(sequences = africa_alignment, names = names(africa_alignment), file.out = "africa_aln.fa")


## From result of downsample_monophyl_cluster.R script to reduce sequences
sub_meta <- read.delim("~/Dropbox/rabies/Africa_Project/sub_meta.tab")

africa_sub <- aln[which(names(aln) %in% sub_meta$Accession)]
write.fasta(sequences = africa_sub, names = names(africa_sub), file.out = "africa_sub_aln.fa")


## Analysis of the sequence database is represented by visualizations shown below

DT::datatable(africa_total, extensions = 'Buttons', options = list(
  autoWidth = TRUE, 
  pageLength = 10,
  dom = 'Bfrtip',
  buttons = c('copy', 'csv', 'excel', 'pdf', 'print')),
  filter = list(
    position = 'top', clear = FALSE))

host <- africa_metadata_hostDog %>% group_by(dog) %>% tally() %>% 
  arrange(desc(n)) %>% 
  slice(1:10)

host$host_simple <- factor(host$dog, levels = 
                             host$dog[order(-host$n)])

host <- ggplot(host, aes(dog,n)) + geom_bar(stat = "identity", fill = "#00879f") +
  labs(title = "Host Species Full Data Set", x = "Host Species", y =  "Count") + theme_bw() +
  theme(axis.text.x = element_text(angle = 70)) 

ggplotly()

country<- africa_total %>% group_by(Country) %>% tally() %>% 
  arrange(n) 

country$Country <- factor(country$Country, levels = 
                            country$Country[order(-country$n)])

country <- ggplot(country, aes(Country,n)) + geom_bar(stat = "identity", fill = "#2DA37F") +
  labs(title = "Country of Isolation Full Data Set", x = "Country of Isolation", y =  "Count") + theme_bw() +
  theme(axis.text.x = element_text(angle = 70)) 

ggplotly()


gene <- africa_total %>% group_by(fragment) %>% tally() %>% 
  arrange(n) 

gene$fragment <- factor(gene$fragment, levels = 
                          gene$fragment[order(-gene$n)])

gene <- ggplot(gene, aes(fragment, n)) + geom_bar(stat = "identity", fill = "#745185") +
  labs(title = "RABV Sequence Fragment (>200bp) Full Data Set", x = "Gene(s)", y =  "Count") +
  theme_bw()

ggplotly()


clade <- africa_total %>% group_by(Clade) %>% tally() %>% 
  arrange(n) 

cladeCclade <- factor(clade$Clade, levels = 
                               clade$Clade[order(-clade$n)])

clade <- ggplot(clade, aes(Clade, n)) + geom_bar(stat = "identity", fill = "#c6f754") +
  labs(title = "RABV Clade Classification Full Data Set", x = "Clade", y =  "Count") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 70)) 

ggplotly()
