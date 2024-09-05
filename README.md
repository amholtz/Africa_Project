# About

This repository contains scripts and data that were used to analyze the African spread of canine-rabies virus. This work is ongoing.



---
# Introduction

The following work was performed with [MAFFT(v7.505)](https://doi.org/10.1093/nar/gkf436), [FastTree(v2.1.11)](https://doi.org/10.1371/journal.pone.0009490), [Goalign(v0.3.5)](https://github.com/evolbioinfo/goalign), [Gotree(v0.4.4)](https://github.com/evolbioinfo/gotree), [TempEst(v1.5.3)](https://doi.org/10.1093/ve/vew007), [HyPhy(v2.5.40)](https://github.com/veg/hyphy), [IQTREE2(v2.2.2.2)](10.1093/molbev/msaa015), [LSD2(v1.8.8)](https://doi.org/10.1093/sysbio/syv068), [BEAST X v10.5.0-beta4](https://github.com/beast-dev/beast-mcmc/releases/tag/v10.5.0-beta4), and [iTol](https://itol.embl.de/tree/1579917420235811657296942#). In addition custom scripts in R and Python were used, which can be found in [R](https://github.com/amholtz/GlobalRabies/tree/main/R) and [Python](https://github.com/amholtz/GlobalRabies/tree/main/python) folders. R version 4.2.1 was used with the following packages: dplyr, tidyverse, ggplot2, plotly, treeio, phangorn, cepiigeogist, countrycode, reshape, data.table, DT, optparse, lubridate, seqinR, readr, taxize, rworldmap, googleVis,rgdal, scales, wesanderson, ape, and Quartet. Python version 3.8 was used with the following packages: numpy, pandas, and pastml   

The intermediate data files can be found in the [data folder](https://github.com/amholtz/Africa_Project/tree/main/data). To reproduce the analyses from scratch follow the instructions below.


### Set up

Data can  be downloaded from the [data folder](https://github.com/amholtz/Africa_Project/tree/main/data).

The initial sequence dataset was downloaded from [NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=Lyssavirus%20rabies,%20taxid:11292) by  searching for taxid:11292.


### Sequence alignment
1.  Global alignment by [MAFFT(v7.505)](https://doi.org/10.1093/nar/gkf436)
  ```
  mafft --reorder --keeplength --maxambiguous 0.05 --addfragments data/allRABV.fasta --auto data/rabv_reference_1988.fasta  > data/with_keeplength_RABV.fasta
  ```
### Filtering Metadata and Creating Concatenated Sequence Alignment

1.  The custom script [DataCleaning.R] also organized sequences by subgenomic region. NC_001542, the reference genome was cut at the positions in the table below [(partition_RABVGenes.txt)](https://github.com/amholtz/Africa_Project/blob/main/data/sequence_alignments/gene_specific_analysis/partition_RABVGenes.txt) which represent start and stop codons for each gene. As an example, sequences categorized as G gene, contain more than 200 nucleotides between start and stop codons and were saved as a new line in a text file. A quality check was conducted to remove sequences that were (1) missing date and country information, (2) older than 1972, (3) identified as vaccine or laboratory strains, (4) with coding regions shorter than 200 nucleotides. As a result, 19,329 sequences were retained for this study.

  | Gene      | Position Start | Position End |
  |-----------|----------------|--------------|
  | N protein | 71             | 1423         |
  | P protein | 1514           | 2407         |
  | M protein | 2496           | 3104         |
  | G protein | 3318           | 4892         |
  | L protein | 5418           | 11846        |

```
 Input:
 meta <- read.csv("Africa_Project/data/ncbi_meta.csv")
 glue <- read.delim("Africa_Project/data/RABVGlue_data.txt")
 aln <- read.fasta("Africa_Project/data/allRABV.fa")
```

4. The sequence alignment was  split into 4 different files, representing the coding regions of the 4 genes defined in [partition_RABVGenes.txt](https://github.com/amholtz/Africa_Project/blob/main/data/sequence_alignments/gene_specific_analysis/partition_RABVGenes.txt).
  ```
  goalign split -i data/with_keeplength_RABV.fasta --partition data/sequence_alignments/gene_specific_analysis/partition_RABVGenes.txt --out-prefix data/sequence_alignments/gene_specific_analysis/cutalign_
  ```

5.  Each gene-specific alignment is then subsected for sequences grouped in step 3 using [Goalign(v0.3.5)](https://github.com/evolbioinfo/goalign) along with WGS  (G gene example)

  ```
  cat data/sequence_alignments/gene_specific_analysis/g.txt data/sequence_alignments/gene_specific_analysis/wgs.txt > data/sequence_alignments/gene_specific_analysis/G_wgs.txt

  goalign subset -i data/sequence_alignments/gene_specific_analysis/cutalign_G.fa -f data/sequence_alignments/gene_specific_analysis/G_wgs.txt --unaligned -o data/sequence_alignments/gene_specific_analysis/cutalign_G_only.fa
  ```

5.  Each gene is then aligned independently by [MAFFT(v7.505)](https://doi.org/10.1093/nar/gkf436) according to the cut reference sequence (Ex: G gene from reference + all RABV sequences that were classified as G gene)
```
mafft --reorder --keeplength --maxambiguous 0.05 --addfragments data/sequence_alignments/gene_specific_analysis/cutalign_G_only.fa --auto data/sequence_alignments/gene_specific_analysis/Ggene_Ref.fa > data/sequence_alignments/gene_specific_analysis/Ggene_aln.fa
```

6. [Goalign(v0.3.5)](https://github.com/evolbioinfo/goalign) concat was then used to concatenate the aligned sequences back together (without noncoding regions)
```
goalign concat -i data/sequence_alignments/gene_specific_analysis/Ngene_aln.fa data/sequence_alignments/gene_specific_analysis/Pgene_aln.fa data/sequence_alignments/gene_specific_analysis/Mgene_aln.fa data/sequence_alignments/gene_specific_analysis/Ggene_aln.fa data/sequence_alignments/gene_specific_analysis/Lgene_aln.fa -o data/concat_seq_genes.fasta
```
Note: Download [concat_seq_genes](https://www.dropbox.com/scl/fi/mjcvf6viux6zocm54uw8n/concat_seq_genes.fasta?rlkey=uo36rvsjayeydflnuzxzqggp7&dl=0)

7. The custom script [DataCleaning.R] creates a subset of sequences that are African in origin. We defined African sequences as those with African origin country in metadata or sequences that are identified by RABVGlue as being in an African clade (Africa 1,2,3,4). Five random sequences from each other clade were conserved to help with rooting and dating of the tree. A subsect of the sequence alignment is produced as a new fasta file with 4,632 sequences (data/africa_aln_July2024.fa).


### Phylogenetic Tree Reconstruction & Dating
1. A global phylogenetic tree was reconstructed using [FastTree(v2.1.11)](https://doi.org/10.1371/journal.pone.0009490)  on African sequences (data/africa_aln_July2024.fa)
```
~/FastTreeMP -gtr -gamma -nt data/africa_aln_July2024.fa > data/Africa_RABV_July2024.nwk
```

2.  Bat clades were removed from the tree using iTol

3. Rooting was accomplished via [TempEst(v1.5.3)](https://doi.org/10.1093/ve/vew007) by residual-mean square function. Date file [tempest_dates.tab](https://github.com/amholtz/Africa_Project/blob/main/data/tempest_dates.tab) is adapted to [Tempest format](https://beast.community/tempest_tutorial)  from our [metadata file](https://github.com/amholtz/Africa_Project/blob/main/data/meta_edited_July2024.tab)

  **Input**: [Africa_RABV_July2024.nwk_collapsed_0.5.nwk](https://github.com/amholtz/Africa_Project/blob/main/data/Africa_RABV_July2024.nwk_collapsed_0.5.nwk), [tempest_dates.tab](https://github.com/amholtz/Africa_Project/blob/main/data/tempest_dates.tab)

  282 outliers were detected after Tempest Analysis. Here, we estimated the best-fitting root in Tempest. All sequences that were 0.03 residuals or more away from the mean were identified as outliers. The following sequences were removed:

  LC682854 LC682856 LC683186 MH507336 MF197744 LC682857 MF197745 MH507337 KP662552 KY681372 GQ918329 EU643570 JN162084 JQ692999 LC682839 LC683187 MZ417923 MZ418084 JQ692996 KX148255 MZ418099 MZ417954 MZ417991 MZ418081 LC683169 MZ417939 MZ417885 MN196561 MZ417969 GQ983487 MZ417951 MZ417889 MZ417877 MZ417907 GQ983502 JN162076 KX148260 MZ418078 MZ417904 MZ417982 LC683163 KF620487 MZ417911 MN196559 MZ418098 MZ417894 MZ417919 MZ417998 MZ417994 GQ983520 KY681359 MZ417992 MZ417902 GU647092 MN196564 MZ417876 GQ983459 MN196568 KC535504 MZ417881 MZ418000 MZ418082 KX148264 HQ450385 MZ418079 MZ417922 MZ417958 MZ417888 MZ417883 MZ417952 AB569299 MZ417901 MZ418088 MZ417900 MZ417993 JX088736 MZ417931 MZ417899 MZ417989 MZ417940 MZ417977 MZ417959 MZ417973 MZ417915 MZ417886 MZ417966 FJ712196 MZ417962 JN162088 MZ417941 MZ417986 MZ417976 MZ417946 MZ417920 MZ418002 MZ418036 MZ418006 MZ418004 MZ417882 MZ417934 KX148245 MZ417963 MZ418097 MZ417926 MZ417927 JQ970485 KX148254 MZ417935 MN196572 MZ418076 MN196570 GU345746 FJ465384 FJ465397 MZ417949 MZ417890 MZ417891 MZ417880 MZ417936 KF154999 MZ417943 LC682830 MZ417953 KX148267 JQ692982 KY451767 KX148252 MZ417933 MZ417896 GU936878 MZ417878 MN196571 MZ417965 MZ417912 MZ417975 EU888773 MG201922 MZ418089 MZ417909 MK689675 KF620489 MZ417980 MZ417950 MZ417893 MZ417887 JQ692992 MZ417918 MZ417968 LC683162 EU643571 MZ417929 MZ417938 MZ417944 MZ418083 MZ417916 FJ712195 MZ418001 JX088737 MZ417990 MZ417910 MZ417917 MZ417932 MF537558 MZ418090 MZ418080 MZ418077 MZ417895 MZ417972 MZ417945 MZ417961 MZ417930 MZ417928 EU888744 MZ417897 MZ417964 LC682828 MZ417925 MZ417985 FJ465393 KC535506 MZ417956 FJ465388 LC682829 MN196563 MN196566 MZ417960 JQ730682 LC380154 MZ417914 MZ417978 MZ417948 MN196565 MZ417971 GQ918323 LC683166 MN196569 JQ692993 EU643589 MN196562 MN726845 MZ417983 MZ418087 MZ417921 GQ983492 MZ417974 MZ417892 MZ417884 MZ417898 MZ417908 MZ417999 MZ417967 MN726838 MZ417913 MZ417905 MN196560 LC682826 EU643590 MZ417937 MZ417987 MZ417906 MZ417947 MF537559 MZ417981 MZ417979 MZ417988 MZ417879 LC682842 MZ418003 MN726843 EU038107 MZ417942 GQ983444 MZ417957 JQ692997 MZ417970 MN196567 DQ786033 MZ418031 MZ418034 MZ418029 MZ418008 MZ418075 MZ418028 MZ418032 MZ418035 MZ418033 MZ418030 MZ418007

  **Output**: [bats_removed_pruned_TempestRooted.tree](https://github.com/amholtz/Africa_Project/blob/main/data/bats_removed_pruned_TempestRooted.tree)

4. Sequences were subsampled by removing similar sequences within the same cluster using custom script by Sam Hong at KU Leuven [downsample_monophyl_cluster.R](https://github.com/amholtz/Africa_Project/blob/main/data/subsampling_SamScript/downsample_monophyl_cluster.R)

  **Input**: [bats_removed_pruned_TempestRooted.tree](https://github.com/amholtz/Africa_Project/blob/main/data/bats_removed_pruned_TempestRooted.tree)

  **Output**: [subsampled_tree.tre](https://github.com/amholtz/Africa_Project/blob/main/data/subsampling_SamScript/subsampled_tree.tre)

  Goalign subset was used to create the new fasta with just the pruned sequences.

  ```
  goalign subset -i data/africa_aln_July2024.fa -f /data/subsampling_SamScript/keep.csv -o sub_aln.fasta
  ```
5. [BEAUTi v10.5.0-beta4](https://github.com/beast-dev/beast-mcmc/releases/tag/v10.5.0-beta4) was used to create the XML for BEAST. First running one BEAST rendition without any metadata to create a large tree file to use as an empirical tree input for the GLM analysis to save time.

  [sub_aln.fasta](https://github.com/amholtz/Africa_Project/blob/main/data/subsampling_SamScript/sub_aln.fasta) imported as partition (636 taxa), [tempest_dates.tab](https://github.com/amholtz/Africa_Project/blob/main/data/tempest_dates.tab) imported for dates.

**BEAUTi Configuration:**

#### Sites
| Page      | Prompt | Selection |
|-----------|----------------|--------------|
| Sites | Substitution model             | GTR         |
| Sites | Base frequencies           | estimated         |
| Sites | Site heterogeneity model           | Gamma (equal weights) + Invariant sites         |
| Sites | Number of Gamma Categories           | 4         |
| Sites | Parition into codon positions           | OFF        |

#### Clocks
| Page   | Prompt                        | Selection                       |
|--------|-------------------------------|---------------------------------|
| Clocks | Clock type                    | Hamiltonian Monte Carlo relaxed Clock |
| Clocks | Relaxed distribution          | Lognormal                       |

#### Trees
| Page  | Prompt                         | Selection                       |
|-------|--------------------------------|---------------------------------|
| Trees | Tree Prior                     | Coalescent: Hamiltonian Monte Carlo SkyGrid |
| Trees | Number of parameters           | 50                              |
| Trees | Time at last transition point  | 673                             |

#### MCMC
| Page | Prompt              | Selection   |
|------|---------------------|-------------|
| MCMC | Length of chain     | 500000000   |
| MCMC | Log parameters every| 50000       |
| MCMC | Checkpoint every    | 1000000     |

6. [BEAST X v10.5.0-beta4](https://github.com/beast-dev/beast-mcmc/releases/tag/v10.5.0-beta4) is used on the HPC using the following command:

```
beast -beagle_cuda -beagle_GPU -beagle_double -save_every 10000 -save_state beast_Africa_RABV_1_checkpoint.state Africa_RABV_sub_aln.xml
```

*Note:* It's best to set up sbatch files to run on the HPC with modules loaded and the commands inside. See [africa_aln_beast.sh] for an example on how to build this file.

# Get in touch

To report any issues, please [open an isssue](https://github.com/amholtz/GlobalRabies/issues).
