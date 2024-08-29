# plot acinetobacter subset with poppunk lineages 
# Load Packages
library(ape)
library(tidyverse)
library(ggtree)
library(phytools)
library(janitor)
library(treeio)
library(ggnewscale)
library(ggtreeExtra)
library(ggstar)

# Read Tree
Acinetobacter_tree <- ape::read.tree("../09_Acinetobacter_Phylogeny_Subset /gtdbtk.bac120.user_msa.fasta_2.treefile")
Acinetobacter_tree_root <- root(Acinetobacter_tree, "VD15_contigs")

# Rooted but drop the root for the display
Acinetobacter_tree_root <- drop.tip(Acinetobacter_tree_root, tip = "VD15_contigs")

# rename tip labels
Acinetobacter_tree_tib <- as_tibble(Acinetobacter_tree_root) %>%
  #separate(col = label, into = c("junk", "label2"), sep = "../02_16s_extracted_genes/", remove = F) %>%
  separate(col = label, into = c("label2", "junk"), sep = "_contigs", remove = F) %>%
  select(label, label2) 

Acinetobacter_tree_rename <- rename_taxa(Acinetobacter_tree_root, Acinetobacter_tree_tib, label, label2)
#MosAIC_tree_rename %>%
#print(n = 800)

dev.new(width=6, height=15)

Acinetobacter_ggtree <- ggtree(Acinetobacter_tree_rename, layout="rectangular", size=0.2) + 
  xlim_tree(xlim = 0.5) + 
  geom_treescale(x = 0.01, y = -1) + 
  ggtree::vexpand(.1, -1)+ 
  geom_tiplab(offset = 0.001, size = 3)

PopPunk <- read_csv("bgmm_K_3/bgmm_K_3_clusters.csv")


## Function to Make a Renamed Tree based on PopPUNK Data 
RenameTreeTipsWithPopPUNK <- function(tree_data, PopPUNK_data){
  # Fix tip Labels
  tree_data_tib <- as_tibble(tree_data) %>%
    mutate(label2 = str_sub(label, 1, 15)) %>%
    separate(col = label2, into = c("label2", "junk"), sep = "_contig", remove = F) %>%
    select(label, label2) 
  
  # Fix PopPUNK Labels
  PopPUNK_data_clean <- PopPUNK_data %>%
    mutate(Taxon = str_sub(Taxon, 1, 15)) %>%
    separate(col = Taxon, into = c("Taxon", "junk"), sep = "_contig", remove = F) %>%
    select(Taxon, Cluster)
  
  # Join the Datasets
  tree_tib_join <- tree_data_tib %>%
    left_join(PopPUNK_data_clean, by = c("label2" = "Taxon")) %>%
    select(label, Cluster)
  
  # Rename Tips with PopPUNK Clusters
  tree_rename <- rename_taxa(tree_data, tree_tib_join, label, Cluster)  
  
  return(tree_rename)
  
}


tree <- RenameTreeTipsWithPopPUNK(Acinetobacter_tree_root, PopPunk)

g1 <- ggtree(tree) + geom_tiplab()



ggsave("Fig_S14_AcinetobacterPhylogenyWithPopPUNKSubset.pdf", height = 10, width = 10)

# alternative with label tips as well 

Acinetobacter_ggtree

PopPUNK_data_clean <- as.data.frame(PopPunk %>%
                                      mutate(Taxon = str_sub(Taxon, 1, 15)) %>%
                                      separate(col = Taxon, into = c("Taxon", "junk"), sep = "_contig", remove = F) %>%
                                      select(Taxon, Cluster) %>% 
                                      mutate(hold = 1))

Acinetobacter_ggtree + geom_fruit(
  data = PopPUNK_data_clean, 
  geom = geom_text, 
  mapping = aes(x = hold, y = Taxon, label = Cluster),
  offset = 0.2, 
  pwidth = 0.1,
  size = 3,
  grid.params=list(
    linetype=8,
    size=2, 
    color = "black", 
    lineend = "square", 
    alpha = 0
  ) 
) + 
  scale_fill_manual(values=c("black","black", "black", "#A6BD9D","#FFB177"),
                    guide=guide_legend(keywidth = 0.3, 
                                       keyheight = 0.3, order=4)) 
ggsave("Fig_S10_AcinetobacterPhylogenyWithPopPUNKSubsetAndNames.pdf", height = 7, width = 10)
