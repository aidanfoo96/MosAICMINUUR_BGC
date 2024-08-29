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
Acinetobacter_tree <- ape::read.tree("gtdbtk.bac120.user_msa.fasta_2.treefile")
Acinetobacter_tree_root <- root(Acinetobacter_tree, "VD15_contigs")

# Rooted but drop the root for the display
Acinetobacter_tree_root <- drop.tip(Acinetobacter_tree_root, tip = "VD15_contigs")

# read poppunk 
PopPunk <- read_csv("../11_Acinetobacter_PopPUNK/bgmm_K_3/bgmm_K_3_clusters.csv")

# Read antismash
MosAIC_MINUUR_antismash <- read_tsv("../05_MosAIC_BGC/030424_MINUUR_MosAIC_classifications_antismash.tsv")

# Read MosAIC & MINUUR metadata 
MosAIC_MINUUR_metadata <- read_tsv("190124_MINUUR_MosAIC_classifications.tsv") 

# Read MosAIC complete metadata 
MosAIC_metadata_complete <- read_csv("Nov_S1_Table_MosAIC_Metadata_Final.csv")

# Read BigSCAPE results 
bigscape_clusters <- read_tsv("../06_BigScapeAnalysis/200524_antismash_hits_gcfs_classifications_MiBIG.tsv")

# rename tip labels
Acinetobacter_tree_tib <- as_tibble(Acinetobacter_tree_root) %>%
  #separate(col = label, into = c("junk", "label2"), sep = "../02_16s_extracted_genes/", remove = F) %>%
  separate(col = label, into = c("label2", "junk"), sep = "_contigs", remove = F) %>%
  select(label, label2) 

Acinetobacter_tree_rename <- rename_taxa(Acinetobacter_tree_root, Acinetobacter_tree_tib, label, label2)
#MosAIC_tree_rename %>%
#print(n = 800)

dev.new(width=6, height=15)

# clean poppunk 
PopPUNK_data_clean <- as.data.frame(PopPunk %>%
  mutate(Taxon = str_sub(Taxon, 1, 15)) %>%
  separate(col = Taxon, into = c("Taxon", "junk"), sep = "_contig", remove = F) %>%
  select(Taxon, Cluster) %>% 
  mutate(hold = 1))

# Join the Datasets
label_PopUNK <- Acinetobacter_tree_rename %>%
  left_join(PopPUNK_data_clean, by = c("label" = "Taxon")) %>%
  select(label, Cluster)

# make tree
Acinetobacter_ggtree <- ggtree(Acinetobacter_tree_rename, layout="rectangular", size=0.2) + 
  xlim_tree(xlim = 0.5) + 
  geom_treescale(x = 0.01, y = -1) + 
  ggtree::vexpand(.1, -1) 
  #geom_tiplab(size = 3, hjust = -0.5)

ggsave("Acinetobacter_tree_with_labels.pdf")

# function to extract column and convert for gheatmap function
extract_column_and_convert <- function(data, column){
  target <- data[, c(1, column)]
  target2 <- target %>%
    remove_rownames %>%
    column_to_rownames(var = "user_genome")
  
  return(target2)
}

MosAIC_species <- extract_column_and_convert(MosAIC_MINUUR_metadata, 7)

# simplify the chart the product category 
MosAIC_MINUUR_antismash %>%
  filter(genus == "Acinetobacter") %>%
  group_by(product) %>%
  summarise(n())

MosAIC_MINUUR_antismash_simplified <- MosAIC_MINUUR_antismash %>%
  mutate(product = if_else(product == "NRP-metallophore / NRPS", "NRP-metallophore", product)) %>%
  mutate(product = if_else(product == "NRP-metallophore / NRPS / redox-cofactor", "NRP-metallophore", product)) %>%
  mutate(product = if_else(product == "NRPS / hserlactone", "hserlactone", product)) %>%
  mutate(product = if_else(product == "NRPS / hserlactone / arylpolyene", "hserlactone", product)) %>%
  mutate(product = if_else(product == "arylpolyene / NRPS / hserlactone", "arylpolyene", product)) %>%
  mutate(product = if_else(product == "redox-cofactor / NRP-metallophore / NRPS", 'NRP-metallophore', product)) 

MosAIC_antismash_summarised = as.data.frame(MosAIC_MINUUR_antismash_simplified %>%
                                              filter(genus == "Acinetobacter") %>%
                                              select(sampleID, product) %>%
                                              group_by(sampleID, product) %>%
                                              summarise(number = n()) %>%
                                              left_join(PopPUNK_data_clean, by = c("sampleID" = "Taxon")))
                                              

# pivot wider
MosAIC_antismash_summarised_wide = MosAIC_antismash_summarised %>% 
  pivot_wider(names_from = product, values_from = number) %>% 
  column_to_rownames(var = "sampleID")

# Wrangle complete MosAIC metadata 
MosAIC_metadata_complete_edit <- as.data.frame(MosAIC_metadata_complete %>%
                                                 mutate(line = 1))

# Wrangle and filter bigscape clusters for Acinetobacter 
bigscape_clusters_wrangle <- bigscape_clusters %>%
  filter(genus == "Acinetobacter") %>%
  filter(str_detect(`Product Prediction`, "siderophore") | str_detect(`Product Prediction`, "NRP-metallophore")) %>%
  filter(family_number != 2805) %>%
  mutate(SiderophoreClass = case_when(
    str_detect(`Product Prediction`, "siderophore") ~ "NI-Siderophore",
    str_detect(`Product Prediction`, "metallophore") ~ "NRP-Siderophore",
    TRUE ~ "Other"
  )) %>%
  select(sampleID, family_number, SiderophoreClass) %>%
  mutate() %>%
  group_by(family_number) %>%
  mutate(gcf_identifier = cur_group_id()) 


bigscape_clusters_wrangle_clean = as.data.frame(bigscape_clusters_wrangle %>%
                                              select(sampleID, gcf_identifier, SiderophoreClass))

bigscape_clusters_wrangle_clean$gcf_identifier <- as.factor(bigscape_clusters_wrangle_clean$gcf_identifier)

P1 <- gheatmap(Acinetobacter_ggtree, MosAIC_antismash_summarised_wide,
         legend_title = "BGCs (n)",
         colnames_position = "bottom",
         colnames_angle = 30,
         font.size = 3, 
         hjust = 1,
         colnames_offset_y = 0, 
         color = "black", 
         width = 1) + 
  scale_fill_gradient2(
    low = "#84a98c",
    mid = "#52796f",
    high = "#ffd166",
    midpoint = 1.5,
    space = "Lab",
    na.value = "white",
    aesthetics = "fill", 
    guide=guide_legend(keywidth = 0.5, 
                       keyheight = 1, 
                       title = "Biosynthetic Gene Clusters (n)")) 

P2 <- P1 + new_scale_fill() + geom_fruit(
  data = MosAIC_metadata_complete_edit, 
  geom = geom_star, 
  mapping = aes(x = line, y = internal_id, starshape=sample_source),
  offset = 1.1, 
  pwidth = 0.1,
  grid.params=list(
    linetype=11,
    size=0.3
  )
) 

P2 + 
  geom_strip('MMO-131', 'HN041', label = "A. pitti", offset = 0.18, offset.text = 0.005) + 
  geom_strip('AA010', 'AA010', label = "A. oleivorans", offset = 0.18, offset.text = 0.005) + 
  geom_strip('MMO-143', 'M51', label = "A. nosocomialis", offset = 0.18, offset.text = 0.005) + 
  geom_strip('HN127', 'HN047', label = "A. modestus", offset = 0.18, offset.text = 0.005) + 
  geom_strip('AS116', 'AS116', label = "A. junii", offset = 0.18, offset.text = 0.005) + 
  geom_strip('AS096', 'AS096', label = "A. beijerinckii", offset = 0.18, offset.text = 0.005) + 
  geom_strip('USHLN143', 'USHLN143', label = "A. fasciculus", offset = 0.18, offset.text = 0.005) + 
  geom_strip('HN020', 'HN020', label = "A. johnsonii", offset = 0.18, offset.text = 0.005) + 
  geom_strip('AS013', 'AS013', label = "A. tandoii", offset = 0.18, offset.text = 0.005) + 
  geom_strip('HN105', 'HN105', label = "A. bereziniae", offset = 0.18, offset.text = 0.005) + 
  geom_strip('AS167', 'AS167', label = "A. spp", offset = 0.18, offset.text = 0.005) +
  geom_strip('USMM068', 'HN121', label = "A. soli", offset = 0.18, offset.text = 0.005) + 
  theme(legend.position = c(0.89, 0.8))

ggsave(filename = "130524_Acinetobacter_BGCSubset_WithAllBGCs_AndNoLabels.pdf")

# simplify even more - just show the siderophore implicated pathways
MosAIC_antismash_summarised_2 = as.data.frame(MosAIC_MINUUR_antismash_simplified %>%
                                              filter(genus == "Acinetobacter") %>%
                                              filter(product == "NI-siderophore" | product == "NRP-metallophore") %>%
                                              select(sampleID, product) %>%
                                              group_by(sampleID, product) %>%
                                              summarise(number = n()))

MosAIC_antismash_summarised_wide = MosAIC_antismash_summarised_2 %>% 
  pivot_wider(names_from = product, values_from = number) %>% 
  column_to_rownames(var = "sampleID")

# Wrangle complete MosAIC metadata 
MosAIC_metadata_complete_edit <- as.data.frame(MosAIC_metadata_complete %>%
                                                 mutate(line = 1))


P1_1 <- gheatmap(Acinetobacter_ggtree, MosAIC_antismash_summarised_wide,
               legend_title = "BGCs (n)",
               colnames_position = "bottom",
               colnames_angle = 30,
               font.size = 3, 
               hjust = 1,
               colnames_offset_y = 0, 
               color = "black", 
               width = 0.5) + 
  scale_fill_gradient2(
    low = "#84a98c",
    mid = "#52796f",
    high = "#ffd166",
    midpoint = 1.5,
    space = "Lab",
    na.value = "white",
    aesthetics = "fill", 
    guide=guide_legend(keywidth = 0.5, 
                       keyheight = 1, 
                       title = "Biosynthetic Gene Clusters (n)")) 

P1_1 <- Acinetobacter_ggtree + geom_fruit(
  data = PopPUNK_data_clean, 
  geom = geom_text, 
  mapping = aes(x = hold, y = Taxon, label = Cluster),
  offset = -0.05, 
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


P2_1 <- P1_1 + geom_fruit(
  data = MosAIC_metadata_complete_edit, 
  geom = geom_star, 
  mapping = aes(x = line, y = internal_id, starshape=sample_source, fill = sample_source),
  offset = 0, 
  pwidth = 0.1,
  size = 2,
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

P3_1 = P2_1 + geom_fruit(data = MosAIC_antismash_summarised_2, geom=geom_tile, 
                                  mapping = aes(y = sampleID, x = product, fill = product, height = 0.9, width = 0.02), 
                                  offset = 0.1, 
                                  pwidth = 0.2, 
                                  axis.params = list(
                                    axis = "x", 
                                    text.angle = 30, 
                                    text.size = 3, 
                                    vjust = 1, 
                                    hjust = 1, 
                                    line.size = 0
                                  ))  +
  scale_fill_manual(values=c("#A6BD9D","#FFB177", "black","black", "black"),
                    guide=guide_legend(keywidth = 0.3, 
                                       keyheight = 0.3, order=4)) 

P4_1 <- P3_1 + new_scale_fill() + 
  geom_fruit(data = bigscape_clusters_wrangle_clean, geom=geom_star, 
             mapping = aes(y = sampleID, x = gcf_identifier, starshape = SiderophoreClass, fill = SiderophoreClass), 
             offset = 0.05, 
             pwidth = 1, 
             grid.params = list(
               linetype=8,
               size=0.5, 
               color = "black", 
               lineend = "square", 
               alpha = 0.6
             ), 
             axis.params = list(
               axis = "x", 
               title = "", 
               text.angle = 0, 
               text.size = 3, 
               line.size = 1, 
               vjust = 1.5
             )) +
  scale_fill_manual(values=c("#A6BD9D","#FFB177", "black","black", "black"),
                    guide=guide_legend(keywidth = 0.1, 
                                       keyheight = 0.1, order=4)) + 
  theme(legend.position = c(0.9, 0.5)) 


P4_1 + 
  geom_strip('MMO-131', 'HN041', label = "A. pitti", offset = 0.24, offset.text = 0.008) + 
  geom_strip('AA010', 'AA010', label = "A. oleivorans", offset = 0.24, offset.text = 0.008) + 
  geom_strip('MMO-143', 'M51', label = "A. nosocomialis", offset = 0.24, offset.text = 0.008) + 
  geom_strip('HN127', 'HN047', label = "A. modestus", offset = 0.24, offset.text = 0.008) + 
  geom_strip('AS116', 'AS116', label = "A. junii", offset = 0.24, offset.text = 0.008) + 
  geom_strip('AS096', 'AS096', label = "A. beijerinckii", offset = 0.24, offset.text = 0.008) + 
  geom_strip('USHLN143', 'USHLN143', label = "A. fasciculus", offset = 0.24, offset.text = 0.008) + 
  geom_strip('HN020', 'HN020', label = "A. johnsonii", offset = 0.24, offset.text = 0.008) + 
  geom_strip('AS013', 'AS013', label = "A. tandoii", offset = 0.24, offset.text = 0.008) + 
  geom_strip('HN105', 'HN105', label = "A. bereziniae", offset = 0.24, offset.text = 0.008) + 
  geom_strip('AS167', 'AS167', label = "A. spp", offset = 0.24, offset.text = 0.008) +
  geom_strip('USMM068', 'HN121', label = "A. soli", offset = 0.24, offset.text = 0.008) + 
  theme(legend.position = c(2, 0.8)) 

ggsave(filename = "260824_Acinetobacter_BGCSubset_With_SiderophoreGCFs_with_legend.pdf", height = 6.5, width = 6)
ggsave(P3_1, filename = "260824_Acinetobacter_BGCSubset_With_SiderophoreGCFs_with_legend.pdf", height = 6.5, width = 10)

# quick exploratory analysis of these acinetobacters 
MosAIC_metadata_complete %>%
  filter(gtdb_genus == "Acinetobacter") %>%
  group_by(internal_id) %>%
  summarise(n()) # how many isolates

MosAIC_metadata_complete %>%
  filter(gtdb_genus == "Acinetobacter") %>%
  group_by(sample_source) %>%
  summarise(n())

MosAIC_MINUUR_antismash_simplified %>% 
  filter(genus == "Acinetobacter") %>% 
  filter(product == "NI-siderophore" | product == "NRP-metallophore") %>%
  group_by(sampleID, product) %>%
  summarise(n())
