# get subsets of target genera and overlay with BGC predicted product

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
MosAIC_tree <- ape::read.tree("021123_MosAIC_16s_aligned_trimmed.afa.bionj")
MosAIC_tree_root <- root(MosAIC_tree, "../02_16s_extracted_genes/GCF_006385615_rrna.fa_16s_rRNA.fa_16S_rRNA")

# Read antismash
MosAIC_antismash <- read_csv("081123_MosAIC_antismash.csv")

# Read MosAIC metadata 
MosAIC_metadata <- read_csv("Nov_S1_Table_MosAIC_Metadata_Final.csv") %>%
  filter(qc_pass == "Pass")

MosAIC_tree_tib <- as_tibble(MosAIC_tree_root) %>%
  separate(col = label, into = c("junk", "label2"), sep = "../02_16s_extracted_genes/", remove = F) %>%
  separate(col = label2, into = c("label2", "junk"), sep = "_contigs", remove = F) %>%
  select(label, label2) 

MosAIC_tree_rename <- rename_taxa(MosAIC_tree_root, MosAIC_tree_tib, label, label2)
MosAIC_tree_rename %>%
  print(n = 800)

### Function to subset tree to desired location 
subset_tree <- function(node, number_levels){
  
  subset <- tree_subset(MosAIC_tree_rename, node, levels_back = number_levels)
  
  return(subset)
}

Acinetobacter_subset <- midpoint_root(subset_tree("HN059", 19))
Acinetobacter_subset_tree <- ggtree(Acinetobacter_subset, layout = "circular", open.angle = 180, options(ignore.negative.edge=TRUE))

Serratia_subset <- midpoint_root(subset_tree("MMO-89", 4))
Serratia_subset_tree <- ggtree(Serratia_subset, layout = "circular") 

Enterobacter_subset <- subset_tree("AS212", 5)
Enterobacter_subset_tree <- ggtree(Enterobacter_subset, layout = "circular", options(ignore.negative.edge=TRUE)) + geom_tiplab(size = 2)

# summarise BGC products 
summarise_BGC_products <- function(target_genera, target_genera_2){
  products = as.data.frame(MosAIC_antismash %>%
                             filter(gtdb_genus == target_genera | gtdb_genus == target_genera_2) %>%
                             group_by(sampleID, product) %>%
                             summarise(number_BGCs = n()))
  
  return(products)
  
}

acinteobacter_products <- summarise_BGC_products("Acinetobacter", "Acinetobacter")
serratia_products <- summarise_BGC_products("Serratia", "Serratia_A")
enterobacter_products <- summarise_BGC_products("Enterobacter", "Enterobacter")

# wrangle species dataframe
MosAIC_species <- MosAIC_antismash %>%
  select(sampleID, gtdb_species) %>%
  distinct(sampleID, .keep_all = T) %>%
  tibble::column_to_rownames(var = "sampleID")


# make species subset trees
# subsetted trees with species metadata and BGC product prediction
plot_tree_subsets <- function(subset_tree, species_BGC_products){
  
  P1 <- subset_tree %>% gheatmap(MosAIC_species, colnames_offset_y = 20,
                                                colnames_angle = 60,
                                                hjust = 0.15,
                                                colnames_position = "top",
                                                colnames = T,
                                                width = 0.05,
                                                offset = 0.0, 
                                                font.size = 5) + 
    guides(fill = guide_legend(ncol = 2, title = "Source")) + 
    scale_fill_manual(values=c("#96761A","#228E3C","#2A4287",
                               "#776739", "#C2A1A1","#966AA5",
                               "#8492BA", "#C6DFBE", "#EFEADD", "#4E1BD6", 
                               "#ECC105", "#FF0000")) +
    new_scale_fill()
  
  P2 <- P1 + 
    geom_fruit(data=species_BGC_products, geom=geom_tile,
               mapping=aes(y=sampleID, x=product, alpha=number_BGCs, fill=product),
               color = "grey50", offset = 0.05, size = 0.03, pwidth = 0.5) +
    scale_alpha_continuous(range=c(0.5, 1),
                           guide=guide_legend(keywidth = 0.3, 
                                              keyheight = 0.3, order=1)) +
    scale_fill_manual(values=c("#96761A","#228E3C","#2A4287",
                               "#776739", "#C2A1A1","#966AA5",
                               "#8492BA", "#C6DFBE", "#EFEADD", "#4E1BD6", 
                               "#ECC105", "#FF0000"),
                      guide=guide_legend(keywidth = 0.3, 
                                         keyheight = 0.3, order=4))+
    geom_treescale(fontsize=1, linesize=0.3, x=0.01, y=0.1) +
    theme(legend.position=c(1.1, 0.5),
          legend.background=element_rect(fill=NA),
          legend.title=element_text(size=6.5),
          legend.text=element_text(size=4.5),
          legend.spacing.y = unit(0.02, "cm"),
    ) 
  
  return(P2)
  
}

plot_tree_subsets(Acinetobacter_subset_tree, acinteobacter_products)
ggsave("101123_acinetobacter_subset.pdf")

plot_tree_subsets(Serratia_subset_tree, serratia_products)
plot_tree_subsets(Enterobacter_subset_tree, enterobacter_products)


MosAIC_antismash2 <- MosAIC_antismash %>%
  select(sampleID, sample_source) %>%
  mutate(host = 1)

AS2 + geom_fruit(
  data = MosAIC_antismash2, 
  geom = geom_star, 
  mapping = aes(x = host, y = sampleID, starshape=sample_source),
  offset = -0.1, 
  grid.params=list(
    linetype=3,
    size=0.5
  )
) + 
  scale_starshape_manual(values = c("mosquito" = "square diamond", "mosquito_larval_water" = "circle", "non_mosquito_Diptera_environment" = "hexagonal star", 
                                    "non_mosquito_Diptera_host" = "regular pentagon", "Valiente Moro, Claire" = "antiparallelogram")) + 
  new_scale_fill()


S1 <- Serratia_subset_tree %>% gheatmap(acinetobacter_species, colnames_offset_y = 20,
                                              colnames_angle = 60,
                                              hjust = 0.15,
                                              colnames_position = "top",
                                              colnames = T,
                                              width = 0.05,
                                              offset = 0.0, 
                                              font.size = 5) + 
  guides(fill = guide_legend(ncol = 2, title = "Source")) + 
  new_scale_fill()

S2 <- S1 + geom_fruit(data=serratia_products, geom=geom_tile,
             mapping=aes(y=sampleID, x=product, alpha=number_BGCs, fill=product),
             color = "grey50", offset = 0.05, size = 0.03, pwidth = 0.5) +
  scale_alpha_continuous(range=c(0.5, 1),
                         guide=guide_legend(keywidth = 0.3, 
                                            keyheight = 0.3, order=1)) +
  scale_fill_manual(values=c("black","#996E3F","#E05040",
                             "#B6ED4F", "#61F7BB","#800080",
                             "#3A8A8C", "#DF8E90", "black", "black", 
                             "black", "black", "black", "black"),
                    guide=guide_legend(keywidth = 0.3, 
                                       keyheight = 0.3, order=4))+
  geom_treescale(fontsize=1, linesize=0.3, x=0.01, y=0.1) +
  theme(legend.position=c(1, 0.5),
        legend.background=element_rect(fill=NA),
        legend.title=element_text(size=6.5),
        legend.text=element_text(size=4.5),
        legend.spacing.y = unit(0.02, "cm"),
  ) 

S2 + geom_fruit(
  data = MosAIC_antismash2, 
  geom = geom_star, 
  mapping = aes(x = host, y = sampleID, starshape=sample_source),
  offset = -0.1, 
  grid.params=list(
    linetype=3,
    size=0.5
  )
) + 
  scale_starshape_manual(values = c("mosquito" = "square diamond", "mosquito_larval_water" = "circle", "non_mosquito_Diptera_environment" = "hexagonal star", 
                                    "non_mosquito_Diptera_host" = "regular pentagon", "Valiente Moro, Claire" = "antiparallelogram")) + 
  new_scale_fill()
