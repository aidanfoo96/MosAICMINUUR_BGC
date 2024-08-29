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
MosAIC_tree <- ape::read.tree("080524_gtdbtk.bac120.user_msa.fasta.treefile")

# outgroup
MosAIC_tree_root <- root(MosAIC_tree, "GCF_006385615.1_ASM638561v1_genomic")

# Read antismash
MosAIC_MINUUR_antismash <- read_tsv("../05_MosAIC_BGC/170524_MINUUR_MosAIC_classifications_antismash_category_counts.tsv")

# Read MosAIC metadata 
MosAIC_MINUUR_metadata <- read_tsv("190124_MINUUR_MosAIC_classifications.tsv") 

# Read MosAIC MINUUR Metadata 

MosAIC_tree_tib <- as_tibble(MosAIC_tree_root) %>%
  #separate(col = label, into = c("junk", "label2"), sep = "../02_16s_extracted_genes/", remove = F) %>%
  separate(col = label, into = c("label2", "junk"), sep = "_contigs", remove = F) %>%
  select(label, label2) 

MosAIC_tree_rename <- rename_taxa(MosAIC_tree_root, MosAIC_tree_tib, label, label2)
#MosAIC_tree_rename %>%
  #print(n = 800)

MosAIC_ggtree <- ggtree(MosAIC_tree_rename, layout="fan", size=0.2, open.angle = 90) #+ geom_tiplab(size = 0.4) + geom_nodelab(size = 0.4)
#MosAIC_ggtree + geom_hilight(data = node_highlight, mapping = aes(node=node),
                               #extendto=1, alpha=0.3, fill="grey", color="grey50",
                               #size=0.05)

ggsave(filename = "170524_tree_with_labels.pdf")

#geom_hilight(data=nodedf, mapping=aes(node=node),
             #extendto=6.8, alpha=0.3, fill="grey", color="grey50",
             #size=0.05) 

# function to extract column and convert for gheatmap function
extract_column_and_convert <- function(data, column){
  target <- data[, c(1, column)]
  target2 <- target %>%
    remove_rownames %>%
    column_to_rownames(var = "user_genome")
  
  return(target2)
}

MosAIC_family <- extract_column_and_convert(MosAIC_MINUUR_metadata, 5)

MosAIC_antismash_summarised = as.data.frame(MosAIC_MINUUR_antismash %>%
  select(sampleID, `BiG-SCAPE class`) %>%
  group_by(sampleID, `BiG-SCAPE class`) %>%
  summarise(number = n()))

MosAIC_MINUUR_antismash_genome_type = as.data.frame(MosAIC_MINUUR_antismash %>%
                                                      select(sampleID, genome))

MosAIC_MINUUR_antismash_simplified <- MosAIC_antismash_summarised %>%
  mutate(product = if_else(product == "NRP-metallophore / NRPS", "NRP-metallophore", product)) %>%
  mutate(product = if_else(product == "NRP-metallophore / NRPS / redox-cofactor", "NRP-metallophore", product)) %>%
  mutate(product = if_else(product == "NRPS / hserlactone", "hserlactone", product)) %>%
  mutate(product = if_else(product == "NRPS / hserlactone / arylpolyene", "hserlactone", product)) %>%
  mutate(product = if_else(product == "arylpolyene / NRPS / hserlactone", "arylpolyene", product)) %>%
  mutate(product = if_else(product == "redox-cofactor / NRP-metallophore / NRPS", 'NRP-metallophore', product)) 

# add extra metadata where no antismash hits were found 
names <- MosAIC_family %>% 
  filter(family == "Anaplasmataceae") %>%
  mutate(genome_type = "MAG") %>%
  select(user_genome, genome_type) %>%
  rename(sampleID = "user_genome")

MosAIC_MINUUR_antismash_genome_type_added <- rbind(MosAIC_MINUUR_antismash_genome_type, names)

MosAIC_antismash_total_BGC = as.data.frame(MosAIC_MINUUR_antismash %>%
  select(sampleID, category) %>%
  group_by(sampleID) %>%
  summarise(number = n()))

MosAIC_order = as.data.frame(MosAIC_MINUUR_metadata %>%
  select(user_genome, order))

MosAIC_source = as.data.frame(MosAIC_metadata %>%
                                select(internal_id, sample_source))


MosAIC_ggtree2 <- MosAIC_ggtree %<+% MosAIC_order + geom_point2(
  mapping=aes(col=order),
  position="identity", 
  size = 0.3) +
  scale_colour_manual(values=c("#005C66", "#BD5129", "#9FBAE7", "#1000D2", "#0070C8", "#52CEFF", 
                               "#009464", "#E1B698", "#06BF34", "#9D49AB", "#94E4FF", "#E3F8FF", 
                               "#A2581B", "#00B9C8", "#9AC7D9", "#50E4BC", "#6E2A4D", "#D582D2", 
                               "#F0AC7A", "#F7D7EF", "#97ADFF", "#7500FF", "#FFF4B1"),
                      guide=guide_legend(keywidth = 0.5, 
                                         keyheight = 0.5, order=1,
                                         override.aes=list(starshape=15)),
                      na.translate=FALSE) +
  scale_starshape_manual(values=c(15, 1),
                         guide=guide_legend(keywidth = 0.5, 
                                            keyheight = 0.5, order=2),
                         na.translate=FALSE) +
  scale_size_continuous(range = c(0.1, 5),
                        guide = guide_legend(keywidth = 0.5, 
                                             keyheight = 0.5, order=3,
                                             override.aes=list(starshape=15)))

MosAIC_ggtree_3 <- MosAIC_ggtree2 + new_scale_fill() +
  geom_fruit(data=MosAIC_antismash_summarised, geom=geom_tile,
             mapping=aes(y=sampleID, x=`BiG-SCAPE class`, alpha=number, fill=`BiG-SCAPE class`),
             color = "grey50", offset = 0.05, size = 0.03, pwidth = 0.5) +
  scale_alpha_continuous(breaks = c(2, 4, 6, 8, 10, 12),
                         guide=guide_legend(keywidth = 0.3, 
                                            keyheight = 0.3)) +
  scale_fill_manual(values=c("#E8DA85","#A8A8A8","#FF8A87",
                             "#A6BD9D", "#C0D291","#2EE651", 
                             "#FFB177", "#FF98C2", "#E7D4FF"),
                    guide=guide_legend(keywidth = 0.3, 
                                       keyheight = 0.3, order=4)) + 
  geom_fruit(data=MosAIC_antismash_total_BGC, geom=geom_bar,
             mapping=aes(y=sampleID, x=number),
             pwidth=0.5,
             offset = 0.08, 
             orientation="y", 
             stat="identity",
             axis.params = list(
               axis = "x", 
               text.size = 1,
               hjust = 0.5, 
               vjust = 1.5, 
             ),
             grid.params = list(size = 0.1)
  ) 

MosAIC_ggtree_4 <- MosAIC_ggtree_3 + new_scale_fill() + geom_fruit(data=MosAIC_MINUUR_antismash_genome_type_added, geom=geom_tile,
           mapping=aes(y=sampleID, fill=genome_type),
           color = "grey50", offset = 0.08, size = 0.03, pwidth = 0.1) + 
  scale_fill_manual(values=c("#188977","orange"),
                    guide=guide_legend(keywidth = 0.3, 
                                       keyheight = 0.3, order=4)) + 
  geom_treescale(fontsize=1, linesize=0.3, x=0.5, y=0.1) 

MosAIC_ggtree_4

ggsave(plot = MosAIC_ggtree_4, filename = "Figure_2A_MosAICMINUUR_antiSMASH_tree.pdf")


  