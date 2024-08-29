# script to investigate BGC diversity in Enterobacter, acinetobacter and serratia 

# packages
library(tidyverse)

# data 
MosAIC_antismash <- read_csv("271123_MosAIC_antismash.csv")

MosAIC_antismash_filt <- MosAIC_antismash %>%
  filter(gtdb_genus == "Enterobacter" | 
           gtdb_genus == "Acinetobacter" | gtdb_genus == "Serratia" | 
           gtdb_genus == "Serratia_A" |
           gtdb_genus == "Pseudomonas_F" |
           gtdb_genus == "Pseudomonas_E")

# how many BGCs in each species? 
MosAIC_antismash_filt %>%
  group_by(internal_id, gtdb_species) %>%
  summarise(number_BGCs = n()) %>%
  ggplot() + 
  aes(x = gtdb_species, y = number_BGCs) + 
  geom_boxplot() + 
  coord_flip() + 
  theme_bw(base_size = 25) + 
  geom_jitter(width = 0.1, alpha = 0.1) + 
  ylim(0, 15) + 
  xlab("Species") + 
  ylab("Number of BGCs")

# how many BGCs of each type, facet by species
MosAIC_antismash_filt %>%
  group_by(internal_id, gtdb_species, product, gtdb_genus) %>%
  summarise(number_BGCs = n()) %>%
  ggplot() + 
  aes(x = product, y = number_BGCs) + 
  geom_boxplot() + 
  coord_flip() + 
  theme_bw(base_size = 25) + 
  geom_jitter(width = 0.1, alpha = 0.1) + 
  facet_wrap(~gtdb_genus, nrow = 1) + 
  ylim(0, 15) + 
  xlab("Species") + 
  ylab("Number of BGCs")  

# how many BGCs of each type in a species 
MosAIC_antismash %>%
  filter(gtdb_genus == "Pseudomonas_F" | gtdb_genus == "Pseudomonas_E") %>%
  group_by(internal_id, gtdb_species, product) %>%
  summarise(number_BGCs = n()) %>%
  ggplot() + 
  aes(x = product, y = internal_id, fill = number_BGCs) + 
  geom_tile() +
  theme_bw(base_size = 15) + 
  theme(axis.text.x = element_text(angle = 15, hjust = 1)) + 
  xlab("BGC Product") + 
  ylab("Species")  
  
