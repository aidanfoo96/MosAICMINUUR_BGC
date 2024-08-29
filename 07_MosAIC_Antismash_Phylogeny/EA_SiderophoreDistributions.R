# script to explore BGC across different isolation sources 
library(tidyverse)
library(janitor)

MosAIC_antismash <- read_csv("271123_MosAIC_antismash.csv")

# boxplot to show BGC diversity facetted by location 
MosAIC_antismash %>% 
  group_by(sampleID, category, sample_source...17) %>%
  summarise(number_BGCs = n()) %>%
  ggplot() + 
  aes(x = sample_source...17, y = number_BGCs) +
  geom_point() + 
  geom_boxplot() + 
  facet_wrap(facets = "category")

# general trend of siderophore and isolation source
MosAIC_antismash %>% 
  filter(product == "NRP-metallophore" | product == "NI-siderophore") %>%
  group_by(product, sample_source...17) %>%
  summarise(number_BGCs = n()) %>%
  ggplot() + 
  aes(x = sample_source...17, y = number_BGCs) +
  geom_bar(stat = "identity") + 
  facet_wrap(facets = "product")
  
# do bacteria have differences
MosAIC_antismash %>% 
  filter(product == "NRP-metallophore" | product == "NI-siderophore") %>%
  group_by(gtdb_genus, product, sample_source...17) %>%
  summarise(number_BGCs = n()) %>%
  ggplot() + 
  aes(x = sample_source...17, y = gtdb_genus, fill = number_BGCs) + 
  geom_tile() + 
  facet_wrap(facets = "product")
  

Serratia_table <- MosAIC_antismash %>%
  filter(gtdb_genus == "Serratia") %>%
  select(sampleID, contignumber, product)

write.table(x = Serratia_table, file = "Serratia_NatrualProduct_list.tsv", sep = "\t", quote = F)
  