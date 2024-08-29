library(tidyverse)

mosaic <- read_csv("S1_Table_MosAIC_Metadata_Final2 (3).csv")
#-----------------------------------------------------------------------------#

# plot Serratia and Acinetobacter
mosaic %>%
  filter(gtdb_genus == "Serratia" | gtdb_genus == "Acinetobacter") %>%
  group_by(sample_source, gtdb_genus) %>%
  summarise(number_isolates = n()) %>%
  ggplot() + 
  aes(x = sample_source, y = number_isolates) + 
  geom_bar(stat = "identity") + 
  theme_bw(base_size = 35) + 
  theme(axis.text.x = element_text(angle = 10, hjust = 1)) + 
  xlab("Isolation source") + 
  ylab("Number of isolates") + 
  facet_wrap(facets = "gtdb_genus", nrow = 3)
  
ggsave("Isolation_Sources.pdf", height = 10, width = 10)

#-----------------------------------------------------------------------------#
# Summaries
mosaic %>%
  filter(gtdb_genus == "Serratia" | gtdb_genus == "Acinetobacter") %>%
  group_by(sample_source, gtdb_genus) %>%
  summarise(number_isolates = n())
