# load packages
library(tidyverse)
library(janitor)

# read network file
OthersNetwork <- read_tsv("2024-02-06_20-17-55_hybrids_glocal/Others/Others_c0.30.network") %>%
  clean_names() 

NRPSNetwork <- read_tsv("060224_NRPS_c0.30.network") %>%
  clean_names()

# read classifications 
MosAIC_MINUUR_classifications <- read_tsv("190124_MINUUR_MosAIC_classifications.tsv") 

# function to clean and join the network file with MosAIC classifications 
clean_and_join_network <- function(network, classifications){
  
  network_clean <- network %>%
    separate(clustername_1, into = c("name_1", "file_ext"), sep = "_contigs", remove = F) %>%
    separate(name_1, into = c("name_1", "file_ext"), sep = ".fa_antismash", remove = F) %>%
    separate(clustername_2, into = c("name_2", "file_ext"), sep = "_contigs", remove = F) %>%
    separate(name_2, into = c("name_2", "file_ext"), sep = ".fa_antismash", remove = F) 
  
  joined_table <- network_clean %>%
    left_join(classifications, by = c("name_1" = "user_genome"))
  
  return(joined_table)
  
}

NRPSNetworkJoined <- clean_and_join_network(OthersNetwork, MosAIC_MINUUR_classifications)
write.table(x = NRPSNetworkJoined, file = "030424_OthersNetwork.tsv", quote = F, sep = "\t", row.names = F)
