# script to join and wrangle BigSCAPE clusters so you can make a nice 
# visualisation of siderophore BGC and their GCF distribution 

# load packages
library(tidyverse)

colnames_gcfs <- c("BGC", "family_number")

# load data 
other_gcfs <- read_tsv(file = "2024-02-01_17-18-50_hybrids_auto/Others/Others_clustering_c0.30.tsv",
                 col_names = colnames_gcfs, 
                 skip = 1)

NRPS_gcfs <- read_tsv(file = "2024-02-01_17-18-50_hybrids_auto/NRPS/NRPS_clustering_c0.30.tsv", 
                      col_names = colnames_gcfs, 
                      skip = 1)

NRPS_PKS_gcfs <- read_tsv(file = "2024-02-01_17-18-50_hybrids_auto/PKS-NRP_Hybrids/PKS-NRP_Hybrids_clustering_c0.30.tsv", 
                          col_names = colnames_gcfs, 
                          skip = 1)

PKSI_gcfs <- read_tsv("2024-02-01_17-18-50_hybrids_auto/PKSI/PKSI_clustering_c0.30.tsv", 
                      col_names = colnames_gcfs, 
                      skip = 1)

PKS_other_gcfs <- read_tsv("2024-02-01_17-18-50_hybrids_auto/PKSother/PKSother_clustering_c0.30.tsv", 
                           col_names = colnames_gcfs, 
                           skip = 1)

RiPPS_gcfs <- read_tsv("2024-02-01_17-18-50_hybrids_auto/RiPPs/RiPPs_clustering_c0.30.tsv", 
                       col_names = colnames_gcfs, 
                       skip = 1)

terpene_gcfs <- read_tsv("2024-02-01_17-18-50_hybrids_auto/Terpene/Terpene_clustering_c0.30.tsv", 
                         col_names = colnames_gcfs, 
                         skip = 1)

gcfs <- rbind(other_gcfs, NRPS_gcfs, NRPS_PKS_gcfs, PKSI_gcfs, PKS_other_gcfs, RiPPS_gcfs, terpene_gcfs)

bigscape_annotations <- read_tsv(file = "2024-02-01_17-18-50_hybrids_auto/Network_Annotations_Full.tsv")

classifications <- read_tsv("../01_MINUURMosAICIntegration/160524_MosAIC_MINUUR_CheckM_clean_classifications.tsv")

# join
antismash_hits_gcfs <- gcfs %>% 
  left_join(bigscape_annotations) %>%
  distinct() %>%
  filter(!startsWith(BGC, prefix = "BGC")) # exlude BigSCAPE BGCs 

antismash_hits_gcfs_MiBIG <- gcfs %>% 
  left_join(bigscape_annotations) %>%
  distinct()
  # include MiBIG BGCs 

antismash_hits_gcfs %>%
  filter(startsWith(BGC, prefix = "BGC")) # exlude BigSCAPE BGCs 

# cleean the ajoined data so you can join with classifications 
antismash_hits_gcfs_clean <- antismash_hits_gcfs %>%
  separate(BGC, into = c("sampleID", "junk"), sep = "_contigs") %>%
  separate(sampleID, into = c("sampleID", "junk2"), sep = ".fa_antismash") %>%
  select(!c(junk, junk2)) 

antismash_hits_gcfs_clean_MiBIG <- antismash_hits_gcfs_MiBIG %>%
  separate(BGC, into = c("sampleID", "junk"), sep = "_contigs") %>%
  separate(sampleID, into = c("sampleID", "junk2"), sep = ".fa_antismash") %>%
  select(!c(junk, junk2)) 
  
# join with classifications 
antismash_hits_gcfs_clean_classifications <- antismash_hits_gcfs_clean %>% 
  left_join(classifications, by = c("sampleID" = "user_genome"))

antismash_hits_gcfs_clean_classifications_MiBIG <- antismash_hits_gcfs_clean_MiBIG %>% 
  left_join(classifications, by = c("sampleID" = "user_genome"))

write.table(x = antismash_hits_gcfs_clean_classifications_MiBIG, 
            file = "200524_antismash_hits_gcfs_classifications_MiBIG.tsv",
            quote = F, sep = "\t", row.names = F)

# total number of unique products (GCFs)
antismash_hits_gcfs_clean_classifications %>%
  group_by(`Product Prediction`, family_number) %>%
  summarise(n()) %>%
  group_by(`Product Prediction`) %>%
  summarise(number_GCF_unique = n()) %>%
  print(n = 500) %>%
  filter(`Product Prediction` != "NA") %>%
  ggplot() + 
  aes(x = reorder(`Product Prediction`, number_GCF_unique), y = number_GCF_unique) + 
  geom_bar(stat = "identity") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

# total number of unique categories (GCFs)
format_labels <- function(x) {
  gsub("_", " ", x)  # Replace underscores with spaces
}


antismash_hits_gcfs_clean_classifications %>%
  group_by(`BiG-SCAPE class`, family_number) %>%
  summarise(n()) %>%
  group_by(`BiG-SCAPE class`) %>%
  summarise(number_GCF_unique = n()) %>%
  print(n = 500) %>%
  filter(`BiG-SCAPE class` != "NA") %>%
  ggplot() + 
  aes(x = reorder(`BiG-SCAPE class`, number_GCF_unique), y = number_GCF_unique) + 
  geom_bar(stat = "identity") + 
  theme_bw(base_size = 30) + 
  theme(axis.text.x = element_text(angle = 20, hjust = 1)) + 
  ylab("Number of GCFs (N)") + 
  xlab("BigSCAPE category") + 
  scale_x_discrete(labels = format_labels)   # Apply custom tick labels
  

ggsave("Figure_S9_TotalGCFs_per_category.pdf", height = 10, width = 10)

# group by 
antismash_hits_gcfs_clean_classifications %>%
  distinct(`Product Prediction`, family_number, family) %>%
  group_by(`Product Prediction`, family) %>%
  summarise(number_GCF_unique = n()) %>%
  print(n = 500) %>%
  filter(`Product Prediction` != "NA") %>%
  #filter(number_GCF_unique > 4) %>%
  ggplot() + 
  aes(x = reorder(`Product Prediction`, number_GCF_unique), y = number_GCF_unique, fill = family) + 
  geom_bar(stat = "identity", position = "stack") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 20, hjust = 1)) + 
  ylab("GCF Count (N)") + 
  xlab("BigSCAPE Product Prediction") + 
  scale_fill_manual(values = c("#005C66", "#BD5129", "#9FBAE7", "#1000D2", "#0070C8", "#52CEFF", 
                               "#009464", "#E1B698", "#06BF34", "#9D49AB", "#94E4FF", "#E3F8FF", 
                               "#A2581B", "#00B9C8", "#9AC7D9", "#50E4BC", "#6E2A4D", "#D582D2", 
                               "#F0AC7A", "#F7D7EF", "#176E93", "#7500FF", "#FFF4B1", "#52E685", 
                               "#FF9146", "#E0FFF1", "#A3A3FF", "#FFAAB7", "#FFE7CC", "black")) 

antismash_hits_gcfs_clean_classifications %>%
  group_by(family_number, `BiG-SCAPE class`, sampleID) %>%
  summarise(number_BGCs = n())

# calculate total number of GCFs per species 
total_GCFs_per_species <- antismash_hits_gcfs_clean_classifications %>%
  group_by(family_number, species) %>%
  summarise(number_GCF = n()) %>%
  group_by(family_number) %>%
  summarise(number_species = n()) %>%
  group_by(number_species) %>%
  summarise(number_different_species = n())
  

total_GCFs_per_species %>%
  mutate(percentage_within_species = number_different_species / sum(number_different_species) * 100)


# plot this 
total_GCFs_per_species %>% 
  ggplot() + 
  aes(x = number_species, y = number_different_species) + 
  geom_bar(stat = "identity") + 
  theme_bw(base_size = 30) + 
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 13, 23)) + 
  xlab("Number of different species") + 
  ylab("GCF count (N)") + 
  ylim(0, 800)

ggsave("Figure_3A_Total_GCFs_per_species.pdf", height = 10, width = 10)

# calculate total number of GCFs per genera 
total_GCFs_per_genus <- antismash_hits_gcfs_clean_classifications %>%
  group_by(family_number, genus) %>%
  summarise(number_GCF = n()) %>%
  group_by(family_number) %>%
  summarise(number_genus = n()) %>%
  group_by(number_genus) %>%
  summarise(number_different_genera = n())
  
# calculate percentages of each GCF 
total_GCFs_per_genus %>%
  mutate(percentage_within_genera = number_different_genera / sum(number_different_genera) * 100)

# plot number GCFs in each genera 
total_GCFs_per_genus %>% 
  ggplot() + 
  aes(x = number_genus, y = number_different_genera) + 
  geom_bar(stat = "identity") + 
  theme_bw(base_size = 30) + 
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 15)) + 
  xlab("Number of different genera") + 
  ylab("GCF count (N)")
  
ggsave("Figure_3B_Total_GCFs_per_genera.pdf", height = 10, width = 10)

# what products are predicted to be encoded by BGCs shared between genera? 
shared_GCF_hits_genera <- antismash_hits_gcfs_clean_classifications %>%
  group_by(family_number, genus, `Product Prediction`) %>%
  summarise(number_GCF = n()) %>%
  group_by(family_number, `Product Prediction`) %>%
  summarise(number_genus = n()) %>%
  filter(number_genus >= 2)  %>%
  print(n = 200)

shared_GCF_hits_genera_families <- as.vector(shared_GCF_hits_genera %>% select(family_number))

# plot genera with at least one hit to >=2 genera 
antismash_hits_gcfs_clean_classifications %>%
  filter(family_number %in% shared_GCF_hits_genera_families$family_number) %>% # filter using vector
  select(family_number, `Product Prediction`, genus, family) %>%
  ggplot() + 
  aes(x = genus, fill = `Product Prediction`) + 
  geom_histogram(stat = "count") + 
  facet_wrap(facets = "family_number", scales = "free", ncol = 4) + 
  ylab("Number of Genomes") + 
  xlab("Genus") +
  coord_flip() + 
  theme_bw() + 
  theme(legend.position = "bottom", 
        axis.text.y = element_text(face = "italic") ) + 
  guides(fill = guide_legend(ncol = 10)) + 
  ylim(0, 50)
  
ggsave("Fig_S9_GCFs_>2_per_genera.pdf", height = 14, width = 10)

# calculate total number of GCFs per family
total_GCFs_per_family <- antismash_hits_gcfs_clean_classifications %>%
  group_by(family_number, family) %>%
  summarise(number_GCF = n()) %>%
  group_by(family_number) %>%
  summarise(number_familys = n()) %>%
  group_by(number_familys) %>%
  summarise(number_different_families = n())

# calculate total number of GCFs per genome 
total_GCFs_per_genome <- antismash_hits_gcfs_clean_classifications %>%
  group_by(family_number, sampleID) %>%
  summarise(number_GCF = n()) %>%
  group_by(family_number) %>%
  summarise(number_genomes = n()) %>%
  group_by(number_genomes) %>%
  summarise(number_different_genome = n())
  
total_GCFs_per_genome %>%
  ggplot() + 
  aes(x = number_genomes, y = number_different_genome) + 
  geom_bar(stat = "identity") + 
  theme_bw(base_size = 30) + 
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65,70, 75)) + 
  xlab("Number of different genomes") + 
  ylab("GCF count (N)") + 
  ylim(0, 800)

ggsave("Figure_3C_Total_GCFs_per_genome.pdf", height = 10, width = 10)

# calculate how many had a match to the MiBIG repository 
antismash_hits_gcfs_clean_classifications_MiBIG %>%
  group_by(sampleID, family_number, species, `Product Prediction`) %>%
  summarise(number_of_GCF = n()) %>%
  group_by(family_number, )

antismash_hits_gcfs_clean_classifications_MiBIG_BGC_values <- antismash_hits_gcfs_clean_classifications_MiBIG %>%
  filter(str_detect(sampleID, "BGC")) %>%
  pull(family_number) # families where there was a match to the MiBIG 

# use above df to filter GCFs hits with a match to the MiBIG repository
antismash_hits_gcfs_filtered_MiBIG <- antismash_hits_gcfs_clean_classifications_MiBIG %>%
  filter(family_number %in% antismash_hits_gcfs_clean_classifications_MiBIG_BGC_values)

# identify family number values and the BGC sampleIDs 
bgc_family_numbers <- antismash_hits_gcfs_filtered_MiBIG %>%
  filter(startsWith(sampleID, prefix = "BGC")) %>% # filter BigSCAPE BGCs
  select(family_number, bgc_sampleID = sampleID) %>%
  distinct()

# join the antismash_hits_gcfs_filtered_MiBIG df to create the new column with the matching BGC in MiBIG 
antismash_hits_gcfs_filtered_MiBIG_joined <- antismash_hits_gcfs_filtered_MiBIG %>%
  left_join(bgc_family_numbers, by = "family_number") %>%
  filter(!startsWith(sampleID, prefix = "BGC")) # exlude BigSCAPE BGCs 

# how many GCFs with a match
antismash_hits_gcfs_filtered_MiBIG_joined %>%
  group_by(family_number) %>%
  summarise(n())

# percentage GCFs with a hit to a BGC in MiBIG
38 / 864 * 100

# plot!
antismash_hits_gcfs_filtered_MiBIG_joined %>%
  ggplot() + 
  aes(x = species, y = bgc_sampleID) + 
  geom_tile() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

fill_colours <- c(
  "NRPS" = "#E8DA85",
  "Others" = "#A8A8A8",
  "PKS-NRP_Hybrids"= "#FF8A87",
  "PKSI" = "#A68D9D",
  "PKSother" = "#C0D291",
  "RiPPS" = "#2EE651",
  "Saccharides" = "#FFB177",
  "Terpene" = "#E7D4FF"
)


antismash_hits_gcfs_filtered_MiBIG_joined %>%
  group_by(bgc_sampleID, family, `BiG-SCAPE class`) %>%
  summarise(number_GCFs = n()) %>%
  ggplot() + 
  aes(x = reorder(bgc_sampleID,desc(number_GCFs)), y = number_GCFs, fill = `BiG-SCAPE class`) + 
  geom_bar(stat = "identity", position = "stack") + 
  #geom_jitter(width = 0.1, length = 0.1) + 
  theme_bw(base_size = 30) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 18), 
        legend.position = c(0.89, 0.82), 
        legend.box = "vertical",
        legend.background = element_rect(fill = "white", colour = "black", size = 0.5), 
        legend.title = element_blank()) + 
  scale_fill_manual(values=fill_colours) + 
  xlab("MiBIG Biosynthetic gene cluster") + 
  ylab("BGC count (N)") 

ggsave("Figure_3D_BGC_match_to_mibig_repository.pdf", height = 9, width = 20)

# filter to get siderophores
siderophore_BGCs <- antismash_hits_gcfs_clean_classifications %>%
  filter(str_detect(`Product Prediction`, "siderophore") | str_detect(`Product Prediction`, "NRP-metallophore")) %>% # filter BigSCAPE BGCs
  mutate(SiderophoreClass = case_when(
    str_detect(`Product Prediction`, "siderophore") ~ "NI-siderophore",
    str_detect(`Product Prediction`, "metallophore") ~ "NRP-metallophore",
    TRUE ~ "Other"
  ))

## how many siderophore BGCs are there in each category? 
siderophore_BGCs %>%
  group_by(SiderophoreClass) %>%
  summarise(n())

## How many siderophore GCFs in each category 
siderophore_BGCs %>%
  group_by(SiderophoreClass, family_number) %>%
  summarise(number_BGC_per_GCF = n()) %>%
  group_by(SiderophoreClass) %>%
  summarise(number_GCF_per_class = n())

# gene cluster families distributed across bacterial order
siderophore_BGCs %>%
  filter(order != "NA") %>%
  distinct(SiderophoreClass, family_number, order) %>%
  group_by(order, SiderophoreClass) %>%
  summarise(number_GCFs = n()) %>%
  ggplot() + 
  aes(x = reorder(SiderophoreClass, number_GCFs), y = number_GCFs, fill = order) + 
  scale_fill_manual(values=c("#B999CC","#1D9A6C","#002ACD","#FFC60B",
                               "#B64521", "#9ACD32","#D15FEE","#FFC0CB",
                               "#EE6A50","#8DEEEE", "#006400","#800000",
                               "#B0171F","#191970")) + 
  geom_bar(stat = "identity", position = "dodge", colour = "black", linewidth = 0.3) + 
  theme_bw(base_size = 30) + 
  xlab("Siderophore class") + 
  ylab("GCF count (N)") + 
  theme(legend.position = "right") + 
  guides(fill = "none")

ggsave("Figure_3_E_siderophore_GCF_distribution_to_order_no_legend.pdf", height = 10, width = 10)
ggsave("Figure_3_E_siderophore_GCF_distribution_to_order_legend.pdf", height = 10, width = 20)


# total number of gene cluster families 
siderophore_BGCs %>%
  filter(family != "NA") %>%
  distinct(SiderophoreClass, family_number, family) %>%
  group_by(family, SiderophoreClass) %>%
  summarise(unique_number_GCFs = n()) %>%
  ggplot() + 
  aes(x = reorder(SiderophoreClass, -unique_number_GCFs), y = unique_number_GCFs, fill = family) + 
  scale_fill_manual(values = c("#005C66", "#BD5129", "#9FBAE7", "#1000D2", "#0070C8", "#52CEFF", 
                               "#009464", "#E1B698", "#06BF34", "#9D49AB", "#94E4FF", "#E3F8FF", 
                               "#A2581B", "#00B9C8", "#9AC7D9", "#50E4BC", "#6E2A4D", "#D582D2", 
                               "#F0AC7A", "#F7D7EF", "#176E93", "#7500FF", "#FFF4B1", "#52E685", 
                               "#FF9146", "#E0FFF1", "#A3A3FF", "#FFAAB7", "#FFE7CC")) + 
  geom_bar(stat = "identity", position = position_dodge(preserve = "single"), colour = "black", linewidth = 0.3, width = 0.6) + 
  theme_bw(base_size = 30) + 
  xlab("Siderophore class") + 
  ylab("GCF count (N)") 

ggsave("180524_siderophore_GCF_distribution_to_family.pdf", height = 10, width = 10)

# total number of BGCs 
siderophore_BGCs %>%
  filter(family != "NA") %>%
  distinct(SiderophoreClass, family_number, family, sampleID) %>%
  group_by(family, SiderophoreClass) %>%
  summarise(total_BGCs = n()) %>%
  ggplot() + 
  aes(x = reorder(SiderophoreClass, -total_BGCs), y = total_BGCs, fill = family) + 
  scale_fill_manual(values = c("#005C66", "#BD5129", "#9FBAE7", "#1000D2", "#0070C8", "#52CEFF", 
                               "#009464", "#E1B698", "#06BF34", "#9D49AB", "#94E4FF", "#E3F8FF", 
                               "#A2581B", "#00B9C8", "#9AC7D9", "#50E4BC", "#6E2A4D", "#D582D2", 
                               "#F0AC7A", "#F7D7EF", "#176E93", "#7500FF", "#FFF4B1", "#52E685", 
                               "#FF9146", "#E0FFF1", "#A3A3FF", "#FFAAB7", "#FFE7CC")) + 
  geom_bar(stat = "identity", position = position_dodge(preserve = "single"), colour = "black", linewidth = 0.3, width = 0.6) + 
  theme_bw(base_size = 30) + 
  xlab("Siderophore class") + 
  ylab("BGC count (N)") + 
  guides(fill = "none")

ggsave("180524_siderophore_BGC_distribution_to_order.pdf", height = 10, width = 10)

# calculate total number of GCFs per genera 
# overall
siderophore_BGCs %>%
  group_by(family_number, genus) %>%
  summarise(number_GCF = n()) %>%
  group_by(family_number) %>%
  summarise(number_genus = n()) %>%
  group_by(number_genus) %>%
  summarise(number_different_genera = n()) %>%
  mutate(percentage = number_different_genera / sum(number_different_genera) * 100)

# split between the two classess
siderophores_per_genus <- siderophore_BGCs %>%
  group_by(family_number, genus, SiderophoreClass) %>%
  summarise(number_GCF = n()) %>%
  group_by(family_number, SiderophoreClass) %>%
  summarise(number_genus = n()) %>%
  group_by(number_genus, SiderophoreClass) %>%
  summarise(number_different_genera = n())

siderophores_per_genus %>%
  ggplot() + 
  aes(x = number_genus, y = number_different_genera) + 
  geom_bar(stat = "identity") + 
  theme_bw(base_size = 30) + 
  scale_x_continuous(breaks = c(1, 2, 3)) + 
  xlab("Number of different genera") + 
  ylab("GCF count (N)") + 
  facet_wrap(facets = "SiderophoreClass") + 
  theme(legend.position = c(1.85, 1.2),  # Specify x, y coordinates of the legend
        legend.justification = c(1, 1),  # Specify justification (1,1) = top right
        legend.margin = margin(t = -10, r = -10, b = -10, l = -10), 
        strip.background = element_rect(fill = "white", colour = NA))

ggsave("Figure_3F_Total_Siderophore_GCFs_per_genera.pdf", height = 10, width = 9)

# wrangle and plot siderophore classes > 2
siderophore_classes_over_two <- as.vector(siderophore_BGCs %>%
  group_by(family_number, genus, SiderophoreClass) %>%
  summarise(number_GCF = n()) %>%
  group_by(family_number, SiderophoreClass) %>%
  summarise(number_genus = n()) %>%
  filter(number_genus >= 2))

siderophore_BGCs %>%
  filter(family_number %in% siderophore_classes_over_two$family_number) %>% # filter using vector
  select(family_number, SiderophoreClass, genus, family) %>%
  ggplot() + 
  aes(x = genus, fill = SiderophoreClass) + 
  geom_histogram(stat = "count") + 
  facet_wrap(facets = "family_number", scales = "free", ncol = 2) + 
  ylab("Number of Genomes") + 
  xlab("Genus") +
  coord_flip() + 
  theme_bw(base_size = 20) + 
  theme(legend.position = "bottom", 
        axis.text.y = element_text(face = "italic") ) + 
  guides(fill = guide_legend(ncol = 10)) + 
  ylim(0, 30) + 
  scale_fill_manual(values=c("#A6BD9D","#FFB177")) 

ggsave("Fig_S11_SiderophoreCounts_split_over_two_genera.pdf", height = 10, width = 10)


# species
siderophore_BGCs %>%
  group_by(family_number, species) %>%
  summarise(number_GCF = n()) %>%
  group_by(family_number) %>%
  summarise(number_species = n()) %>%
  group_by(number_species) %>%
  summarise(number_different_species = n()) %>%
  mutate(percentage = number_different_species / sum(number_different_species) * 100)

siderophores_per_species <- siderophore_BGCs %>%
  group_by(family_number, species, SiderophoreClass) %>%
  summarise(number_GCF = n()) %>%
  group_by(family_number, SiderophoreClass) %>%
  summarise(number_species = n()) %>%
  group_by(number_species, SiderophoreClass) %>%
  summarise(number_different_species = n()) %>%
  mutate(percentage = number_different_species / sum(number_different_species) * 100)

siderophores_per_species %>%
  ggplot() + 
  aes(x = number_species, y = number_different_species) + 
  geom_bar(stat = "identity") + 
  theme_bw(base_size = 30) + 
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 10, 15, 20, 25, 30, 35, 40, 45, 50)) + 
  xlab("Number of different species") + 
  ylab("GCF count (N)") + 
  facet_wrap(facets = "SiderophoreClass") + 
  theme(legend.position = c(1.85, 1.2),  # Specify x, y coordinates of the legend
        legend.justification = c(1, 1),  # Specify justification (1,1) = top right
        legend.margin = margin(t = -10, r = -10, b = -10, l = -10), 
        strip.background = element_rect(fill = "white", colour = NA))

ggsave("Figure_3_G_Total_Siderophore_GCFs_per_species.pdf", height = 10, width = 9)


siderophore_classes_over_two_sepcies <- as.vector(siderophore_BGCs %>%
                                            group_by(family_number, species, SiderophoreClass) %>%
                                            summarise(number_GCF = n()) %>%
                                            group_by(family_number, SiderophoreClass) %>%
                                            summarise(number_species = n()) %>%
                                            filter(number_species >= 3))

siderophore_BGCs %>%
  filter(family_number %in% siderophore_classes_over_two_sepcies$family_number) %>% # filter using vector
  select(family_number, SiderophoreClass, species, family) %>%
  ggplot() + 
  aes(x = species, fill = SiderophoreClass) + 
  geom_histogram(stat = "count") + 
  facet_wrap(facets = "family_number", scales = "free", ncol = 2) + 
  ylab("Number of Genomes") + 
  xlab("Genus") +
  coord_flip() + 
  theme_bw() + 
  theme(legend.position = "bottom", 
        axis.text.y = element_text(face = "italic") ) + 
  guides(fill = guide_legend(ncol = 10)) + 
  ylim(0, 30) + 
  scale_fill_manual(values=c("#A6BD9D","#FFB177")) + 
  scale_x_discrete(labels = format_labels)   # Apply custom tick labels


ggsave("Fig_S12_SiderophoreCounts_split_over_two_species.pdf", height = 12, width = 10)

# how many siderophore BGCs with a match to the MiBIG
siderophore_BGCs_with_MiBIG_hit <- antismash_hits_gcfs_filtered_MiBIG_joined %>%
  filter(str_detect(`Product Prediction`, "siderophore") | str_detect(`Product Prediction`, "NRP-metallophore")) %>% # filter BigSCAPE BGCs
  mutate(SiderophoreClass = case_when(
    str_detect(`Product Prediction`, "siderophore") ~ "NI-Siderophore",
    str_detect(`Product Prediction`, "metallophore") ~ "NRP-Siderophore",
    TRUE ~ "Other"
  ))

siderophore_BGCs_with_MiBIG_hit %>%
  group_by(bgc_sampleID, family_number) %>%
  summarise(n())
# percentage from all other MiBIG hits

15 / 38 * 100

# figure with hits to MiBIG
siderophore_BGCs_with_MiBIG_hit %>%
  ggplot() + 
  aes(x = species, fill = SiderophoreClass) + 
  geom_histogram(stat = "count") + 
  facet_wrap(facets = "bgc_sampleID", scales = "free", ncol = 2) + 
  ylab("Number of Genomes") + 
  xlab("Species") +
  coord_flip() + 
  theme_bw() + 
  theme(legend.position = "bottom", 
        axis.text.y = element_text(face = "italic") ) + 
  guides(fill = guide_legend(ncol = 10)) + 
  ylim(0, 20) + 
  scale_fill_manual(values=c("#A6BD9D","#FFB177")) 

ggsave("Fig_S10_SiderophoreCounts_with_MIBIG_representative.pdf", height = 14, width = 10)

# how many specific to acinetobacter? 
siderophore_BGCs_with_MiBIG_hit %>%
  filter(genus == "Acinetobacter") %>%
  group_by(bgc_sampleID, sampleID, species, SiderophoreClass) %>%
  summarise(n())

# BGC0000294 = acinetobactin
# BGC0000295 = acinetoferrin

## statistics 
### How many GCFs are there in total? 
antismash_hits_gcfs_clean_classifications %>%
  group_by(family_number) %>%
  summarise(n())
  
