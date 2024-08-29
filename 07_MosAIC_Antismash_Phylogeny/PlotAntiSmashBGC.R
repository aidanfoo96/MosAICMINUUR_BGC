# load packages 
library(tidyverse)

MosAIC_antismash <- read_csv("081123_MosAIC_antismash.csv")

# summarise how many BGCs are within a family 
BGC_summary_family <- MosAIC_antismash %>%
  select(sampleID, gtdb_family, gtdb_species) %>%
  group_by(sampleID, gtdb_family, gtdb_species) %>%
  summarise(number = n())

# Plot this in a boxplot 
BGC_summary_family %>%
  filter(gtdb_family != "NA") %>%
  ggplot() + 
  aes(x = gtdb_family, y = number) + 
  geom_boxplot() + 
  geom_jitter(width = 0.1, alpha = 0.2) + 
  coord_flip() + 
  theme_bw(base_size = 30) + 
  ylab("Number of BGCs") + 
  xlab("Family") + 
  ylim(0, 35)

ggsave(filename = "021123_BGC_BoxPlots_Family.pdf", height = 10, width = 10)

# summarise number of BGCs and genome sizes of each sample 
BGC_summary_family_GS <- MosAIC_antismash %>%
  select(sampleID, gtdb_family, gtdb_species, genome_size_bp) %>%
  group_by(sampleID, gtdb_family, gtdb_species, genome_size_bp) %>%
  summarise(number = n())

# relationship between genomesize and number of predicted BGCs
BGC_summary_family_GS %>%
  ggplot() + 
  aes(x = genome_size_bp, y = number) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  theme_bw(base_size = 25) + 
  xlab("Genome Size (bp)") + 
  ylab("Number of BGCs")

ggsave(filename = "021123_BGC_GenomeSize_Scatter.pdf", height = 10, width = 10)

MosAIC_antismash %>%
  filter(category == "NRPS") %>%
  select(sampleID, gtdb_family, gtdb_species) %>%
  group_by(sampleID, gtdb_family, gtdb_species) %>%
  summarise(number = n()) %>%
  ggplot() + 
  aes(x = gtdb_family, y = number) + 
  geom_boxplot() + 
  geom_jitter(width = 0.1, alpha = 0.2) + 
  coord_flip() + 
  theme_bw(base_size = 30) + 
  ylab("Number of BGCs") + 
  xlab("Family")

# function to filter for BGC category and make boxplot
filter_and_plt_box <- function(data, target){
  
  plot <- data %>%
    filter(gtdb_family != "NA") %>%
    filter(category == target) %>%
    select(sampleID, gtdb_family, gtdb_species) %>%
    group_by(sampleID, gtdb_family, gtdb_species) %>%
    summarise(number = n()) %>%
    ggplot() + 
    aes(x = gtdb_family, y = number) + 
    geom_boxplot() + 
    geom_jitter(width = 0.1, alpha = 0.2) + 
    coord_flip() + 
    theme_bw(base_size = 30) + 
    ylab("Number of BGCs") + 
    xlab("Family") + 
    ylim(0, 15)
  
  return(plot)
  
}

NRPS_boxplot <- filter_and_plt_box(MosAIC_antismash, "NRPS")  
ggsave(filename = "081123_BGC_NRPS_boxplot.pdf", height = 10, width = 10)

other_boxplot <- filter_and_plt_box(MosAIC_antismash, "other")  
ggsave(filename = "081123_BGC_other_boxplot.pdf", height = 10, width = 10)

terpene_boxplot <- filter_and_plt_box(MosAIC_antismash, "terpene")  
ggsave(filename = "081123_BGC_terpene_boxplot.pdf", height = 10, width = 10)

PKS_boxplot <- filter_and_plt_box(MosAIC_antismash, "PKS")  
ggsave(filename = "081123_BGC_PKS_boxplot.pdf", height = 10, width = 10)

RiPP_boxplot <- filter_and_plt_box(MosAIC_antismash, "RiPP")  
ggsave(filename = "081123_BGC_RiPP_boxplot.pdf", height = 10, width = 10)

MosAIC_antismash %>%
  select(sampleID, gtdb_family, gtdb_species, category) %>%
  group_by(sampleID, gtdb_family, gtdb_species, category) %>%
  summarise(number = n()) %>%
  filter(gtdb_family != "NA" & category != "NA" & category != "saccharide") %>%
  ggplot() + 
  aes(x = gtdb_family, y = number) + 
  geom_boxplot() + 
  geom_jitter(width = 0.1, alpha = 0.2) + 
  coord_flip() + 
  theme_bw(base_size = 30) + 
  ylab("Number of BGCs") + 
  xlab("Family") + 
  facet_wrap(~category, nrow = 1) + 
  theme(axis.text.x = element_text(angle = 50, hjust = 1))

ggsave(filename = "021123_BGC_boxplot_facet.pdf", height = 10, width = 11)

# differences in number of BGCs between different sources? 
MosAIC_antismash %>%
  filter(gtdb_family != "NA") %>%
  filter(category == "other") %>%
  select(sampleID, category, sample_source) %>%
  group_by(sampleID, sample_source, category) %>%
  summarise(number = n()) %>%
  ggplot() + 
  aes(x = sample_source, y = number) + 
  geom_boxplot() + 
  geom_jitter(width = 0.1, alpha = 0.2) + 
  coord_flip() + 
  theme_bw(base_size = 30) + 
  ylab("Number of BGCs") + 
  xlab("Family") + 
  facet_wrap(~category, nrow = 1) + 
  theme(axis.text.x = element_text(angle = 50, hjust = 1))

ggsave(filename = "021123_BGC_boxplot_sample_sourcet.pdf", height = 10, width = 11)
