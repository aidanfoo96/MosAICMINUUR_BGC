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
library(phylogram)
library(spiralize)
library(gplots)
library(heatmaply)
library(dendextend)
library(gtools)

# Read Tree
Acinetobacter_tree <- ape::read.tree("gtdbtk.bac120.user_msa.fasta_2.treefile")
Acinetobacter_tree_root <- root(Acinetobacter_tree, "VD15_contigs")

# read inhibition scores 
Acinetobacer_inhibition_scores <- read_tsv("140524_AcinetobacterInhibitionScoresSummarised.tsv")

# rename to match tree tip labels
Acinetobacer_inhibition_scores_renamed <- Acinetobacer_inhibition_scores %>%
  mutate(target_strain_right = if_else(target_strain_right == "MMO128", "MMO-128", target_strain_right)) %>%
  mutate(competing_strain_left = if_else(competing_strain_left == "MMO128", "MMO-128", competing_strain_left)) %>%
  mutate(target_strain_right = if_else(target_strain_right == "MMO163", "MMO-163", target_strain_right)) %>%
  mutate(competing_strain_left = if_else(competing_strain_left == "MMO163", "MMO-163", competing_strain_left)) %>%
  mutate(target_strain_right = if_else(target_strain_right == "USMMO68", "USMM068", target_strain_right)) %>%
  mutate(competing_strain_left = if_else(competing_strain_left == "USMMO68", "USMM068", competing_strain_left))

  
# collpase the tree using PopPUNK defined lineages
Acinetobacter_tree_root_collapsed <- Acinetobacter_tree_root %>% 
  drop.tip(tip = "HN041_contigs") %>%
  drop.tip(tip = "HN067_contigs") %>%
  drop.tip(tip = "HN068_contigs") %>%
  drop.tip(tip = "HN068_contigs") %>%
  drop.tip(tip = "HN071_contigs") %>%
  drop.tip(tip = "HN093_contigs") %>%
  drop.tip(tip = "HN099_contigs") %>%
  drop.tip(tip = "HN122_contigs") %>%
  drop.tip(tip = "HN123_contigs") %>%
  drop.tip(tip = "HN126_contigs") %>%
  drop.tip(tip = "VD14_contigs") %>%
  drop.tip(tip = "USMM066_contigs") %>%
  drop.tip(tip = "MMO-130_contigs") %>%
  drop.tip(tip = "MMO-131_contigs") %>%
  drop.tip(tip = "MMO-143_contigs") 

# Rooted but drop the root for the display
Acinetobacter_tree_root_collapsed <- drop.tip(Acinetobacter_tree_root_collapsed, tip = "VD15_contigs")

# rename tip labels
Acinetobacter_tree_tib <- as_tibble(Acinetobacter_tree_root_collapsed) %>%
  #separate(col = label, into = c("junk", "label2"), sep = "../02_16s_extracted_genes/", remove = F) %>%
  separate(col = label, into = c("label2", "junk"), sep = "_contigs", remove = F) %>%
  select(label, label2) 


Acinetobacter_tree_rename <- rename_taxa(Acinetobacter_tree_root_collapsed, Acinetobacter_tree_tib, label, label2)
#MosAIC_tree_rename %>%
#print(n = 800)

# convert to a dendrogram object
Acinetobacter_dendo <- chronos(Acinetobacter_tree_rename)

# tree with collapsed tips 
Acinetobacter_ggtree <- ggtree(Acinetobacter_tree_rename, layout="rectangular", size=0.2) + 
  xlim_tree(xlim = 0.5) + 
  geom_treescale(x = 0.001) + 
  ggtree::vexpand(.1, -1) + 
  geom_tiplab()

ggsave("Acinetobacter_Collapsed_Tree.pdf")

# function to conver to appropirate format and filter condition 

convert_and_filter_matrix <- function(target_condition){
  
  Acinetobacer_inhibition_scores_WideSummary <- Acinetobacer_inhibition_scores_renamed %>% 
    filter(condition == target_condition) %>%
    filter(inhibition_score != "NA") %>%
    select(target_strain_right, competing_strain_left, inhibition_score) %>%
    group_by(target_strain_right, competing_strain_left) %>%
    summarise(mean_inhibition_score = mean(inhibition_score)) %>%
    pivot_wider(names_from = target_strain_right, values_from = mean_inhibition_score) %>%
    column_to_rownames(var = "competing_strain_left")
  
  return(Acinetobacer_inhibition_scores_WideSummary)
  
}

AcinetobacterInhibitionBip <- convert_and_filter_matrix("2,2′-Bipyridyl")
AcinetobacterInhibitionFe <- convert_and_filter_matrix("FeCl3")


gheatmap(Acinetobacter_ggtree, Acinetobacer_inhibition_scores_WideSummary,
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
                       title = "Inhibition Score")) 

AcinetobacterInhibitionBipMatrix <- as.matrix(AcinetobacterInhibitionBip)
AcinetobacterInhibitionFeMatrix <- as.matrix(AcinetobacterInhibitionFe)

AcinetobacterInhibitionBipMatrix <- AcinetobacterInhibitionBipMatrix[Acinetobacter_tree_rename$tip.label, Acinetobacter_tree_rename$tip.label]
AcinetobacterInhibitionFeMatrix <- AcinetobacterInhibitionFeMatrix[Acinetobacter_tree_rename$tip.label, Acinetobacter_tree_rename$tip.label]

dendrogram <- as.dendrogram.phylo(Acinetobacter_tree_rename)
conv <- chronos(Acinetobacter_tree_rename)
hc <- as.hclust.phylo(conv)
dendro <- as.dendrogram(hc) 

dendro <- dendro %>%
  set("branches_lwd", 0.3) 


hmcols <- rev(redgreen(2750));

colours <- c("#ffd166", "#52796f")


heatmaply(AcinetobacterInhibitionFeMatrix, 
                 Rowv = dendro, 
                 Colv = dendro, 
                 linewidth = 0.1, 
                 na.value = "black", 
                 grid_size = 0.1, 
                 col = viridis(100),
                 #subplot_widths = c(0.2, 0.8), 
                 #subplot_heights = c(0.3, 0.5), 
                 scale_fill_gradient_fun = scale_fill_gradient(low = "dark green", high = "dark orange", limits = c(0, 100)), 
                 margins = c(40, 130), 
                 showticklabels = c(T, T), 
                 show_dendrogram = c(T, T),
          branches_lwd = 0.2, 
          plot_method = "ggplot", 
          margins(c(50, 50)), 
          grid_gap = 0.2, 
          heatmap_layers = list(theme_minimal(base_size = 15), 
                                theme(axis.text.x = element_text(angle = 40, hjust = 1), 
                                      panel.grid.major = element_blank(), 
                                      panel.grid.minor = element_blank(), 
                                      legend.position = "bottom", 
                                      panel.border = element_rect(fill = NA, linewidth = 2, linetype = "solid"), 
                                      axis.ticks = element_line(color = "black", size = 0.5))), 
          zlim = c(0, 100), 
          file = "Figure_6_D_Fe_Heatmap.pdf")




heatmaply(AcinetobacterMatrixOrd_Bip, 
          Rowv = dendro, 
          Colv = dendro, 
          linewidth = 0.1, 
          na.value = "black", 
          grid_size = 0.1, 
          col = viridis(100),
          #subplot_widths = c(0.2, 0.8), 
          #subplot_heights = c(0.3, 0.5), 
          scale_fill_gradient_fun = scale_fill_gradient(low = "dark green", high = "dark orange", limits = c(0, 100)), 
          showticklabels = c(T, T), 
          show_dendrogram = c(T, T),
          branches_lwd = 0.2, 
          plot_method = "ggplot", 
          margins(c(50, 50)), 
          grid_gap = 0.2, 
          heatmap_layers = list(theme_minimal(base_size = 15), 
                                theme(axis.text.x = element_text(angle = 40, hjust = 1), 
                                      panel.grid.major = element_blank(), 
                                      panel.grid.minor = element_blank(), 
                                      legend.position = "bottom", 
                                      panel.border = element_rect(fill = NA, linewidth = 2, linetype = "solid"), 
                                      axis.ticks = element_line(color = "black", size = 0.5))), 
          file = "Figure_6_E_Bip_Heatmap.pdf")


# make co-inhibition heatmap with MiBIG predicted biosynthetic gene clusters 


# extract some stats from this heatmap 

## boxplot with average inhibition scores
Acinetobacer_inhibition_scores_renamed %>%
  filter(inhibition_score != "NA") %>%
  group_by(competing_strain_left, target_strain_right, condition) %>%
  summarise(average_inhibition_score = mean(inhibition_score)) %>%
  ggplot() + 
  aes(x = condition, y = average_inhibition_score, fill = condition) + 
  scale_fill_manual(values = c("FeCl3" = "#52796f", "2,2′-Bipyridyl" = "#ffd166")) +   # Specify fill colors for groups
  #geom_bar(stat = "identity", position = "dodge") + 
  geom_violin() + 
  geom_jitter(width = 0.1, height = 2, alpha = 0.2) + 
  theme_bw(base_size = 30) + 
  ylab("Average Inhibition Score") + 
  xlab("Condition") + 
  ylim(0, 100) + 
  labs(fill = "Condition") + 
  theme(legend.position = "bottom") 

ggsave(filename = "140524_InhibitionDistribution.pdf", height = 10, width = 10)

# how many interactions were scored 
# 305 distinct interactions screened 
Acinetobacer_inhibition_scores_renamed %>%
  group_by(target_strain_right, competing_strain_left) %>%
  summarise(n())

# 609 from both FeCl3 and 2'2'-Bip
Acinetobacer_inhibition_scores_renamed %>%
  group_by(target_strain_right, competing_strain_left, condition) %>%
  summarise(n())

# What was the difference between distribution grey values in both conditions 
InhibitionSummary <- Acinetobacer_inhibition_scores_renamed %>%
  filter(colony_grey_value != "NA") %>%
  group_by(condition) %>%
  summarise(mean(colony_grey_value), 
            range(colony_grey_value), 
            sd(colony_grey_value)) 

# How many have inhibition scores == 1000% 
Acinetobacer_inhibition_scores_renamed %>%
  filter(condition == "2,2′-Bipyridyl") %>%
  filter(inhibition_score == 100) %>%
  #group_by(target_strain_right, competing_strain_left) %>%
  summarise(n())

Acinetobacer_inhibition_scores_renamed %>%
  filter(condition == "2,2′-Bipyridyl") %>%
  filter(inhibition_score <=99 & inhibition_score >= 90) %>%
  group_by(target_strain_right, competing_strain_left) %>%
  summarise(n())

Acinetobacer_inhibition_scores_renamed %>%
  filter(condition == "2,2′-Bipyridyl") %>%
  filter(inhibition_score <90 & inhibition_score >= 10) %>%
  group_by(target_strain_right, competing_strain_left) %>%
  summarise(n())

Acinetobacer_inhibition_scores_renamed %>%
  filter(condition == "2,2′-Bipyridyl") %>%
  filter(inhibition_score <=10) %>%
  group_by(target_strain_right, competing_strain_left) %>%
  summarise(n())

# percentage inhibition? 
Acinetobacer_inhibition_scores_renamed %>%
  filter(condition == "2,2′-Bipyridyl") %>%
  group_by(target_strain_right, competing_strain_left, condition) %>%
  summarise(mean_inhibition_score = mean(inhibition_score)) %>%
  group_by(condition, mean_inhibition_score) %>%
  summarise(number_of_interactions = n()) %>%
  mutate(percentage_strains = number_of_interactions / sum(number_of_interactions) * 100) %>%
  print(n = 30)

Acinetobacer_inhibition_scores_renamed %>%
  filter(condition == "FeCl3") %>%
  group_by(target_strain_right, competing_strain_left, condition) %>%
  summarise(mean_inhibition_score = mean(inhibition_score)) %>%
  group_by(condition, mean_inhibition_score) %>%
  summarise(number_of_interactions = n()) %>%
  mutate(percentage_strains = number_of_interactions / sum(number_of_interactions) * 100) %>%
  print(n = 30)
  

# make boxplots comparing USMMO68 inhibition and HN121 inhibition 
Acinetobacer_inhibition_scores_renamed %>%
  filter(competing_strain_left == "USMM068" | competing_strain_left == "HN121") %>%
  ggplot() + 
  aes(x = competing_strain_left, y = inhibition_score, fill = condition) + 
  scale_fill_manual(values = c("FeCl3" = "#52796f", "2,2′-Bipyridyl" = "#ffd166")) +   # Specify fill colors for groups
  #geom_bar(stat = "identity", position = "dodge") + 
  geom_violin(outlier.alpha = 0) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.2, jitter.height = 1, dodge.width = 0.9), alpha=0.2, size = 3) +
  theme_bw(base_size = 30) + 
  ylab("Inhibition score (%)") + 
  xlab("Competing isolate") + 
  ylim(0, 100) + 
  labs(fill = "Condition") + 
  theme(legend.position = "bottom") 

ggsave("Fig_S15_HN121USMM068_InhibitionScores.pdf", height = 10, width = 10)

# plot inverse interactions of USMMO68
Acinetobacer_inhibition_scores_renamed %>%
  filter(colony_grey_value != "NA") %>%
  filter(target_strain_right == "USMM068" | target_strain_right == "HN121") %>%
  group_by(target_strain_right, competing_strain_left, condition) %>%
  summarise(mean_inhibition = mean(inhibition_score)) %>%
  ggplot() + 
  aes(x = reorder(competing_strain_left, mean_inhibition), y = mean_inhibition, fill = condition) + 
  geom_bar(stat = "identity", position = "dodge", width = 0.5) + 
  ylim(0, 100) + 
  theme_bw(base_size = 30) + 
  xlab("Target strain") + 
  ylab("Inhibition score (%)") + 
  labs(fill = "Condition") + 
  theme(axis.text.x = element_text(angle = 50, hjust = 1), 
        legend.position = "bottom", 
        strip.background = element_rect(fill = "white", colour = NA)) + 
  facet_wrap(facets = "target_strain_right", ncol = 1) + 
  scale_fill_manual(values = c("FeCl3" = "#52796f", "2,2′-Bipyridyl" = "#ffd166"))    # Specify fill colors for groups
  
ggsave("Fig_S16_TargetInhibitionScores.pdf", height = 10, width = 10)


# get some stats from this 
Acinetobacer_inhibition_scores_renamed %>%
  filter(competing_strain_left == "USMM068" | competing_strain_left == "HN121") %>%
  group_by(condition, competing_strain_left, colony_grey_value) %>%
  summarise(n()) %>%
  print(n = 400 )

Acinetobacer_inhibition_scores_renamed %>%
  filter(competing_strain_left == "USMM068" | competing_strain_left == "HN121") %>%
  select(condition, competing_strain_left, colony_grey_value) %>%
  print(n = 400 )

group1 <- Acinetobacer_inhibition_scores_renamed %>%  
  filter(competing_strain_left == "USMM068" | competing_strain_left == "HN121") %>%
  select(condition, competing_strain_left, colony_grey_value) %>%
  filter(condition == "2,2′-Bipyridyl" & competing_strain_left == "USMM068") %>%
  select(colony_grey_value) %>%
  as.vector()

group2 <- Acinetobacer_inhibition_scores_renamed %>%  
  filter(competing_strain_left == "USMM068" | competing_strain_left == "HN121") %>%
  select(condition, competing_strain_left, colony_grey_value) %>%
  filter(condition == "2,2′-Bipyridyl" & competing_strain_left == "HN121") %>%
  select(colony_grey_value) %>%
  as.vector()

group1$colony_grey_value

wilcox.test(group1$colony_grey_value, group2$colony_grey_value)

# overall inhibition scores between these two strains 
Acinetobacer_inhibition_scores_renamed %>%
  filter(colony_grey_value != "NA") %>%
  filter(competing_strain_left == "USMM068" | competing_strain_left == "HN121") %>%
  group_by(condition, competing_strain_left) %>%
  summarise(mean(inhibition_score), 
            range(inhibition_score), 
            sd(inhibition_score)) 

# fold reduction of average inhibition score by USMMO68
Acinetobacer_inhibition_scores_renamed %>%
  filter(colony_grey_value != "NA") %>%
  filter(competing_strain_left == "USMM068") %>%
  group_by(condition, competing_strain_left) %>%
  summarise(mean(inhibition_score))

foldchange(93.4, 4.41)
  
