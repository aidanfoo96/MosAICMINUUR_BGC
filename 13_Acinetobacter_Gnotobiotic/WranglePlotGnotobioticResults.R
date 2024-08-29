library(tidyverse)
library(car)
library(ggsignif)
library(ggh4x)
library(ggforce)
library(gtools)

cfu_data <- read_csv("060525_GnotoResults.csv")

cfu_data_sorted <- cfu_data %>%
  arrange(lifestage, plateID)

write.table(file = "010524_PCR_results_sorted.tsv", cfu_data_sorted, sep = "\t", quote = F, row.names = F)

sample_order = c("a_soli_cfu", "a_johnsonii_cfu")
condition_order = c("mono", "co")
life_stage_order = c("L4", "pupae", "adult")

cfu_data_long <- cfu_data %>%
  pivot_longer(cols = c("a_soli_cfu", "a_johnsonii_cfu"), 
               names_to = "acinetobacter_spp", 
               values_to = "cfu_ml_2") %>%
  mutate(acinetobacter_spp = factor(acinetobacter_spp, levels = sample_order)) %>%
  mutate(condition = factor(condition, levels = condition_order)) %>%
  mutate(lifestage = factor(lifestage, levels = life_stage_order)) 

cfu_data_long %>%
  #filter(lifestage != "adult") %>%
  filter(lifestage != "NA") %>%
  ggplot() + 
  aes(x = condition, y = cfu_ml_2, fill = acinetobacter_spp) +
  geom_boxplot(outlier.alpha = 0, alpha = 0.7) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.25), alpha=0.4, size = 2) +
  facet_wrap(facets = "lifestage") + 
  theme_bw(base_size = 30) + 
  ylab("CFU / mL") + 
  xlab("Condition") + 
  labs(fill = "", 
       labels = c("fsdf", "fsldif")) + 
  scale_fill_manual(values = c("#ECAE92", "#A4CDD9"), 
                    labels = c("A. soli", "A. johnsonii")) + 
  guides(fill = guide_legend(ncol = 2)) +# Change the number of columns
  theme(legend.position = c(1.3, 1.2),  # Specify x, y coordinates of the legend
        legend.justification = c(1, 1),  # Specify justification (1,1) = top right
        legend.margin = margin(t = -10, r = -10, b = -10, l = -10), 
        plot.margin = margin(1, 8, 1, 0.5, "cm"), 
        strip.background = element_rect(fill = "white", colour = NA)) + # Set margins for the whole plot) + 
  ylim(0, 3000000)

ggsave(filename = "200524_A_soli_A_johnsonii.pdf", height = 10, width = 12)

cfu_data_long %>%
  filter(lifestage == "adult") %>%
  filter(lifestage != "NA") %>%
  ggplot() + 
  aes(x = condition, y = cfu_ml_2, fill = acinetobacter_spp) +
  geom_boxplot(outlier.alpha = 0, alpha = 0.7) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.25), alpha=0.4, size = 2) +
  facet_wrap(facets = "lifestage") + 
  theme_bw(base_size = 30) + 
  ylab("CFU / mL") + 
  xlab("Condition") + 
  labs(fill = "") + 
  guides(fill = "none") +# Change the number of columns
  scale_fill_manual(values = c("#ECAE92", "#A4CDD9"), 
                    labels = c("A. soli", "A. johnsonii")) + 
  theme(legend.position = c(1.85, 1.2),  # Specify x, y coordinates of the legend
        legend.justification = c(1, 1),  # Specify justification (1,1) = top right
        legend.margin = margin(t = -10, r = -10, b = -10, l = -10), 
        plot.margin = margin(1, 8, 1, 0.5, "cm"), 
        strip.background = element_rect(fill = "white", colour = NA)) + # Set margins for the whole plot) + 
  ylim(0, 70000)

ggsave(filename = "200524_A_soli_A_johnsonii_adult.pdf", height = 10, width = 7)


### plot relative abundance 
cfu_data_long %>%
  filter(condition == "co") %>%
  ggplot() +
  aes(x = lifestage, y = cfu_ml_2, fill = acinetobacter_spp) + 
  geom_bar(stat = "identity", position = "fill", alpha = 0.8) + 
  theme_bw(base_size = 20) + 
  labs(fill = "") + 
  ylab("Co-inoculated abundance (%)") + 
  xlab("Lifestage") + 
  scale_fill_manual(values = c("dark orange", "dark green"), 
                    labels = c("A. soli", "A. johnsonii")) 


# stats 

# function to perform normality tests
normality_tests <- cfu_data_long %>%
  select(condition, lifestage, acinetobacter_spp, cfu_ml_2) %>%
  group_by(condition, lifestage, acinetobacter_spp) %>%
  summarise_all(.funs = funs(statistic = shapiro.test(cfu_ml_2)$statistic,
                             p.value = shapiro.test(cfu_ml_2)$p.value))

# function to extract target group 
extract_target_group <- function(lifestage_target, species_target) {
  
  cfu_data_long %>%
    filter(lifestage == lifestage_target) %>%
    filter(acinetobacter_spp == species_target) %>%
    filter(cfu_ml_2 != "NA")
  
}

pupae_ajon <- extract_target_group("pupae", "a_johnsonii_cfu")
pupae_asol <- extract_target_group("pupae", "a_soli_cfu")
larvae_ajon <- extract_target_group("L4", "a_johnsonii_cfu")
larvae_asol <- extract_target_group("L4", "a_soli_cfu")
adult_ajon <- extract_target_group("adult", "a_johnsonii_cfu")
adult_asol <- extract_target_group("adult", "a_soli_cfu")

leveneTest(cfu_ml_2 ~ condition, pupae_asol) # variance across samples is equal 
leveneTest(cfu_ml_2 ~ condition, pupae_ajon) # variance across samples are not equal
leveneTest(cfu_ml_2 ~ condition, larvae_ajon) # variance across samples are equal
leveneTest(cfu_ml_2 ~ condition, larvae_asol) # variance across samples are not equal
leveneTest(cfu_ml_2 ~ condition, adult_ajon) # variance across samples are not equal
leveneTest(cfu_ml_2 ~ condition, adult_asol) # variance across samples are not equal


### uneven distributions, and variance across samples uneven, therefore do non-parametric 
#### peform Mann-Whitney U test 

# significance tests 
# tests between a johnsonii in mono and co inoculatations in different lifestages
wilcox.test(pupae_ajon$cfu_ml_2 ~ condition, pupae_ajon) # sig difference between cfu in co and mono pupae
wilcox.test(larvae_ajon$cfu_ml_2 ~ condition, larvae_ajon) # sig difference between cfu in co and mono pupae
wilcox.test(adult_ajon$cfu_ml_2 ~ condition, adult_ajon) # no sig difference between cfu in co and mono adults

# tests between a soli in mono and co inoculatations in different lifestages
wilcox.test(larvae_asol$cfu_ml_2 ~ condition, larvae_asol) # sig difference between cfu in co and mono pupae
wilcox.test(pupae_asol$cfu_ml_2 ~ condition, pupae_asol) # sig difference between cfu in co and mono pupae
wilcox.test(adult_asol$cfu_ml_2 ~ condition, adult_asol) # no sig difference between cfu in co and mono adults

# fold changes 
larvae_ajon %>%
  group_by(condition) %>%
  summarise(mean(cfu_ml_2))

foldchange(1558175, 88833)

pupae_ajon %>%
  group_by(condition) %>%
  summarise(mean(cfu_ml_2))

foldchange(326500, 35520)

adult_ajon %>%
  group_by(condition) %>%
  summarise(mean(cfu_ml_2))

foldchange(8334, 719)


# group cfu by lifestage and condition to compare a soli and a johnsonii

extract_target_group_condition <- function(lifestage_target, target_condition) {
  
  cfu_data_long %>%
    filter(lifestage == lifestage_target) %>%
    filter(condition == target_condition) %>%
    filter(cfu_ml_2 != "NA")
  
}

pupae_mono <- extract_target_group_condition("pupae", "mono")
pupae_co <- extract_target_group_condition("pupae", "co")
larvae_mono <- extract_target_group_condition("L4", "mono")
larvae_co <- extract_target_group_condition("L4", "co")
adult_mono <- extract_target_group_condition("adult", "mono")
adult_co <- extract_target_group_condition("adult", "co")

# test differences between a soli and johnsonii in mono-inoculated groups
wilcox.test(larvae_mono$cfu_ml_2 ~ acinetobacter_spp, larvae_mono) # non-sig
wilcox.test(pupae_mono$cfu_ml_2 ~ acinetobacter_spp, pupae_mono) # non-sig
wilcox.test(adult_mono$cfu_ml_2 ~ acinetobacter_spp, adult_mono) # non-sig

# test differences between a soli and johnsonii in co-inoculated groups
wilcox.test(larvae_co$cfu_ml_2 ~ acinetobacter_spp, larvae_co) # sig
wilcox.test(pupae_co$cfu_ml_2 ~ acinetobacter_spp, pupae_co) # sig
wilcox.test(adult_co$cfu_ml_2 ~ acinetobacter_spp, adult_co) # sig

# fold changes between groups  
pupae_mono %>%
  group_by(acinetobacter_spp) %>%
  summarise(mean(cfu_ml_2))

foldchange(326500, 158700)

larvae_co %>%
  group_by(acinetobacter_spp) %>%
  summarise(mean(cfu_ml_2))

foldchange(696629, 88833)

pupae_co %>%
  group_by(acinetobacter_spp) %>%
  summarise(mean(cfu_ml_2))

foldchange(1162832, 35520)

adult_co %>%
  group_by(acinetobacter_spp) %>%
  summarise(mean(cfu_ml_2))

foldchange(62070, 719)

# replot with signficance values 
annotation_df <- data.frame(
  lifestage = c("a L4", "b pupae"),
  acinetobacter_spp = c(""),
  start = c("mono_a_johnsonii_cfu", "mono_a_johnsonii_cfu"), 
  end = c("co_a_johnsonii_cfu", "co_a_johnsonii_cfu"), 
  y = c(3.25e+06, 3.25e+06), 
  label = c("*", "***")
)

annotation_df_2 <- data.frame(
  lifestage = c("a L4", "b pupae", "c adult"),
  acinetobacter_spp = c("", "", ""),
  start = c("co_a_soli_cfu", "co_a_soli_cfu", "co_a_soli_cfu"), 
  end = c("co_a_johnsonii_cfu", "co_a_johnsonii_cfu", "co_a_johnsonii_cfu"), 
  y = c(3e+06, 3e+06, 3e+06), 
  label = c("**", "***", "***")
)

annotation_df_3 <- data.frame(
  lifestage = c("b pupae"),
  acinetobacter_spp = c("", "", ""),
  start = c("mono_a_soli_cfu"), 
  end = c("mono_a_johnsonii_cfu"), 
  y = c(3.5e+06),
  label = c("*")
)

# pivot wider by joining condition and acinetobacter spp 
condition_species_order <- c("mono_a_soli_cfu", "mono_a_johnsonii_cfu", "co_a_soli_cfu", "co_a_johnsonii_cfu")

cfu_data_long_2 <- cfu_data_long %>%
  mutate(condition_and_species = paste(str_c(condition, acinetobacter_spp, sep = "_"))) %>%
  mutate(condition_and_species = factor(condition_and_species, levels = condition_species_order)) %>%
  mutate(lifestage = if_else(lifestage == "adult", "c adult",lifestage)) %>%
  mutate(lifestage = if_else(lifestage == "L4", "a L4",lifestage)) %>%                                    
  mutate(lifestage = if_else(lifestage == "pupae", "b pupae",lifestage)) 

  
cfu_data_long_2 %>%
  #filter(lifestage != "adult") %>%
  filter(lifestage != "NA") %>%
  ggplot() + 
  aes(x = condition_and_species, "&", condition, y = cfu_ml_2, fill = acinetobacter_spp) +
  geom_boxplot(outlier.alpha = 0, alpha = 0.7) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.4), alpha=0.3, size = 2) + 
  theme_bw(base_size = 30) + 
  facet_wrap(facets = "lifestage") +
  ylab("CFU / mL") + 
  xlab("Condition") + 
  scale_fill_manual(values = c("#ECAE92", "#A4CDD9", "#ECAE92", "#A4CDD9"), 
                    labels = c("A. soli", "A. johnsonii")) + 
  theme(legend.position = c(1.7, 1.2),  # Specify x, y coordinates of the legend
        legend.justification = c(1, 1),  # Specify justification (1,1) = top right
        legend.margin = margin(t = -10, r = -10, b = -10, l = -10), 
        plot.margin = margin(1, 8, 1, 0.5, "cm"), 
        strip.background = element_rect(fill = "white", colour = NA), 
        axis.text.x = element_text(angle = 25, hjust = 1, size = 20)) + # Set margins for the whole plot) + 
  ylim(0, 4000000) + 
  geom_signif(
    data = annotation_df_2, 
    aes(xmin = start, xmax = end, annotations = label, y_position = y), 
    textsize = 10, tip_length = .01, vjust = 0.5,
    manual = TRUE
  ) + 
  geom_signif(
    data = annotation_df, 
    aes(xmin = start, xmax = end, annotations = label, y_position = y), 
    textsize = 10, tip_length = .01, vjust = 0.5,
    manual = TRUE
  ) + 
  geom_signif(
    data = annotation_df_3, 
    aes(xmin = start, xmax = end, annotations = label, y_position = y), 
    textsize = 10, tip_length = .01, vjust = 0.5,
    manual = TRUE
  ) + #guides(fill = guide_legend(ncol = 2)) # Change the number of columns
  scale_x_discrete(labels = c("mono_a_soli_cfu" = "A. soli", 
                              "mono_a_johnsonii_cfu" = "A. johnsonii", 
                              "co_a_soli_cfu" = "A. soli", 
                              "co_a_johnsonii_cfu" = "A. johnsonii")) 

ggsave(filename = "090624_A_soli_A_johnsonii_sigbars.pdf", height = 10, width = 12)
