library(tidyverse)
library(ggsignif)

premadecultures <- read_csv("070524_5050_PremadeCultures.csv") %>%
  filter(plate != "NA")

premadecultures_long <- premadecultures %>%
  pivot_longer(cols = c("a_soli_cfu", "a_johnsonii_cfu"), 
               names_to = "acinetobacter_spp", 
               values_to = "cfu_ml_2") 


facet_labels <- c("25_75" = "75:25", "50_50" = "50:50", "75_25" = "25:75")


premadecultures_long_5050 <- premadecultures_long %>% filter(condition == "50_50")
premadecultures_long_2575 <- premadecultures_long %>% filter(condition == "25_75")
premadecultures_long_7525 <- premadecultures_long %>% filter(condition == "75_25")

wilcox.test(premadecultures_long_5050$cfu_ml_2 ~ acinetobacter_spp, premadecultures_long_5050) # non-sig
wilcox.test(premadecultures_long_2575$cfu_ml_2 ~ acinetobacter_spp, premadecultures_long_2575) # sig
wilcox.test(premadecultures_long_7525$cfu_ml_2 ~ acinetobacter_spp, premadecultures_long_7525) # sig


annotation_df_art <- data.frame(
  condition = "25_75",
  acinetobacter_spp = c(""),
  start = "a_johnsonii_cfu", 
  end = "a_soli_cfu",
  y = 8.5e+08,
  label = "*"
)

annotation_df_art_2 <- data.frame(
  condition = "75_25",
  acinetobacter_spp = c(""),
  start = "a_johnsonii_cfu", 
  end = "a_soli_cfu",
  y = 8.5e+08,
  label = "*"
)

annotation_df_art_3 <- data.frame(
  condition = "50_50",
  acinetobacter_spp = c(""),
  start = "a_johnsonii_cfu", 
  end = "a_soli_cfu",
  y = 8.5e+08,
  label = "n.s"
)

# plot premade cultures
premadecultures_long %>%
  ggplot() +
  aes(x = acinetobacter_spp, y = cfu_ml_2, fill = acinetobacter_spp) + 
  geom_boxplot(outlier.alpha = 0) + 
  geom_point(alpha = 0.3) + 
  ylim(0, 9e+08) + 
  theme_bw(base_size = 30) + 
  xlab("Species") + 
  ylab("CFU / mL") + 
  facet_grid(. ~ condition, labeller = labeller(condition = facet_labels)) + 
  scale_fill_manual(values = c("#ECAE92", "#A4CDD9", "#ECAE92", "#A4CDD9")) +
  scale_x_discrete(labels = c("a_soli_cfu" = "A. soli", 
                              "a_johnsonii_cfu" = "A. johnsonii", 
                              "co_a_soli_cfu" = "A. soli", 
                              "co_a_johnsonii_cfu" = "A. johnsonii")) + 
  theme(axis.text.x = element_text(face = "italic", angle = 10, hjust = 1), 
        legend.position = "none") + 
  geom_signif(
    data = annotation_df_art, 
    aes(xmin = start, xmax = end, annotations = label, y_position = y), 
    textsize = 7, tip_length = .01, vjust = 0.3,
    manual = TRUE
  ) + 
  geom_signif(
    data = annotation_df_art_2, 
    aes(xmin = start, xmax = end, annotations = label, y_position = y), 
    textsize = 7, tip_length = .01, vjust = 0.3,
    manual = TRUE
  ) + 
  geom_signif(
    data = annotation_df_art_3, 
    aes(xmin = start, xmax = end, annotations = label, y_position = y), 
    textsize = 7, tip_length = .01, vjust = -0.3,
    manual = TRUE
  ) 

ggsave("260824_SuppFigX_PremadeCultures.pdf", height = 10, width = 11)

# 
premadecultures %>%
  pivot_longer(cols = c("a_soli_cfu", "a_johnsonii_cfu"), 
               names_to = "acinetobacter_spp", 
               values_to = "cfu_ml_2") %>%
  ggplot() +
  aes(x = condition, y = cfu_ml_2, fill = acinetobacter_spp) + 
  geom_bar(stat = "identity", position = "fill", alpha = 0.8) + 
  theme_bw(base_size = 25) + 
  xlab("Species") + 
  ylab("cfu / mL") + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + 
  scale_fill_manual(values = c("dark orange", "dark green"), 
                    labels = c("A. soli", "A. johnsonii")) 



