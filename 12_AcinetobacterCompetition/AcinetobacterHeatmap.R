library(tidyverse)

tbl <- read_csv("140524_Acinetobacter_PixelIntensity.csv")

# rename spelling error 
tbl <- tbl %>%
  mutate(condition = ifelse(condition == "2'2-Bipyrydyl", "2,2′-Bipyridyl", condition))

colours <- c("#ffd166", "#52796f")

## summarise pixel instensity in iron rich and iron deficient media 
tbl %>%
  mutate(rounded_values = round(colony_grey_value/100) * 100) %>%
  group_by(rounded_values, condition) %>%
  #filter(rounded_values != 0) %>%
  ggplot() + # plot into histogram
  aes(x = rounded_values, fill = condition) + 
  #geom_bar(stat = "identity") + 
  geom_histogram(bins = 10, binwidth = 800, na.rm = T, alpha = 0.8) + 
  scale_fill_manual(values = colours) + 
  theme_bw(base_size = 25) + 
  theme(legend.position = "bottom") + 
  xlab("Mean grey value") + 
  ylab("Number of samples") + 
  labs(fill = "Condition")

tbl_summary <- tbl %>%
  mutate(rounded_values = round(colony_grey_value/100) * 100) %>%
  group_by(rounded_values, condition) 

ggplot(tbl_summary, aes(x = rounded_values, fill = condition)) + 
  geom_histogram(alpha = 0.5, position = "identity", bins = 100) +  # Set a default alpha value for both groups
  geom_histogram(data = subset(tbl_summary, condition == "FeCl3"), alpha = 0.5, position = "identity", bins = 100) +   # Set a different alpha value for Group 1
  geom_histogram(data = subset(tbl_summary, condition == "2,2′-Bipyridyl"), alpha = 0.5, position = "identity", bins = 100) +    # Set a different alpha value for Group 1
  scale_fill_manual(values = c("FeCl3" = "#52796f", "2,2′-Bipyridyl" = "#ffd166")) +   # Specify fill colors for groups
  theme_bw(base_size = 30) + 
  theme(legend.position = "bottom") + 
  xlab("Mean grey value") + 
  ylab("Interactions (N)") + 
  labs(fill = "Condition")

ggsave(filename = "Figure_6B_InhibitionHistogram.pdf", height = 10, width = 10)

group1 <- tbl_summary %>%
  select(condition, rounded_values) %>%
  filter(condition == "2,2′-Bipyridyl") %>%
  as.vector()

group2 <- tbl_summary %>%
  select(condition, rounded_values) %>%
  filter(condition == "FeCl3") %>%
  as.vector()

t.test(group1$rounded_values, group2$rounded_values)

# do a t-test on these two conditions 


# Using this distribution, reassign values above a given grey value to an inhibition score 
## >10000, no inhibition e.g. 0
## 7500 -->10000, 25
## 5000 --> 7500, 50
## 2500 --> 5000, 75
## 1000 --> 2500, 90
## 50 --> 1000, 99
## <10, 1000 

tbl_summary_inhibition_scores <- tbl_summary %>%
  mutate(inhibition_score = ifelse(colony_grey_value > 10000, 0,
                                   ifelse(colony_grey_value <= 10000 & colony_grey_value > 7500, 25,
                                          ifelse(colony_grey_value <= 7500 & colony_grey_value > 5000, 50, 
                                                 ifelse(colony_grey_value <= 5000 & colony_grey_value > 2500, 75, 
                                                        ifelse(colony_grey_value <= 2500 & colony_grey_value > 1500, 90, 
                                                               ifelse(colony_grey_value <= 1500 & colony_grey_value > 1000, 97, 
                                                                      ifelse(colony_grey_value <= 1000 & colony_grey_value > 30, 99, 
                                                                             ifelse(colony_grey_value <= 30, 100, colony_grey_value)))))))))
  
tbl_summary_inhibition_scores$condition <- factor(tbl_summary_inhibition_scores$condition, levels = c("FeCl3", "2,2′-Bipyridyl"))

write.table(x = tbl_summary_inhibition_scores, file = "140524_AcinetobacterInhibitionScoresSummarised.tsv", quote = F, sep = "\t", row.names = F)

tbl_summary_inhibition_scores %>%
  filter(target_strain_right != "NA") %>%
  filter(competing_strain_left != "NA") %>%
  filter(condition != "NA") %>%
  ggplot() + 
  aes(x = target_strain_right, y = competing_strain_left, fill = inhibition_score) + 
  geom_tile() + 
  theme_minimal(base_size = 15) + 
  facet_wrap(facets = "condition") + 
  scale_fill_gradient(low = "dark green", high = "dark orange") + 
  theme(axis.text.x = element_text(angle = 40, hjust = 1), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "bottom", 
        panel.border = element_rect(fill = NA, linewidth = 2, linetype = "solid"), 
        axis.ticks = element_line(color = "black", size = 0.5)) + 
  ylab("Competing strain") + 
  xlab("Target strain") + 
  labs(fill = "Inhibition Score")

ggsave(filename = "170424_Inhibition_Heatmap.pdf", height = 6, width = 10)

### summary statistics 
#### how many interactions did we test? 
tbl_summary_inhibition_scores 

