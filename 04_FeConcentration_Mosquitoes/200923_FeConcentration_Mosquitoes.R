# packages
library(tidyverse)

# load data 
feconc = read_tsv("200923_FeConc_LifeStages.txt")
sample_order = c("LW Day 1", "LW Day 7", "Larvae (L2/3)", "Larvae (L4)", "Pupae", "Adult Male", "Adult Female", "Adult Female 24hr PBM")

# rename the groups 
feconc_clean = feconc %>%
  mutate(group = ifelse(group == "d1_larval_water", "LW Day 1", group), 
         group = ifelse(group == "d7_larval_water", "LW Day 7", group), 
         group = ifelse(group == "L2_3", "Larvae (L2/3)", group), 
         group = ifelse(group == "L4", "Larvae (L4)", group), 
         group = ifelse(group == "Pupae", "Pupae", group), 
         group = ifelse(group == "adult_male", "Adult Male", group), 
         group = ifelse(group == "adult_female", "Adult Female", group), 
         group = ifelse(group == "adult_female_BF", "Adult Female 24hr PBM", group))

# plot fe concentration 
feconc_clean %>%
  mutate(group = factor(group, levels = sample_order)) %>%
  ggplot() + 
  aes(x = group, y = OD_Norm) + 
  geom_boxplot() + 
  geom_point() + 
  theme_bw(base_size = 25) + 
  xlab("Life stage") + 
  ylab("OD593") + 
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

ggsave(filename = "250923_FeConcentration_Aedesaegypti_LVP_Lifestages.pdf", height = 10, width = 10)

# split tibbles into experimental groups 
#
 # Make a lm using the standard values
#

larvae_pupae_LW_group = feconc_clean %>%
  filter(group == "Larvae (L2/3)" | group == "Larvae (L4)" |
           group == "Pupae" | group == "LW Day 1" | group == "LW Day 7") %>%
  select(c(PC0, PC2, PC4, PC6, PC8, PC10)) %>%
  slice(1) %>%
  pivot_longer(cols = c(PC0, PC2, PC4, PC6, PC8, PC10), names_to = "conc", values_to = "OD") %>%
  separate(conc, into = c("letter", "conc"), sep = "PC") %>%
  mutate(conc = as.numeric(conc)) %>%
  select(!letter)

lw_pupae_larvae_values = feconc_clean %>%
  filter(group == "Larvae (L2/3)" | group == "Larvae (L4)" |
           group == "Pupae" | group == "LW Day 1" | group == "LW Day 7") %>%
  select(group, OD_Norm)

adults_NBF = feconc_clean %>%
  filter(group == "Adult Male" | group == "Adult Female") %>%
  select(c(PC0, PC2, PC4, PC6, PC8, PC10)) %>%
  slice(1) %>%
  pivot_longer(cols = c(PC0, PC2, PC4, PC6, PC8, PC10), names_to = "conc", values_to = "OD") %>%
  separate(conc, into = c("letter", "conc"), sep = "PC") %>%
  mutate(conc = as.numeric(conc)) %>%
  select(!letter)

adults_NBF_values = feconc_clean %>%
  filter(group == "Adult Male" | group == "Adult Female") %>%
  select(group, OD_Norm)

adults_BF = feconc_clean %>%
  filter(group == "Adult Female 24hr PBM") %>%
  select(c(PC0, PC2, PC4, PC6, PC8, PC10)) %>%
  slice(1) %>%
  pivot_longer(cols = c(PC0, PC2, PC4, PC6, PC8, PC10), names_to = "conc", values_to = "OD") %>%
  separate(conc, into = c("letter", "conc"), sep = "PC") %>%
  mutate(conc = as.numeric(conc)) %>%
  select(!letter)

adults_BF_values = feconc_clean %>%
  filter(group == "Adult Female 24hr PBM") %>%
  select(group, OD_Norm)

  # function to predict Fe Concentration using OD values 

PredictFeConc <- function(filtered_table_samples, filtered_table_PCs){
  # filtered_table_samples = sample names
  # filtered_table_PCs = positive control values  
  
  model <- lm(conc ~ OD, data = filtered_table_PCs)
  
  # make empty DF
  predicted_values_df <- data.frame(OD = numeric(0), PredictedValue = numeric(0))
  
  # for loop to predict iron concentration using standard curve
  for (i in 1:nrow(filtered_table_samples)){
    current_od <- filtered_table_samples$OD_Norm[i]
    
    # predict the value based on current OD valus 
    predicted_value <- predict(model, newdata = data.frame(OD = current_od))
    
    # add the result to the predicted values in the dataframe 
    predicted_values_df <- rbind(predicted_values_df, data.frame(OD = current_od, PredictedValue = predicted_value))
  }
  
  predicted_iron_concentration <- filtered_table_samples %>%
    left_join(predicted_values_df, by = c("OD_Norm" = "OD" )) %>%
    mutate(PredictedValueConv = ifelse(PredictedValue < 0, 0, PredictedValue))
  
  return(predicted_iron_concentration)
  
}

LW_L2_L3_Pupae_Fe = PredictFeConc(lw_pupae_larvae_values, larvae_pupae_LW_group)
adult_NBF = PredictFeConc(adults_NBF_values, adults_NBF)
adult_BF = PredictFeConc(adults_BF_values, adults_BF)

final_table = rbind(LW_L2_L3_Pupae_Fe, adult_NBF, adult_BF)

final_table %>%
  mutate(group = factor(group, levels = sample_order)) %>%
  ggplot() + 
  aes(x = group, y = PredictedValueConv) + 
  geom_boxplot() +
  geom_point() + 
  theme_bw(base_size = 25) + 
  xlab("Life stage") + 
  ylab("Fe (nmole/sample)") + 
  ylim(0, 4) + 
  theme(axis.text.x = element_text(angle = 15, hjust = 1))

ggsave("250923_FeConcentrationLifestages.pdf", width = 10)

#### Plot Pools 
feconc_clean %>%
  mutate(group = factor(group, levels = sample_order)) %>%
  mutate(pool = ifelse(group == "LW Day 1" | group == "LW Day 7", "Larval Water", pool)) %>%
  mutate(pool = ifelse(group == "Adult Female 24hr PBM", "Blood-Fed Female", pool)) %>%
  ggplot() + 
  aes(x = group, y = OD_Norm) + 
  geom_boxplot() + 
  geom_point() + 
  facet_wrap(facets = "pool", scales = "free_x", nrow = 1) + 
  theme_bw(base_size = 20) + 
  theme(axis.text.x = element_text(angle = 40, hjust = 1)) + 
  ylab("Normalised OD") + 
  xlab("Life Stage")

ggsave("250923_FeConcentrationPools.pdf", width = 10)

# get statistics 
feconc_clean %>%
  group_by(group, pool) %>%
  summarise(average_OD = mean(OD))

feconc_clean %>%
  group_by(group, pool) %>%
  summarise(average_OD = mean(OD_Norm))

feconc_clean %>%
  group_by(group, NC) %>%
  summarise(average_OD = mean(NC))

feconc_clean %>%
  select(group, OD_Norm) %>%
  print(n = 400)

# anova - test sig differences between groups
fe_conc_anova = kruskal.test(PredictedValueConv ~ group, data = final_table)
summary(fe_conc_anova)

