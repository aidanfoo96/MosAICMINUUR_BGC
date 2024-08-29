## Load packages
library(tidyverse)

## load data 
feconc <- read_tsv("200923_FeConc_LifeStages.txt") %>%
  mutate(split_group = ifelse(grepl("^PC", group), group, NA)) %>%
  separate(split_group, into = c("prefix", "fe_conc"), sep = "_", remove = FALSE, convert = TRUE) %>%
  select(-split_group) 
  
conc_order = c("PC_0", "PC_2", "PC_4", "PC_6", "PC_8", "PC_10")
sample_order = c("d1_larval_water", "d7_larval_water", "L2_3", "L4", "Pupae", "adult_male", "adult_female")

# wrangled DF 
feconc2 <- feconc %>%
  select(!fe_conc) %>%
  mutate(average_value = rowMeans(across(where(is.numeric)), na.rm = TRUE)) %>%
  mutate(true_OD = ifelse(group %in% c("d1_larval_water", "d7_larval_water"), group, average_value - 0.159))

# alternative where all values are modified 
feconc3 <- feconc %>%
  mutate_at(vars(dup1, dup2, dup3, dup4), ~ . - 0.159)


# Standard curve plot
feconc %>%
  filter(grepl("^PC", group)) %>%
  pivot_longer(cols = !c(group, pool, prefix, fe_conc), names_to = "duplicate", values_to = "OD") %>%
  #mutate(group = factor(group, ordered = T, levels = conc_order)) %>%
  ggplot() + 
  aes(x = fe_conc, y = OD) + 
  geom_point() + 
  geom_smooth(method = "lm", se = F) + 
  theme_bw(base_size = 25) + 
  xlab("Iron (nmole/well)") + 
  ylab("Absorbance (OD593nm)") # line of best fit

ggsave(filename = "180923_PosControl_Controlsnt.pdf", height = 10, width = 10)


# Values from the standard curve
fe_conc_summary <- feconc %>%
  filter(grepl("^PC", group)) %>%
  pivot_longer(cols = !c(group, pool, prefix, fe_conc), names_to = "duplicate", values_to = "OD") %>%
  select(fe_conc, OD)

# linear model on value from standard curve 
model <- lm(fe_conc ~ OD, data = fe_conc_summary)
summary(model)

# make df with OD values from sample groups 
OD_values <- feconc3 %>%
  mutate(dup1 = ifelse(pool == 2.0, dup1 / 2, dup1), 
         dup2 = ifelse(pool == 2.0, dup2 / 2, dup2), 
         dup3 = ifelse(pool == 2.0, dup3 / 2, dup3), 
         dup4 = ifelse(pool == 2.0, dup4 / 2, dup4)) %>%
  mutate(group = factor(group, levels = sample_order)) %>%
  filter(!grepl("^PC", group)) %>%
  filter(group != "NC" & group != "NC_Bead") %>%
  filter(group != "d1_larval_water" & group != "d7_larval_water") %>%
  pivot_longer(cols = !c(group, pool, prefix, fe_conc), names_to = "duplicate", values_to = "OD") %>%
  select(group, pool, OD) %>%
  mutate(pool_fac = as.factor(pool)) 

# facetted plot with OD values and bools 
OD_values_2 <- feconc3 %>%
  filter(!grepl("^PC", group)) %>%
  filter(group != "NC" & group != "NC_Bead") %>%
  pivot_longer(cols = !c(group, pool, prefix, fe_conc), names_to = "duplicate", values_to = "OD") %>%
  select(group, pool, OD) %>%
  mutate(pool_fac = as.factor(pool)) 

OD_values_3 <- feconc %>%
  filter(!grepl("^PC", group)) %>%
  filter(group == "NC" | group == "NC_Bead") %>%
  pivot_longer(cols = !c(group, pool, prefix, fe_conc), names_to = "duplicate", values_to = "OD") %>%
  select(group, pool, OD) %>%
  mutate(pool_fac = as.factor(pool)) 

OD_values_4 <- feconc2 %>%
  filter(!grepl("^PC", group)) %>%
  filter(group == "d1_larval_water" | group == "d7_larval_water") %>%
  mutate(group = factor(group, levels = sample_order)) %>%
  pivot_longer(cols = !c(group, pool, prefix, average_value, true_OD), names_to = "duplicate", values_to = "OD") %>%
  select(group, pool, OD) %>%
  mutate(pool_fac = as.factor(pool)) 

OD_values_5 = rbind(OD_values_4, OD_values)

# Plot OD values 
OD_values_5 %>%
  ggplot() + 
  aes(x = group, y = OD) + 
  geom_point() + 
  geom_boxplot()  + 
  theme_bw(base_size = 25) + 
  ylim(0, 0.3) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + 
  #facet_wrap(facets = "group", nrow = 3) + 
  xlab("Life stage") + 
  ylab("OD593")

ggsave(filename = "180923_OD593_LarvalWater.pdf", height = 10, width = 10)


# create blank df 
predicted_values_df <- data.frame(OD = numeric(0), PredictedValue = numeric(0))

# for loop to predict iron concentration using standard curve
for (i in 1:nrow(OD_values)){
  current_od <- OD_values$OD[i]
  
  # predict the value based on current OD valus 
  predicted_value <- predict(model, newdata = data.frame(OD = current_od))
  
  # add the result to the predicted values in the dataframe 
  predicted_values_df <- rbind(predicted_values_df, data.frame(OD = current_od, PredictedValue = predicted_value))
}

predicted_values_df

# rejoin dataframe, appending the metadata from the original table
predicted_mosquito_iron_concentration <- OD_values %>%
  left_join(predicted_values_df, by = c("OD" = "OD" )) %>%
  mutate(value = ifelse(value == 0, 1, value)) %>%
  mutate(IronConcPerMos = (PredictedValue / value)) %>%
  mutate(SampleID = paste(group, value, sep = "_"))

# plot iron concentration again pooled samples
predicted_mosquito_iron_concentration %>%
  mutate(SampleID = factor(SampleID, levels = larvae_order)) %>%
  ggplot() + 
  aes(x = SampleID, y = PredictedValue) + 
  geom_boxplot() +
  geom_point() + 
  theme_bw(base_size = 25) + 
  xlab("Pooled Larvae Number") + 
  ylab("Iron (nmole/well)") + 
  ylim(0, 10)

ggsave(filename = "040923_FeConcentration_Larvae_Well.pdf", height = 10, width = 10)


predicted_mosquito_iron_concentration %>%
  mutate(SampleID = factor(SampleID, levels = larvae_order)) %>%
  ggplot() + 
  aes(x = SampleID, y = IronConcPerMos) + 
  geom_boxplot() +
  geom_point() + 
  theme_bw(base_size = 25) + 
  xlab("Pooled Larvae Number") + 
  ylab("Iron (nmole/larvae)") + 
  ylim(0, 10)

ggsave(filename = "040923_FeConcentration_Per_Larvae.pdf", height = 10, width = 10)


