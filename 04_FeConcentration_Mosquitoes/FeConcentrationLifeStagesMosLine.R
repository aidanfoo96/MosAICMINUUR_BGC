# packages
library(tidyverse)
library(ggsignif)
library(DHARMa)
library(gtools)

# load data 
feconc = read_tsv("191023_FeConc_LifeStages_MosLines.txt")
sample_order = c("LW Day 1", "LW Day 7", "Larvae (L1)", "Larvae (L2/3)", "Larvae (L4)", "Pupae", "Adult Male", "Adult Female", "Adult Female 24hr PBM")

# rename the groups 
feconc_clean = feconc %>%
  mutate(group = ifelse(group == "d1_larval_water", "LW Day 1", group), 
         group = ifelse(group == "d7_larval_water", "LW Day 7", group), 
         group = ifelse(group == "L1", "Larvae (L1)", group), 
         group = ifelse(group == "L2_3", "Larvae (L2/3)", group), 
         group = ifelse(group == "L4", "Larvae (L4)", group), 
         group = ifelse(group == "Pupae", "Pupae", group), 
         group = ifelse(group == "adult_male", "Adult Male", group), 
         group = ifelse(group == "adult_female", "Adult Female", group), 
         group = ifelse(group == "adult_female_bf", "Adult Female 24hr PBM", group)) %>%
  arrange(line, OD_Norm) %>%
  mutate(row_id = row_number()) 

# plot fe concentration 
feconc_clean %>%
  mutate(group = factor(group, levels = sample_order)) %>%
  ggplot() + 
  aes(x = group, y = OD_Norm) + 
  geom_boxplot(outlier.alpha = 0) + 
  geom_jitter(color="black", size=0.8, alpha=0.8, width = 0.1) +
  theme_bw(base_size = 25) + 
  xlab("Life stage") + 
  ylab("OD593") + 
  theme(axis.text.x = element_text(angle = 20, hjust = 1)) + 
  facet_wrap(facets = "line", scales = "free_x", nrow = 1) 

ggsave(filename = "250923_FeConcentration_Aedesaegypti_LVP_Lifestages.pdf", height = 10, width = 10)

# split tibbles into experimental groups 
#
# Make a lm using the standard values
#

L1Group = feconc_clean %>%
  filter(group == "Larvae (L1)") %>%
  select(c(PC0, PC2, PC4, PC6, PC8, PC10)) %>%
  slice(1) %>%
  pivot_longer(cols = c(PC0, PC2, PC4, PC6, PC8, PC10), names_to = "conc", values_to = "OD") %>%
  separate(conc, into = c("letter", "conc"), sep = "PC") %>%
  mutate(conc = as.numeric(conc)) %>%
  select(!letter)

L1Group_Values = feconc_clean %>%
  filter(group == "Larvae (L1)") %>%
  select(group, OD_Norm, line, row_id)

L2_L3 = feconc_clean %>%
  filter(group == "Larvae (L2/3)") %>%
  select(c(PC0, PC2, PC4, PC6, PC8, PC10)) %>%
  slice(1) %>%
  pivot_longer(cols = c(PC0, PC2, PC4, PC6, PC8, PC10), names_to = "conc", values_to = "OD") %>%
  separate(conc, into = c("letter", "conc"), sep = "PC") %>%
  mutate(conc = as.numeric(conc)) %>%
  select(!letter)

L2_L3_Values = feconc_clean %>%
  filter(group == "Larvae (L2/3)") %>%
  select(group, OD_Norm, line, row_id)

L4 = feconc_clean %>%
  filter(group == "Larvae (L4)") %>%
  select(c(PC0, PC2, PC4, PC6, PC8, PC10)) %>%
  slice(1) %>%
  pivot_longer(cols = c(PC0, PC2, PC4, PC6, PC8, PC10), names_to = "conc", values_to = "OD") %>%
  separate(conc, into = c("letter", "conc"), sep = "PC") %>%
  mutate(conc = as.numeric(conc)) %>%
  select(!letter)

L4_Values = feconc_clean %>%
  filter(group == "Larvae (L4)") %>%
  select(group, OD_Norm, line, row_id)

Pupae = feconc_clean %>%
  filter(group == "Pupae") %>%
  select(c(PC0, PC2, PC4, PC6, PC8, PC10)) %>%
  slice(1) %>%
  pivot_longer(cols = c(PC0, PC2, PC4, PC6, PC8, PC10), names_to = "conc", values_to = "OD") %>%
  separate(conc, into = c("letter", "conc"), sep = "PC") %>%
  mutate(conc = as.numeric(conc)) %>%
  select(!letter)

Pupae_Values = feconc_clean %>%
  filter(group == "Pupae") %>%
  select(group, OD_Norm, line, row_id)

LarvalWater = feconc_clean %>%
  filter(group == "LW Day 1" | group == "LW Day 7" ) %>%
  select(c(PC0, PC2, PC4, PC6, PC8, PC10)) %>%
  slice(1) %>%
  pivot_longer(cols = c(PC0, PC2, PC4, PC6, PC8, PC10), names_to = "conc", values_to = "OD") %>%
  separate(conc, into = c("letter", "conc"), sep = "PC") %>%
  mutate(conc = as.numeric(conc)) %>%
  select(!letter)

LarvalWater_Values = feconc_clean %>%
  filter(group == "LW Day 1" | group == "LW Day 7" ) %>%
  select(group, OD_Norm, line, row_id)

Adults = feconc_clean %>%
  filter(group == "Adult Male" | group == "Adult Female" ) %>%
  select(c(PC0, PC2, PC4, PC6, PC8, PC10)) %>%
  slice(1) %>%
  pivot_longer(cols = c(PC0, PC2, PC4, PC6, PC8, PC10), names_to = "conc", values_to = "OD") %>%
  separate(conc, into = c("letter", "conc"), sep = "PC") %>%
  mutate(conc = as.numeric(conc)) %>%
  select(!letter)

Adults_Values = feconc_clean %>%
  filter(group == "Adult Male" | group == "Adult Female" ) %>%
  select(group, OD_Norm, line, row_id)

bf_adults = feconc_clean %>%
  filter(group == "Adult Female 24hr PBM") %>%
  select(c(PC0, PC2, PC4, PC6, PC8, PC10)) %>%
  slice(1) %>%
  pivot_longer(cols = c(PC0, PC2, PC4, PC6, PC8, PC10), names_to = "conc", values_to = "OD") %>%
  separate(conc, into = c("letter", "conc"), sep = "PC") %>%
  mutate(conc = as.numeric(conc)) %>%
  select(!letter)

bf_adults_values = feconc_clean %>%
  filter(group == "Adult Female 24hr PBM") %>%
  select(group, OD_Norm, line, row_id)

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
    current_row <- filtered_table_samples$row_id[i]
    # predict the value based on current OD valus 
    predicted_value <- predict(model, newdata = data.frame(OD = current_od))
    
    # add the result to the predicted values in the dataframe 
    predicted_values_df <- rbind(predicted_values_df, 
                                 data.frame(OD = current_od, 
                                            PredictedValue = predicted_value, 
                                            row_id = current_row))
  }
  
  predicted_iron_concentration <- filtered_table_samples %>%
    left_join(predicted_values_df, by = c("row_id" = "row_id")) %>%
    mutate(PredictedValueConv = ifelse(PredictedValue < 0, 0, PredictedValue))
  
  predicted_iron_concentration_uM <- predicted_iron_concentration %>%
    mutate(uM_Fe = (PredictedValueConv / 50) * 4) # dilution factor
  
  predicted_iron_concentration_nM <- predicted_iron_concentration_uM %>%
    mutate(nM_Fe = uM_Fe * 1000)
  
  predicted_iron_concentration_ng <- predicted_iron_concentration_nM %>%
    mutate(ng_Fe = (nM_Fe * 55.845)) # use molecular weight of iron to calculate ng/Sample
  
  predicted_iron_concentration_mg <- predicted_iron_concentration_uM %>%
    mutate(ug_Fe = (uM_Fe * 55.845))
  
  predicted_iron_concentration_mg_distinct <- predicted_iron_concentration_ng %>% 
    distinct(OD_Norm, .keep_all = TRUE)
  
  return(predicted_iron_concentration_ng)
  
}

L1Group_Fe = PredictFeConc(L1Group_Values, L1Group)
L2_L3_Fe = PredictFeConc(L2_L3_Values, L2_L3)
L4_Fe = PredictFeConc(L4_Values, L4)
Pupae_Fe = PredictFeConc(Pupae_Values, Pupae)
LarvalWater_Fe = PredictFeConc(LarvalWater_Values, LarvalWater)
Adults_Fe = PredictFeConc(Adults_Values, Adults)
bf_adults_Fe = PredictFeConc(bf_adults_values, bf_adults)

final_table = rbind(L1Group_Fe, L2_L3_Fe, L4_Fe, Pupae_Fe, LarvalWater_Fe, Adults_Fe, bf_adults_Fe)

final_table %>%
  mutate(group = factor(group, levels = sample_order)) %>%
  ggplot() + 
  aes(x = group, y = nM_Fe) + 
  geom_boxplot(outlier.alpha = 0) + 
  geom_jitter(alpha=0.3, width = 0.1) +
  theme_bw(base_size = 20) + 
  xlab("Life stage") + 
  ylab("Fe (nM)") + 
  theme(axis.text.x = element_text(angle = 15, hjust = 1)) + 
  facet_wrap(facets = "line", scales = "fixed", nrow = 3)
  
#facet_grid(rows = "line")

ggsave("201023_FeConcentrationLifestages.pdf", width = 10, height =10)

# get statistics 
final_table %>%
  filter(group == "Adult Female" | group == "Larvae (L4)") %>%
  group_by(group, line) %>%
  summarise(avg_ng_Fe = mean(ng_Fe))

feconc_clean %>%
  group_by(group, NC) %>%
  summarise(average_OD = mean(NC))

feconc_clean %>%
  select(group, OD_Norm) %>%
  print(n = 400)

# anova - test sig differences between groups
fe_conc_anova = kruskal.test(PredictedValueConv ~ group, data = final_table)
summary(fe_conc_anova)

final_table %>%
  mutate(group = factor(group, levels = sample_order)) %>%
  ggplot() + 
  aes(x = group, y = PredictedValueConv) + 
  geom_boxplot(outlier.alpha = 0) + 
  geom_jitter(alpha=0.6, width = 0.1) +
  theme_bw(base_size = 25) + 
  xlab("Life stage") + 
  ylab("ng Fe / sample") + 
  theme(axis.text.x = element_text(angle = 15, hjust = 1)) 
  #facet_wrap(facets = "line", scales = "free_y", nrow = 1) 
  #facet_grid(rows = "line")


# function to fit glm and compare L4 values to all other values
glm_L4_to_Pupae <- function(data, target){
  
  sample_order2 = c("Larvae (L4)", "LW Day 1", "Larvae (L1)", "Larvae (L2/3)", "Pupae", "Adult Male", "Adult Female", "Adult Female 24hr PBM")
  
  data_filt = data %>%
    mutate(group = factor(group, levels = sample_order2)) %>%
    filter(line == target)
  
  glm = glm(PredictedValueConv ~ group, data = data_filt)
  
  return(glm)
  
}

LVP_glm = glm_L4_to_Pupae(final_table, "LVP")
print(summary.glm(LVP_glm), show.residuals = T)
simulationLVP <- simulateResiduals(LVP_glm)
plot(simulationLVP)


galveston_glm = glm_L4_to_Pupae(final_table, "Galveston")
print(summary.glm(galveston_glm), show.residuals = T)
simulationGal <- simulateResiduals(galveston_glm)
plot(simulationGal)


Juchitan_glm = glm_L4_to_Pupae(final_table, "Juchitan")
print(summary.glm(Juchitan_glm), show.residuals = T)
simulationJuc <- simulateResiduals(Juchitan_glm)
plot(simulationJuc)

Thai1_glm = glm_L4_to_Pupae(final_table, "Thai1")
print(summary.glm(Thai1_glm), show.residuals = T)
simulationT1 <- simulateResiduals(Thai1_glm)
plot(simulationT1)

Thai2_glm = glm_L4_to_Pupae(final_table, "Thai2")
print(summary.glm(Thai2_glm), show.residuals = T)
simulationT2 <- simulateResiduals(Thai2_glm)
plot(simulationT2)

### Interpreting the GLM model 
#### Estimate Std - for every unit of change, that is the unit increase or decrease 
#### Std error - how much fuzzyness around the slope
#### t-value - estimate / by the std error :
  #### smaller std error = less uncertainty 
  #### t-values that are larger, more confidence in 

#### Null-model : what are the results where we look at the ntercept to nothing 
#### Residual deviennce : take the predictors into account

### values are sig but looks like the data doesn't meet the assumptions of normality 


## Run non-parametric tests 

WilCox_Test <- function(data, target, group1, group2){
  
  sample_order2 = c("Larvae (L4)", "LW Day 1", "Larvae (L1)", "Larvae (L2/3)", "Pupae", "Adult Male", "Adult Female", "Adult Female 24hr PBM")
  
  data_filt = data %>%
    mutate(group = factor(group, levels = sample_order2)) %>%
    filter(line == target) %>%
    filter(group == group1 | group == group2) 
  
  wilcox = wilcox.test(PredictedValueConv ~ group, 
                       data = data_filt, 
                       alternative = "two.sided")
  
  return(wilcox)
  
}

WilCox_L4_to_Pupae(final_table, "LVP", "Larvae (L4)", "Pupae")
WilCox_L4_to_Pupae(final_table, "LVP", "Larvae (L4)", "Adult Male")
WilCox_L4_to_Pupae(final_table, "LVP", "Larvae (L4)", "Adult Female")

WilCox_L4_to_Pupae(final_table, "Galveston", "Larvae (L4)", "Pupae")
WilCox_L4_to_Pupae(final_table, "Galveston", "Larvae (L4)", "Adult Male")
WilCox_L4_to_Pupae(final_table, "Galveston", "Larvae (L4)", "Adult Female")

WilCox_L4_to_Pupae(final_table, "Juchitan", "Larvae (L4)", "Pupae")
WilCox_L4_to_Pupae(final_table, "Juchitan", "Larvae (L4)", "Adult Male")
WilCox_L4_to_Pupae(final_table, "Juchitan", "Larvae (L4)", "Adult Female")

WilCox_L4_to_Pupae(final_table, "Thai1", "Larvae (L4)", "Pupae")
WilCox_L4_to_Pupae(final_table, "Thai1", "Larvae (L4)", "Adult Male")
WilCox_L4_to_Pupae(final_table, "Thai1", "Larvae (L4)", "Adult Female")

WilCox_L4_to_Pupae(final_table, "Thai2", "Larvae (L4)", "Pupae")
WilCox_L4_to_Pupae(final_table, "Thai2", "Larvae (L4)", "Adult Male")
WilCox_L4_to_Pupae(final_table, "Thai2", "Larvae (L4)", "Adult Female")

plot_Fe = function(data, target, sigval1, sigval2, sigval3){
  p = data %>%
    filter(group != "LW Day 1") %>%
    filter(group != "LW Day 7") %>%
    filter(group != "Adult Female 24hr PBM") %>%
    filter(line == target) %>%
    mutate(group = factor(group, levels = sample_order)) %>%
    ggplot() + 
    aes(x = group, y = nM_Fe) + 
    geom_boxplot(outlier.alpha = 0) + 
    geom_jitter(alpha=0.3, width = 0.1, size = 4) +
    theme_bw(base_size = 30) + 
    xlab("Life stage") + 
    ylab("Total Fe (nM)") + 
    theme(axis.text.x = element_text(angle = 20, hjust = 1)) + 
    ylim(0, 320) + 
    geom_signif(comparisons = list(c("Larvae (L4)", "Adult Female")),
                annotations = c(sigval1), 
                y_position = 290, 
                textsize = 10, 
                vjust = 0.5) +
    geom_signif(comparisons = list(c("Larvae (L4)", "Adult Male")),
                annotations = c(sigval2), 
                y_position = 270, 
                textsize = 10, 
                vjust = 0.5) + 
    geom_signif(comparisons = list(c("Larvae (L4)", "Pupae")),
                annotations = c(sigval3), 
                y_position = 250, 
                textsize = 10, 
                vjust = 0.5) 
    
    return(p)
}

Galveston_Line = plot_Fe(final_table, "Galveston", "***", "***", "**")
ggsave("Figure_4B_GalvestonFeConcentrationLines.pdf", width = 10, height = 10)

LVP_Fe = plot_Fe(final_table, "LVP", "***", "***", "*")
ggsave("Figure_4B_LVPFeConcentrationLines.pdf", width = 10, height = 10)

Juchitan_Fe = plot_Fe(final_table, "Juchitan", "***", "***", "n.s")
ggsave("Figure_4C_JuchitanFeConcentrationLines.pdf", width = 10, height = 10)

Thai1_Fe = plot_Fe(final_table, "Thai1", "***", "***", "***")
ggsave("Figure_4D_Thai1FeConcentrationLines.pdf", width = 10, height = 10)

Thai2_Fe = plot_Fe(final_table, "Thai2", "***", "**", "n.s")
ggsave("Figure_4E_Thai2FeConcentrationLines.pdf", width = 10, height = 10)

# Measure fish food!
fish_food = feconc_clean %>%
  filter(group == "fish_food")

fish_food = feconc_clean %>%
  filter(group == "fish_food") %>%
  select(c(PC0, PC2, PC4, PC6, PC8, PC10)) %>%
  slice(1) %>%
  pivot_longer(cols = c(PC0, PC2, PC4, PC6, PC8, PC10), names_to = "conc", values_to = "OD") %>%
  separate(conc, into = c("letter", "conc"), sep = "PC") %>%
  mutate(conc = as.numeric(conc)) %>%
  select(!letter)

fish_food_values = feconc_clean %>%
  filter(group == "fish_food") %>%
  select(group, OD_Norm, line)

fish_food_fe = PredictFeConc(fish_food_values, fish_food)

fish_food_fe = fish_food_fe %>%
  mutate(ng_Fe_per_g = ((ng_Fe * 2) * 100))

format_labels <- function(x) {
  x <- gsub("_", "\n", x)  # Replace underscores with newlines
  x
}

fish_food_fe %>%
  mutate(group = "Tetramin 10%") %>%
  ggplot() + 
  aes(x = group, y = ng_Fe_per_g) + 
  geom_boxplot(outlier.alpha = 0) + 
  geom_jitter(alpha=0.6, width = 0.1) +
  theme_bw(base_size = 30) + 
  xlab("Diet") + 
  ylab("ng Fe / mL tetramin") + 
  ylim(0, 160000) + 
  scale_x_discrete(labels = function(x) gsub("Tetramin 10%", "Tetramin\n10%", x))  # Apply multi-line label
  

ggsave("Figure_4H_FishFoodFeConcentrationLines.pdf", width = 3.3, height = 13)

# overall distribution with males
final_table %>%
  filter(group != "Adult Female" & group != "LW Day 1" & group != "LW Day 7") %>%
  mutate(group = factor(group, levels = sample_order)) %>%
  rename("Group" = group) %>%
  ggplot() + 
  aes(x = Group, y =  ng_Fe, group = 1) + 
  geom_smooth(aes(group = line, col = line), alpha = 0.2, method = "loess", level = 0.95) + 
  theme_bw(base_size = 30) + 
  xlab("Sample") + 
  ylab("ng Fe / sample") + 
  theme(axis.text.x = element_text(angle = 15, hjust = 1)) + 
  scale_color_manual(values=c("#5C5934", "#996E3F", "#E05040", "#B6ED4F", "#61F7BB"))

ggsave("FeConcentrationMosLinesOverallMales.pdf", width = 12, height = 10)

# overall distribution with females
final_table %>%
  filter(group != "Adult Male" & group != "LW Day 1" & group != "LW Day 7") %>%
  mutate(group = factor(group, levels = sample_order)) %>%
  rename("Group" = group) %>%
  ggplot() + 
  aes(x = Group, y =  ng_Fe, group = 1) + 
  geom_smooth(aes(group = line, col = line), alpha = 0.2, method = "loess", level = 0.95) + 
  theme_bw(base_size = 30) + 
  xlab("Sample") + 
  ylab("ng Fe / sample") + 
  theme(axis.text.x = element_text(angle = 15, hjust = 1)) + 
  scale_color_manual(values=c("#5C5934", "#996E3F", "#E05040", "#B6ED4F", "#61F7BB"))

ggsave("FeConcentrationMosLinesOverallFemales.pdf", width = 12, height = 10)

# plot larval water 
final_table %>%
  filter(group == "LW Day 1" | group == "LW Day 7") %>%
  mutate(group = factor(group, levels = sample_order)) %>%
  ggplot() + 
  aes(x = group, y =  nM_Fe) + 
  geom_boxplot() +
  geom_jitter(alpha=0.6, width = 0.1) +
  theme_bw(base_size = 30) + 
  xlab("Larval \nwater") + 
  ylab("Total Fe (nM)") + 
  scale_x_discrete(labels = c("D1", "D7"))   # Apply multi-line label


ggsave("Figure_4G_FeConcentrationMosLinesLarvalWater.pdf", width = 3.3, height = 13)

# plot differences between males and females 
final_table %>%
  filter(group == "Adult Female" | group == "Adult Male") %>%
  mutate(group = factor(group, levels = sample_order)) %>%
  ggplot() + 
  aes(x = group, y = nM_Fe) + 
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(alpha=0.6, width = 0.1) +
  theme_bw(base_size = 30) + 
  xlab("Sex") + 
  ylab("Total Fe (nM)") + 
  scale_x_discrete(labels = c("M", "F")) 


ggsave("Figure_4_F_FeConcentrationMosAdultsMaleFemales.pdf", width = 3.3, height = 10)

# Stats!
## what was the change between larvae L1 to Larvae L4 
final_table %>%
  filter(group == "Larvae (L1)" | group == "Larvae (L4)") %>%
  select(group, nM_Fe) %>%
  group_by(group) %>%
  summarise(average_FeConc = mean(nM_Fe), 
            sd(nM_Fe)) 

foldchange(117, 11.2)

# concentration in mosquito lines of larvae 
final_table %>%
  filter(group == "Larvae (L1)" | group == "Larvae (L4)") %>%
  select(group, line, nM_Fe) %>%
  group_by(group, line) %>%
  summarise(average_FeConc = mean(nM_Fe), 
            sd(nM_Fe)) 

# concentration between L1 to L2/3 
final_table %>%
  filter(group == "Larvae (L1)" | group == "Larvae (L2/3)") %>%
  select(group, line, nM_Fe) %>%
  group_by(group, line) %>%
  summarise(average_FeConc = mean(nM_Fe), 
            sd(nM_Fe)) 

# concentration after emergence overall
final_table %>% 
  filter(group == "Larvae (L4)" | group == "Pupae") %>%
  select(group, nM_Fe) %>%
  group_by(group) %>%
  summarise(average_FeConc = mean(nM_Fe), 
            sd(nM_Fe)) 

foldchange(117, 62)

# concentration after emergence overall
final_table %>% 
  filter(group == "Larvae (L4)" | group == "Adult Female" | group == "Adult Male") %>%
  select(group, nM_Fe) %>%
  group_by(group) %>%
  summarise(average_FeConc = mean(nM_Fe), 
            sd(nM_Fe)) 

foldchange(117, 20.2) # fold change females 
foldchange(117, 5.10) # fold change males 


# concentration after emergence per line
final_table %>% 
  filter(group == "Larvae (L4)" | group == "Pupae") %>%
  select(group, line, nM_Fe) %>%
  group_by(group, line) %>%
  summarise(average_FeConc = mean(nM_Fe), 
            sd(nM_Fe)) 

# comparing males and females 
final_table %>% 
  filter(group == "Adult Female" | group == "Adult Male") %>%
  select(group, nM_Fe, line) %>%
  group_by(group, line) %>%
  summarise(average_FeConc = mean(nM_Fe), 
            sd(nM_Fe)) 

foldchange(20.2, 5.10)

# compare larval water
final_table %>% 
  filter(group == "LW Day 1" | group == "LW Day 7") %>%
  select(group, nM_Fe) %>%
  group_by(group) %>%
  summarise(average_FeConc = mean(nM_Fe), 
            sd(nM_Fe)) 

foldchange(1.28, 43.4)

fish_food_fe %>%
  summarise(mean(ng_Fe_per_g))

