library(tidyverse)

iron_data <- read_tsv("maldonado_adapted_data.txt")

iron_data_long <- iron_data %>%
  pivot_longer(cols = !c(sampleID, strain), names_to = "replicate", values_to = "iron_conc" )

iron_data_long$sampleID <- factor(iron_data_long$sampleID, levels = c("egg", "L1", "L2", 
                                                                      "L3", "L4", "pupae", 
                                                                      "male", "female"))
iron_data_long %>%
  ggplot() + 
  aes(x = sampleID, y = iron_conc) + 
  geom_boxplot() + 
  theme_bw(base_size = 20) + 
  xlab("Lifestage") + 
  ylab("mg Fe / g dry weight") + 
  facet_wrap(facets = "strain") + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
  
ggsave(filename = "maldonado_Fe_Conc_adapted.png", height = 8, width = 10)
