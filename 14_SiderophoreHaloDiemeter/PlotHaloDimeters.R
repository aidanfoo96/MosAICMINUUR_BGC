library(tidyverse)

halo <- read_csv("270324_halo_perimeter_measurements.csv")

# plot  
halo %>%
  ggplot() + 
  aes(x = reorder(strain, colony_halo_ratio, FUN = median, decreasing = TRUE), y = colony_halo_ratio) +
  geom_boxplot() + 
  geom_point(alpha = 0.5) + 
  theme_bw(base_size = 25) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1)) + 
  xlab("Strain") + 
  ylab("Culture/halo ratio")

ggsave(filename = "260824_siderophore_halo_diemeter_measurements.pdf", heigh = 5, width = 22)
