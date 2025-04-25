# Load required packages
library(tidyverse)
library(reshape2)
library(ggtext)  
library(scales)  

#  data 
dat_cor=readRDS("data/corr_data_output.rds")

# Convert  the data from wide to long format  and Add a tiered significance column based on p-value thresholds
dat_cor_long <- dat_cor %>%
  select(haplotype, pv_pf_rho, api_rho, pv_pf_p, api_p) %>%
  melt(id.vars = c("haplotype", "pv_pf_p", "api_p")) %>%
  mutate(
    signif = case_when(
      variable == "pv_pf_rho" & pv_pf_p < 0.001 ~ "***",
      variable == "pv_pf_rho" & pv_pf_p < 0.01 ~ "**",
      variable == "pv_pf_rho" & pv_pf_p < 0.05 ~ "*",
      variable == "api_rho" & api_p < 0.001 ~ "***",
      variable == "api_rho" & api_p < 0.01 ~ "**",
      variable == "api_rho" & api_p < 0.05 ~ "*",
      TRUE ~ ""
    ),
    variable = factor(variable, 
                      levels = c("pv_pf_rho", "api_rho"),
                      labels = c("Pv proportion", "API"))
  )


#  heatmap plot
corr_plot <- ggplot(dat_cor_long, aes(x = variable, y = haplotype, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = signif), color = "black", size = 5) +
  scale_fill_gradient2(
    low = "#E15759", 
    mid = "white", 
    high = "#59A14F",
    limits = c(-1, 1), 
    name = expression(~rho)
  ) +
  labs(
    x = "", 
    y = "Haplotype",
    caption = "Significance: *** p < 0.001; ** p < 0.01; * p < 0.05"  
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.caption = element_text(hjust = 0.5, size = 10),  
    plot.caption.position = "panel"  
  )


corr_plot
# Save the plot
ggsave("corr_plot.tiff", corr_plot,
       width = 18, height = 12, units = "cm", dpi = 600, compression = "lzw")


