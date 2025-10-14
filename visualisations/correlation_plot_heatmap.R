
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#required packages
library(tidyverse)
library(reshape2)
library(ggtext)  
library(scales)  

#  data 
dat_cor=readRDS("data/corr_data_output.rds")

# Convert  the data from wide to long format  and Add a tiered significance column based on p-value thresholds
dat_cor_long <- dat_cor %>%
  select(haplotype, pv_pf_rho, api_rho, pv_pf_p, api_p,hap_freq) %>%
  melt(id.vars = c("haplotype", "pv_pf_p", "api_p","hap_freq")) %>%
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


dat_cor_long$variable <- factor(dat_cor_long$variable, levels = c("API", "Pv proportion"))
dat_cor_long$haplotype <- factor(dat_cor_long$haplotype, levels = rev(c(
  "AIRNI",
  "AICNI",
  "ISGEAA",
  "ISAKAA (WT)",
  "AIRNI / ISGEAA",
  "AICNI / ISGEAA",
  "AICNI / ISAKAA",
  "CVIET",
  "CVMNK (WT)",
  "NFSND"
)))


#  heatmap plot

corr_plot <- ggplot(dat_cor_long, aes(x = variable, y = haplotype)) +
  geom_point(aes(size = hap_freq, fill = value), shape = 21, color = "black", stroke = 0.5) +
  geom_text(aes(label = signif), vjust = -0.9, size = 4, color = "black") +  # add significance stars above the bubbles
  scale_fill_gradient2(
    low = "#D73027",   # red
    mid = "white",   # yellow/white
    high = "#1A9850",  # green
    midpoint = 0,
    limits = c(-1, 1),
    name = expression(~rho)
  ) +
  scale_size_continuous(
    range = c(3, 10),   # size of bubbles
    name = "Haplotype frequency (%)",
    breaks =c(12,25,50,75)
  ) +
  labs(
    x = "",
    y = "Haplotype",
    caption = "Bubble size: Frequency (%)\nFill color: Spearman's ρ\nAsterisks: Significance"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    plot.caption = element_text(hjust = 0.5, size = 10)
  )


corr_plot
# Save the plot
ggsave("corr_plot_corrected.pdf", corr_plot,
       width = 17.5, height = 17, units = "cm", dpi = 600)


