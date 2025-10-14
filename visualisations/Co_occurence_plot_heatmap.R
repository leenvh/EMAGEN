
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#required packages
library(tidyverse)
library(reshape2)
library(ggtext)  
library(scales)  
library(readxl)

#  data 
dat_API=read.csv("data/merged_results_api_cnd_annot_f.csv")
dat_geo=read.csv("data/merged_results_geography_cnd_annot_LA_f.csv")
dat_crt=read.csv("data/PfcrtMut_stratified_Cooccurance_result_cnd_annot_f.csv")


dat_API$cat <- "API"
colnames(dat_API)[1] <- "strata"
dat_API <- dat_API%>%
  select(,c(1:9,11))%>%
  mutate(strata = ifelse(strata == "Low Transmission (0-10)", "Low", strata))%>%
  mutate(strata = ifelse(strata == "Moderate Transmission (10-50)", "Moderate", strata))%>%
  mutate(strata = ifelse(strata == "High Transmission (>=50)", "High", strata))
         
dat_geo$cat <- "Geo cluster"
colnames(dat_geo)[1] <- "strata"
dat_geo <- dat_geo[,c(1:9,12)]
dat_crt$cat <- "PfCRT variant prevalence"
colnames(dat_crt)[1] <- "strata"
dat_crt <- dat_crt[,c(1:9,13)]

summary <- rbind(dat_API,dat_geo,dat_crt) 
write.csv(summary,"supp_table5.csv")
summary<-summary %>%
  mutate(Variant_pair = paste0(Mutation1,"_",Mutation2)) %>%
  filter(Mutation1 != "Pfk13_P441L" & Mutation2 != "Pfk13_P441L") %>%
  filter(Mutation1 != "Pfk13_A675V" & Mutation2 != "Pfk13_A675V") %>%
  filter(Mutation1 != "Pfk13_R622I" | Mutation2 != "Pfmdr1_N86Y") %>%
  filter(Mutation1 != "Pfcrt_I356T" | Mutation2 != "Pfmdr1_N86Y") %>%
  filter(Mutation1 != "Pfcrt_I356T" | Mutation2 != "Pfmdr1_Y184F") %>%
  mutate(Variant_pair = ifelse(Variant_pair == "Pfmdr1_Y184F_Pfk13_R622I", 
                               "Pfk13_R622I_Pfmdr1_Y184F", 
                               Variant_pair))%>%
  mutate(Mutation1=ifelse(Variant_pair == "Pfk13_R622I_Pfmdr1_Y184F", 
                          "Pfk13_R622I", 
                          Mutation1)) %>%
  mutate(Mutation2=ifelse(Variant_pair == "Pfk13_R622I_Pfmdr1_Y184F", 
                          "Pfmdr1_Y184F", 
                          Mutation2))

# Convert  the data from wide to long format  and Add a tiered significance column based on p-value thresholds
#summary_long <- summary %>%
#select(strata, Mutation1, Mutation2, Co.occurrence.Count, P.value,N) %>%
#melt(id.vars = c("strata","Mutation1", "Mutation2")) 


variant_order <- c(
  # Group 1: Pfk13_R622I pairs
  "Pfcrt_M74I_Pfk13_R622I",
  "Pfcrt_N75E_Pfk13_R622I",
  "Pfcrt_K76T_Pfk13_R622I",
  "Pfcrt_A220S_Pfk13_R622I",
  "Pfcrt_Q271E_Pfk13_R622I",
  "Pfcrt_N326S_Pfk13_R622I",
  "Pfcrt_I356T_Pfk13_R622I",
  "Pfcrt_R371I_Pfk13_R622I",
  
  # Group 2: Pfmdr1_Y184F pairs
  "Pfk13_R622I_Pfmdr1_Y184F",
  "Pfcrt_M74I_Pfmdr1_Y184F",
  "Pfcrt_N75E_Pfmdr1_Y184F",
  "Pfcrt_K76T_Pfmdr1_Y184F",
  "Pfcrt_A220S_Pfmdr1_Y184F",
  "Pfcrt_Q271E_Pfmdr1_Y184F",
  "Pfcrt_N326S_Pfmdr1_Y184F",
  "Pfcrt_I356T_Pfmdr1_Y184F",
  "Pfcrt_R371I_Pfmdr1_Y184F",
  
  # Group 3: Pfmdr1_N86Y pairs
  "Pfk13_R622I_Pfmdr1_N86Y",
  "Pfcrt_M74I_Pfmdr1_N86Y",
  "Pfcrt_N75E_Pfmdr1_N86Y",
  "Pfcrt_K76T_Pfmdr1_N86Y",
  "Pfcrt_A220S_Pfmdr1_N86Y",
  "Pfcrt_Q271E_Pfmdr1_N86Y",
  "Pfcrt_N326S_Pfmdr1_N86Y",
  "Pfcrt_I356T_Pfmdr1_N86Y",
  "Pfcrt_R371I_Pfmdr1_N86Y"
)

summary$Variant_pair <- factor(summary$Variant_pair, levels = rev(variant_order))
summary <- summary %>%
  filter(!is.na(Variant_pair))

strata_order <- c(
  "Low",
  "Moderate ",
  "Moderate",
  "High",
  "1",
  "2",
  "3")

summary$strata <- factor(summary$strata, levels = strata_order)



#  heatmap plot

# Create formatted labels showing only Mutation1 with italicized gene names
summary <- summary %>%
  mutate(Mutation1_label = case_when(
    # Pfcrt variants
    grepl("^Pfcrt_", Mutation1) ~ paste0("*Pf*CRT ", gsub("Pfcrt_", "", Mutation1)),
    # Pfmdr1 variants
    grepl("^Pfmdr1_", Mutation1) ~ paste0("*Pf*MDR1 ", gsub("Pfmdr1_", "", Mutation1)),
    # Pfk13 variants
    grepl("^Pfk13_", Mutation1) ~ paste0("*Pf*K13 ", gsub("Pfk13_", "", Mutation1)),
    TRUE ~ Mutation1
  )) %>%
  mutate(Mutation2_label = case_when(
    # Pfmdr1 variants
    grepl("^Pfmdr1_", Mutation2) ~ paste0("*Pf*MDR1 ", gsub("Pfmdr1_", "", Mutation2)),
    # Pfk13 variants
    grepl("^Pfk13_", Mutation2) ~ paste0("*Pf*K13 ", gsub("Pfk13_", "", Mutation2)),
    TRUE ~ Mutation2
  ))%>%
  mutate(cat_label = case_when(
    grepl("^Pf", cat) ~ paste0("*Pf*", gsub("Pf", "", cat)),
    TRUE ~ cat 
  ))


cat_order <- c(
  "Geo cluster",
  "*Pf*CRT variant prevalence",
  "API"
)

summary$cat_label <- factor(summary$cat_label, levels = cat_order)

Mutation1_label_order <- c(
  "*Pf*K13 R622I",
  "*Pf*CRT M74I",
  "*Pf*CRT N75E",
  "*Pf*CRT K76T",
  "*Pf*CRT A220S",
  "*Pf*CRT Q271E",
  "*Pf*CRT N326S",
  "*Pf*CRT I356T",
  "*Pf*CRT R371I")


summary$Mutation1_label <- factor(summary$Mutation1_label, 
                                  levels = rev(Mutation1_label_order))

summary<-summary[-which(summary$Co.occurrence.Count %in% c("0","1")),]

# Plot
cooccur_plot <- ggplot(summary, aes(x = strata, y = Mutation1_label)) +
  geom_point(aes(size = Co.occurrence.Count, fill = P.value), shape = 21, color = "black", stroke = 0.5) +
  facet_grid(Mutation2_label ~ cat_label, scales = "free") +
  scale_fill_gradientn(
    colors = c("#B2182B", "#D73027", "grey90", "grey90"),
    values = c(0, 0.4, 0.76, 1),  
    trans = "log10",
    limits = c(0.00001, 1), 
    breaks = c(0.0001, 0.001, 0.01, 0.05, 0.5),
    labels = c("0.0001", "0.001", "0.01", "0.05", "0.5"),
    name = "p-value"
  ) +
  scale_size_continuous(
    range = c(1, 10),
    name = "Number of co-occurences"
  ) +
  labs(
    x = "",
    y = ""
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_markdown(),     
    strip.text.y = element_markdown(),     
    strip.text.x = element_markdown(),  
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    plot.caption = element_text(hjust = 0.5, size = 10)
  )



cooccur_plot



# Save the plot
ggsave("co_occur_plot_stratified.pdf", cooccur_plot,
       width = 28, height = 22, units = "cm", dpi = 600)




# Heatmap
cooccur_plot_heatmap <- ggplot(summary, aes(x = strata, y = Mutation1_label)) +
  geom_tile(aes(fill = P.value), color = "white", linewidth = 0.5) +  # creates heatmap cells
  geom_text(aes(label = Co.occurrence.Count), size = 3) +  # adds count numbers
  facet_grid(Mutation2_label ~ cat_label, scales = "free") +
  scale_fill_gradientn(
    colors = c("#B2182B", "#D73027", "grey90", "grey90"),
    values = c(0, 0.4, 0.76, 1),  
    trans = "log10",
    limits = c(0.00001, 1), 
    breaks = c(0.0001, 0.001, 0.01, 0.05, 0.5),
    labels = c("0.0001", "0.001", "0.01", "0.05", "0.5"),
    name = "p-value"
  ) +
  labs(
    x = "",
    y = "",
    caption = "Numbers indicate co-occurrence counts"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_markdown(),     
    strip.text.y = element_markdown(),     
    strip.text.x = element_markdown(),  
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    plot.caption = element_text(hjust = 0.5, size = 10)
  )

cooccur_plot_heatmap

ggsave("co_occur_heatmap_stratified.pdf", cooccur_plot_heatmap,
       width = 28, height = 22, units = "cm", dpi = 600)


