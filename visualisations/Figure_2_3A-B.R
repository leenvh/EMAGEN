       

   # Load necessary libraries
# Install remotes package if not already installed
library(sf)
library(stringr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggforce)
library(purrr)
library(grid)
library(tmap)
library(cowplot)
library(scatterpie)
library(gridExtra)
library(ggtext)
library(waffle)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyr)
library(readxl)
library(ggforce) 

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load and Read the data
Africa <-st_read("../Data/afric.shp")
border <-st_read("../Data/border.shp")
Et_region <- st_read("../Data/Et_region.shp")
Haplotype <- read_excel("../Data/All_haplotype.xlsx")
Haplotype_all <- st_read("../Data/Haplotype_all_corrected.shp")
District <- st_read("../Data/District.shp")

# pfcrt combined plot ##########
# --- 1. Load and Prepare the Data --- 
pfcrt_data_long <- Haplotype %>%
  filter(Gene == "Pfcrt") %>%
  select(District, haplotype, prevalence) %>%
  group_by(District, haplotype) %>%
  summarise(Prevalence = sum(prevalence), .groups = "drop") %>%
  group_by(District) %>%
  mutate(Percentage = Prevalence / sum(Prevalence) * 100,
        ymax = cumsum(Percentage),
        ymin = c(0, head(ymax, n = -1)))
pfcrt_data_long$haplotype <- factor(pfcrt_data_long$haplotype, levels = c("CVMNK", "CVIET"))

# --- 2. Define Colors for Pfcrt Haplotypes ---
#colors_pfcrt <- c("CVIET" = "#963696", "CVMNK" = "#e4b4e4")
colors_pfcrt <- c(
  "CVMNK" = "#9ed1ed",   
  "CVIET" = "#2471a3" 
)

# --- 3. Create Pie Chart Plot for Pfcrt ---
pie_chart_plot_pfcrt <- ggplot(pfcrt_data_long, aes(x = "", y = Percentage, fill = haplotype)) +
     geom_bar(stat = "identity", width = 1, color = "white", size = 0.3) +
     coord_polar("y", start = 0) +
     facet_wrap(~ District, ncol = 5) +
     theme_void() +
     theme(
       legend.position = "none",
       strip.text = element_text(size = 8, face = "bold"),
       plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
       panel.spacing = unit(0.8, "lines"),
       plot.background = element_rect(fill = "transparent", color = NA),
       plot.title = element_text(size = 10)  # Minimized title size
     ) +
     labs(title = "District-Level Prevalence of Pfcrt Haplotypes") +
     scale_fill_manual(values = colors_pfcrt,breaks = c("CVMNK", "CVIET"))
   
# --- 4. Create Waffle Chart Data for Pfcrt ---
overall_data_pfcrt <- c(
  CVMNK = 191,
  CVIET = 301
)
   
# --- 5. Create Waffle Chart (Pfcrt) ---
waffle_chart_pfcrt <- waffle(
 parts = overall_data_pfcrt,
 rows = 10,  # You can increase rows for more squares, or adjust for desired visual
 colors = colors_pfcrt,
 title = "Overall Prevalence of Pfcrt Haplotypes (N=492)",
 xlab = "1 square = 1 sample",
 legend_pos = "bottom",
 size = 0.5  # Adjust size for square size
) +
 theme(
   plot.background = element_rect(fill = "transparent", color = NA),
   plot.title = element_text(face = "bold", size = 16),           # Bold and larger title
   legend.title = element_text(face = "bold", size = 14),         # Bold and larger legend title
   legend.text = element_text(face = "bold", size = 12),          # Bold and larger legend labels
   axis.title.x = element_text(face = "bold", size = 12),         # Bold and larger x-axis label
   panel.spacing = unit(0, "lines")                               # Remove spacing between squares
 )
   
# --- 6. Combine Pie Chart and Waffle Chart ---
combined_plot_pfcrt <- grid.arrange(waffle_chart_pfcrt, pie_chart_plot_pfcrt, ncol = 1, heights = c(1, 1.5))  # Adjust height ratio as needed

# --- 7. Save Combined Plot ---
ggsave("../output/combined_pfcrt_visualization.png", combined_plot_pfcrt, bg = "transparent", width = 15, height = 12, dpi = 600)

# Display the combined plot
print(combined_plot_pfcrt)



# pfdhfr combined plot ##########
district_data_pfdhfr <- data.frame(
  District = c("Amibara", "Batu", "Dasenech", "Dera", "Fedis",
  "Fentale Metehara", "Gambella zuria", "Gondar", "Guba",
  "K.Humera", "Metema", "Mizan", "Pawe", "Selamago",
  "Tanqua Abergelle"),
  AICNI = c(13, 1, 7, 26, 4, 10, 23, 19, 9, 26, 23, 2, 5, 10, 16),
  AIRNI = c(18, 34, 26, 1, 39, 32, 18, 2, 2, 14, 7, 39, 5, 5, 3),
  ANRNI = c(0, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  ANCSI = c(0, 0, 2, 0, 0, 0, 0, 0, 22, 0, 3, 0, 16, 3, 0)
)
   
# Reshape district data for pie charts
data_long_pfdhfr <- district_data_pfdhfr %>%
  pivot_longer(cols = -District, names_to = "Haplotype", values_to = "Count") %>%
  group_by(District) %>%
  mutate(Percentage = Count / sum(Count) * 100,
         ymax = cumsum(Percentage),
         ymin = c(0, head(ymax, n = -1)))
   
# --- 2. Overall Haplotype Data for Waffle Chart (Pfdhfr) ---
overall_data_pfdhfr <- c(
  ANCSI = 46,
  ANRNI = 5,
  AICNI = 194,
  AIRNI = 245
)
   
# --- 3. Consistent Colors (Pfdhfr) ---
#colors_pfdhfr <- c("AICNI" = "#844da3", "AIRNI" = "#896deb", "ANCSI" = "#e4b4e4", "ANRNI" = "#ff6f4b")
colors_pfdhfr <- c(
  "ANCSI" = "#9ed1ed",  
  "ANRNI" = "#2471a3",   
  "AICNI" = "#cd5e08",
  "AIRNI" = "#844da3"
)
   
# --- 4. Create Pie Chart Plot (Pfdhfr) ---
pie_chart_plot_pfdhfr <- ggplot(data_long_pfdhfr, aes(x = "", y = Percentage, fill = Haplotype)) +
  geom_bar(stat = "identity", width = 1, color = "white", size = 0.3) +
  coord_polar("y", start = 0) +
  facet_wrap(~ District, ncol = 5) +
  theme_void() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 8, face = "bold"),
    plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
    panel.spacing = unit(0.8, "lines"),
    plot.background = element_rect(fill = "transparent", color = NA),
    plot.title = element_text(size = 8)  # Minimized title size
  ) +
  scale_fill_manual(values = colors_pfdhfr)
   

# --- 5. Create Waffle Chart (Pfdhfr) ---
waffle_chart_pfdhfr <- waffle(
  parts = overall_data_pfdhfr,
  rows = 10,  # Keep the number of rows as is
  colors = colors_pfdhfr,
  title = "Overall Prevalence of Pfdhfr Haplotypes (N=490)",
  xlab = "1 square = 1 sample",
  legend_pos = "bottom",
  size = 0.3  # Reduced square size to stretch the chart horizontally
) +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    plot.title = element_text(face = "bold", size = 18),           # Increased title size
    legend.title = element_text(face = "bold", size = 14),         # Bold and larger legend title
    legend.text = element_text(face = "bold", size = 12),          # Bold and larger legend labels
    axis.title.x = element_text(face = "bold", size = 12),         # Bold and larger x-axis label
    panel.spacing = unit(0, "lines")                               # Remove spacing between squares
  )

# --- 6. Combine Plots (Pfdhfr) ---
combined_plot_pfdhfr <- grid.arrange(waffle_chart_pfdhfr, pie_chart_plot_pfdhfr, ncol = 1, heights = c(1, 1))  # Equal size for both plots
   
# --- 7. Save Combined Plot (Pfdhfr) ---
ggsave("../output/combined_pfdhfr_visualization.png", combined_plot_pfdhfr, bg = "transparent", width = 15, height = 12, dpi = 600)
   
# Display the combined plot
print(combined_plot_pfdhfr)




# pfdhps combined plot ##########
# --- 1. Load Data ---
# --- 1. Data Processing for Pfdhps ---
pfdhps_data_long <- Haplotype %>%
  filter(Gene == "Pfdhps") %>%
  select(District, haplotype, prevalence) %>%
  group_by(District, haplotype) %>%
  summarise(Prevalence = sum(prevalence), .groups = "drop") %>%
  group_by(District) %>%
  mutate(Percentage = Prevalence / sum(Prevalence) * 100,
         ymax = cumsum(Percentage),
         ymin = c(0, head(ymax, n = -1)))
   
# --- 2. Define Colors for Pfdhps Haplotypes ---
colors_pfdhps <- c(
 "SAKAA" = "#9ed1ed",   # Light blue
 "SGEAA" = "#2471a3",   # Orange
 "SGEGA" = "#cd5e08"    # Dark blue
)
 
# --- 3. Create Pie Chart Plot for Pfdhps ---
pie_chart_plot_pfdhps <- ggplot(pfdhps_data_long, aes(x = "", y = Percentage, fill = haplotype)) +
 geom_bar(stat = "identity", width = 1, color = "white", size = 0.3) +
 coord_polar("y", start = 0) +
 facet_wrap(~ District, ncol = 5) +
 theme_void() +
 theme(
   legend.position = "right",
   strip.text = element_text(size = 8, face = "bold"),
   plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
   panel.spacing = unit(0.8, "lines"),
   plot.background = element_rect(fill = "transparent", color = NA)
 ) +
 labs(title = "District-Level Prevalence") +
 scale_fill_manual(values = colors_pfdhps)
   
# --- 4. Create Waffle Chart Data for Pfdhps ---
overall_data_pfdhps <- c(
 SAKAA = 155,
 SGEAA = 330,
 SGEGA = 17
)

# --- 5. Create Waffle Chart (Pfdhps) ---
waffle_chart_pfdhps <- waffle(
 parts = overall_data_pfdhps,
 rows = 8,  # Keeping the same number of rows
 colors = colors_pfdhps,
 title = "Overall Prevalence of Pfdhps Haplotypes (N=502)",
 legend_pos = "bottom",
 size = 0.15  # Further reduced square size to make the chart longer
) +
 theme(
   plot.background = element_rect(fill = "transparent", color = NA),
   plot.title = element_text(face = "bold", size = 12),         # Increased title size
   legend.title = element_text(face = "bold", size = 10),       # Bold and larger legend title
   legend.text = element_text(face = "bold", size = 8),         # Bold and larger legend labels
   legend.spacing = unit(0, "cm"),                              # Remove spacing between legend items
   legend.margin = margin(t = 0, r = 0, b = -5, l = 0, unit = "pt"), # Adjust margin to reduce gap
   panel.spacing = unit(0, "lines")                             # Remove spacing between squares
 )

# --- 6. Combine Pie Chart and Waffle Chart ---
combined_plot_pfdhps <- grid.arrange(waffle_chart_pfdhps, pie_chart_plot_pfdhps, ncol = 1, heights = c(1, 2))

# --- 7. Save Combined Plot ---
ggsave("../output/combined_pfdhps_visualization.png", combined_plot_pfdhps, bg = "transparent", width = 15, height = 12, dpi = 600)

# Display the combined plot
print(combined_plot_pfdhps)
   




# pfmdr1 combined plot ##########
# --- 1. Data Processing for Pfmdr1 ---
pfmdr1_data_long <- Haplotype %>%
  filter(Gene == "Pfmdr1") %>%
  select(District, haplotype, prevalence) %>%
  group_by(District, haplotype) %>%
  summarise(Prevalence = sum(prevalence), .groups = "drop") %>%
  group_by(District) %>%
  mutate(Percentage = Prevalence / sum(Prevalence) * 100,
         ymax = cumsum(Percentage),
         ymin = c(0, head(ymax, n = -1)))

# --- 2. Define Colors for Pfmdr1 Haplotypes ---
#colors_pfmdr1 <- c("NFSND" ="#0086ad", "NYSND" = "#a6d0f0", "YFSND" = "#48CAE4", "YYSND" = "#6da3ff")
colors_pfmdr1 <- c(
  "NYSND" = "#9ed1ed",  
  "NFSND" = "#2471a3",   
  "YYSND" = "#cd5e08",
  "YFSND" = "#844da3"
)

# --- 3. Create Pie Chart Plot for Pfmdr1 ---
pie_chart_plot_pfmdr1 <- ggplot(pfmdr1_data_long, aes(x = "", y = Percentage, fill = haplotype)) +
  geom_bar(stat = "identity", width = 1, color = "white", size = 0.3) +
  coord_polar("y", start = 0) +
  facet_wrap(~ District, ncol = 5) +
  theme_void() +
  theme(
    legend.position = "right",
    strip.text = element_text(size = 8, face = "bold"),
    plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
    panel.spacing = unit(0.8, "lines"),
    plot.background = element_rect(fill = "transparent", color = NA)
  ) +
  labs(title = "") +
  scale_fill_manual(values = colors_pfmdr1)

# --- 4. Create Waffle Chart Data for Pfmdr1 ---
overall_data_pfmdr1 <- c(
  NFSND = 449,
  NYSND = 26,
  YFSND = 2,
  YYSND = 6
)

# --- 5. Create Waffle Chart (Pfmdr1) ---
waffle_chart_pfmdr1 <- waffle(
  parts = overall_data_pfmdr1,
  rows = 10,  # Keep the number of rows constant
  colors = colors_pfmdr1,
  title = "Overall Prevalence of Pfmdr1 Haplotypes (N=483)",
  xlab = "1 square = 1 sample",
  legend_pos = "bottom",
  size = 0.15  # Further reduced square size to make the chart even longer
) +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    plot.title = element_text(face = "bold", size = 12),          # Increased title size
    legend.title = element_text(face = "bold", size = 10),        # Bold and larger legend title
    legend.text = element_text(face = "bold", size = 8),         # Bold and larger legend labels
    axis.title.x = element_text(face = "bold", size = 8),        # Bold and larger x-axis label
    panel.spacing = unit(0, "lines")                              # Remove spacing between squares
  )

# Adjust the aspect ratio for tighter columns
waffle_chart_pfmdr1 <- waffle_chart_pfmdr1 + coord_fixed(ratio = 1.2)

# --- 6. Combine Pie Chart and Waffle Chart ---
combined_plot_pfmdr1 <- grid.arrange(waffle_chart_pfmdr1, pie_chart_plot_pfmdr1, ncol = 1, heights = c(1, 2))

# --- 7. Save Combined Plot ---
ggsave("../output/combined_pfmdr1_visualization.png", combined_plot_pfmdr1, bg = "transparent", width = 15, height = 12, dpi = 600)

# Display the combined plot
print(combined_plot_pfmdr1)








# pfk13 combined plot  ###############
# --- 1. Data Processing for Pfk13 ---
pfk13_data_long <- Haplotype %>%
  filter(Gene == "Pfk13") %>%
  select(District, haplotype, prevalence) %>%
  group_by(District, haplotype) %>%
  summarise(Prevalence = sum(prevalence), .groups = "drop") %>%
  group_by(District) %>%
  mutate(Percentage = Prevalence / sum(Prevalence) * 100,
         ymax = cumsum(Percentage),
         ymin = c(0, head(ymax, n = -1)))

pfk13_data_long$haplotype <- factor(
  pfk13_data_long$haplotype,
  levels = c("R622I",  "A675V","P441L", "A578S", "Wildtype")
)

# --- 2. Define Colors for Pfk13 Haplotypes ---
colors_pfk13 <- c("R622I" = "#993366", "P441L" = "#730099", "A675V" = "#6b66c6", "A578S" = "#ffa345", "Wildtype" = "#e4b4e4")

# --- 3. Create Pie Chart Plot for Pfk13 ---
ggplot(pfk13_data_long, aes(x = "", y = Percentage, fill = haplotype)) +
  geom_bar(stat = "identity", width = 1, color = "white", size = 0.3) +
  coord_polar("y", start = 0) +
  facet_wrap(~ District, ncol = 5) +
  theme_void() +
  theme(
    legend.position = "right",
    strip.text = element_text(size = 8, face = "bold"),
    plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
    panel.spacing = unit(0.8, "lines"),
    plot.background = element_rect(fill = "transparent", color = NA)
  ) +
  labs(title = "") +
  scale_fill_manual(values = colors_pfk13, name = expression(italic("Pf")* "K13 variant"))

ggsave("../output/pies_k13_visualization.png",bg = "transparent", width = 15, height = 12, dpi = 600)





 
# Base Map Pv #############

# --- Define Limits for Zoom ---
min_longitude <- 30.35396  # xmin
max_longitude <- 49.28536  # xmax
min_latitude <- 3          # ymin
max_latitude <- 15         # ymax

# --- Create Base Map with World Outline ---
ggplot() +
  # Plot world map with a light gray fill and gray outline
  geom_sf(data = Africa, fill = "gray90", color = "white", alpha = 0.3) +
  
  # Overlay Ethiopian regions in white
  geom_sf(data = Et_region, fill = "white", color = "gray80", alpha = 0.9) +
  
  # Plot Pv.proportion data with colors based on 'pvpropdf' field
  geom_sf(data = Haplotype_all, aes(fill = pvpropdf), color = NA, size = 0.2) +
  
  # Custom theme adjustments
  scale_fill_viridis_c(option = "D", name = "") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "gray75", color = NA),   # Set background color of the panel
    panel.grid.major = element_blank(),                             # Remove major grid lines
    panel.grid.minor = element_blank(),                             # Remove minor grid lines
    plot.background = element_rect(fill = "transparent", color = NA),  # Transparent plot background
    legend.position = "right",
    panel.border = element_rect(color = "gray", fill = NA, size = 1)   # Set border color without filling
  ) +
  labs(
    title = "Pv.Proportion Map"
  ) +
  
  # Zoom in on specific limits
  coord_sf(xlim = c(min_longitude, max_longitude), ylim = c(min_latitude, max_latitude), expand = FALSE)












# pfmdr1 Map ##############
# Define specific districts and colors for Pfmdr1 haplotypes
#colors_pfmdr1 <- c("NFSND" ="#0086ad", "NYSND" = "#a6d0f0", "YFSND" = "#48CAE4", "YYSND" = "#6da3ff")

colors_pfmdr1 <- c(
  "NYSND" = "#9ed1ed",  
  "NFSND" = "#2471a3",   
  "YYSND" = "#cd5e08",
  "YFSND" = "#844da3"
)



# Define zoom limits for the map
min_longitude <- 32.35396
max_longitude <- 48.28536
min_latitude <- 3
max_latitude <- 15

# Classify `pvpropdf` into categories
Haplotype_all <- Haplotype_all %>%
  mutate(pv_class = cut(
    pvpropdf, 
    breaks = c(-Inf, 0, 10, 20, 30, Inf),
    labels = c("0", "<10", "10-20", "20-30", "30-48"),
    right = FALSE
  ))

# Pre-compute centroids for each district (no need to calculate this multiple times)
Haplotype_all$centroid <- st_centroid(st_geometry(Haplotype_all))

# Define offsets for each district
offsets <- list(
  "Abergele" = c(0.5, 0.5),
  "Amibara" = c(0.1, 0.8),
  "Batu" = c(0, 0.8),
  "Dasenech" = c(0.8, 0.1),
  "Dera" = c(0.8, -0.4),
  "Fedis" = c(0.2, 0.7),
  "Fentale" = c(0.1, -0.8),
  "Gambella zuria" = c(-0.7, 0.2),
  "Gondar" = c(0.8, 0.1),
  "Guba" = c(-0.7, -0.2),
  "K. Humera" = c(-0.6, 0.1),
  "Metema" = c(-0.8, 0.1),
  "Mizan" = c(0.6, 0.4),
  "Pawe" = c(0.1, -0.7),
  "Selamago" = c(-0.7, 0.1)
)

# Create the base map
# Calculate centroids for each district
library(sf)
district_centroids1 <- st_centroid(Haplotype_all)

# Define offsets for each district
offsets1 <- list(
  "Abergele" = c(0.7, -0.1),
  "Amibara" = c(0.46, 0.22),
  "Batu" = c(-0.3, -0.3),
  "Dasenech" = c(2.1, 0.1),
  "Dera" = c(0.5, 0.3),
  "Fedis" = c(0.5, 0.1),
  "Fentale" = c(0.8, -0.1),
  "Gambella zuria" = c(1.5, 0.2),
  "Gondar" = c(-0.7, -0.2),
  "Guba" = c(0.1, 0.5),
  "K. Humera" = c(0.8, 0.2),
  "Metema" = c(0.9, 0.3),
  "Mizan" = c(-0.4, 0.3),
  "Pawe" = c(0.1, 0.3),
  "Selamago" = c(0.8, 0.1)
)

# Create a data frame for label positions and offsets
label_positions <- data.frame(
  District = district_centroids1$District,  # Assuming 'District' is the column name for districts
  x = st_coordinates(district_centroids1)[, 1],  # X-coordinates of centroids
  y = st_coordinates(district_centroids1)[, 2]   # Y-coordinates of centroids
)

# Add offset columns based on the offsets list
label_positions <- label_positions %>%
  rowwise() %>%
  mutate(
    x_offset = offsets1[[District]][1],
    y_offset = offsets1[[District]][2],
    x_adjusted = x + x_offset,
    y_adjusted = y + y_offset
  )

# Plot with adjusted label positions
# Add the border layer and label layer to the base map
# Ensure all data layers have the same CRS
common_crs <- st_crs(border) # Using border's CRS as the base
Africa <- st_transform(Africa, common_crs)
Et_region <- st_transform(Et_region, common_crs)
Haplotype_all <- st_transform(Haplotype_all, common_crs)

# Adjust the map plot
# Updated map with pvpropdf and external legend
# Updated map with vertical legend title
# Updated map with consistent font sizes
base_map <- ggplot() +
  geom_sf(data = Africa, fill = "gray80", color = "white", alpha = 0.3) +
  geom_sf(data = Et_region, fill = "white", color = "gray80", alpha = 0.9) +
  geom_sf(data = border, fill = NA, color = "NA", size = 0.5) +
  
  # Fill by pvpropdf
  geom_sf(data = Haplotype_all, aes(fill = pvpropdf), color = NA, size = 0.2) +
  
  # Continuous color scale with reversed colors and vertical title
  scale_fill_viridis_c(
    option = "D",
    direction = -1,  # Reverse the color scale
    name = expression(italic("P. vivax")~"proportion (%)"),
    na.value = "gray80",
    guide = guide_colorbar(
      barwidth = 1,    # Width of the bar
      barheight = 15,  # Length of the bar
      title.position = "right",  # Position the title next to the bar
      title.theme = element_text(
        angle = 90,     # Vertical title
        size = 14, 
        face = "bold"
      )
    )
  ) +
  
  # Title and theme adjustments
  ggtitle('') +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    #legend.position = "right",
    #legend.key.size = unit(0.5, "lines"),
    #legend.title = element_text(size = 14, face = "bold",family = "Calibri"),
    #legend.text = element_text(size = 14),
    legend.position = "none",
    panel.background = element_rect(fill = "white", color = NA),  # Make water white
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.background = element_blank()
  ) +
  
  # Set coordinate limits
  coord_sf(
    xlim = c(min_longitude, max_longitude), 
    ylim = c(min_latitude, max_latitude),   
    expand = FALSE
  ) +
  
  # Add district labels
  geom_text(
    data = label_positions,
    aes(x = x_adjusted, y = y_adjusted, label = District),
    size = 4.5, color = "gray20", family = "Calibri"
  ) +
  
  # Add border labels
  geom_sf_text(
    data = border,
    aes(label = Name_1),
    size = 3.5,
    color = "gray20",
    family = "Calibri"
  )

# Display the map
print(base_map)

# Loop through each district to create pie charts for Pfcrt prevalence
for (i in 1:nrow(Haplotype_all)) {
  district <- Haplotype_all$District[i]
  
  # Get the centroid coordinates for the current district
  centroid <- st_coordinates(Haplotype_all$centroid[i, ])
  
  # Get the offset for the district from the offsets list
  offset <- offsets[[district]]
  x_offset <- ifelse(!is.null(offset), offset[1], 0)
  y_offset <- ifelse(!is.null(offset), offset[2], 0)
  
  # Apply offset to the centroid coordinates
  adjusted_x <- centroid[1] + x_offset
  adjusted_y <- centroid[2] + y_offset
  
  # Correct the filtering to use `district` variable
  pie_df <- Haplotype_all %>%
    filter(District == district & haplotype %in% c("NFSND", "NYSND", "YFSND", "YYSND")) %>%
    select(haplotype, prevalence)
  
  # Check if pie_df is not empty and has valid prevalence values
  if (nrow(pie_df) > 0) {
    pie_chart <- ggplot(pie_df, aes(x = "", y = prevalence, fill = haplotype)) +
      geom_bar(stat = "identity", width = 1) +
      coord_polar(theta = "y") +
      scale_fill_manual(values = colors_pfmdr1) +
      theme_void() +
      theme(legend.position = "none")
    
    # Convert pie chart to graphical object (grob)
    pie_grob <- ggplotGrob(pie_chart)
    
    # Use fixed pie size for all pie charts (adjust as needed)
    pie_size <- 0.6  # Adjust for uniform pie size across districts
    
    # Add pie chart with the recalculated centroids and fixed pie size
    base_map <- base_map +
      annotation_custom(
        grob = pie_grob,
        xmin = adjusted_x - pie_size, xmax = adjusted_x + pie_size,
        ymin = adjusted_y - pie_size, ymax = adjusted_y + pie_size
      )
  }
}

# Save the final map as a TIFF file
output_path <- "../output/Pfmdr1_Map.tiff"
ggsave(
  filename = output_path,
  plot = base_map,
  device = "tiff",
  dpi = 300,  # High-resolution output
  width = 900 / 96,  # Convert pixels to inches (divide by 96 DPI)
  height = 577 / 96  # Convert pixels to inches (divide by 96 DPI)
)










# pfdhfr Map ###########################################################################################################################
# Define specific districts and colors for Pfdhfr haplotypes
colors_pfdhfr <- c(
  "ANCSI" = "#9ed1ed",  
  "ANRNI" = "#2471a3",   
  "AICNI" = "#cd5e08",
  "AIRNI" = "#844da3"
)


# Define zoom limits for the map
min_longitude <- 32.35396
max_longitude <- 48.28536
min_latitude <- 3
max_latitude <- 15

# Classify `pvpropdf` into categories
Haplotype_all <- Haplotype_all %>%
  mutate(pv_class = cut(
    pvpropdf, 
    breaks = c(-Inf, 0, 10, 20, 30, Inf),
    labels = c("0", "<10", "10-20", "20-30", "30-48"),
    right = FALSE
  ))

# Precompute centroids for each district (no need to calculate this multiple times)
Haplotype_all$centroid <- st_centroid(st_geometry(Haplotype_all))

# Define offsets for each district
offsets <- list(
  "Abergele" = c(0.5, 0.5),
  "Amibara" = c(0.1, 0.8),
  "Batu" = c(0, 0.8),
  "Dasenech" = c(0.8, 0.1),
  "Dera" = c(0.8, -0.4),
  "Fedis" = c(0.2, 0.7),
  "Fentale" = c(0.1, -0.8),
  "Gambella zuria" = c(-0.7, 0.2),
  "Gondar" = c(0.8, 0.1),
  "Guba" = c(-0.7, -0.2),
  "K. Humera" = c(-0.6, 0.1),
  "Metema" = c(-0.8, 0.1),
  "Mizan" = c(0.6, 0.4),
  "Pawe" = c(0.1, -0.7),
  "Selamago" = c(-0.7, 0.1)
)

# Create the base map
# Calculate centroids for each district
library(sf)
district_centroids1 <- st_centroid(Haplotype_all)

# Define offsets for each district
offsets1 <- list(
  "Abergele" = c(0.7, -0.1),
  "Amibara" = c(0.46, 0.22),
  "Batu" = c(-0.3, -0.3),
  "Dasenech" = c(2.1, 0.1),
  "Dera" = c(0.5, 0.3),
  "Fedis" = c(0.5, 0.1),
  "Fentale" = c(0.8, -0.1),
  "Gambella zuria" = c(1.5, 0.2),
  "Gondar" = c(-0.7, -0.2),
  "Guba" = c(0.1, 0.5),
  "K. Humera" = c(0.8, 0.2),
  "Metema" = c(0.9, 0.3),
  "Mizan" = c(-0.4, 0.3),
  "Pawe" = c(0.1, 0.3),
  "Selamago" = c(0.8, 0.1)
)

# Create a data frame for label positions and offsets
label_positions <- data.frame(
  District = district_centroids1$District,  # Assuming 'District' is the column name for districts
  x = st_coordinates(district_centroids1)[, 1],  # X-coordinates of centroids
  y = st_coordinates(district_centroids1)[, 2]   # Y-coordinates of centroids
)

# Add offset columns based on the offsets list
label_positions <- label_positions %>%
  rowwise() %>%
  mutate(
    x_offset = offsets1[[District]][1],
    y_offset = offsets1[[District]][2],
    x_adjusted = x + x_offset,
    y_adjusted = y + y_offset
  )

# Plot with adjusted label positions
# Add the border layer and label layer to the base map
# Ensure all data layers have the same CRS
common_crs <- st_crs(border) # Using border's CRS as the base
Africa <- st_transform(Africa, common_crs)
Et_region <- st_transform(Et_region, common_crs)
Haplotype_all <- st_transform(Haplotype_all, common_crs)

base_map <- ggplot() +
  geom_sf(data = Africa, fill = "gray80", color = "white", alpha = 0.3) +
  geom_sf(data = Et_region, fill = "white", color = "gray80", alpha = 0.9) +
  geom_sf(data = border, fill = NA, color = "NA", alpha = 0.9) +
  geom_sf(data = Haplotype_all, aes(fill = pv_class), color = NA, size = 0.2) +
  scale_fill_manual(
    values = c("0" = "#ffcccc", "<10" = "#ffcccc", "10-20" = "#ffcccc", "20-30" = "#ffcccc", "30-48" = "#ffcccc"),
    name = expression(italic("P.vivax")~"proportion")
  ) +
  ggtitle('') +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "none",  # Remove the legend
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),  # Make water white
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  coord_sf(
    xlim = c(min_longitude, max_longitude), 
    ylim = c(min_latitude, max_latitude),   
    expand = FALSE
  ) +
  geom_text(
    data = label_positions,
    aes(x = x_adjusted, y = y_adjusted, label = District),
    size = 4.5, color = "gray20", family = "Calibri"
  ) +
  # Add border labels using the "Name_1" field
  geom_sf_text(
    data = border,
    aes(label = Name_1),
    size = 3.5,
    color = "gray20",
    family = "Calibri"
  )

# Display the map
print(base_map)


# Loop through each district to create pie charts for Pfcrt prevalence
for (i in 1:nrow(Haplotype_all)) {
  district <- Haplotype_all$District[i]
  
  # Get the centroid coordinates for the current district
  centroid <- st_coordinates(Haplotype_all$centroid[i, ])
  
  # Get the offset for the district from the offsets list
  offset <- offsets[[district]]
  x_offset <- ifelse(!is.null(offset), offset[1], 0)
  y_offset <- ifelse(!is.null(offset), offset[2], 0)
  
  # Apply offset to the centroid coordinates
  adjusted_x <- centroid[1] + x_offset
  adjusted_y <- centroid[2] + y_offset
  
  # Correct the filtering to use `district` variable
  pie_df <- Haplotype_all %>%
    filter(District == district & haplotype %in% c("AICNI", "AIRNI", "ANCSI", "ANRNI")) %>%
    select(haplotype, prevalence)
  
  # Check if pie_df is not empty and has valid prevalence values
  if (nrow(pie_df) > 0) {
    pie_chart <- ggplot(pie_df, aes(x = "", y = prevalence, fill = haplotype)) +
      geom_bar(stat = "identity", width = 1) +
      coord_polar(theta = "y") +
      scale_fill_manual(values = colors_pfdhfr) +
      theme_void() +
      theme(legend.position = "none")
    
    # Convert pie chart to graphical object (grob)
    pie_grob <- ggplotGrob(pie_chart)
    
    # Use fixed pie size for all pie charts (adjust as needed)
    pie_size <- 0.6  # Adjust for uniform pie size across districts
    
    # Add pie chart with the recalculated centroids and fixed pie size
    base_map <- base_map +
      annotation_custom(
        grob = pie_grob,
        xmin = adjusted_x - pie_size, xmax = adjusted_x + pie_size,
        ymin = adjusted_y - pie_size, ymax = adjusted_y + pie_size
      )
  }
}

# Save the final map as a TIFF file
output_path <- "../output/Pfdhfr_Map.tiff"
ggsave(
  filename = output_path,
  plot = base_map,
  device = "tiff",
  dpi = 300,  # High-resolution output
  width = 900 / 96,  # Convert pixels to inches (divide by 96 DPI)
  height = 577 / 96  # Convert pixels to inches (divide by 96 DPI)
)









# pfdhps Map #########
# Define specific districts and colors for Pfdhps haplotypes
colors_pfdhps <- c(
  "SAKAA" = "#9ed1ed",   # Light blue
  "SGEAA" = "#2471a3", # Orange
  "SGEGA" =  "#cd5e08"   # Dark blue
)

# Define zoom limits for the map
min_longitude <- 32.35396
max_longitude <- 48.28536
min_latitude <- 3
max_latitude <- 15

# Classify `pvpropdf` into categories
Haplotype_all <- Haplotype_all %>%
  mutate(pv_class = cut(
    pvpropdf, 
    breaks = c(-Inf, 0, 10, 20, 30, Inf),
    labels = c("0", "<10", "10-20", "20-30", "30-48"),
    right = FALSE
  ))

# Precompute centroids for each district (no need to calculate this multiple times)
Haplotype_all$centroid <- st_centroid(st_geometry(Haplotype_all))

# Define offsets for each district
offsets <- list(
  "Abergele" = c(0.5, 0.5),
  "Amibara" = c(0.1, 0.8),
  "Batu" = c(0, 0.8),
  "Dasenech" = c(0.8, 0.1),
  "Dera" = c(0.8, -0.4),
  "Fedis" = c(0.2, 0.7),
  "Fentale" = c(0.1, -0.8),
  "Gambella zuria" = c(-0.7, 0.2),
  "Gondar" = c(0.8, 0.1),
  "Guba" = c(-0.7, -0.2),
  "K. Humera" = c(-0.6, 0.1),
  "Metema" = c(-0.8, 0.1),
  "Mizan" = c(0.6, 0.4),
  "Pawe" = c(0.1, -0.7),
  "Selamago" = c(-0.7, 0.1)
)

# Create the base map
# Calculate centroids for each district
library(sf)
district_centroids1 <- st_centroid(Haplotype_all)

# Define offsets for each district
offsets1 <- list(
  "Abergele" = c(0.7, -0.1),
  "Amibara" = c(0.46, 0.22),
  "Batu" = c(-0.3, -0.3),
  "Dasenech" = c(2.1, 0.1),
  "Dera" = c(0.5, 0.3),
  "Fedis" = c(0.5, 0.1),
  "Fentale" = c(0.8, -0.1),
  "Gambella zuria" = c(1.5, 0.2),
  "Gondar" = c(-0.7, -0.2),
  "Guba" = c(0.1, 0.5),
  "K. Humera" = c(0.8, 0.2),
  "Metema" = c(0.9, 0.3),
  "Mizan" = c(-0.4, 0.3),
  "Pawe" = c(0.1, 0.3),
  "Selamago" = c(0.8, 0.1)
)

# Create a data frame for label positions and offsets
label_positions <- data.frame(
  District = district_centroids1$District,  # Assuming 'District' is the column name for districts
  x = st_coordinates(district_centroids1)[, 1],  # X-coordinates of centroids
  y = st_coordinates(district_centroids1)[, 2]   # Y-coordinates of centroids
)

# Add offset columns based on the offsets list
label_positions <- label_positions %>%
  rowwise() %>%
  mutate(
    x_offset = offsets1[[District]][1],
    y_offset = offsets1[[District]][2],
    x_adjusted = x + x_offset,
    y_adjusted = y + y_offset
  )

# Plot with adjusted label positions
# Add the border layer and label layer to the base map
# Ensure all data layers have the same CRS
common_crs <- st_crs(border) # Using border's CRS as the base
Africa <- st_transform(Africa, common_crs)
Et_region <- st_transform(Et_region, common_crs)
Haplotype_all <- st_transform(Haplotype_all, common_crs)

# Create the base map with labels for border using "Name_1"
base_map <- ggplot() +
  geom_sf(data = Africa, fill = "gray80", color = "white", alpha = 0.3) +
  geom_sf(data = Et_region, fill = "white", color = "gray80", alpha = 0.9) +
  geom_sf(data = border, fill = NA, color = "NA", alpha = 0.9) +
  geom_sf(data = Haplotype_all, aes(fill = pv_class), color = NA, size = 0.2) +
  scale_fill_manual(
    values = c("0" = "#ffcccc", "<10" = "#ffcccc", "10-20" = "#ffcccc", "20-30" = "#ffcccc", "30-48" = "#ffcccc"),
    name = expression(italic("P.vivax")~"proportion")
  ) +
  ggtitle('') +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "none",  # Remove the legend
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),  # Make water white
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text=element_blank(),
    axis.ticks = element_blank()
  ) +
  coord_sf(
    xlim = c(min_longitude, max_longitude), 
    ylim = c(min_latitude, max_latitude),   
    expand = FALSE
  ) +
  geom_text(
    data = label_positions,
    aes(x = x_adjusted, y = y_adjusted, label = District),
    size = 4.5, color = "gray20", family = "Calibri"
  ) +
  # Add border labels using the "Name_1" field
  geom_sf_text(
    data = border,
    aes(label = Name_1),
    size = 3.5,
    color = "gray20",
    family = "Calibri"
  )


# Loop through each district to create pie charts for Pfcrt prevalence
for (i in 1:nrow(Haplotype_all)) {
  district <- Haplotype_all$District[i]
  
  # Get the centroid coordinates for the current district
  centroid <- st_coordinates(Haplotype_all$centroid[i, ])
  
  # Get the offset for the district from the offsets list
  offset <- offsets[[district]]
  x_offset <- ifelse(!is.null(offset), offset[1], 0)
  y_offset <- ifelse(!is.null(offset), offset[2], 0)
  
  # Apply offset to the centroid coordinates
  adjusted_x <- centroid[1] + x_offset
  adjusted_y <- centroid[2] + y_offset
  
  # Correct the filtering to use `district` variable
  pie_df <- Haplotype_all %>%
    filter(District == district & haplotype %in% c("SAKAA", "SGEAA","SGEGA")) %>%
    select(haplotype, prevalence)
  
  # Check if pie_df is not empty and has valid prevalence values
  if (nrow(pie_df) > 0) {
    pie_chart <- ggplot(pie_df, aes(x = "", y = prevalence, fill = haplotype)) +
      geom_bar(stat = "identity", width = 1) +
      coord_polar(theta = "y") +
      scale_fill_manual(values = colors_pfdhps) +
      theme_void() +
      theme(legend.position = "none")
    
    # Convert pie chart to graphical object (grob)
    pie_grob <- ggplotGrob(pie_chart)
    
    # Use fixed pie size for all pie charts (adjust as needed)
    pie_size <- 0.6  # Adjust for uniform pie size across districts
    
    # Add pie chart with the recalculated centroids and fixed pie size
    base_map <- base_map +
      annotation_custom(
        grob = pie_grob,
        xmin = adjusted_x - pie_size, xmax = adjusted_x + pie_size,
        ymin = adjusted_y - pie_size, ymax = adjusted_y + pie_size
      )
  }
}

# Save the final map as a TIFF file
output_path <- "../output/Pfdhps_Map.tiff"
ggsave(
  filename = output_path,
  plot = base_map,
  device = "tiff",
  dpi = 300,  # High-resolution output
  width = 900 / 96,  # Convert pixels to inches (divide by 96 DPI)
  height = 577 / 96  # Convert pixels to inches (divide by 96 DPI)
)








# pfcrt Map ###########################################################################################################################
# Define specific districts and colors for Pfcrt haplotypes
#colors_pfcrt <- c("CVIET" = "#963696", "CVMNK" = "#e4b4e4")
colors_pfcrt <- c(
  "CVMNK" = "#9ed1ed",   
  "CVIET" = "#2471a3" 
)


# Define zoom limits for the map
min_longitude <- 32.35396
max_longitude <- 48.28536
min_latitude <- 3
max_latitude <- 15

# Classify `pvpropdf` into categories
Haplotype_all <- Haplotype_all %>%
  mutate(pv_class = cut(
    pvpropdf, 
    breaks = c(-Inf, 0, 10, 20, 30, Inf),
    labels = c("0", "<10", "10-20", "20-30", "30-48"),
    right = FALSE
  ))

# Precompute centroids for each district (no need to calculate this multiple times)
Haplotype_all$centroid <- st_centroid(st_geometry(Haplotype_all))

# Define offsets for each district
offsets <- list(
  "Abergele" = c(0.5, 0.5),
  "Amibara" = c(0.1, 0.8),
  "Batu" = c(0, 0.8),
  "Dasenech" = c(0.8, 0.1),
  "Dera" = c(0.8, -0.4),
  "Fedis" = c(0.2, 0.7),
  "Fentale" = c(0.1, -0.8),
  "Gambella zuria" = c(-0.7, 0.2),
  "Gondar" = c(0.8, 0.1),
  "Guba" = c(-0.7, -0.2),
  "K. Humera" = c(-0.6, 0.1),
  "Metema" = c(-0.8, 0.1),
  "Mizan" = c(0.6, 0.4),
  "Pawe" = c(0.1, -0.7),
  "Selamago" = c(-0.7, 0.1)
)

# Create the base map
# Calculate centroids for each district
library(sf)
district_centroids1 <- st_centroid(Haplotype_all)

# Define offsets for each district
offsets1 <- list(
  "Abergele" = c(0.7, -0.1),
  "Amibara" = c(0.46, 0.22),
  "Batu" = c(-0.3, -0.3),
  "Dasenech" = c(2.1, 0.1),
  "Dera" = c(0.5, 0.3),
  "Fedis" = c(0.5, 0.1),
  "Fentale" = c(0.8, -0.1),
  "Gambella zuria" = c(1.5, 0.2),
  "Gondar" = c(-0.7, -0.2),
  "Guba" = c(0.1, 0.5),
  "K. Humera" = c(0.8, 0.2),
  "Metema" = c(0.9, 0.3),
  "Mizan" = c(-0.4, 0.3),
  "Pawe" = c(0.1, 0.3),
  "Selamago" = c(0.8, 0.1)
)

# Create a data frame for label positions and offsets
label_positions <- data.frame(
  District = district_centroids1$District,  # Assuming 'District' is the column name for districts
  x = st_coordinates(district_centroids1)[, 1],  # X-coordinates of centroids
  y = st_coordinates(district_centroids1)[, 2]   # Y-coordinates of centroids
)

# Add offset columns based on the offsets list
label_positions <- label_positions %>%
  rowwise() %>%
  mutate(
    x_offset = offsets1[[District]][1],
    y_offset = offsets1[[District]][2],
    x_adjusted = x + x_offset,
    y_adjusted = y + y_offset
  )

# Plot with adjusted label positions
# Add the border layer and label layer to the base map
# Ensure all data layers have the same CRS
common_crs <- st_crs(border) # Using border's CRS as the base
Africa <- st_transform(Africa, common_crs)
Et_region <- st_transform(Et_region, common_crs)
Haplotype_all <- st_transform(Haplotype_all, common_crs)

# Updated map with consistent font sizes
base_map <- ggplot() +
  geom_sf(data = Africa, fill = "gray80", color = "white", alpha = 0.3) +
  geom_sf(data = Et_region, fill = "white", color = "gray80", alpha = 0.9) +
  geom_sf(data = border, fill = NA, color = "NA", size = 0.5) +
  
  # Fill by pvpropdf
  geom_sf(data = Haplotype_all, aes(fill = pvpropdf), color = NA, size = 0.2) +
  
  # Continuous color scale with reversed colors and vertical title
  scale_fill_viridis_c(
    option = "D",
    direction = -1,  # Reverse the color scale
    name = expression(italic("P. vivax")~"Proportion (%)"),
    na.value = "gray80",
    guide = guide_colorbar(
      barwidth = 1,    # Width of the bar
      barheight = 15,  # Length of the bar
      title.position = "right",  # Position the title next to the bar
      title.theme = element_text(
        angle = 90,     # Vertical title
        size = 14, 
        face = "bold"
      )
    )
  ) +
  
  # Title and theme adjustments
  ggtitle('') +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "none",  # Remove the legend
    legend.key.size = unit(0.5, "lines"),
    legend.title = element_text(size = 14, face = "bold",family = "Calibri"),
    legend.text = element_text(size = 14),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),  # Make water white
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.background = element_blank()
  ) +
  
  # Set coordinate limits
  coord_sf(
    xlim = c(min_longitude, max_longitude), 
    ylim = c(min_latitude, max_latitude),   
    expand = FALSE
  ) +
  
  # Add district labels
  geom_text(
    data = label_positions,
    aes(x = x_adjusted, y = y_adjusted, label = District),
    size = 4.5, color = "gray20", family = "Calibri"
  ) +
  
  # Add border labels
  geom_sf_text(
    data = border,
    aes(label = Name_1),
    size = 3.5,
    color = "gray20",
    family = "Calibri"
  )

# Display the map
print(base_map)
# Loop through each district to create pie charts for Pfcrt prevalence
for (i in 1:nrow(Haplotype_all)) {
  district <- Haplotype_all$District[i]
  
  # Get the centroid coordinates for the current district
  centroid <- st_coordinates(Haplotype_all$centroid[i, ])
  
  # Get the offset for the district from the offsets list
  offset <- offsets[[district]]
  x_offset <- ifelse(!is.null(offset), offset[1], 0)
  y_offset <- ifelse(!is.null(offset), offset[2], 0)
  
  # Apply offset to the centroid coordinates
  adjusted_x <- centroid[1] + x_offset
  adjusted_y <- centroid[2] + y_offset
  
  # Correct the filtering to use `district` variable
  pie_df <- Haplotype_all %>%
    filter(District == district & haplotype %in% c("CVIET", "CVMNK")) %>%
    select(haplotype, prevalence)
  
  # Check if pie_df is not empty and has valid prevalence values
  if (nrow(pie_df) > 0) {
    pie_chart <- ggplot(pie_df, aes(x = "", y = prevalence, fill = haplotype)) +
      geom_bar(stat = "identity", width = 1) +
      coord_polar(theta = "y") +
      scale_fill_manual(values = colors_pfcrt) +
      theme_void() +
      theme(legend.position = "none")
    
    # Convert pie chart to graphical object (grob)
    pie_grob <- ggplotGrob(pie_chart)
    
    # Use fixed pie size for all pie charts (adjust as needed)
    pie_size <- 0.6  # Adjust for uniform pie size across districts
    
    # Add pie chart with the recalculated centroids and fixed pie size
    base_map <- base_map +
      annotation_custom(
        grob = pie_grob,
        xmin = adjusted_x - pie_size, xmax = adjusted_x + pie_size,
        ymin = adjusted_y - pie_size, ymax = adjusted_y + pie_size
      )
  }
}

# Save the final map as a TIFF file
output_path <- "../output/Pfcrt_Map.tiff"
ggsave(
  filename = output_path,
  plot = base_map,
  device = "tiff",
  dpi = 300,  # High-resolution output
  width = 900 / 96,  # Convert pixels to inches (divide by 96 DPI)
  height = 577 / 96  # Convert pixels to inches (divide by 96 DPI)
)










# pfk13 Map ###########################################################################################################################

# Define specific districts and colors for PfK13 mutations
colors_pfk13 <- c("R622I" = "#993366", "P441L" = "#730099", "A675V" = "#6b66c6", "A578S" = "#ffa345", "Wildtype" = "#e4b4e4")

# Load district shapefile and filter for specific districts
District <- st_read("../Data_used/District.shp")

# Define zoom limits for the map
min_longitude <- 32.35396
max_longitude <- 48.28536
min_latitude <- 3
max_latitude <- 15

# Classify `pvpropdf` into categories
Haplotype_all <- Haplotype_all %>%
  mutate(pv_class = cut(
    pvpropdf, 
    breaks = c(-Inf, 0, 10, 20, 30, Inf),
    labels = c("0", "<10", "10-20", "20-30", "30-48"),
    right = FALSE
  ))

# Precompute centroids for each district (no need to calculate this multiple times)
Haplotype_all$centroid <- st_centroid(st_geometry(Haplotype_all))

# Define offsets for each district
offsets <- list(
  "Abergele" = c(0.5, 0.5),
  "Amibara" = c(0.1, 0.8),
  "Batu" = c(0, 0.8),
  "Dasenech" = c(0.8, 0.1),
  "Dera" = c(0.8, -0.4),
  "Fedis" = c(0.2, 0.7),
  "Fentale" = c(0.1, -0.8),
  "Gambella zuria" = c(-0.7, 0.2),
  "Gondar" = c(0.8, 0.1),
  "Guba" = c(-0.7, -0.2),
  "K. Humera" = c(-0.6, 0.1),
  "Metema" = c(-0.8, 0.1),
  "Mizan" = c(0.6, 0.4),
  "Pawe" = c(0.1, -0.7),
  "Selamago" = c(-0.7, 0.1)
)

# Create the base map
# Calculate centroids for each district
library(sf)
district_centroids1 <- st_centroid(Haplotype_all)

# Define offsets for each district
offsets1 <- list(
  "Abergele" = c(0.7, -0.1),
  "Amibara" = c(0.46, 0.22),
  "Batu" = c(-0.3, -0.3),
  "Dasenech" = c(2.1, 0.1),
  "Dera" = c(0.5, 0.3),
  "Fedis" = c(0.5, 0.1),
  "Fentale" = c(0.8, -0.1),
  "Gambella zuria" = c(1.5, 0.2),
  "Gondar" = c(-0.7, -0.2),
  "Guba" = c(0.1, 0.5),
  "K. Humera" = c(0.8, 0.2),
  "Metema" = c(0.9, 0.3),
  "Mizan" = c(-0.4, 0.3),
  "Pawe" = c(0.1, 0.3),
  "Selamago" = c(0.8, 0.1)
)

# Create a data frame for label positions and offsets
label_positions <- data.frame(
  District = district_centroids1$District,  # Assuming 'District' is the column name for districts
  x = st_coordinates(district_centroids1)[, 1],  # X-coordinates of centroids
  y = st_coordinates(district_centroids1)[, 2]   # Y-coordinates of centroids
)

# Add offset columns based on the offsets list
label_positions <- label_positions %>%
  rowwise() %>%
  mutate(
    x_offset = offsets1[[District]][1],
    y_offset = offsets1[[District]][2],
    x_adjusted = x + x_offset,
    y_adjusted = y + y_offset
  )

# Add the border layer and label layer to the base map
# Ensure all data layers have the same CRS
common_crs <- st_crs(border) # Using border's CRS as the base
Africa <- st_transform(Africa, common_crs)
Et_region <- st_transform(Et_region, common_crs)
Haplotype_all <- st_transform(Haplotype_all, common_crs)
District <- st_transform(District, common_crs)

# Define specific colors for Pfk13 mutations
colors_pfk13 <- c("R622I" = "#993366", "P441L" = "#730099", 
                  "A675V" = "#6b66c6", "A578S" = "#ffa345", 
                  "Wildtype" = "#e4b4e4")

# Create the base map with labels for border using "Name_1"
base_map <- ggplot() +
  geom_sf(data = Africa, fill = "gray80", color = "white", alpha = 0.3) +
  geom_sf(data = Et_region, fill = "white", color = "gray80", alpha = 0.9) +
  geom_sf(data = border, fill = NA, color = "NA", alpha = 0.9) +
  geom_sf(data = Haplotype_all, aes(fill = pv_class), color = NA, size = 0.2) +
  scale_fill_manual(
    values = c("0" = "#ffcccc", "<10" = "#ffcccc", "10-20" = "#ffcccc", "20-30" = "#ffcccc", "30-48" = "#ffcccc"),
    name = expression(italic("P.vivax")~"proportion")
  ) +
  ggtitle('') +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "none",  # Remove the legend
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  coord_sf(
    xlim = c(min_longitude, max_longitude), 
    ylim = c(min_latitude, max_latitude),   
    expand = FALSE
  ) +
  geom_text(
    data = label_positions,
    aes(x = x_adjusted, y = y_adjusted, label = District),
    size = 4.5, color = "gray20", family = "Calibri"
  ) +
  # Add border labels using the "Name_1" field
  geom_sf_text(
    data = border,
    aes(label = Name_1),
    size = 3.5,
    color = "gray20",
    family = "Calibri"
  )
# Loop through each district to add pie charts for Pfk13 prevalence (no changes needed here)
for (i in 1:nrow(Haplotype_all)) {
  district <- Haplotype_all$District[i]
  
  # Get the centroid coordinates for the current district
  centroid <- st_coordinates(Haplotype_all$centroid[i, ])
  
  # Get the offset for the district from the offsets list
  offset <- offsets[[district]]
  x_offset <- ifelse(!is.null(offset), offset[1], 0)
  y_offset <- ifelse(!is.null(offset), offset[2], 0)
  
  # Apply offset to the centroid coordinates
  adjusted_x <- centroid[1] + x_offset
  adjusted_y <- centroid[2] + y_offset
  
  # Correct the filtering to use `district` variable
  pie_df <- Haplotype_all %>%
    filter(District == district & haplotype %in% c("R622I", "P441L", "A675V", "A578S", "Wildtype")) %>%
    select(haplotype, prevalence)
  
  # Check if pie_df is not empty and has valid prevalence values
  if (nrow(pie_df) > 0) {
    pie_chart <- ggplot(pie_df, aes(x = "", y = prevalence, fill = haplotype)) +
      geom_bar(stat = "identity", width = 1) +
      coord_polar(theta = "y") +
      scale_fill_manual(values = colors_pfk13) +
      theme_void() +
      theme(legend.position = "none")
    
    # Convert pie chart to graphical object (grob)
    pie_grob <- ggplotGrob(pie_chart)
    
    # Use fixed pie size for all pie charts (adjust as needed)
    pie_size <- 0.6  # Adjust for uniform pie size across districts
    
    # Add pie chart with the recalculated centroids and fixed pie size
    base_map <- base_map +
      annotation_custom(
        grob = pie_grob,
        xmin = adjusted_x - pie_size, xmax = adjusted_x + pie_size,
        ymin = adjusted_y - pie_size, ymax = adjusted_y + pie_size
      )
  }
}

# Save the final map as a TIFF file
output_path <- "../output/Pfk13_Map.tiff"
ggsave(
  filename = output_path,
  plot = base_map,
  device = "tiff",
  dpi = 300,  # High-resolution output
  width = 900 / 96,  # Convert pixels to inches (divide by 96 DPI)
  height = 577 / 96  # Convert pixels to inches (divide by 96 DPI)
)







# waffle plots  ######

# --- 1. Load and Prepare Data ---
Haplotype <- read_excel("../Data_used/All_haplotype.xlsx")

# --- 3. Create Waffle Chart with White Background and Tight Legend ---
waffle_chart_pfcrt <- waffle(
  parts = overall_data_pfcrt,
  rows = 8,  # Rows for grid layout
  colors = colors_pfcrt,
  legend_pos = "bottom",
  size = 0.15  # Adjust square size
) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),  # White background
    panel.background = element_rect(fill = "white", color = NA), # Ensure panel background is white
    plot.title = element_text(face = "bold", size = 16),         # Bold title
    legend.title = element_text(face = "bold", size = 14),       # Bold legend title
    legend.text = element_text(face = "bold", size = 12),        # Bold legend text
    axis.title.x = element_text(face = "bold", size = 12),       # Bold x-axis label
    legend.spacing.y = unit(-0.5, "cm"),                         # Negative vertical spacing to remove gap
    legend.spacing.x = unit(0, "cm"),                            # No horizontal spacing
    legend.box.margin = margin(-10, 0, -10, 0, "pt"),            # Remove margin between legend and plot
    legend.key.height = unit(0.4, "cm"),                         # Adjust legend key height
    legend.key.width = unit(0.4, "cm")                           # Adjust legend key width
  )

# --- 2. Save the Waffle Chart Separately ---
# Check if the directory exists, and create it if necessary
save_dir <- "../output/Waffle"
if (!dir.exists(save_dir)) {
  dir.create(save_dir, recursive = TRUE)  # Create the directory and all necessary parent directories
}

# Save the Waffle Chart
ggsave(
  filename = file.path(save_dir, "Waffle_Pfcrt_update.tiff"),  # Construct the full file path
  plot = waffle_chart_pfcrt,                          # Specify the waffle chart to save
  width = 5,                                          # Width of the saved image in inches
  height = 4,                                         # Height of the saved image in inches
  dpi = 300,                                          # Resolution of the saved image
  bg = "white"                                        # Set background to white
)

# Display Waffle Chart (optional)
print(waffle_chart_pfcrt)








# --- 2. Create Waffle Chart (Pfdhfr) ---
waffle_chart_pfdhfr <- waffle(
  parts = overall_data_pfdhfr,
  rows = 8,  # Number of rows for the waffle chart
  colors = colors_pfdhfr,
  legend_pos = "bottom",
  size = 0.3  # Adjusted square size
) +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),  # Transparent background
    panel.background = element_rect(fill = "transparent", color = NA), # No panel background
    plot.title = element_text(face = "bold", size = 16),               # Bold and larger title
    legend.title = element_text(face = "bold", size = 14),             # Bold legend title
    legend.text = element_text(face = "bold", size = 12),              # Bold legend labels
    axis.title.x = element_text(face = "bold", size = 12),             # Bold x-axis label
    legend.spacing.y = unit(-0.5, "cm"),                               # Reduce gap between legend and waffle
    legend.box.margin = margin(-10, 0, -10, 0, "pt"),                  # Remove margin between legend and plot
    legend.key.height = unit(0.4, "cm"),                               # Adjust legend key height
    legend.key.width = unit(0.4, "cm"),                                # Adjust legend key width
    panel.spacing = unit(0, "lines")                                   # Remove panel spacing
  )

# --- 6. Save Waffle Chart Separately ---
# Define save directory
save_dir <- "../output/Waffle"
if (!dir.exists(save_dir)) {
  dir.create(save_dir, recursive = TRUE)  # Create the directory and all necessary parent directories
}

# Save Waffle Chart
ggsave(
  filename = file.path(save_dir, "Waffle_Pfdhfr.png"),  # Save as Waffle_Pfdhfr.png
  plot = waffle_chart_pfdhfr,                          # Specify the waffle chart to save
  width = 6,                                           # Adjust width for better visualization
  height = 4,                                          # Adjust height
  dpi = 300,                                           # High resolution
  bg = "transparent"                                   # Transparent background
)

# Display Waffle Chart (optional)
print(waffle_chart_pfdhfr)








# --- 3. Create Waffle Chart Data for Pfdhps ---
overall_data_pfdhps <- c(
  SAKAA = 155,
  SGEAA = 330,
  SGEGA = 17
)

# --- 5. Create Waffle Chart (Pfdhps) ---
waffle_chart_pfdhps <- waffle(
  parts = overall_data_pfdhps,
  rows = 8,  # Keeping the same number of rows
  colors = colors_pfdhps,
  legend_pos = "bottom",  # Position legend below the chart
  size = 0.15             # Adjust square size for visual clarity
) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),       # White plot background
    panel.background = element_rect(fill = "white", color = NA),      # White panel background
    legend.title = element_text(face = "bold", size = 10),            # Bold and larger legend title
    legend.text = element_text(face = "bold", size = 12),             # Bold legend labels
    legend.spacing.y = unit(-0.5, "cm"),                              # Reduce gap between legend and waffle chart
    legend.box.margin = margin(-10, 0, -10, 0, "pt"),                 # Tighten spacing around the legend box
    legend.key.height = unit(0.4, "cm"),                              # Adjust legend key height
    legend.key.width = unit(0.4, "cm"),                               # Adjust legend key width
    panel.spacing = unit(0, "lines")                                  # Remove spacing between chart elements
  )

# --- 6. Save Waffle Chart Separately ---
# Define save directory
save_dir <- "../output/Waffle"
if (!dir.exists(save_dir)) {
  dir.create(save_dir, recursive = TRUE)  # Create the directory and all necessary parent directories
}

# Save Waffle Chart
ggsave(
  filename = file.path(save_dir, "Waffle_Pfdhps.png"),  # Save file path
  plot = waffle_chart_pfdhps,                           # Specify the waffle chart to save
  width = 6,                                            # Width of the saved image in inches
  height = 4,                                           # Height of the saved image in inches
  dpi = 300,                                            # High resolution for clarity
  bg = "white"                                          # White background
)

# Display Waffle Chart (optional)
print(waffle_chart_pfdhps)











# --- 4. Create Waffle Chart Data for Pfmdr1 ---
overall_data_pfmdr1 <- c(
  NFSND = 449,
  NYSND = 26,
  YFSND = 2,
  YYSND = 6
)

# --- 5. Create Waffle Chart (Pfmdr1) ---
waffle_chart_pfmdr1 <- waffle(
  parts = overall_data_pfmdr1,
  rows = 8,  # Number of rows for grid layout
  colors = colors_pfmdr1,
  legend_pos = "bottom",  # Position legend below the chart
  size = 0.15             # Adjust square size for clarity
) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),  # White plot background
    panel.background = element_rect(fill = "white", color = NA), # White panel background
    legend.title = element_text(face = "bold", size = 14),       # Bold and larger legend title
    legend.text = element_text(face = "bold", size = 12),        # Bold legend labels
    legend.spacing.y = unit(-0.5, "cm"),                         # Reduce gap between legend and waffle chart
    legend.box.margin = margin(-10, 0, -10, 0, "pt"),            # Tighten spacing around the legend box
    legend.key.height = unit(0.4, "cm"),                         # Adjust legend key height
    legend.key.width = unit(0.4, "cm"),                          # Adjust legend key width
    panel.spacing = unit(0, "lines")                             # Remove spacing between chart elements
  )

# --- 6. Save Waffle Chart Separately ---
# Define save directory
# Check if the directory exists, and create it if necessary
save_dir <- "../output/Waffle"
if (!dir.exists(save_dir)) {
  dir.create(save_dir, recursive = TRUE)  # Create the directory and all necessary parent directories
}

# Save Waffle Chart
ggsave(
  filename = file.path(save_dir, "Waffle_Pfmdr1_upd.tiff"),  # Save file path
  plot = waffle_chart_pfmdr1,                          # Specify the waffle chart to save
  width = 5,                                           # Width of the saved image in inches
  height = 4,                                          # Height of the saved image in inches
  dpi = 300,                                           # High resolution for clarity
  bg = "white"                                         # White background
)

# Display Waffle Chart (optional)
print(waffle_chart_pfmdr1)









# Bar graphs #######
# --- 1. Load Required Libraries ---
library(ggplot2)
library(dplyr)
# --- 2. Load Data ---
Bar_dat <- read.csv("../Data_used/All bar data.csv")

# --- 3. Filter Data for Pfcrt ---
pfcrt_data <- Bar_dat %>%
  filter(Gene == "Pfcrt") %>%
  mutate(Mutation.Position = recode(Mutation.Position,
                                    `72` = "C72S",
                                    `73` = "V73K",
                                    `74` = "M74I",
                                    `75` = "N75E",
                                    `76` = "K76T"),
         Mutation.Position = factor(Mutation.Position, levels = c("C72S","V73K","M74I","N75E","K76T")))

# --- 4. Define Colors ---
colors_pfcrt <- c("Wildtype" = "white", "Mutant" = "grey")

# --- 5. Create Stacked Bar Plot ---
stacked_bar_plot_pfcrt<-ggplot(pfcrt_data, aes(y = Overall.Prevalence...., 
                                                 x = as.factor(Mutation.Position),
                                                 fill = Variant_Type)) +
  geom_bar(stat = "identity", width = 0.9,color="grey") +
  scale_fill_manual(values = colors_pfcrt) +
  #scale_pattern_manual()+
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))+
  labs(x = "variant", y = "frequency (%)") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.ticks = element_line(color = "black"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    plot.background = element_rect(fill = "transparent", color = NA)
  )

# --- 6. Define Save Directory ---
save_dir <- "../output/Bar"

# Ensure Save Directory Exists
if (!dir.exists(save_dir)) {
  dir.create(save_dir, recursive = TRUE)
}

# --- 7. Save the Stacked Bar Plot ---
ggsave(
  filename = file.path(save_dir, "Stacked_Bar_Pfcrt_NoGrid.png"),  # Save file path
  plot = stacked_bar_plot_pfcrt,                                         # Specify the bar plot to save
  width = 4,                                                       # Width of the saved image in inches
  height = 4,                                                      # Height of the saved image in inches
  dpi = 300,                                                       # High resolution for clarity
  bg = "transparent"                                               # Transparent background
)

# --- 8. Display the Stacked Bar Plot ---
print(stacked_bar_plot_pfcrt)









 # pfdhfr Bar graph
# --- 1. Load Data ---
Bar_dat <- read.csv("../Data_used/All bar data.csv")

# --- 2. Filter Data for Pfdhfr ---
pfdhfr_data <- Bar_dat %>%
  filter(Gene == "Pfdhfr") %>%
  mutate(Mutation.Position = recode(Mutation.Position,
                                    `16` = "A16V",
                                    `51` = "N51I",
                                    `59` = "C59R",
                                    `108` = "S108N",
                                    `164` = "I164L"),
         Mutation.Position = factor(Mutation.Position, levels = c("A16V","N51I","C59R","S108N","I164L")))


# --- 3. Define Colors for Wildtype and Mutant ---
colors_pfdhfr <- c("Wildtype" = "white", "Mutant" = "grey")

# --- 4. Create Stacked Bar Plot for Pfdhfr ---
stacked_bar_plot_pfdhfr<-ggplot(pfdhfr_data, aes(y = Overall.Prevalence...., 
                                                 x = as.factor(Mutation.Position),
                                                 fill = Variant_Type)) +
  geom_bar(stat = "identity", width = 0.9,color="grey") +
  scale_fill_manual(values = colors_pfdhfr) +
  #scale_pattern_manual()+
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))+
  labs(x = "variant", y = "frequency (%)") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.ticks = element_line(color = "black"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    plot.background = element_rect(fill = "transparent", color = NA)
  )

# --- 5. Save Stacked Bar Plot for Pfdhfr ---
# Define save directory
save_dir <- "../output/Bar"
if (!dir.exists(save_dir)) {
  dir.create(save_dir, recursive = TRUE)  # Create the directory if it doesn't exist
}

# Save the bar plot for Pfdhfr
ggsave(
  filename = file.path(save_dir, "Stacked_Bar_Pfdhfr.png"),  # Save file path
  plot = stacked_bar_plot_pfdhfr,                            # Specify the bar plot to save
  width = 4,                                                # Width of the saved image in inches
  height = 4,                                               # Height of the saved image in inches
  dpi = 300,                                                # High resolution for clarity
  bg = "transparent"                                        # Transparent background
)

# Display the Stacked Bar Plot for Pfdhfr (optional)
print(stacked_bar_plot_pfdhfr)




# Pfdhps Bar graph
# --- 1. Load Data ---
Bar_dat <- read.csv("../Data_used/All bar data.csv")

# --- 2. Filter Data for Pfdhps ---
pfdhps_data <- Bar_dat %>%
  filter(Gene == "Pfdhps") %>%
  mutate(Mutation.Position = recode(Mutation.Position,
                                    `431` = "I431V",
                                    `436` = "S436A",
                                    `437` = "A437G",
                                    `540` = "K540E",
                                    `581` = "A581G",
                                    `613` = "A613S"),
         Mutation.Position = factor(Mutation.Position, levels = c("I431V","S436A","A437G","K540E","A581G","A613S")))

# --- 3. Define Colors for Wildtype and Mutant ---
colors_pfdhps <- c("Wildtype" = "white", "Mutant" = "grey")  

# --- 4. Create Stacked Bar Plot for Pfdhps ---
stacked_bar_plot_pfdhps<-ggplot(pfdhps_data, aes(y = Overall.Prevalence...., 
                                                 x = as.factor(Mutation.Position),
                                                 fill = Variant_Type)) +
  geom_bar(stat = "identity", width = 0.9,color="grey") +
  scale_fill_manual(values = colors_pfmdr1) +
  #scale_pattern_manual()+
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))+
  labs(x = "variant", y = "frequency (%)") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 11),
    axis.ticks = element_line(color = "black"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    plot.background = element_rect(fill = "transparent", color = NA)
  )

# --- 5. Save Stacked Bar Plot for Pfdhps ---
# Define save directory
save_dir <- "../output/Bar"
if (!dir.exists(save_dir)) {
  dir.create(save_dir, recursive = TRUE)  # Create the directory if it doesn't exist
}

# Save the bar plot for Pfdhps
ggsave(
  filename = file.path(save_dir, "Stacked_Bar_Pfdhps.png"),  # Save file path
  plot = stacked_bar_plot_pfdhps,                            # Specify the bar plot to save
  width = 4,                                                # Width of the saved image in inches
  height = 4,                                               # Height of the saved image in inches
  dpi = 300,                                                # High resolution for clarity
  bg = "white"                                              # White background
)

# Display the Stacked Bar Plot for Pfdhps (optional)
print(stacked_bar_plot_pfdhps)




   # Bar fo Pfdmr1 bar graph 
# --- 1. Load Data ---
Bar_dat <- read.csv("../Data_used/All bar data.csv")

# --- 2. Filter Data for Pfmdr1 ---
pfmdr1_data <- Bar_dat %>%
  filter(Gene == "Pfmdr1") %>%
  mutate(Mutation.Position = recode(Mutation.Position,
                                    `86` = "N86Y",
                                    `184` = "Y184F",
                                    `1034` = "S1034C",
                                    `1042` = "N1042D",
                                    `1246` = "D1246Y"),
         Mutation.Position = factor(Mutation.Position, levels = c("N86Y","Y184F","S1034C","N1042D","D1246Y")))

# --- 3. Define Colors for Wildtype and Mutant ---
colors_pfmdr1 <- c("Wildtype" = "white", "Mutant" = "darkgrey")

# --- 4. Create Stacked Bar Plot for Pfmdr1 ---
stacked_bar_plot_pfmdr1<-ggplot(pfmdr1_data, aes(y = Overall.Prevalence...., 
                        x = as.factor(Mutation.Position),
                        fill = Variant_Type)) +
  geom_bar(stat = "identity", width = 0.9,color="grey") +
  scale_fill_manual(values = colors_pfmdr1) +
  #scale_pattern_manual()+
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))+
  labs(x = "variant", y = "frequency (%)") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 11),
    axis.ticks = element_line(color = "black"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    plot.background = element_rect(fill = "transparent", color = NA)
  )



#stacked_bar_plot_pfmdr1 <- ggplot(pfmdr1_data, aes(x = as.factor(Mutation.Position), y = Proportion...., fill = Variant_Type)) +
#geom_bar(stat = "identity", position = "stack", width = 0.9) +  # Narrow bars, remove gap between them
#scale_fill_manual(values = colors_pfmdr1) +  # Apply custom colors for Pfmdr1
#labs(x = "Mutation Position", y = "Proportion (%)") +  # Axis labels
#theme_minimal() +  # Minimal theme for clean look
#theme(
#axis.title.x = element_text(face = "bold", size = 14),  # Bold x-axis title
#axis.title.y = element_text(face = "bold", size = 14),  # Bold y-axis title
#axis.text.x = element_text(size = 12),                 # Adjust x-axis text size
#axis.text.y = element_text(size = 12),                 # Adjust y-axis text size
#legend.title = element_blank(),                        # Remove legend title
#legend.text = element_text(size = 12),                 # Adjust legend text size
#legend.position = "bottom",                             # Position legend at the bottom
#legend.justification = c(0.5, 0),                       # Center the legend at the bottom
#legend.direction = "horizontal",                       # Set legend to be horizontal
#legend.box.spacing = unit(0, "cm"),                    # Remove gap between legend and plot
#plot.title = element_text(hjust = 0, face = "bold", size = 14),  # Title aligned to top left
#plot.background = element_rect(fill = "transparent", color = NA),  # Transparent background
#plot.margin = unit(c(1, 1, 3, 1), "cm")  # Add extra space at the bottom for the legend
#)

# --- 5. Save Stacked Bar Plot for Pfmdr1 ---
# Define save directory
save_dir <- "../output/Bar"
if (!dir.exists(save_dir)) {
  dir.create(save_dir, recursive = TRUE)  # Create the directory if it doesn't exist
}

# Save the bar plot for Pfmdr1
ggsave(
  filename = file.path(save_dir, "Stacked_Bar_Pfmdr1.png"),  # Save file path
  plot = stacked_bar_plot_pfmdr1,                            # Specify the bar plot to save
  width = 4,                                                # Width of the saved image in inches
  height = 4,                                               # Height of the saved image in inches
  dpi = 300,                                                # High resolution for clarity
  bg = "transparent"                                        # Transparent background
)

# Display the Stacked Bar Plot for Pfmdr1 (optional)
print(stacked_bar_plot_pfmdr1)

# compile stacked_bar_plot_pfmdr1 and base_mapmdr1 use centriod for bar 45.827383,13.307708






# Pfk13 Mutations Prevalence bar chart
# --- 1. Load Required Libraries ---
# library(ggplot2)
# library(dplyr)
# 
# Pfk13 <- read.csv("../Data_used/PFK13.csv")
# 
# 
# # --- 3. Define Colors Based on SNP_REPORT ---
# color_labels <- c("WHO reported Validated" = "#993366", 
#                   "Others" = "#6b66c6", 
#                   "WHO reported Candidate" = "#ffa345")
# 
# # Reorder levels so that "Others" comes first
# Pfk13$SNP_REPORT <- factor(Pfk13$SNP_REPORT, levels = c("Others", "WHO reported Validated", "WHO reported Candidate"))
# 
# # --- 4. Create Bar Plot ---
# pfk13_bar_plot <- ggplot(Pfk13, aes(x = reorder(Mutation, -Overall.Prevalence....), 
#                                     y = Overall.Prevalence...., fill = SNP_REPORT)) +
#   geom_bar(stat = "identity", width = 0.8) +
#   scale_fill_manual(values = color_labels) +  # Updated legend without title
#   geom_text(aes(label = `Total.isolate.with.Mutation..n.`, color = SNP_REPORT), 
#             vjust = -0.5, size = 4, fontface = "bold") +  # Make label text bold
#   scale_color_manual(values = color_labels) +  # Ensure labels match bar colors
#   labs(x = "Mutation", y = "Overall Prevalence (%)", title = "Pfk13 Mutations Prevalence") +
#   theme_minimal() +
#   theme(
#     panel.grid = element_blank(), # Remove all grid lines
#     axis.line.x = element_line(color = "black"), 
#     axis.line.y = element_line(color = "black"), 
#     axis.ticks = element_line(color = "black"),
#     axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
#     axis.text.y = element_text(size = 11),
#     axis.title.x = element_text(face = "bold", size = 14),
#     axis.title.y = element_text(face = "bold", size = 14),
#     legend.position = c(0.5, 0.5),  # Center the legend
#     legend.justification = c(0.5, 0.8),  # Center the legend both horizontally and vertically
#     legend.title = element_blank(),  # Remove legend title
#     legend.text = element_text(size = 12),  # Adjust legend text size
#     legend.key.size = unit(0.8, "lines"),  # Adjust legend key size
#     plot.background = element_rect(fill = "white", color = NA),  # Set plot background to white
#     panel.background = element_rect(fill = "white", color = NA)  # Set panel background to white
#   ) +
#   guides(fill = guide_legend(direction = "vertical"))  # Vertical legend layout
# 
# # --- 5. Save the Plot ---
# output_path <- "../output/Pfk13_Bar_Plot.png"
# ggsave(output_path, pfk13_bar_plot, width = 5.5, height = 4.5, dpi = 300)  







# Updated Pfk13 Mutations Prevalence bar chart

# --- 1. Load Required Libraries ---
library(ggplot2)
library(dplyr)

# --- 2. Manually Create Data ---
Pfk13 <- data.frame(
  Mutation = c("R622I", "P441L", "A675V", "E433D", "I352T", 
               "K108E", "K189N", "K189T", "L143V", "L258M", "R255K"),
  Prevalence_Percent = c(10.0, 1.1, 1.7, 1.7, 1.9, 3.0, 2.8, 30.8, 1.7, 3.8, 7.0),
  Count_n = c(57, 6,  10, 10, 11, 17, 16, 177, 10, 22, 40), 
  SNP_REPORT = c("WHO validated ART-R marker", "WHO candidate ART-R marker", "WHO validated ART-R marker", "Unlikely associated with ART-R","Unlikely associated with ART-R",
                 "Unlikely associated with ART-R", "Not associated with ART-R", "Not associated with ART-R", "Unlikely associated with ART-R", "Not associated with ART-R", "Unlikely associated with ART-R")
)

# --- 3. Define Custom Colors ---
color_labels <- c(
  "WHO validated ART-R marker" = "#F8766D", 
  "Unlikely associated with ART-R" = "#6b66c6", 
  "WHO candidate ART-R marker" = "#ffa345",
  "Not associated with ART-R" = "#55a66a"
)

# Reorder factor levels for consistent bar order
Pfk13$SNP_REPORT <- factor(
  Pfk13$SNP_REPORT, 
  levels = c("WHO validated ART-R marker", 
             "WHO candidate ART-R marker", 
             "Not associated with ART-R",
             "Unlikely associated with ART-R"
             )
)

# --- 4. Create Bar Plot ---
pfk13_bar_plot <- ggplot(Pfk13, aes(x = reorder(Mutation, -Prevalence_Percent), 
                                    y = Prevalence_Percent, fill = SNP_REPORT)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = color_labels) +
  geom_text(aes(label = paste0("(", Count_n, ")")), 
            vjust = -0.5, size = 2, fontface = "plain", color = "black") +  # Removed bold
  labs(
    x = expression(italic("Pf")*"K13 variant"),
    y = "Prevalence (%)",
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 9),
    axis.text.y = element_text(size = 11),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(size = 12),
    legend.position = c(0.7, 0.7),
    legend.justification = c(0.5, 0.8),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.8, "lines"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  guides(fill = guide_legend(direction = "vertical"))

# --- 5. Save the Plot ---
#output_png_path <- "../output/Pfk13_Bar_Plot.png"
output_pdf_path <- "../output/Pfk13_Bar_Plot.pdf"

# Save as PNG
#ggsave(output_png_path, pfk13_bar_plot, width = 4, height = 3, dpi = 300)

# Save as PDF
ggsave(output_pdf_path, pfk13_bar_plot, width = 4.5, height = 3.5)

















