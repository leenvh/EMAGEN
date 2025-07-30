

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# ==========================
# Site Maps: Figure 1A–1C
# ==========================

# Load required libraries
library(sf)
library(ggplot2)
library(ggrepel)

# ---------------------------------
# Figure 1A: Study Area Base Map
# ---------------------------------

# --- Load shapefiles (Update paths as needed) ---
Africa <- st_read("data/afric.shp")
border <- st_read("data/border.shp")  # Bordering Countries of Ethiopia
Et_region <- st_read("data/et_region.shp")
HF <- st_read("data/Study_Health_Facilties.shp")

# --- Map extent (Ethiopia) ---
min_longitude <- 30.35
max_longitude <- 49.29
min_latitude <- 3
max_latitude <- 15

# --- Create the base map ---
fig1A <- ggplot() +
  geom_sf(data = Africa, fill = "gray90", color = "white", alpha = 0.3) +
  geom_sf(data = Et_region, fill = "white", color = "gray80", alpha = 0.9) +
  geom_sf(data = border, fill = NA, color = "gray50", size = 0.5) +
  geom_text_repel(
    data = border,
    aes(geometry = geometry, label = Name_1),
    stat = "sf_coordinates",
    size = 3,
    color = "gray30"
  ) +
  geom_sf(data = HF, color = "black", size = 2) +
  geom_text_repel(
    data = HF,
    aes(geometry = geometry, label = Woreda),
    stat = "sf_coordinates",
    size = 3,
    color = "black",
    fontface = "bold"
  ) +
  coord_sf(
    xlim = c(min_longitude, max_longitude),
    ylim = c(min_latitude, max_latitude),
    expand = FALSE
  ) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "gray75"),
    panel.grid = element_blank(),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.border = element_rect(color = "gray", fill = NA, size = 1),
    legend.position = "none"
  ) +
  labs(title = "Study Area Map")

# --- Save Figure 1A ---
ggsave("outputs/Fig1A_Study_Area_Map.tiff", fig1A, width = 10, height = 8, dpi = 600, bg = "transparent")



# ---------------------------------
# Figure 1B: API Strata Map
# ---------------------------------

# --- Load API layer ---
API <- st_read("data/API2020_2022.shp")

# --- Harmonize CRS ---
if (st_crs(HF) != st_crs(API)) {
  HF <- st_transform(HF, st_crs(API))
}

# --- Define custom strata colors ---
category_colors <- c(
  "High" = "#d73027",
  "Moderate" = "#cd5e08",
  "Low" = "#2471a3",
  "Very Low" = "#9ed1ed",
  "Free/NA" = "white"
)

API$Strata <- factor(API$Strata, levels = c("High", "Moderate", "Low", "Very Low", "Free/NA"))

# --- Create API map ---
fig1B <- ggplot() +
  geom_sf(data = API, aes(fill = Strata), color = "gray70", size = 0.05) +
  scale_fill_manual(
    values = category_colors,
    na.value = "white",
    name = "Average Annual Parasite\nIncidence (API), 2020–2022",
    labels = c("High (>50)", "Moderate (10–50)", "Low (5–10)", "Very Low (0–5)", "Free/NA (0 or NA)")
  ) +
  geom_sf(data = HF, aes(shape = "Study Health Facilities"), color = "black", size = 3) +
  #geom_sf_text(data = HF, aes(label = HF_Name), size = 6, color = "black", fontface = "bold") +
  #scale_shape_manual(values = c("Study Health Facilities" = 16), guide = guide_legend(order = 2, title = NULL)) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.title = element_blank(),
    legend.position = "right",
    legend.justification = "right",
    legend.box = "vertical",
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

# --- Save Figure 1B ---
ggsave("outputs/Fig1B_API_Strata_Map.png", fig1B, width = 10, height = 7.5, dpi = 300)



# ---------------------------------
# Figure 1C: P. vivax Proportion Map
# ---------------------------------

# --- Load Pv proportion layer ---
Pv_proportion <- st_read("data/Pv_proportions.shp")
Pv_proportion$Pvprop <- ifelse(is.na(Pv_proportion$Pvprop) | Pv_proportion$Pvprop == 0, NA, Pv_proportion$Pvprop)

# Optional: Add 'free' overlay layer if exists
# free <- st_read("data/shapefiles/free_layer.shp")

# --- Create Pv proportion map ---
fig1C <- ggplot() +
  geom_sf(data = Pv_proportion, aes(fill = Pvprop), color = "gray70", size = 0.05) +
  scale_fill_gradientn(
    colors = c("#2471a3", "#9ed1ed", "#cd5e08"),
    na.value = "white"
  ) +
  # Optional white overlay for "free" layer
  # geom_sf(data = free, fill = "white", color = "gray70") +
  geom_sf(data = HF, aes(shape = "Study Health Facilities"), color = "black", size = 2.5) +
  #geom_sf_text(data = HF, aes(label = HF_Name), size = 5.5, color = "black", fontface = "bold") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_blank(),
    legend.position = "right",
    legend.justification = "right",
    axis.line = element_line(color = "black"),
    axis.line.y.right = element_blank(),
    axis.line.x.top = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  labs(fill = expression(italic("P. vivax") ~ "" ~ "proportion (%)")) +
  guides(shape = "none")

# --- Save Figure 1C ---
ggsave("outputs/Fig1C_Pv_Proportion_Map.png", fig1C, width = 8, height = 6, dpi = 300)

