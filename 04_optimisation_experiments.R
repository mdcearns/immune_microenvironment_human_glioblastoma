library(tidyverse)
library(dplyr)
library(readxl)
library(ggplot2)
library(scales)

df <- read_excel("GBM_spatial_analysis/data/optimisation_experiments.xlsx")

df <- df %>%
  mutate(across(ends_with('ET'), as.numeric))

df <- df %>%
  mutate(across(ends_with('dilution'), as.numeric))

#pivot longer so for each observation of antibody-experiment you have the corresponding ET and dilution values.
df_long <- df %>%
  pivot_longer(
    cols = -Antibody,
    names_to = c("Experiment", ".value"),
    names_sep = "_"
  )


# Store original antibody order from df BEFORE reshaping
antibody_order <- df$Antibody  

# Apply correct order to long_df
df_long$Antibody <- factor(df_long$Antibody, levels = rev(antibody_order), ordered = TRUE)

# Ensure experiment order correct for final plot
df_long$Experiment <- factor(df_long$Experiment, levels = unique(df_long$Experiment), ordered = TRUE)


# Plot
oe2 <- ggplot(df_long, aes(x = Experiment, y = Antibody)) +
  # Bottom layer: Black borders
  geom_point(aes(size = ET), color = "black", stroke = 1.5, shape = 21, fill = NA, na.rm = TRUE) +
  
  # Top layer: Color-mapped points
  geom_point(aes(color = dilution, size = ET), stroke = 0.5, na.rm = TRUE) +
  
  # Color scale
  scale_color_gradientn(colors = c("red", "orange", "yellow", "seagreen", "deepskyblue", "dodgerblue4", "royalblue4", "midnightblue"),
                        trans = "log", 
                        limits = c(25, 50000),
                        breaks = c(25, 50, 100, 200, 400, 1000, 3000, 12000, 50000),
                        guide = guide_colorbar(
                          trans = "log", 
                          reverse = TRUE,
                          ticks = TRUE,
                          frame.colour = "black",
                          barwidth = 4,
                          barheight = 20,
                          label.position = "right",
                          title.position = "top"
                        ),
                        name = "Dilution factor") +
  
  # Size scale
  scale_size_continuous(range = c(1.5, 7),
                        breaks = c(70, 100, 200, 300),
                        name = "Exposure Time (ms)") +
  
  
  # Labels and theme
  labs(title = "Summary of cyclic immunofluorescence panel optimisation experiments",
       x = "Experiment", y = "Antibody") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    panel.background = element_rect(fill = "grey90", color = NA),  # Set panel background to light grey
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "grey75"),
    panel.grid.minor = element_line(color = "grey75"))


ggsave("GBM_spatial_analysis/outputs/oe2.png", plot = oe2, height = 500, width = 350, units = "mm")

