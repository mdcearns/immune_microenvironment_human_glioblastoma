# Load libraries ####

packages <- c("openxlsx", "dplyr", "tibble", "stringr", "tidyr", "purrr", "ggplot2", "scales",
             "cowplot", "forcats", "compositions", "lme4", "lmerTest", "ComplexHeatmap",
             "pheatmap", "viridisLite", "circlize", "SPIAT", "emmeans", "EnhancedVolcano",
             "ggpubr", "RANN", "ClusterR", "interactions", "ggforce", "patchwork",
             "survival", "survminer")

sapply(packages, require, character.only = T)

set.seed(2025)

saveRDS(keyTumour, "GBM_spatial_analysis/outputs/keyTumour.RDS")

# Get thesis data ####

load("GBM_spatial_analysis/outputs/thesis_spatial_environment.RData")
cellData <- readRDS("GBM_spatial_analysis/outputs/cellData.RDS")

keyData <- patientData %>%
  filter(Cohort %in% c("LTS", "STS"))

keyCells <- cellData %>%
  mutate(Patient_ID = str_extract(Annotation_ID, "^\\d{3}")) %>%
  select(Cell_ID, Annotation_ID, Patient_ID, Cohort, Timepoint,
         Annotation_Type, Core_ID, TMA, Classifier,
         Cell.Type, X_coord, Y_coord) %>%
  filter(Cohort %in% c("LTS", "STS"))

keyCells$Cohort <- factor(keyCells$Cohort,
                          levels = c("LTS", "STS"))
keyCells$Timepoint <- factor(keyCells$Timepoint,
                             levels = c("Primary", "Recurrence1", "Recurrence2"))
keyCells$Annotation_Type <- factor(keyCells$Annotation_Type,
                                   levels = c("Tumour", "PBZ", "PPN", "Necrosis"))

# Plot cell type proportions across dataset ####

cellTypes <- c("GBM stem cell", "Astrocyte/GBM cell", "Neuron", # Note adapted later to avoid "/"
               "Endothelial cell", "Vascular smooth muscle cell",
               "Fibroblast", "Non-immune (other)",  "Macrophage",
               "Dendritic cell", "Microglia", "Neutrophil", "NK cell", "B cell",
               "Plasma cell", "CD8+ T cell", "Th cell", "Treg",
               "T cell (other)", "Immune (other)")

keyCells$Cell.Type <- factor(keyCells$Cell.Type, levels = cellTypes)

my_cols <- c(
  "grey50", 
  "#FFFDD0", 
  "#8B4513",
  "#00FA9A",    
  "#006400",   
  "#32CD32",   
  "black",
  "#E60026",
  "#C41E3A", 
  "#FF6F61",
  "#D100D1", 
  "#1E90FF", 
  "#0000CD",    
  "#87CEFA",  
  "#FFF200",
  "#D35400",
  "#FFC300",
  "#FF8C00",
  "#FFB6C1"
)

props <- keyCells %>%
  group_by(Annotation_ID, Cell.Type) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(Annotation_ID) %>%
  mutate(total = sum(count),
         prop = (count / total) * 100) %>%
  ungroup()

cd45_neg <- c("GBM stem cell", "Astrocyte/GBM cell", "Neuron", "Endothelial cell",
                      "Vascular smooth muscle cell", "Fibroblast", "Non-immune (other)")
non_immune <- props %>%
  filter(Cell.Type %in% cd45_neg)

immune <- props %>%
  filter(!Cell.Type %in% cd45_neg)

h1 <- 
  ggplot(data = non_immune, 
         mapping = aes(x = Cell.Type, y = prop, fill = Cell.Type)) + 
  geom_boxplot(outlier.color = "black", outlier.shape = 21, outlier.size = 1.3, color = "black") + 
  labs(y = "Proportion of all cells (%)", title = "Non-immune cell types") +
  theme_classic() +
  theme(
    plot.title = element_text(color = "black", size = 14, hjust = 0.5, face = "bold"),
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white"),
    panel.grid.major = element_line(color = "gray", size = 0.2),
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "black", size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title = element_text(color = "black", size = 12),
    axis.title.x = element_blank(),
    legend.position = "none"
  ) +
  scale_fill_manual(values = my_cols)

non_immune_box_data <- ggplot_build(h1)$data[[1]]
write.xlsx(non_immune_box_data, "GBM_spatial_analysis/outputs/LTS_analysis/non_immune_box_data.xlsx")

ggsave("GBM_spatial_analysis/outputs/LTS_analysis/non_immune_box.png",
       h1,
       height = 4.2, width = 4)

h2 <- 
  ggplot(data = immune, 
         mapping = aes(x = Cell.Type, y = prop, fill = Cell.Type)) + 
  geom_boxplot(outlier.color = "black", outlier.shape = 21, outlier.size = 1.3, color = "black") + 
  labs(y = "Proportion of all cells (%)", title = "Immune cell types") +
  coord_cartesian(ylim = c(0, 18.5)) +  
  theme_classic() +
  theme(
    plot.title = element_text(color = "black", size = 14, hjust = 0.5, face = "bold"),
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white"),
    panel.grid.major = element_line(color = "gray", size = 0.2),
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "black", size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title = element_text(color = "black", size = 12),
    axis.title.x = element_blank(),
    legend.position = "none"
  ) +
  scale_fill_manual(values = my_cols[8:length(my_cols)])

immune_box_data <- ggplot_build(h2)$data[[1]]
write.xlsx(immune_box_data, "GBM_spatial_analysis/outputs/LTS_analysis/immune_box_data.xlsx")

ggsave("GBM_spatial_analysis/outputs/LTS_analysis/immune_box.png",
       h2,
       height = 3.5, width = 6)

rm(props, non_immune, immune, h1, h2, immune_box_data)

# Check number of CD44-expressing cancer cells by annotation

numCells <- keyCells %>%
  filter(Annotation_Type == "Tumour",
         Cell.Type %in% c("GBM stem cell", "Astrocyte/GBM cell")) %>%
  mutate(CD44_pos = ifelse(grepl("CD44", Classifier), "Pos", "Neg")) %>%
  group_by(Annotation_ID, CD44_pos) %>%
  summarise(Count = n(), .groups = "drop") %>%
  pivot_wider(names_from = CD44_pos, values_from = Count, values_fill = 0)

  


# Plot cell type proportions by sample ####

comp <- keyCells %>%
  dplyr::count(Annotation_ID, Patient_ID, Cohort, Timepoint, Annotation_Type, Cell.Type, name = "n_cell.type") %>%
  group_by(Annotation_ID) %>%
  mutate(pct = n_cell.type / sum(n_cell.type) * 100) %>%
  mutate(stemcell_pct = pct[Cell.Type == "GBM stem cell"][1],
         astro_pct = pct[Cell.Type == "Astrocyte/GBM cell"][1]) %>%
  ungroup() %>%
  arrange(Cohort, Timepoint, Annotation_Type, desc(stemcell_pct), desc(astro_pct))

comp$Annotation_ID <- factor(comp$Annotation_ID, levels = unique(comp$Annotation_ID))

subtypesPlot <- ggplot(
  comp, aes(x = Annotation_ID, y = pct, fill = Cell.Type)
) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.1) +
  labs(x = "Sample", fill = "Cell Type") +
  scale_fill_manual(values = my_cols,
                    labels = cellTypes) +
  theme_minimal() +
  theme(axis.text.x =  element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        strip.text = element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        legend.position = "bottom")

ggsave("GBM_spatial_analysis/outputs/LTS_analysis/compositional.png",
       subtypesPlot,
       height = 10, width = 14)

rm(subtypesPlot)

# Cohort: CLR-transformation and LMM compositional analysis ####

# Counts per annotation

perAnn <- keyCells %>%
  filter(Timepoint == "Primary") %>%
  group_by(Annotation_ID, Cell.Type) %>%
  reframe(Patient_ID = unique(Patient_ID),
          Cohort = unique(Cohort),
          Annotation_Type = unique(Annotation_Type),
          n_cell.type = n()) %>%
  group_by(Annotation_ID) %>%
  mutate(Total = sum(n_cell.type)) %>%
  ungroup() %>%
  mutate(pct = n_cell.type / Total) %>%
  select(-c(n_cell.type, Total)) %>%
  arrange(Cohort, Patient_ID, Annotation_Type)
          
# CLR transformation

min_nonzero <- min(perAnn$pct[perAnn$pct > 0])
pseudocount <- min_nonzero / 100

clr_adj <- perAnn %>%
  pivot_wider(names_from = "Cell.Type", values_from = "pct")

clr_adj[, 5:ncol(clr_adj)][is.na(clr_adj[, 5:ncol(clr_adj)])] <- pseudocount
clr_adj[, 5:ncol(clr_adj)] <- clr_adj[, 5:ncol(clr_adj)] / rowSums(clr_adj[, 5:ncol(clr_adj)])

clr_trf <- clr_adj
clr_trf[, 5:ncol(clr_adj)] <- clr(clr_trf[, 5:ncol(clr_adj)])

sub_df_cohort <- data.frame(
  CellType = cellTypes,
  cohort_est = NA,
  cohort_SE = NA,
  cohort_p_value = NA,
  cohort_Upper = NA,
  cohort_Lower = NA,
  pbz_est = NA,
  ppn_est = NA,
  zn_est = NA,
  stringsAsFactors = FALSE
)

for (i in seq_along(cellTypes)) {
  var <- cellTypes[i]
  model_formula <- as.formula(paste0("`", var, "` ~ Cohort + Annotation_Type + (1 | Patient_ID)"))
  model <- lmer(model_formula, data = clr_trf)
  
  coefs <- summary(model)$coefficients
  
  cohort_est     <- coefs["CohortSTS", "Estimate"]
  cohort_SE      <- coefs["CohortSTS", "Std. Error"]
  cohort_p_value <- ifelse(coefs["CohortSTS", "Pr(>|t|)"] < 0.0001, "< 0.0001",
                           paste0("= ", format(signif(coefs["CohortSTS", "Pr(>|t|)"], 2), scientific = FALSE)))
  
  pbz_est     <- coefs["Annotation_TypePBZ", "Estimate"]
  pbz_SE      <- coefs["Annotation_TypePBZ", "Std. Error"]
  pbz_p_value <- ifelse(coefs["Annotation_TypePBZ", "Pr(>|t|)"] < 0.0001, "< 0.0001",
                        paste0("= ", format(signif(coefs["Annotation_TypePBZ", "Pr(>|t|)"], 2), scientific = FALSE)))
  
  ppn_est     <- coefs["Annotation_TypePPN", "Estimate"]
  ppn_SE      <- coefs["Annotation_TypePPN", "Std. Error"]
  ppn_p_value <- ifelse(coefs["Annotation_TypePPN", "Pr(>|t|)"] < 0.0001, "< 0.0001",
                        paste0("= ", format(signif(coefs["Annotation_TypePPN", "Pr(>|t|)"], 2), scientific = FALSE)))
  
  zn_est     <- coefs["Annotation_TypeNecrosis", "Estimate"]
  zn_SE      <- coefs["Annotation_TypeNecrosis", "Std. Error"]
  zn_p_value <- ifelse(coefs["Annotation_TypeNecrosis", "Pr(>|t|)"] < 0.0001, "< 0.0001",
                       paste0("= ", format(signif(coefs["Annotation_TypeNecrosis", "Pr(>|t|)"], 2), scientific = FALSE)))
  
  sub_df_cohort[i, "cohort_est"] <- cohort_est
  sub_df_cohort[i, "cohort_SE"] <- cohort_SE
  sub_df_cohort[i, "cohort_p_value"] <- cohort_p_value
  sub_df_cohort[i, "cohort_Upper"] <- cohort_est + 1.96 * cohort_SE
  sub_df_cohort[i, "cohort_Lower"] <- cohort_est - 1.96 * cohort_SE
  sub_df_cohort[i, "pbz_est"] <- pbz_est
  sub_df_cohort[i, "ppn_est"] <- ppn_est
  sub_df_cohort[i, "zn_est"] <- zn_est
  
}

sub_df_cohort <- sub_df_cohort %>%
  mutate(CellType = factor(CellType, levels = rev(unique(CellType))))

write.xlsx(sub_df_cohort, "GBM_spatial_analysis/outputs/LTS_analysis/primary_cohort_diff_model.xlsx")

# Plot compositional analyses by cohort ####

cohortDiff <- ggplot(sub_df_cohort, aes(x = CellType, y = cohort_est)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1) +
  geom_point(color = my_cols, size = 5) +
  geom_linerange(aes(ymin = cohort_Lower, ymax = cohort_Upper),
                 color = my_cols, linewidth = 1) +
  geom_text(
    aes(label = paste0("p ", cohort_p_value), y = max(cohort_Upper) + 0.5, hjust = 0),
    color = "black",
    size = 5
  ) +
  coord_flip() +
  labs(y = "Centered log-ratio change", x = NULL) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15),
    axis.title.x = element_text(size = 19, hjust = 0.4),
    plot.background = element_rect(fill = "grey90")) +
  scale_y_continuous(limits = c(min(sub_df_cohort$cohort_Lower), max(sub_df_cohort$cohort_Upper) + 1.57),
                     breaks = seq(-3, 2, 1),
                     labels = seq(-3, 2, 1)
  )

ggsave("GBM_spatial_analysis/outputs/LTS_analysis/cohort_diff.png",
       cohortDiff,
       height = 5, width = 8)

# Region: CLR-transformation and LMM compositional analysis ####

sub_df_region <- data.frame(
  CellType = cellTypes,
  pbz_est = NA,
  pbz_SE = NA,
  pbz_p_value = NA,
  pbz_Upper = NA,
  pbz_Lower = NA,
  ppn_est = NA,
  ppn_SE = NA,
  ppn_p_value = NA,
  ppn_Upper = NA,
  ppn_Lower = NA,
  zn_est = NA,
  zn_SE = NA,
  zn_p_value = NA,
  zn_Upper = NA,
  zn_Lower = NA,
  Donor = NA,
  stringsAsFactors = FALSE
)

for (i in seq_along(cellTypes)) {
  var <- cellTypes[i]
  model_formula <- as.formula(paste0("`", var, "` ~ Annotation_Type + (1 | Patient_ID)"))
  model <- lmer(model_formula, data = clr_trf)
  
  pbz_est <- fixef(model)["Annotation_TypePBZ"]
  pbz_SE <- summary(model)$coefficients["Annotation_TypePBZ", "Std. Error"]
  pbz_p_value <- ifelse(summary(model)$coefficients["Annotation_TypePBZ", "Pr(>|t|)"] < 0.0001, "< 0.0001",
                        paste0("= ", format(signif(summary(model)$coefficients["Annotation_TypePBZ", "Pr(>|t|)"], 2), scientific = F)))
  ppn_est <- fixef(model)["Annotation_TypePPN"]
  ppn_SE <- summary(model)$coefficients["Annotation_TypePPN", "Std. Error"]
  ppn_p_value <- ifelse(summary(model)$coefficients["Annotation_TypePPN", "Pr(>|t|)"] < 0.0001, "< 0.0001",
                        paste0("= ", format(signif(summary(model)$coefficients["Annotation_TypePPN", "Pr(>|t|)"], 2), scientific = F)))
  zn_est <- fixef(model)["Annotation_TypeNecrosis"]
  zn_SE <- summary(model)$coefficients["Annotation_TypeNecrosis", "Std. Error"]
  zn_p_value <- ifelse(summary(model)$coefficients["Annotation_TypeNecrosis", "Pr(>|t|)"] < 0.0001, "< 0.0001",
                       paste0("= ", format(signif(summary(model)$coefficients["Annotation_TypeNecrosis", "Pr(>|t|)"], 2), scientific = F)))
  Donor <- attr(summary(model)$varcor$Patient_ID, "stddev")
  
  sub_df_region[i, "pbz_est"] <- pbz_est
  sub_df_region[i, "pbz_SE"] <- pbz_SE
  sub_df_region[i, "pbz_p_value"] <- pbz_p_value
  sub_df_region[i, "pbz_Upper"] <- pbz_est + 1.96 * pbz_SE
  sub_df_region[i, "pbz_Lower"] <- pbz_est - 1.96 * pbz_SE
  
  sub_df_region[i, "ppn_est"] <- ppn_est
  sub_df_region[i, "ppn_SE"] <- ppn_SE
  sub_df_region[i, "ppn_p_value"] <- ppn_p_value
  sub_df_region[i, "ppn_Upper"] <- ppn_est + 1.96 * ppn_SE
  sub_df_region[i, "ppn_Lower"] <- ppn_est - 1.96 * ppn_SE
  
  sub_df_region[i, "zn_est"] <- zn_est
  sub_df_region[i, "zn_SE"] <- zn_SE
  sub_df_region[i, "zn_p_value"] <- zn_p_value
  sub_df_region[i, "zn_Upper"] <- zn_est + 1.96 * zn_SE
  sub_df_region[i, "zn_Lower"] <- zn_est - 1.96 * zn_SE
  
  sub_df_region[i, "Donor"] <- Donor
}

sub_df_region <- sub_df_region %>%
  mutate(CellType = factor(CellType, levels = rev(unique(CellType))))

write.xlsx(sub_df_region, "GBM_spatial_analysis/outputs/LTS_analysis/primary_region_diff_models.xlsx")

# Plot PBZ ####

pbzDiff <- ggplot(sub_df_region, aes(x = CellType, y = pbz_est)) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  geom_point(color = my_cols, size = 5) +
  geom_linerange(aes(ymin = pbz_Lower, ymax = pbz_Upper),
                 color = my_cols, linewidth = 1) +
  geom_text(
    aes(label = paste0("p ", pbz_p_value), y = max(pbz_Upper) + 0.3, hjust = 0),
    color = "black",
    size = 5
  ) +
  coord_flip() +
  labs(y = "Centered log-ratio change", x = NULL) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15),
    axis.title.x = element_text(size = 19, hjust = 0.28),
    plot.background = element_rect(fill = "grey90")) +
  scale_y_continuous(limits = c(min(sub_df_region$pbz_Lower), max(sub_df_region$pbz_Upper) + 1.57),
                     breaks = seq(-3, 4, 1),
                     labels = seq(-3, 4, 1)
  )

ggsave("GBM_spatial_analysis/outputs/LTS_analysis/pbz_diff.png",
       pbzDiff,
       height = 5, width = 8)

# Plot PPN ####

ppnDiff <- ggplot(sub_df_region, aes(x = CellType, y = ppn_est)) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  geom_point(color = my_cols, size = 5) +
  geom_linerange(aes(ymin = ppn_Lower, ymax = ppn_Upper),
                 color = my_cols, linewidth = 1) +
  geom_text(
    aes(label = paste0("p ", ppn_p_value), y = max(ppn_Upper) + 0.5, hjust = 0),
    color = "black",
    size = 5
  ) +
  coord_flip() +
  labs(y = "Centered log-ratio change", x = NULL) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15),
    axis.title.x = element_text(size = 19, hjust = 0.54),
    plot.background = element_rect(fill = "grey90")) +
  scale_y_continuous(
    limits = c(min(sub_df_region$ppn_Lower), max(sub_df_region$ppn_Upper) + 1.57),
    breaks = seq(-7, 5, 1),
    labels = seq(-7, 5, 1)
  )

ggsave("GBM_spatial_analysis/outputs/LTS_analysis/ppn_diff.png",
       ppnDiff,
       height = 5, width = 8)

# Plot ZN ####

znDiff <- ggplot(sub_df_region, aes(x = CellType, y = zn_est)) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  geom_point(color = my_cols, size = 5) +
  geom_linerange(aes(ymin = zn_Lower, ymax = zn_Upper),
                 color = my_cols, linewidth = 1) +
  geom_text(
    aes(label = paste0("p ", zn_p_value), y = max(zn_Upper) + 0.3, hjust = 0),
    color = "black",
    size = 5
  ) +
  coord_flip() +
  labs(y = "Centered log-ratio change", x = NULL) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15),
    axis.title.x = element_text(size = 19, hjust = 0.48),
    plot.background = element_rect(fill = "grey90")) +
  scale_y_continuous(limits = c(min(sub_df_region$zn_Lower), max(sub_df_region$zn_Upper) + 1.57),
                     breaks = seq(-6, 5, 1),
                     labels = seq(-6, 5, 1)
  )

ggsave("GBM_spatial_analysis/outputs/LTS_analysis/zn_diff.png",
       znDiff,
       height = 5, width = 8)

# Recurrence: CLR-transformation and LMM compositional analysis ####

clr_input <- comp %>%
  mutate(pct = pct / 100)

min_nonzero <- min(clr_input$pct[clr_input$pct > 0])
pseudocount <- min_nonzero / 100

clr_adj <- clr_input %>%
  select(-c(n_cell.type, stemcell_pct, astro_pct)) %>%
  pivot_wider(names_from = "Cell.Type", values_from = "pct")

clr_adj[, 6:24][is.na(clr_adj[, 6:24])] <- pseudocount
clr_adj[, 6:24] <- clr_adj[, 6:24] / rowSums(clr_adj[, 6:24])

clr_trf <- clr_adj
clr_trf[, 6:24] <- clr(clr_trf[, 6:24])

clr_trf$Timepoint <- as.character(clr_trf$Timepoint)
clr_trf$Timepoint[clr_trf$Timepoint %in% c("Recurrence1", "Recurrence2")] <- "Recurrence"
clr_trf$Timepoint <- factor(clr_trf$Timepoint,
                            levels = c("Primary", "Recurrence"))

sub_df <- data.frame(
  CellType = cellTypes,
  timepoint_est = NA,
  timepoint_SE = NA,
  timepoint_p_value = NA,
  timepoint_Upper = NA,
  timepoint_Lower = NA,
  Donor = NA,
  stringsAsFactors = FALSE
)

for (i in seq_along(cellTypes)) {
  var <- cellTypes[i]
  model_formula <- as.formula(paste0("`", var, "` ~ Timepoint + Annotation_Type + (1 | Patient_ID)"))
  model <- lmer(model_formula, data = clr_trf)
  
  timepoint_est <- fixef(model)["TimepointRecurrence"]
  timepoint_SE <- summary(model)$coefficients["TimepointRecurrence", "Std. Error"]
  timepoint_p_value <- ifelse(summary(model)$coefficients["TimepointRecurrence", "Pr(>|t|)"] < 0.0001, "< 0.0001",
                           paste0("= ", format(signif(summary(model)$coefficients["TimepointRecurrence", "Pr(>|t|)"], 2), scientific = F)))
  
  Donor <- attr(summary(model)$varcor$Patient_ID, "stddev")
  
  sub_df[i, "timepoint_est"] <- timepoint_est
  sub_df[i, "timepoint_SE"] <- timepoint_SE
  sub_df[i, "timepoint_p_value"] <- timepoint_p_value
  sub_df[i, "timepoint_Upper"] <- timepoint_est + 1.96 * timepoint_SE
  sub_df[i, "timepoint_Lower"] <- timepoint_est - 1.96 * timepoint_SE
  
  sub_df[i, "Donor"] <- Donor
}

sub_df <- sub_df %>%
  mutate(CellType = factor(CellType, levels = rev(unique(CellType))))

write.xlsx(sub_df, "GBM_spatial_analysis/outputs/LTS_analysis/secondary_diff_model.xlsx")

# Plot recurrence ####

timepointDiff <- ggplot(sub_df, aes(x = CellType, y = timepoint_est)) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  geom_point(color = my_cols, size = 5) +
  geom_linerange(aes(ymin = timepoint_Lower, ymax = timepoint_Upper),
                 color = my_cols, linewidth = 1) +
  geom_text(
    aes(label = paste0("p ", timepoint_p_value), y = max(timepoint_Upper) + 0.5, hjust = 0),
    color = "black",
    size = 5
  ) +
  coord_flip() +
  labs(y = "Centered log-ratio change", x = NULL) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15),
    axis.title.x = element_text(size = 19, hjust = 0.37),
    plot.background = element_rect(fill = "grey90")) +
  scale_y_continuous(limits = c(min(sub_df$timepoint_Lower), max(sub_df$timepoint_Upper) + 1.57),
                     breaks = seq(-3, 2, 1),
                     labels = seq(-3, 2, 1)
  )

ggsave("GBM_spatial_analysis/outputs/LTS_analysis/timepoint_diff.png",
       timepointDiff,
       height = 5, width = 8)

# Heatmap of marker expression ####

keyPrim <- keyCells %>%
  filter(Timepoint == "Primary")

markerList <- c(
  "SOX2", "GFAP", "Ki67", "HIF1A",  "MHC Class I",  "CD44", "NeuN", "CD31", "SMA", "Periostin", "CXADR", "CLEC2D", "CD45", 
  "CD68", "CD163", "CD11c", "TMEM119", "CD66b", "CD177", "CD33",  "MHC Class II", "PD-L1", "PD-L2", "CD56", "CD20", 
  "CD38", "CD3", "CD8", "CD4", "FOXP3", "CD45RO", "Granzyme K", "Granzyme A", "Granzyme B", 
  "CD69", "CD103",  "TCF1", "BCL6", "CXCR5", "CXCL13", "PD-1", "CTLA-4",  "LAG-3", "TOX"
)

pctExp <- data.frame(
  matrix(
    NA, nrow = length(cellTypes), ncol = length(markerList)
  )
)

colnames(pctExp) <- markerList
rownames(pctExp) <- cellTypes

IQR_lower <- data.frame(
  matrix(
    NA, nrow = length(cellTypes), ncol = length(markerList)
  )
)

colnames(IQR_lower) <- markerList
rownames(IQR_lower) <- cellTypes

IQR_upper <- data.frame(
  matrix(
    NA, nrow = length(cellTypes), ncol = length(markerList)
  )
)

colnames(IQR_upper) <- markerList
rownames(IQR_upper) <- cellTypes

for (ct in cellTypes){
  
  cat("Processing", ct, "\n")
  
  data <- keyPrim %>%
    filter(Cell.Type == ct) %>%
    group_by(Annotation_ID)
  
  for (marker in markerList){

    result <- data %>%
      summarise(exp = sum(grepl(paste0(marker, "($|:)"), Classifier)),
                total = n(),
                pct = exp / total)
    
    pct_median <- median(result$pct, na.rm = TRUE)
    pct_iqr <- IQR(result$pct, na.rm = TRUE)
    iqr_lower <- quantile(result$pct, 0.25, na.rm = TRUE)
    iqr_upper <- quantile(result$pct, 0.75, na.rm = TRUE)
    
    pctExp[ct, marker] <- pct_median
    IQR_lower[ct, marker] <- iqr_lower
    IQR_upper[ct, marker] <- iqr_upper
      }
}

palette <- inferno(100)[30:100]
col_fun <- colorRamp2(
  seq(0, 1, length.out = length(palette)),
  palette
)

ht <- Heatmap(
  pctExp,
  name = "Proportion", 
  col = col_fun,                   
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 8.5),
  column_names_gp = gpar(fontsize = 8.5),
  column_names_rot = 45,       
  border = FALSE,   
  heatmap_legend_param = list(
    at = c(0, 0.25, 0.5, 0.75, 1),
    labels = c("0%", "25%", "50%", "75%", "100%"),
    title = "Proportion"
  ),
  column_title = "Proportion of marker expression by cell type"
)
ht
png("GBM_spatial_analysis/outputs/LTS_analysis/marker_heatmap.png",
    width = 2800, height = 1400, res = 300) 
draw(ht, padding = unit(c(5, 5, 5, 10), "mm")) # add buffer
dev.off()

pctExp100 <- pctExp * 100
IQR_lower100 <- IQR_lower * 100
IQR_upper100 <- IQR_upper * 100
IQR_diff100 <- IQR_diff * 100

write.xlsx(pctExp100, "GBM_spatial_analysis/outputs/LTS_analysis/pct_exp.xlsx")
write.xlsx(IQR_upper100, "GBM_spatial_analysis/outputs/LTS_analysis/pct_exp_iqr_upper.xlsx")
write.xlsx(IQR_diff100, "GBM_spatial_analysis/outputs/LTS_analysis/pct_exp_iqr_lower.xlsx")

x <- "Th cell"
y <- "TOX"
med <- signif(pctExp100[x, y], 3)
lower <- signif(IQR_lower100[x,y], 3)
upper <- signif(IQR_upper100[x,y], 3)
paste0(x, " ", y, " median ", med, "% (IQR: ", lower, "% – ", upper, "%)")

# Dotplot comparing LTS and STS ####

annots <- unique(keyPrim$Annotation_ID)
results <- data.frame()

for (ct in cellTypes){
  
  for (marker in markerList){
    
    cat("Processing", marker, "expression in cell type:", ct, "\n")
    
    for (a in annots){
        
        annData <- filter(keyPrim,
                          Annotation_ID == a)
        output <- data.frame(
          Annotation_ID = a,
          Patient_ID = unique(annData$Patient_ID),
          Cohort = unique(annData$Cohort),
          Timepoint = unique(annData$Timepoint),
          Annotation_Type = unique(annData$Annotation_Type),
          Cell.Type = ct,
          Marker = marker,
          Pos = sum(annData$Cell.Type == ct & grepl(paste0(marker, "($|:)"), annData$Classifier)),
          Total = sum(annData$Cell.Type == ct)
        )
        results <- rbind(results, output)
    }
  }
  }
  
# Gather results for plotting

results <- results %>%
  filter(Total >= 10)

pctDiff_est <- data.frame(
  matrix(
    NA, nrow = length(cellTypes), ncol = length(markerList)
  )
)
colnames(pctDiff_est) <- markerList
rownames(pctDiff_est) <- cellTypes

pctDiff_p <- data.frame(
  matrix(
    NA, nrow = length(cellTypes), ncol = length(markerList)
  )
)
colnames(pctDiff_p) <- markerList
rownames(pctDiff_p) <- cellTypes

for (ct in cellTypes){
    
  for (marker in markerList){
    
    cat("Processing", marker, "expression in cell type:", ct, "\n")
    
    input <- results %>%
      filter(Cell.Type == ct,
             Marker == marker)
    
    # Check levels
    if(length(unique(input$Cohort)) < 2 ||
       length(unique(input$Annotation_Type)) < 2 ||
       length(unique(input$Pos / input$Total)) <= 1) { # excludes those with no variance e.g. absolute requirements for lineage assignment
      cat("Skipping", ct, marker, "due to insufficient factor levels\n")
      pctDiff_est[ct, marker] <- NA
      pctDiff_p[ct, marker] <- NA
      next
    }
    
    model <- glmer(cbind(Pos, Total - Pos) ~ Cohort + Annotation_Type + (1 | Patient_ID), family = binomial, data = input)
    coefs <- summary(model)$coefficients
    
    if ("CohortSTS" %in% rownames(coefs)) {
      pctDiff_est[ct, marker] <- coefs["CohortSTS", "Estimate"]
      pctDiff_p[ct, marker] <- coefs["CohortSTS", "Pr(>|z|)"]
    } else {
      pctDiff_est[ct, marker] <- NA
      pctDiff_p[ct, marker] <- NA
      }
    }
}

est <- pctDiff_est %>%
  mutate(CellType = rownames(pctDiff_est)) %>%
  pivot_longer(-CellType, names_to = "Marker", values_to = "Logodds")

pval <- pctDiff_p %>%
  mutate(CellType = rownames(pctDiff_p)) %>%
  pivot_longer(-CellType, names_to = "Marker", values_to = "pvalue") %>%
  mutate(padj = p.adjust(pvalue, method = "BH"))

plot_df <- est %>%
  left_join(pval, by = c("CellType", "Marker")) %>%
  filter(pvalue < 0.05) %>%
  mutate(dotsize = -log10(pvalue))

plot_df$CellType <- factor(plot_df$CellType,
                           levels = rev(cellTypes))
plot_df$Marker <- factor(plot_df$Marker,
                           levels = markerList)

palette <- colorRampPalette(
  c("blue",    
    "#E6E0EB",
    "firebrick3" 
  )
)(100)


diff <- ggplot(plot_df, aes(x = Marker, y = CellType)) +
  geom_point(aes(size = dotsize, color = Logodds)) +
  scale_color_gradient2(
    low = palette[1], mid = palette[round(length(palette)/2)], high = palette[length(palette)],
    midpoint = 0
  ) +
  labs(color = "log(OR)", size = "-log10(p)", title = "Differential proportions (LTS vs. STS)") +
  scale_y_discrete(position = "right") +
  scale_size_continuous(range = c(3, 8)) +  # Adjust sizes as you like
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 13),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 20, unit = "pt")
  )

diff
ggsave("GBM_spatial_analysis/outputs/LTS_analysis/diff_dotplot.png",
       diff,
       height = 5, width = 9)


write.xlsx(plot_df, "GBM_spatial_analysis/outputs/LTS_analysis/pct_diff_lts_sts.xlsx")



# Read data into SPE object ####

keySpatial <- cellData %>%
  mutate(Patient_ID = str_extract(Annotation_ID, "^\\d{3}")) %>%
  select(Cell_ID, Annotation_ID, Patient_ID, Cohort, Timepoint,
         Annotation_Type, Core_ID, TMA, Classifier,
         Cell.Type, X_coord, Y_coord, everything()) %>%
  filter(Cohort %in% c("LTS", "STS"),
         Timepoint == "Primary")

keySpatial$Cell.Type[keySpatial$Cell.Type == "Astrocyte/GBM cell"] <- "AstrocyteGBM cell"

cellTypes <- c("GBM stem cell", "AstrocyteGBM cell", "Neuron",
               "Endothelial cell", "Vascular smooth muscle cell",
               "Fibroblast", "Non-immune (other)",  "Macrophage",
               "Dendritic cell", "Microglia", "Neutrophil", "NK cell", "B cell",
               "Plasma cell", "CD8+ T cell", "Th cell", "Treg",
               "T cell (other)", "Immune (other)")

options(future.globals.maxSize = 6.5 * 1024^3)
# Average minimum distance ####

av_min_dist <- data.frame()

for (a in annots){
  
  cat("Calculating for annotation", which(annots == a), "/", length(annots), "\n")
  
  annData <- filter(keySpatial,
                    Annotation_ID == a)
  rownames(annData) <- annData$Cell_ID
  annData$Cell_ID <- NULL
  im <- t(annData[, 14:ncol(annData)]) # may vary
  spe <- format_image_to_spe(format = "general", 
                             intensity_matrix = im,
                             coord_x = annData$X_coord,
                             coord_y = annData$Y_coord,
                             phenotype = NA)
  
  meta_data <- colnames(annData[, 1:9])
  remove <- c("Cell.ID", "Phenotype", "sample_id")  
  colData(spe) <- cbind(colData(spe), annData[, meta_data])
  colData(spe)[, remove] <- NULL
  rownames(spatialCoords(spe)) <- colnames(spe)
  
  # extract coordinates & metadata
  coords <- spatialCoords(spe)       
  meta <- as.data.frame(colData(spe))  
  
 # cellTypes <- unique(meta$Cell.Type) 
  
  # Compute full distance matrix
  dist_mat <- as.matrix(dist(coords))
  cells <- rownames(dist_mat)
  
  distances <- vector()
  
  for (c in cells) {
    
    ref_c <- which(cells == c)
    d <- min(dist_mat[ref_c, -ref_c])
    distances <- c(distances, d)
  }
  
  output <- data.frame(
    Annotation_ID = a,
    Patient_ID = unique(annData$Patient_ID),
    Annotation_Type = unique(annData$Annotation_Type),
    Cohort = unique(annData$Cohort),
    Distance = mean(distances),
    num_cells = length(cells)
  )
  
  av_min_dist <- rbind(av_min_dist, output)
}

write.xlsx(av_min_dist, "GBM_spatial_analysis/outputs/LTS_analysis/distances/av_min_dist.xlsx")

av_min_dist <- read.xlsx(
  "GBM_spatial_analysis/outputs/LTS_analysis/distances/av_min_dist.xlsx"
) %>%
  filter(num_cells >= 100) %>%
  mutate(Group = paste(Cohort, Annotation_Type, sep = " "))

av_min_dist$Group <- gsub("Necrosis", "ZN", av_min_dist$Group)

av_min_dist$Group <- factor(
  av_min_dist$Group,
  levels = c("LTS Tumour", "LTS PBZ", "LTS ZN",   "LTS PPN",
             "STS Tumour", "STS PBZ", "STS ZN", "STS PPN")
)

av_min_dist <- read.xlsx("GBM_spatial_analysis/outputs/LTS_analysis/distances/av_min_dist.xlsx")

results <- av_min_dist %>%
  group_by(Annotation_Type) %>%
  summarise(Median = median(Distance),
            Lower = quantile(Distance, 0.25),
            Upper = quantile(Distance, 0.75))

results

lts_t <-
  t.test(av_min_dist$Distance[av_min_dist$Group == "LTS Tumour"],
         av_min_dist$Distance[av_min_dist$Group == "LTS PBZ"],
         alternative = "two.sided")

lts_t <-
  wilcox.test(av_min_dist$Distance[av_min_dist$Group == "LTS Tumour"],
         av_min_dist$Distance[av_min_dist$Group == "LTS PBZ"])

lts_pvalues <- data.frame(
  group1 = "LTS Tumour",
  group2 = "LTS PBZ",
  y.position = 13.5,
  p = "***" 
)

sts_t <-
  wilcox.test(av_min_dist$Distance[av_min_dist$Group == "STS Tumour"],
         av_min_dist$Distance[av_min_dist$Group == "STS PBZ"])

sts_pvalues <- data.frame(
  group1 = "STS Tumour",
  group2 = "STS PBZ",
  y.position = 13.5,
  p = "****" 
)

min_dist_bp <- ggplot(data = av_min_dist, 
                      mapping = aes(x = Group, y = Distance, fill = Group)) + 
  geom_boxplot(outlier.color = "black", outlier.shape = 21, outlier.size = 1.3, color = "black") + 
  labs(y = bquote("Average minimum distance (" * ~mu * "m)"), title = "Average minimum intercellular distance") +
  theme_classic() +
  theme(
    plot.title = element_text(color = "black", size = 14, hjust = 0.5, face = "bold"),
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white"),
    panel.grid.major = element_line(color = "gray", size = 0.2),
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "black", size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title = element_text(color = "black", size = 12),
    axis.title.x = element_blank(),
    legend.position = "none"
  ) +
  stat_pvalue_manual(data = lts_pvalues, 
                     label = "p",
                     tip.length = 0.01, inherit.aes = F) +
  stat_pvalue_manual(data = sts_pvalues, 
                     label = "p",
                     tip.length = 0.01, inherit.aes = F)
  
min_dist_bp
ggsave("GBM_spatial_analysis/outputs/LTS_analysis/distances/min_dist_bp.png",
       min_dist_bp,
       height = 5, width = 6)

# Extract colours
plot_build <- ggplot_build(min_dist_bp)
colours_used <- unique(plot_build$data[[1]]$fill)

# Pairwise distances ####
output <- data.frame()
annots <- unique(keySpatial$Annotation_ID)

for (a in annots){
  
  cat("Processing annotation:", which(annots == a), "of", length(annots), "\n")
  
  annData <- filter(keySpatial,
                    Annotation_ID == a)
  rownames(annData) <- annData$Cell_ID
  annData$Cell_ID <- NULL
  im <- t(annData[, 14:ncol(annData)]) # may vary
  spe <- format_image_to_spe(format = "general", 
                             intensity_matrix = im,
                             coord_x = annData$X_coord,
                             coord_y = annData$Y_coord,
                             phenotype = NA)
  
  meta_data <- colnames(annData[, 1:9])
  remove <- c("Cell.ID", "Phenotype", "sample_id")  
  colData(spe) <- cbind(colData(spe), annData[, meta_data])
  colData(spe)[, remove] <- NULL
  rownames(spatialCoords(spe)) <- colnames(spe)
  
  # extract coordinates & metadata
  coords <- spatialCoords(spe)       
  meta <- as.data.frame(colData(spe))  
  
  #cellTypes <- unique(meta$Cell.Type) 
  
  # Compute full distance matrix
  dist_mat <- as.matrix(dist(coords))
  cells <- rownames(dist_mat)
  
  for (ref in cellTypes) {
    ref_inds <- which(annData$Cell.Type == ref)
    num_ref <- sum(annData$Cell.Type == ref)
    
    cat("Calculating distances from", ref, "\n")
    
    for (target in cellTypes) {
      
      target_inds <- which(annData$Cell.Type == target)
      
      if (length(ref_inds) > 0 && length(target_inds) > 0) {
        
        # Matrix slice: all ref cells × all target cells
        dist_sub <- dist_mat[ref_inds, target_inds, drop = FALSE]
        
        # Avoid comparisons between the same cell
        if (ref == target) {
          diag(dist_sub) <- NA
        }
        
        min_dist <- min(dist_sub, na.rm = TRUE)
        mean_dist <- mean(dist_sub, na.rm = TRUE)
        
        result <- data.frame(
          Reference = ref,
          num_ref_cells = num_ref,
          Target = target,
          Pair = paste(ref, target, sep = "_"),
          Min = min_dist,
          Mean = mean_dist,
          Annotation_ID = a,
          Patient_ID = unique(annData$Patient_ID),
          Annotation_Type = unique(annData$Annotation_Type),
          Cohort = unique(annData$Cohort)
        )
        
        output <- rbind(output, result)
      }
    }
  }
}

write.xlsx(output, "GBM_spatial_analysis/outputs/LTS_analysis/distances/distances_outputs.xlsx")

anns <- c("Tumour", "PBZ")

dists <- read.xlsx("GBM_spatial_analysis/outputs/LTS_analysis/distances/distances_outputs.xlsx") %>%
  filter(Annotation_Type %in% anns,
         num_ref_cells >= 5) %>%
  mutate(
    Annotation_Type = factor(Annotation_Type, levels = anns),
    Cohort = factor(Cohort),  # Important if not already a factor
    Pair = factor(Pair),      # Ensure Pair is a factor too
    Patient_ID = factor(Patient_ID)  # For safety
    )
model <- lmer(Mean ~ Cohort * Annotation_Type +
                (1 | Patient_ID) +
                (1 | Pair),
              data = dists)
summary(model)

em <- emmeans(model, ~ Cohort | Annotation_Type)
contrast(em, method = "pairwise")

# Run models and extract data

ann_types <- c("Tumour", "PBZ")

for (ann in ann_types){
  
  cat("Calculating for", ann, "samples \n")
  
  filtered_output <- read.xlsx(
    "GBM_spatial_analysis/outputs/LTS_analysis/distances/distances_outputs.xlsx"
  ) %>%
    filter(num_ref_cells >= 5, # pair data only from annotations with min 5 ref cells
           Annotation_Type == ann)
  
  pairs <- filtered_output %>%
    group_by(Pair, Cohort) %>%
    summarise(num_pts = n_distinct(Patient_ID)) %>%
    pivot_wider(
      names_from = Cohort,
      values_from = num_pts
    ) %>%
    filter(LTS >= 5, STS >= 5) %>% # pair data only where min 5 donors per cohort
    pull(Pair)
  
  min_fit_data <- data.frame()
  mean_fit_data <- data.frame()
  
  for (p in pairs){
    df <- filtered_output %>%
      filter(Pair == p)
    
    safe_lmer <- function(formula, data) {
      warning_flag <- FALSE
      fit <- withCallingHandlers(
        tryCatch(
          lmer(formula, data = data),
          error = function(e) {
            message("Model error: ", e$message)
            return(NULL)
          }
        ),
        warning = function(w) {
          message("Model warning: ", conditionMessage(w))
          # If you want to treat some warnings as errors, set flag:
          warning_flag <<- TRUE
          invokeRestart("muffleWarning") # suppress warning output here
        }
      )
      if (warning_flag) return(NULL)
      return(fit)
    }
    
    # Model minimum distances
    
    min_fit <- safe_lmer(Min ~ Cohort + (1 | Patient_ID), data = df)
    if (is.null(min_fit)) {
      cat("Skipping pair", p, "due to model fit issues\n")
      next
    }
    
    emm <- emmeans(min_fit, ~ Cohort)
    contrast <- contrast(emm, method = "pairwise")
    contrast <- as.data.frame(contrast)
    min_fit_output <- data.frame(
      Pair = p,
      est = contrast$estimate,
      pval = contrast$p.value
    )
    min_fit_data <- rbind(min_fit_data, min_fit_output)
    
    # Model mean distances
    
    mean_fit <- safe_lmer(Mean ~ Cohort + (1 | Patient_ID), data = df)
    if (is.null(mean_fit)) {
      cat("Skipping pair", p, "due to model fit issues\n")
      next
    }
    
    emm <- emmeans(mean_fit, ~ Cohort)
    contrast <- contrast(emm, method = "pairwise")
    contrast <- as.data.frame(contrast)
    mean_fit_output <- data.frame(
      Pair = p,
      est = contrast$estimate,
      pval = contrast$p.value
    )
    mean_fit_data <- rbind(mean_fit_data, mean_fit_output)
  }
  
  # Plot minimum distance
  
  min_fit_data$Pair <- gsub("Fibroblast", "Fibro", min_fit_data$Pair)
  min_fit_data$Pair <- gsub("AstrocyteGBM cell", "Astro", min_fit_data$Pair)
  min_fit_data$Pair <- gsub("Th cell", "Th", min_fit_data$Pair)
  min_fit_data$Pair <- gsub("CD8+ T cell", "CD8", min_fit_data$Pair, fixed = T)
  min_fit_data$Pair <- gsub("Macrophage", "MDM", min_fit_data$Pair)
  min_fit_data$Pair <- gsub("Microglia", "MG", min_fit_data$Pair)
  min_fit_data$Pair <- gsub("Dendritic cell", "DC", min_fit_data$Pair)
  min_fit_data$Pair <- gsub("Immune (other)", "Immune", min_fit_data$Pair, fixed = T)
  min_fit_data$Pair <- gsub("Non-immune (other)", "Non-immune", min_fit_data$Pair, fixed = T)
  min_fit_data$Pair <- gsub("GBM stem cell", "GSC", min_fit_data$Pair)
  min_fit_data$Pair <- gsub("NK cell", "NK", min_fit_data$Pair)
  min_fit_data$Pair <- gsub("Endothelial cell", "EC", min_fit_data$Pair)
  min_fit_data$Pair <- gsub("Vascular smooth muscle cell", "VSMC", min_fit_data$Pair)
  
  min_sig_data <- min_fit_data[min_fit_data$pval <= 0.05, ]
  #min_sig_data <- min_sig_data[!grepl("(?<!Non-)Immune", min_sig_data$Pair, ignore.case = TRUE, perl = TRUE), ]
  min_ordered <- min_sig_data[order(abs(min_sig_data$est), decreasing = TRUE), ]
  min_top_labels <- min_ordered$Pair
  
  EnhancedVolcano(min_fit_data,
                  lab = min_fit_data$Pair,
                  selectLab = min_top_labels,
                  x = 'est',
                  y = 'pval',
                  xlab = expression(paste("Difference in average minimum distance (", mu, "m; LTS - STS)")),
                  pCutoff = 0.05,
                  FCcutoff = 0)
  
  p <- EnhancedVolcano(min_fit_data,
                       lab = min_fit_data$Pair,
                       selectLab = min_top_labels,
                       x = 'est',
                       y = 'pval',
                       xlab = expression(paste("Difference in average minimum distance (", mu, "m; LTS - STS)")),
                       pCutoff = 0.05,
                       FCcutoff = 0,
                       pointSize = 2,
                       labSize = 3.5,
                       max.overlaps = Inf,
                       title = paste0("LTS ", ann, " ←→ STS ", ann),
                       subtitle = NULL,
                       boxedLabels = TRUE,
                       legendPosition = "none",
                       drawConnectors = T,
                       colConnectors = "black") +
    theme(plot.title = element_text(hjust = 0.34))
  
  ggsave(paste0("GBM_spatial_analysis/outputs/LTS_analysis/distances/celltype_min_distance_", ann, ".png"),
         p,
         height = 6, width = 8.5)
  
  write.xlsx(min_ordered, paste0("GBM_spatial_analysis/outputs/LTS_analysis/distances/top_hits_min_distance_", ann, ".xlsx"))
  
  # Plot mean distance
  
  mean_fit_data$Pair <- gsub("Fibroblast", "Fibro", mean_fit_data$Pair)
  mean_fit_data$Pair <- gsub("AstrocyteGBM cell", "Astro", mean_fit_data$Pair)
  mean_fit_data$Pair <- gsub("Th cell", "Th", mean_fit_data$Pair)
  mean_fit_data$Pair <- gsub("CD8+ T cell", "CD8", mean_fit_data$Pair, fixed = T)
  mean_fit_data$Pair <- gsub("Macrophage", "MDM", mean_fit_data$Pair)
  mean_fit_data$Pair <- gsub("Microglia", "MG", mean_fit_data$Pair)
  mean_fit_data$Pair <- gsub("Dendritic cell", "DC", mean_fit_data$Pair)
  mean_fit_data$Pair <- gsub("Immune (other)", "Immune", mean_fit_data$Pair, fixed = T)
  mean_fit_data$Pair <- gsub("Non-immune (other)", "Non-immune", mean_fit_data$Pair, fixed = T)
  mean_fit_data$Pair <- gsub("GBM stem cell", "GSC", mean_fit_data$Pair)
  mean_fit_data$Pair <- gsub("NK cell", "NK", mean_fit_data$Pair)
  mean_fit_data$Pair <- gsub("Endothelial cell", "EC", mean_fit_data$Pair)
  mean_fit_data$Pair <- gsub("Vascular smooth muscle cell", "VSMC", mean_fit_data$Pair)
  
  mean_sig_data <- mean_fit_data[mean_fit_data$pval <= 0.05, ]
 # mean_sig_data <- mean_sig_data[!grepl("(?<!Non-)Immune", mean_sig_data$Pair, ignore.case = TRUE, perl = TRUE), ]
  mean_ordered <- mean_sig_data[order(abs(mean_sig_data$est), decreasing = TRUE), ]
  mean_top_labels <- mean_ordered$Pair[1:20]
  
  EnhancedVolcano(mean_fit_data,
                  lab = mean_fit_data$Pair,
                  selectLab = mean_top_labels,
                  x = 'est',
                  y = 'pval',
                  xlab = expression(paste("Difference in mean distance (", mu, "m; LTS - STS)")),
                  pCutoff = 0.05,
                  FCcutoff = 0)
  
  q <- EnhancedVolcano(mean_fit_data,
                       lab = mean_fit_data$Pair,
                       selectLab = mean_top_labels,
                       x = 'est',
                       y = 'pval',
                       xlab = expression(paste("Difference in mean distance (", mu, "m; LTS - STS)")),
                       pCutoff = 0.05,
                       FCcutoff = 0,
                       pointSize = 2,
                       labSize = 3.5,
                       max.overlaps = Inf,
                       title = paste0("LTS ", ann, " ←→ STS ", ann),
                       subtitle = NULL,
                       boxedLabels = TRUE,
                       legendPosition = "none",
                       drawConnectors = T,
                       colConnectors = "black") +
    theme(plot.title = element_text(hjust = 0))
  
  
  ggsave(paste0("GBM_spatial_analysis/outputs/LTS_analysis/distances/celltype_mean_distance_", ann, ".png"),
         q,
         height = 6, width = 8.5)
  
  write.xlsx(mean_ordered, paste0("GBM_spatial_analysis/outputs/LTS_analysis/distances/top_hits_mean_distance_", ann, ".xlsx"))
  
}

# (Check minimum distances of top hits) ####

topHits <- c("Th cell/Neutrophil", "Th cell/Fibroblast", "Neutrophil/CD8+ T cell", "Neutrophil/Dendritic cell")
checkDist <- output %>%
  filter(num_ref_cells <= 5,
         Annotation_Type %in% c("Tumour", "PBZ"),
         Pair %in% topHits) %>%
  group_by(Cohort, Annotation_Type, Pair) %>%
  summarise(n = n(),
            Avg_min = mean(Mean)) %>%
  arrange(Cohort, Pair)
  
  
# Mixing score ####

mix_df <- data.frame()

for (a in annots){
  
  cat("Processing annotation", which(annots == a), "/", length(annots), "\n")
  
  annData <- filter(keySpatial,
                    Annotation_ID == a)
  rownames(annData) <- annData$Cell_ID
  annData$Cell_ID <- NULL
  im <- t(annData[, 14:ncol(annData)]) # may vary
  spe <- format_image_to_spe(format = "general", 
                             intensity_matrix = im,
                             coord_x = annData$X_coord,
                             coord_y = annData$Y_coord,
                             phenotype = NA)
  
  meta_data <- colnames(annData[, 1:9])
  remove <- c("Cell.ID", "Phenotype", "sample_id")  
  colData(spe) <- cbind(colData(spe), annData[, meta_data])
  colData(spe)[, remove] <- NULL
  rownames(spatialCoords(spe)) <- colnames(spe)
  
  for (ref in cellTypes){
    
    for (target in cellTypes){
      
      # Skip if reference and target are the same cell type
      if (ref == target) {
        next
      }
      
      # Check that both cell types are present in this spe
      ref_present <- any(colData(spe)$Cell.Type == ref)
      target_present <- any(colData(spe)$Cell.Type == target)
      
      if (!ref_present || !target_present) {
        cat("Skipping", ref, "-", target, "because one or both not present in this spe.\n")
        next
      }
      
      m <- mixing_score_summary(spe_object = spe, reference_celltype = ref, 
                                target_celltype = target, radius=100, feature_colname ="Cell.Type")
      
      output <- data.frame(
        Reference = m$Reference,
        Target = m$Target,
        Annotation_ID = a,
        Patient_ID = unique(annData$Patient_ID),
        Annotation_Type = unique(annData$Annotation_Type),
        Cohort = unique(annData$Cohort),
        Mix_score = m$Mixing_score,
        Norm_mix_score = m$Normalised_mixing_score,
        num_ref_cells = m$Number_of_reference_cells
      ) %>%
        mutate(Pair = paste(Reference, Target, sep = "_"))
      
      mix_df <- rbind(mix_df, output)
    }
  }
}

write.xlsx(mix_df, "GBM_spatial_analysis/outputs/LTS_analysis/mixing_scores.xlsx")

mix_df <- read.xlsx("GBM_spatial_analysis/outputs/LTS_analysis/mixing_scores.xlsx") %>%
  mutate(Group = paste(Cohort, Annotation_Type, sep = " ")) %>%
  filter(num_ref_cells >= 5,
         Annotation_Type %in% c("Tumour", "PBZ"))

pairs <- unique(mix_df$Pair)
tumour_sig <- data.frame()
pbz_sig <- data.frame()

for (pair in pairs){
  data <- mix_df %>%
    filter(Pair == pair)
  
  n1 <- sum(data$Group == "LTS Tumour")
  n2 <- sum(data$Group == "STS Tumour")
  
  if (n1 >=5 && n2 >= 5){
    tumour_t <-
      t.test(data$Norm_mix_score[data$Group == "LTS Tumour"],
             data$Norm_mix_score[data$Group == "STS Tumour"],
             alternative = "two.sided")
    
    tumour_output <- data.frame(
      Pair = pair,
      pval = tumour_t$p.value
    )
    tumour_sig <- rbind(tumour_sig, tumour_output)
  } else {
    cat("Skipping comparison: not enough observations (n1 =", n1, ", n2 =", n2, ")\n")
  }
  
  n3 <- sum(data$Group == "LTS PBZ")
  n4 <- sum(data$Group == "STS PBZ")
  
  if (n3 >= 5 && n4 >= 5){
    pbz_t <-
      t.test(data$Norm_mix_score[data$Group == "LTS PBZ"],
             data$Norm_mix_score[data$Group == "STS PBZ"],
             alternative = "two.sided")
    
    pbz_output <- data.frame(
      Pair = pair,
      pval = pbz_t$p.value
    )
    
    pbz_sig <- rbind(pbz_sig, pbz_output)
  } else {
    cat("Skipping comparison: not enough observations (n3 =", n3, ", n4 =", n4, ")\n")
  }
}

# Plot for tumour

tumour_sig <- filter(tumour_sig, pval <= 0.05)
tumour_plot_data <- mix_df %>%
  filter(Pair %in% tumour_sig$Pair,
         Annotation_Type == "Tumour")

tumour_plot_data$Pair <- gsub("Macrophage", "MDM", tumour_plot_data$Pair)
tumour_plot_data$Pair <- gsub("Dendritic cell", "DC", tumour_plot_data$Pair)
tumour_plot_data$Pair <- gsub("Vascular smooth muscle cell", "VSMC", tumour_plot_data$Pair)
tumour_plot_data$Pair <- gsub("GBM stem cell", "GSC", tumour_plot_data$Pair)
tumour_plot_data$Pair <- gsub("AstrocyteGBM cell", "Astro", tumour_plot_data$Pair)
tumour_plot_data$Pair <- gsub("Non-immune (other)", "Non-immune", tumour_plot_data$Pair, fixed = T)

tumour_pairs <- c("GSC_Astro", "GSC_Non-immune", "MDM_GSC", "MDM_Astro", "Neutrophil_Non-immune",  "VSMC_MDM", "VSMC_DC", "Neutrophil_MDM")
tumour_plot_data <- filter(tumour_plot_data,
                           Pair %in% tumour_pairs)

tumour_plot_data$Pair <- factor(tumour_plot_data$Pair, levels = tumour_pairs)

pval_df <- data.frame(
  Pair = tumour_pairs,        # matches facet levels exactly
  group1 = c("LTS Tumour"),
  group2 = c("STS Tumour"),
  y.position = c(1.1, 1, 1.2, 1.2, 1.55, 1.2, 0.8, 1.6),
  p.signif = c("*", "**", "*", "*", "**", "*", "*", "*")
)
pval_df$Pair <- factor(pval_df$Pair, levels = tumour_pairs)

tumour_mix_bp <- ggplot(data = tumour_plot_data, 
                        mapping = aes(x = Group, y = Norm_mix_score, fill = Group)) + 
  geom_boxplot(outlier.color = "black", outlier.shape = 21, outlier.size = 1.3, color = "black") + 
  facet_wrap(~ Pair, scales = "free_y", nrow = 2, dir = "v") +
  labs(y = "Normalised mixing score", title = "Cell type mixing scores") +
  theme_classic() +
  theme(
    plot.title = element_text(color = "black", size = 14, hjust = 0.5, face = "bold"),
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white"),
    panel.grid.major = element_line(color = "gray", size = 0.2),
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "black", size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title = element_text(color = "black", size = 14),
    axis.title.x = element_blank(),
    strip.text = element_text(size = 11),
    legend.position = "none"
  ) +
  stat_pvalue_manual(data = pval_df, 
                     label = "p.signif",
                     tip.length = 0.01, inherit.aes = F)
tumour_mix_bp
ggsave("GBM_spatial_analysis/outputs/LTS_analysis/tumour_mix.png",
       tumour_mix_bp,
       height = 5, width = 9)

# Plot for PBZ

#dc_neut, gsc_MG, non-immune_cd8, fibro_cd8

pbz_sig <- filter(pbz_sig, pval <= 0.05)
pbz_plot_data <- mix_df %>%
  filter(Pair %in% pbz_sig$Pair,
         Annotation_Type == "PBZ")

pbz_plot_data$Pair <- gsub("Dendritic cell", "DC", pbz_plot_data$Pair)
pbz_plot_data$Pair <- gsub("GBM stem cell", "GSC", pbz_plot_data$Pair)
pbz_plot_data$Pair <- gsub("Fibroblast", "Fibro", pbz_plot_data$Pair)
pbz_plot_data$Pair <- gsub("Microglia", "MG", pbz_plot_data$Pair)
pbz_plot_data$Pair <- gsub("CD8+ T cell", "CD8", pbz_plot_data$Pair, fixed = T)
pbz_plot_data$Pair <- gsub("Non-immune (other)", "Non-immune", pbz_plot_data$Pair, fixed = T)

pbz_pairs <- c("GSC_MG", "Non-immune_CD8", "Fibro_CD8", "DC_Neutrophil")
pbz_plot_data <- filter(pbz_plot_data,
                        Pair %in% pbz_pairs)

pbz_plot_data$Pair <- factor(pbz_plot_data$Pair, levels = pbz_pairs)

pval_df_pbz <- data.frame(
  Pair = pbz_pairs,        # matches facet levels exactly
  group1 = c("LTS PBZ"),
  group2 = c("STS PBZ"),
  y.position = c(1.6, 2.4, 1.2, 3.8),
  p.signif = c("**", "***", "**", "*")
)
pval_df_pbz$Pair <- factor(pval_df_pbz$Pair, levels = pbz_pairs)

pbz_mix_bp <- ggplot(data = pbz_plot_data, 
                     mapping = aes(x = Group, y = Norm_mix_score, fill = Group)) + 
  geom_boxplot(outlier.color = "black", outlier.shape = 21, outlier.size = 1.3, color = "black") + 
  facet_wrap(~ Pair, scales = "free_y") +
  labs(y = "Normalised mixing score", title = "PBZ mixing scores") +
  theme_classic() +
  theme(
    plot.title = element_text(color = "black", size = 14, hjust = 0.5, face = "bold"),
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white"),
    panel.grid.major = element_line(color = "gray", size = 0.2),
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "black", size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title = element_text(color = "black", size = 14),
    axis.title.x = element_blank(),
    strip.text = element_text(size = 11),
    legend.position = "none"
  ) +
  scale_fill_manual(values = c(colours_used[2], colours_used[6])) +
  stat_pvalue_manual(data = pval_df_pbz, 
                     label = "p.signif",
                     tip.length = 0.01, inherit.aes = F)

ggsave("GBM_spatial_analysis/outputs/LTS_analysis/pbz_mix.png",
       pbz_mix_bp,
       height = 5, width = 5)

# Within radius ####

radius <- c(12, 30, 100)
annots <- unique(keySpatial$Annotation_ID)
output <- data.frame()

for (r in radius){
  
  for (a in annots){
    
    cat("Calculating % cell types within radius =", r, "um ------------- Annotation", which(annots == a), "/", length(annots), "\n")
    annData <- filter(keySpatial, Annotation_ID == a)
    rownames(annData) <- annData$Cell_ID
    annData$Cell_ID <- NULL
    
    # create temp spe object
    
    im <- t(annData[, 14:ncol(annData)])  # intensity matrix
    spe <- format_image_to_spe(format = "general", 
                               intensity_matrix = im,
                               coord_x = annData$X_coord,
                               coord_y = annData$Y_coord,
                               phenotype = NA)
    
    meta_data <- colnames(annData[, 1:9])
    remove <- c("Cell.ID", "Phenotype", "sample_id")  
    colData(spe) <- cbind(colData(spe), annData[, meta_data])
    colData(spe)[, remove] <- NULL
    rownames(spatialCoords(spe)) <- colnames(spe)
    
    # extract coordinates & metadata
    coords <- spatialCoords(spe)       
    meta <- as.data.frame(colData(spe))  
    radius <- 30
    n_cells <- nrow(coords)
    
   # cellTypes <- unique(meta$Cell.Type) 
    
    # Compute full distance matrix
    dist_mat <- as.matrix(dist(coords))
    
    neighbor_pct_df <- data.frame(
      Cell_ID = rownames(meta),
      Reference = meta$Cell.Type,
      matrix(0, nrow = n_cells, ncol = length(cellTypes),
             dimnames = list(NULL, cellTypes)),
      check.names = FALSE  # prevents name modification
    )
    
    # Loop through cells
    
    for (i in seq_len(n_cells)) {
      neighbors <- which(dist_mat[i, ] <= radius & seq_len(n_cells) != i)
      
      if (length(neighbors) > 0) {
        nbr_types <- meta$Cell.Type[neighbors]
        nbr_counts <- table(factor(nbr_types, levels = cellTypes))
        nbr_pct <- nbr_counts / sum(nbr_counts) * 100
        
        # Assign percentages directly by column names
        for (ct in names(nbr_pct)) {
          neighbor_pct_df[i, ct] <- nbr_pct[ct]
        }
      }
    }
    
    results <- neighbor_pct_df %>%
      group_by(Reference) %>%
      summarise(
        num_ref_cells = n(),
        across(
          .cols = -c(Cell_ID),  # all columns except Cell_ID and Cell_Type
          .fns = mean,
          na.rm = TRUE
        )
      ) %>%
      pivot_longer(
        cols = -c(Reference, num_ref_cells),
        names_to = "Target",
        values_to = "Percent"
      ) %>%
      mutate(
        Pair = paste(Reference, Target, sep = "_"),
        Annotation_ID = a,
        Patient_ID = unique(annData$Patient_ID),
        Annotation_Type = unique(annData$Annotation_Type),
        Cohort = unique(annData$Cohort)
      )
    output <- rbind(output, results)
  }
  write.xlsx(output, paste0("GBM_spatial_analysis/outputs/LTS_analysis/celltype_within_", r, "_radius.xlsx"))
}

rm(meta_data, spe, distances, im, input)

# Extract data for plotting

radius <- c(12, 30, 100)
ann_types <- c("Tumour", "PBZ")

for (r in radius){
  
  for (type in ann_types){
    
    cat("Processing", type, "results for radius", r, "\n")
    
    filtered_output <- read.xlsx(paste0(
      "GBM_spatial_analysis/outputs/LTS_analysis/within_radius/celltype_within_", r, "_radius.xlsx")
    ) %>%
      filter(num_ref_cells >= 5,
             !is.na(Percent),
             Annotation_Type == type) %>%
      mutate(Prop = Percent / 100,
             Prop_adj = pmin(pmax(Prop, 1e-5), 1 - 1e-5),
             Logit = log(Prop_adj / (1 - Prop_adj)))
    
    pairs <- filtered_output %>%
      group_by(Pair, Cohort) %>%
      summarise(num_pts = n_distinct(Patient_ID)) %>%
      pivot_wider(
        names_from = Cohort,
        values_from = num_pts
      ) %>%
      filter(LTS >= 5, STS >= 5) %>% # pair data only where min 5 donors per cohort
      pull(Pair)
    
    fit_data <- data.frame()
    
    for (p in pairs){
      df <- filtered_output %>%
        filter(Pair == p)
      
      safe_lmer <- function(formula, data) {
        warning_flag <- FALSE
        fit <- withCallingHandlers(
          tryCatch(
            lmer(formula, data = data),
            error = function(e) {
              message("Model error: ", e$message)
              return(NULL)
            }
          ),
          warning = function(w) {
            message("Model warning: ", conditionMessage(w))
            # If you want to treat some warnings as errors, set flag:
            warning_flag <<- TRUE
            invokeRestart("muffleWarning") # suppress warning output here
          }
        )
        if (warning_flag) return(NULL)
        return(fit)
      }
      
      fit <- safe_lmer(Logit ~ Cohort + (1 | Patient_ID), data = df)
      if (is.null(fit)) {
        cat("Skipping pair", p, "due to model fit issues\n")
        next
      }
      
      emm <- emmeans(fit, ~ Cohort)
      contrast <- contrast(emm, method = "pairwise")
      contrast <- as.data.frame(contrast)
      
      emm_df <- emm %>%
        as.data.frame() %>%
        mutate(Pct = plogis(emmean) * 100)
      
      diff_pct <- emm_df$Pct[emm_df$Cohort == "STS"] - emm_df$Pct[emm_df$Cohort == "LTS"]
      
      output <- data.frame(
        Pair = p,
        est = diff_pct,
        pval = contrast$p.value
      )
      fit_data <- rbind(fit_data, output)
      
    }
    
    fit_data$Pair <- gsub("Fibroblast", "Fibro", fit_data$Pair)
    fit_data$Pair <- gsub("AstrocyteGBM cell", "Astro", fit_data$Pair)
    fit_data$Pair <- gsub("Th cell", "Th", fit_data$Pair)
    fit_data$Pair <- gsub("CD8+ T cell", "CD8", fit_data$Pair, fixed = T)
    fit_data$Pair <- gsub("Macrophage", "MDM", fit_data$Pair)
    fit_data$Pair <- gsub("Microglia", "MG", fit_data$Pair)
    fit_data$Pair <- gsub("Dendritic cell", "DC", fit_data$Pair)
    fit_data$Pair <- gsub("Immune (other)", "Immune", fit_data$Pair, fixed = T)
    fit_data$Pair <- gsub("Non-immune (other)", "Non-immune", fit_data$Pair, fixed = T)
    fit_data$Pair <- gsub("GBM stem cell", "GSC", fit_data$Pair)
    fit_data$Pair <- gsub("NK cell", "NK", fit_data$Pair)
    fit_data$Pair <- gsub("Endothelial cell", "EC", fit_data$Pair)
    fit_data$Pair <- gsub("Vascular smooth muscle cell", "VSMC", fit_data$Pair)
    
    sig_data <- fit_data[fit_data$pval <= 0.05, ]
    specific_celltypes <- sig_data[!grepl("(?<!Non-)Immune", sig_data$Pair, ignore.case = TRUE, perl = TRUE), ]
    ordered <- specific_celltypes[order(abs(specific_celltypes$est), decreasing = TRUE), ]
    top_labels <- ordered$Pair[1:20]
    
    write.xlsx(ordered, paste0(
      "GBM_spatial_analysis/outputs/LTS_analysis/within_radius/top_hits_within_radius_", r, "_", type, ".xlsx"))
    
    EnhancedVolcano(fit_data,
                    lab = fit_data$Pair,
                    selectLab = top_labels,
                    x = 'est',
                    y = 'pval',
                    xlab = expression(paste("Difference in % target cell contribution to reference cell neighbourhood")),
                    pCutoff = 0.05,
                    FCcutoff = 0.1)
    
    p <- EnhancedVolcano(fit_data,
                         lab = fit_data$Pair,
                         selectLab = top_labels,
                         x = 'est',
                         y = 'pval',
                         xlab = bquote("Differential neighbourhood enrichment, " * .(r) * ~mu * "m radius, % (STS - LTS)"),
                         pCutoff = 0.05,
                         FCcutoff = 0,
                         pointSize = 2,
                         labSize = 3.5,
                         max.overlaps = 30,
                         title = paste0("LTS ", type, " ←→ STS ", type),
                         subtitle = NULL,
                         boxedLabels = TRUE,
                         legendPosition = "none",
                         drawConnectors = T,
                         colConnectors = "black") +
      theme(plot.title = element_text(hjust = 0.61))
    
    ggsave(paste0(
      "GBM_spatial_analysis/outputs/LTS_analysis/within_radius/within_radius_", r, "_", type, ".png"),
      p,
      height = 6, width = 8.6)
    
  }
}

# Cellular neighbourhood analysis ####

keyTumour <- keySpatial %>%
  filter(Annotation_Type == "Tumour")

comp_matrix <- list()
annots <- unique(keyTumour$Annotation_ID)

for (a in annots){
  
  cat("Assigning CNs for annotation", which(annots == a), "/", length(annots), "\n")
  
  annData <- filter(keyTumour, Annotation_ID == a)
  
  if (nrow(annData) < 50) {
    cat("Skipping annotation", a, "due to insufficient cell numbers \n")
    next }
  
  rownames(annData) <- annData$Cell_ID
  annData$Cell_ID <- NULL
  
  # create temp spe object
  
  im <- t(annData[, 14:ncol(annData)])  # intensity matrix
  spe <- format_image_to_spe(format = "general", 
                             intensity_matrix = im,
                             coord_x = annData$X_coord,
                             coord_y = annData$Y_coord,
                             phenotype = NA)
  
  meta_data <- colnames(annData[, 1:9])
  remove <- c("Cell.ID", "Phenotype", "sample_id")  
  colData(spe) <- cbind(colData(spe), annData[, meta_data])
  colData(spe)[, remove] <- NULL
  rownames(spatialCoords(spe)) <- colnames(spe)
  
  # extract coordinates
  coords <- spatialCoords(spe) 
  cells <- rownames(coords)

  # Set number of NN
  
  k <- 10
  
  # Compute cell types of k-NN
  nn <- RANN::nn2(coords, k = 1 + k)
  nn_idx <- nn$nn.idx[, -1]
  rownames(nn_idx) <- cells
  neighb <- apply(nn_idx, 2, function(col){annData$Cell.Type[col]})
  rownames(neighb) <- cells
  
  # Adapt into composition matrix of all cell types
  
  comp_matrix[[a]] <- t(apply(neighb, 1, function(row) { # t() because apply() natively returns results as columns
    counts <- table(factor(row, levels = cellTypes))  # count all 19 types in right order
    pct <- counts / length(row) * 100               # as percentages
    as.numeric(pct)                           # return numeric vector
  }))
  colnames(comp_matrix[[a]]) <- cellTypes
}

# Unlist final composition matrix across all annotations
comp_matrix <- do.call(rbind, comp_matrix)

# Set number of CNs for clustering
K <- 9

# Cluster
set.seed(2025)
mbk <- MiniBatchKmeans(data = comp_matrix, clusters = K, batch_size = 1024, num_init = 3)
cluster_assignments <- predict_MBatchKMeans(comp_matrix, mbk$centroids)

CN_k10_K9_df <- data.frame(
  Cell_ID = rownames(comp_matrix),
  CN_k10_K9 = cluster_assignments
)

# Add cluster assignments to metadata
keyTumour <- keyTumour %>%
  left_join(CN_k10_K9_df, by = "Cell_ID") %>%
  select(1:10, CN_k10_K9, everything())

keyTumour$CN_k10_K9 <- paste0("CN", keyTumour$CN_k10_K9)
keyTumour$CN_k10_K9 <- factor(keyTumour$CN_k10_K9,
                              levels = sort(unique(keyTumour$CN_k10_K9)))
keyTumour$Cell.Type <- factor(keyTumour$Cell.Type,
                              levels = cellTypes)

# Plot heatmap of cell types per CN ####

hm_data <- keyTumour %>%
  filter(CN_k10_K9 != "CNNA") %>%
  mutate(Overall_total = n()) %>%
  group_by(Cell.Type) %>%
  mutate(Overall_count = n()) %>%
  ungroup() %>%
  group_by(CN_k10_K9) %>%
  mutate(CN_total = n()) %>%
  group_by(CN_k10_K9, Cell.Type) %>%
  mutate(CN_count = n()) %>%
  ungroup() %>%
  mutate(Overall_pct = Overall_count / Overall_total * 100,
         CN_pct = CN_count / CN_total * 100) %>%
  select(CN_k10_K9, Cell.Type, CN_pct, Overall_pct) %>%
  distinct()

neighbourhoods <- sort(unique(hm_data$CN_k10_K9))
pct_CN <- matrix(NA, nrow = length(neighbourhoods), ncol = length(cellTypes),
                 dimnames = list(neighbourhoods, cellTypes))

for (cn in neighbourhoods){
  
  cat("Calculating composition of", cn, "\n")
  
  cn_data <- hm_data %>%
    filter(CN_k10_K9 == cn)
  
  for (ct in cellTypes){
    
    if (ct %in% unique(cn_data$Cell.Type)){
      
      ct_data <- cn_data %>%
        filter(Cell.Type == ct)
      
      pct_CN[cn, ct] <-
        log2(ct_data$CN_pct / ct_data$Overall_pct)
    } else {
      next
    }
  }
}

pct_CN[is.na(pct_CN)] <- 0

palette <- colorRampPalette(
  c("blue",    
    "#E6E0EB",
    "firebrick3" 
  )
)(100)

ht <- Heatmap(
  pct_CN,
  name = "Enrichment score",
  col = palette,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 8.5),
  column_names_gp = gpar(fontsize = 8.5),
  column_names_rot = 45,
  border = FALSE,
  heatmap_legend_param = list(
    at = c(-6, -3, 0, +3, +6),
    labels = c("-6", "-3", "0", "+3", "+6"),
    title = "Enrichment score"
  ),
  column_title = "Cellular neighbourhood composition"
)

draw(ht)

png("GBM_spatial_analysis/outputs/LTS_analysis/neighbourhoods/CN_k9_K10_heatmap.png",
    width = 1800, height = 1200, res = 300) 
draw(ht, padding = unit(c(10, 15, 15, 15), "mm")) # add buffer
dev.off()

# Extract percentages, enrichment and significance

plotdata <- keyTumour %>%
  filter(CN_k10_K9 != "CNNA") %>%
  group_by(Annotation_ID, Patient_ID, Cohort, CN_k10_K9) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Annotation_ID, Patient_ID, Cohort) %>%
  mutate(Total = sum(Count),
         pct = Count / Total * 100)

cn4_props <- plotdata %>%
  filter(CN_k10_K9 == "CN4")

neighbourhoods <- keyTumour %>%
  filter(CN_k10_K9 != "CNNA") %>%
  pull(CN_k10_K9) %>%
  as.character() %>%
  unique()

neighbourhoods <- neighbourhoods[order(as.numeric(gsub("\\D", "", neighbourhoods)))]
neighbourhoods <- factor(neighbourhoods,
                         levels = neighbourhoods)

cn_model_output <- list()
for (cn in neighbourhoods){
  cn_data <- filter(plotdata, CN_k10_K9 == cn)
  
  cn_model <- glmer(cbind(Count, Total - Count) ~ Cohort + (1 | Patient_ID),
                    family = binomial, data = cn_data)
  coefs <- summary(cn_model)$coefficients
  
  enriched_cells <- colnames(pct_CN)[pct_CN[cn, ] >= 1]
  enriched_cells_string <- paste(enriched_cells, collapse = "_")
  
  df <- as.data.frame(coefs["CohortSTS", c("Estimate", "Std. Error", "Pr(>|z|)")])
  cn_model_output[[cn]] <- t(df) %>%
    as.data.frame() %>%
    mutate(Enriched_cells = enriched_cells_string)
  rownames(cn_model_output[[cn]]) <- cn
}

cn_model_output <- do.call(rbind, cn_model_output)

write.xlsx(cn, "GBM_spatial_analysis/outputs/LTS_analysis/neighbourhoods/CN_k10_K9_models_enriched_cells.xlsx", rowNames = T)
cn <- cn_model_output %>%
  select(Estimate, `Pr(>|z|)`)

# Plot percentages

pval <- data.frame(
  CN_k10_K9 = "CN4", 
  Pair = pbz_pairs,        # matches facet levels exactly
  group1 = "LTS",
  group2 = "STS",
  y.position = 35,
  p.signif = "****"
)

cn_bp <- ggplot(data = plotdata,
                mapping = aes(x = Cohort, y = pct, fill = Cohort)) + 
  geom_boxplot(outlier.color = "black", outlier.shape = 21, outlier.size = 1.3, color = "black") + 
  facet_wrap(~ CN_k10_K9, scales = "free_x", nrow = 1) +
  labs(y = "CN as percentage of total cells (%)") +
  theme_classic() +
  theme(
    plot.title = element_blank(),
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white"),
    panel.grid.major = element_line(color = "gray", size = 0.2),
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "black", size = 11),
    axis.text.x = element_blank(),
    axis.title = element_text(color = "black", size = 12),
    axis.title.x = element_blank(),
    legend.position = "right"
  ) +
  stat_pvalue_manual(data = pval, 
                     label = "p.signif",
                     tip.length = 0.01, inherit.aes = F)

ggsave("GBM_spatial_analysis/outputs/LTS_analysis/neighbourhoods/cn_k10_K9_props.png",
       cn_bp,
       height = 4, width = 8)

# Pairwise distances by CN ####
keyTumour$temp_CN <- ifelse(keyTumour$CN_k10_K9 == "CN4", "CN4", "Other")
keyTumour <- select(keyTumour, 1:10, temp_CN, everything())

output <- data.frame()
annots <- unique(keyTumour$Annotation_ID)

for (a in annots){
  
  cat("Processing annotation:", which(annots == a), "of", length(annots), "\n")
  
  annData <- filter(keyTumour,
                    Annotation_ID == a)
  
  for (cn in c("CN4", "Other")) {
    
    ref_cells <- annData %>%
      filter(temp_CN == cn)
    
    for (ref in cellTypes) {
      
      ref_inds <- ref_cells$Cell_ID[ref_cells$Cell.Type == ref]
      num_ref <- length(ref_inds)
      
      for (target in cellTypes) {
        
        cat("Pairwise distances between", cn, ref, "and all", target, "\n")
        
        target_inds <- annData$Cell_ID[annData$Cell.Type == target]
        
        if (length(ref_inds) == 0 | length(target_inds) == 0) { next }
        
        subset_inds <- c(ref_inds, target_inds)
        cnData <- annData %>%
          filter(Cell_ID %in% subset_inds)
        unique(cnData$Cell.Type)
        
        rownames(cnData) <- cnData$Cell_ID
        cnData$Cell_ID <- NULL
        
        im <- t(cnData[, 14:ncol(cnData)])  # adjust if needed
        spe <- format_image_to_spe(format = "general",
                                   intensity_matrix = im,
                                   coord_x = cnData$X_coord,
                                   coord_y = cnData$Y_coord,
                                   phenotype = NA)
        
        meta_data <- colnames(cnData[, 1:11])
        remove <- c("Cell.ID", "Phenotype", "sample_id")
        colData(spe) <- cbind(colData(spe), cnData[, meta_data])
        colData(spe)[, remove] <- NULL
        rownames(spatialCoords(spe)) <- colnames(spe)
        
        # extract coordinates & metadata
        coords <- spatialCoords(spe)       
        meta <- as.data.frame(colData(spe))  
        
        # Compute full distance matrix
        dist_mat <- as.matrix(dist(coords))
        cells <- rownames(dist_mat)
        
        # Avoid comparisons between the same cell
        if (ref == target) {
          diag(dist_mat) <- NA
        }
        
        min_dist <- min(dist_mat, na.rm = TRUE)
        mean_dist <- mean(dist_mat, na.rm = TRUE)
        
        result <- data.frame(
          CN = cn,
          Reference = ref,
          num_ref_cells = num_ref,
          Target = target,
          Pair = paste(ref, target, sep = "_"),
          Min = min_dist,
          Mean = mean_dist,
          Annotation_ID = a,
          Patient_ID = unique(annData$Patient_ID),
          Annotation_Type = unique(annData$Annotation_Type),
          Cohort = unique(annData$Cohort)
        )
        
        output <- rbind(output, result)
      }
    }
  }
}

write.xlsx(output, "GBM_spatial_analysis/outputs/LTS_analysis/distances/cn_distances_outputs.xlsx")

dists <- read.xlsx("GBM_spatial_analysis/outputs/LTS_analysis/distances/cn_distances_outputs.xlsx") %>%
  filter(num_ref_cells >= 5) %>%
  mutate(
    Annotation_Type = factor(Annotation_Type, levels = anns),
    Cohort = factor(Cohort),  # Important if not already a factor
    Pair = factor(Pair),      # Ensure Pair is a factor too
    Patient_ID = factor(Patient_ID),  # For safety
    CN = factor(CN)
  )
model <- lmer(Mean ~ CN +
                (1 | Patient_ID) +
                (1 | Pair),
              data = dists)
summary(model)

em <- emmeans(model, ~ CN)
contrast(em, method = "pairwise")

# Run models and extract data

# Filter to pairs with at least 5 patients in both CNs
pairs <- dists %>%
  group_by(Pair, CN) %>%
  summarise(num_pts = n_distinct(Patient_ID), .groups = "drop") %>%
  pivot_wider(
    names_from = CN,
    values_from = num_pts
  ) %>%
  filter(CN4 >= 5, Other >= 5) %>%
  pull(Pair)

mean_fit_data <- data.frame()

safe_lmer <- function(formula, data) {
  warning_flag <- FALSE
  fit <- withCallingHandlers(
    tryCatch(
      lmer(formula, data = data),
      error = function(e) {
        message("Model error: ", e$message)
        return(NULL)
      }
    ),
    warning = function(w) {
      message("Model warning: ", conditionMessage(w))
      warning_flag <<- TRUE
      invokeRestart("muffleWarning")
    }
  )
  if (warning_flag) return(NULL)
  return(fit)
}

for (p in pairs) {
  cat("==== Pair:", p, "====\n")
  
  df <- dists %>%
    filter(Pair == p) %>%
    drop_na(Min, Mean)
  
  # Check unique CNs
  cn_levels <- unique(df$CN)
  if (length(cn_levels) < 2) {
    cat("Skipping pair", p, "- fewer than 2 CN levels\n")
    next
  }
  
  # Check patients per CN
  pat_count <- df %>%
    group_by(CN) %>%
    summarise(n = n_distinct(Patient_ID), .groups = "drop")
  
  if (any(pat_count$n < 2)) {
    cat("Skipping pair", p, "- one CN has <2 patients\n")
    print(pat_count)
    next
  }
  
  ## Model mean distances
  mean_fit <- safe_lmer(Mean ~ CN + (1 | Patient_ID), data = df)
  
  if (is.null(mean_fit)) {
    cat("Retrying with lm()...\n")
    try_lm <- try(lm(Mean ~ CN, data = df), silent = TRUE)
    if (inherits(try_lm, "try-error")) {
      cat("Skipping pair", p, "- lm also failed\n")
      next
    }
    emm <- emmeans(try_lm, ~ CN)
  } else {
    emm <- emmeans(mean_fit, ~ CN)
  }
  
  contrast <- contrast(emm, method = "pairwise") %>% as.data.frame()
  mean_fit_data <- rbind(mean_fit_data, data.frame(
    Pair = p,
    est = contrast$estimate,
    pval = contrast$p.value
  ))
}

# Plot mean distance

mean_fit_data$Pair <- gsub("Fibroblast", "Fibro", mean_fit_data$Pair)
mean_fit_data$Pair <- gsub("AstrocyteGBM cell", "Astro", mean_fit_data$Pair)
mean_fit_data$Pair <- gsub("Th cell", "Th", mean_fit_data$Pair)
mean_fit_data$Pair <- gsub("CD8+ T cell", "CD8", mean_fit_data$Pair, fixed = T)
mean_fit_data$Pair <- gsub("Macrophage", "MDM", mean_fit_data$Pair)
mean_fit_data$Pair <- gsub("Microglia", "MG", mean_fit_data$Pair)
mean_fit_data$Pair <- gsub("Dendritic cell", "DC", mean_fit_data$Pair)
mean_fit_data$Pair <- gsub("Immune (other)", "Immune", mean_fit_data$Pair, fixed = T)
mean_fit_data$Pair <- gsub("Non-immune (other)", "Non-immune", mean_fit_data$Pair, fixed = T)
mean_fit_data$Pair <- gsub("GBM stem cell", "GSC", mean_fit_data$Pair)
mean_fit_data$Pair <- gsub("NK cell", "NK", mean_fit_data$Pair)
mean_fit_data$Pair <- gsub("Endothelial cell", "EC", mean_fit_data$Pair)
mean_fit_data$Pair <- gsub("Vascular smooth muscle cell", "VSMC", mean_fit_data$Pair)

mean_sig_data <- mean_fit_data[mean_fit_data$pval <= 0.05, ]
mean_sig_data <- mean_sig_data[!grepl("(?<!Non-)Immune", mean_sig_data$Pair, ignore.case = TRUE, perl = TRUE), ]
mean_ordered <- mean_sig_data[order(abs(mean_sig_data$est), decreasing = TRUE), ]
mean_top_labels <- mean_ordered$Pair[1:20]

EnhancedVolcano(mean_fit_data,
                lab = mean_fit_data$Pair,
                selectLab = mean_top_labels,
                x = 'est',
                y = 'pval',
                xlab = expression(paste("Difference in mean distance (", mu, "m; CN4 - Other)")),
                pCutoff = 0.05,
                FCcutoff = 0)

q <- EnhancedVolcano(mean_fit_data,
                     lab = mean_fit_data$Pair,
                     selectLab = mean_top_labels,
                     x = 'est',
                     y = 'pval',
                     xlab = expression(paste("Difference in mean distance (", mu, "m; CN4 - Other)")),
                     pCutoff = 0.05,
                     FCcutoff = 0,
                     pointSize = 2,
                     labSize = 3.5,
                     max.overlaps = Inf,
                     title = paste0("CN4 ←→ Other"),
                     subtitle = NULL,
                     boxedLabels = TRUE,
                     legendPosition = "none",
                     drawConnectors = T,
                     colConnectors = "black") +
  theme(plot.title = element_text(hjust = 0.89))


ggsave(paste0("GBM_spatial_analysis/outputs/LTS_analysis/distances/cn_celltype_mean_distance.png"),
       q,
       height = 6, width = 8.5)

write.xlsx(mean_ordered, paste0("GBM_spatial_analysis/outputs/LTS_analysis/distances/cn_top_hits_mean_distance.xlsx"))
# Check marker expression ####

cancer_cells <- c("GBM stem cell",
                  "AstrocyteGBM cell",
                  "Non-immune (other)",
                  "Fibroblast")
                  
checkData <- keyTumour %>%
  select(1:22,
         -c(11:21)) %>%
  filter(Cell.Type %in% cancer_cells)

checkData$Cell.Type <- factor(checkData$Cell.Type,
                              levels = cancer_cells)

markerList <- c(
  "SOX2", "GFAP", "Ki67", "HIF1A",  "MHC Class I",  "CD44", "NeuN", "CD31", "SMA", "Periostin", "CXADR", "CLEC2D", "CD45", 
  "CD68", "CD163", "CD11c", "TMEM119", "CD66b", "CD177", "CD33",  "MHC Class II", "PD-L1", "PD-L2", "CD56", "CD20", 
  "CD38", "CD3", "CD8", "CD4", "FOXP3", "CD45RO", "Granzyme K", "Granzyme A", "Granzyme B", 
  "CD69", "CD103",  "TCF1", "BCL6", "CXCR5", "CXCL13", "PD-1", "CTLA-4",  "LAG-3", "TOX"
)

markers <- c("MHC Class I", "Periostin")
  
for (m in markers) {
  
  bpData <- checkData %>%
    mutate(temp_celltype = case_when( 
      grepl(paste0(m, "($|:)"), Classifier) ~ paste0(m, "+ ", Cell.Type),
      !grepl(paste0(m, "($|:)"), Classifier) ~ paste0(m, "- ", Cell.Type),
      TRUE ~ Cell.Type  # leave other cell types unchanged
    )) %>%
    mutate(temp_CN = case_when(
      CN_k10_K9 == "CN4" ~ "CN4",
      CN_k10_K9 != "CN4" ~ "Other",
      TRUE ~ CN_k10_K9
    )) %>%
    group_by(Annotation_ID, temp_CN, Cell.Type) %>%
    mutate(Total = n()) %>%
    group_by(Annotation_ID, temp_CN, temp_celltype) %>%
    reframe(Count = n(),
            Total = Total,
            temp_CN = temp_CN,
            temp_celltype = temp_celltype) %>%
    mutate(pct = Count / Total * 100) %>%
    distinct() %>%
    filter(!grepl(paste0(m, "-"), temp_celltype))
  
  bpData$temp_CN <- factor(bpData$temp_CN,
                           levels = unique(bpData$temp_CN))
  bpData$temp_celltype <- factor(bpData$temp_celltype,
                           levels = c(
                             paste0(m, "+ ", cancer_cells[1]),
                             paste0(m, "+ ", cancer_cells[2]),
                             paste0(m, "+ ", cancer_cells[3]),
                             paste0(m, "+ ", cancer_cells[4])
                             )
  )

  
  if (length(bpData$temp_celltype) == 0) {
    next } else {
      
      pvalues <- data.frame(
        group1 = "CN4",
        group2 = "Other",
        y.position = 73,
        p = "****",
        temp_celltype = c(
          paste0(m, "+ ", cancer_cells[1]),
          paste0(m, "+ ", cancer_cells[2])
        )
      )
      bp <- ggplot(data = bpData,
                   mapping = aes(x = temp_CN, y = pct, fill = temp_CN)) + 
        geom_boxplot(outlier.color = "black", outlier.shape = 21, outlier.size = 1.3, color = "black") + 
        facet_wrap(~ factor(temp_celltype, levels = c(
          paste0(m, "+ ", cancer_cells[1]),
          paste0(m, "+ ", cancer_cells[2])
        )), scales = "free_x", nrow = 1) +
        labs(y = "Percentage of cell type expressing marker") +
        theme_classic() +
        theme(
          plot.title = element_blank(),
          panel.background = element_rect(fill = "white", color = "white"),
          plot.background = element_rect(fill = "white", color = "white"),
          panel.grid.major = element_line(color = "gray", size = 0.2),
          panel.grid.minor = element_blank(),
          axis.text = element_text(color = "black", size = 11),
          axis.text.x = element_blank(),
          axis.title = element_text(color = "black", size = 12),
          axis.title.x = element_blank(),
          legend.position = "right",
          legend.title = element_blank(),
          strip.text = element_text(size = 9)
        ) +
        scale_fill_manual(values = c("#e7298a", "grey80")) +
        stat_pvalue_manual(data = pvalues, 
                           label = "p",
                           tip.length = 0.01, inherit.aes = F)
        
    }

  ggsave(paste0("GBM_spatial_analysis/outputs/LTS_analysis/neighbourhoods/markers/MHC_POSTN/", m, ".png"),
         bp,
         height = 4, width = 5.7)
}

# Test differences

astroData <- bpData %>%
  filter(temp_celltype == "Periostin+ AstrocyteGBM cell") %>%
  mutate(Donor = str_sub(Annotation_ID, 1, 3))

astroData$temp_CN <- factor(astroData$temp_CN,
                            levels = c("Other", "CN4"))
model <- glmer(cbind(Count, Total - Count) ~ temp_CN + (1 | Donor), family = binomial, data = astroData)
coefs <- summary(model)$coefficients
summary(model)
paste0("OR = ", signif(exp(coefs["temp_CNOther", "Estimate"]), 3),
       ", 95% CI: ", signif(exp(coefs["temp_CNOther", "Estimate"] - 1.96 * coefs["temp_CNOther", "Std. Error"]), 3),
                         " - ", signif(exp(coefs["temp_CNOther", "Estimate"] + 1.96 * coefs["temp_CNOther", "Std. Error"]), 3),
                                    ", p = ", signif(coefs["temp_CNOther", "Pr(>|z|)"], 3))
# Force CNs #### 

########## note - not looped

comp_matrix <- list()
annots <- unique(keyTumour$Annotation_ID)

for (a in annots){
  
  cat("Assigning CNs for annotation", which(annots == a), "/", length(annots), "\n")
  
  annData <- filter(keyTumour, Annotation_ID == a)
  
  if (nrow(annData) < 50) {
    cat("Skipping annotation", a, "due to insufficient cell numbers \n")
    next }
  
  rownames(annData) <- annData$Cell_ID
  annData$Cell_ID <- NULL
  
  # create temp spe object
  
  im <- t(annData[, marker_start:marker_end])  # intensity matrix
  spe <- format_image_to_spe(format = "general", 
                             intensity_matrix = im,
                             coord_x = annData$X_coord,
                             coord_y = annData$Y_coord,
                             phenotype = NA)
  
  meta_data <- colnames(annData[, 1:which(colnames(annData) == "Y_coord")])
  remove <- c("Cell.ID", "Phenotype", "sample_id")  
  colData(spe) <- cbind(colData(spe), annData[, meta_data])
  colData(spe)[, remove] <- NULL
  rownames(spatialCoords(spe)) <- colnames(spe)
  
  # extract coordinates
  coords <- spatialCoords(spe) 
  cells <- rownames(coords)
  
  # Set number of NN
  
  k <- 10
  
  # Compute cell types of k-NN
  nn <- RANN::nn2(coords, k = 1 + k)
  nn_idx <- nn$nn.idx[, -1]
  rownames(nn_idx) <- cells
  neighb <- apply(nn_idx, 2, function(col){annData$Cell.Type[col]})
  rownames(neighb) <- cells
  
  # Adapt into composition matrix of all cell types
  
  comp_matrix[[a]] <- t(apply(neighb, 1, function(row) { # t() because apply() natively returns results as columns
    counts <- table(factor(row, levels = cellTypes))  # counts all 19 types in right order
    pct <- counts / length(row) * 100               # as percentages
    as.numeric(pct)                           # return numeric vector
  }))
  colnames(comp_matrix[[a]]) <- cellTypes
}

# Unlist final composition matrix across all annotations
comp_matrix <- do.call(rbind, comp_matrix)

# Set number of CNs for clustering
K <- 27

# Cluster
set.seed(2025)
mbk <- MiniBatchKmeans(data = comp_matrix, clusters = K, batch_size = 1024, num_init = 3)
cluster_assignments <- predict_MBatchKMeans(comp_matrix, mbk$centroids)

df <- data.frame(
  Cell_ID = rownames(comp_matrix),
  CN_k10_K27 = cluster_assignments
)

# Add cluster assignments to metadata
keyTumour <- keyTumour %>%
  left_join(df, by = "Cell_ID") %>%
  select(1:10, CN_k10_K27, everything())

keyTumour$CN_k10_K27 <- paste0("CN", keyTumour$CN_k10_K27)
keyTumour$CN_k10_K27 <- factor(keyTumour$CN_k10_K27,
                               levels = sort(unique(keyTumour$CN_k10_K27)))
keyTumour$Cell.Type <- factor(keyTumour$Cell.Type,
                              levels = cellTypes)

# Plot heatmap of cell types per CN

hm_data <- keyTumour %>%
  filter(CN_k10_K27 != "CNNA") %>%
  mutate(Overall_total = n()) %>%
  group_by(Cell.Type) %>%
  mutate(Overall_count = n()) %>%
  ungroup() %>%
  group_by(CN_k10_K27) %>%
  mutate(CN_total = n()) %>%
  group_by(CN_k10_K27, Cell.Type) %>%
  mutate(CN_count = n()) %>%
  ungroup() %>%
  mutate(Overall_pct = Overall_count / Overall_total * 100,
         CN_pct = CN_count / CN_total * 100) %>%
  select(CN_k10_K27, Cell.Type, CN_pct, Overall_pct) %>%
  distinct()

neighbourhoods <- unique(hm_data$CN_k10_K27)[order(as.numeric(gsub("\\D", "", unique(hm_data$CN_k10_K27))))]
pct_CN <- matrix(NA, nrow = length(neighbourhoods), ncol = length(cellTypes),
                 dimnames = list(neighbourhoods, cellTypes))

for (cn in neighbourhoods){
  
  cat("Calculating composition of", cn, "\n")
  
  cn_data <- hm_data %>%
    filter(CN_k10_K27 == cn)
  
  for (ct in cellTypes){
    
    if (ct %in% unique(cn_data$Cell.Type)){
      
      ct_data <- cn_data %>%
        filter(Cell.Type == ct)
      
      pct_CN[cn, ct] <-
        log2(ct_data$CN_pct / ct_data$Overall_pct)
    } else {
      next
    }
  }
}

pct_CN[is.na(pct_CN)] <- 0

write.xlsx(pct_CN, "GBM_spatial_analysis/outputs/LTS_analysis/neighbourhoods/CN_k10_K27_pct.xlsx", rowNames = TRUE)

palette <- colorRampPalette(
  c("blue",    
    "#E6E0EB",
    "firebrick3" 
  )
)(100)

ht <- Heatmap(
  pct_CN,
  name = "Enrichment score",
  col = palette,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 8.5),
  column_names_gp = gpar(fontsize = 8.5),
  column_names_rot = 45,
  border = FALSE,
  heatmap_legend_param = list(
    at = c(-6, -3, 0, +3, +6),
    labels = c("-6", "-3", "0", "+3", "+6"),
    title = "Enrichment score"
  ),
  column_title = "Cellular neighbourhood composition"
)

draw(ht)

png("GBM_spatial_analysis/outputs/LTS_analysis/neighbourhoods/CN_k9_K27_heatmap.png",
    width = 1800, height = 1900, res = 300) 
draw(ht, padding = unit(c(10, 15, 15, 15), "mm")) # add buffer
dev.off()


# Extract CN percentages, enrichment and significance

plotdata <- keyTumour %>%
  filter(CN_k10_K27 != "CNNA") %>%
  group_by(Annotation_ID, Patient_ID, Cohort, CN_k10_K27) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Annotation_ID, Patient_ID, Cohort) %>%
  mutate(Total = sum(Count),
         pct = Count / Total * 100)

cn_model_output <- list()
for (cn in neighbourhoods){
  cn_data <- filter(plotdata, CN_k10_K27 == cn)
  
  cn_model <- glmer(cbind(Count, Total - Count) ~ Cohort + (1 | Patient_ID),
                    family = binomial, data = cn_data)
  coefs <- summary(cn_model)$coefficients
  
  enriched_cells <- colnames(pct_CN)[pct_CN[cn, ] >= 1]
  enriched_cells_string <- paste(enriched_cells, collapse = "_")
  
  df <- as.data.frame(coefs["CohortSTS", c("Estimate", "Std. Error", "Pr(>|z|)")])
  cn_model_output[[cn]] <- t(df) %>%
    as.data.frame() %>%
    mutate(Enriched_cells = enriched_cells_string)
  rownames(cn_model_output[[cn]]) <- cn
}

cn_model_output <- do.call(rbind, cn_model_output)

write.xlsx(cn_model_output, "GBM_spatial_analysis/outputs/LTS_analysis/neighbourhoods/CN_k10_K27_models_enriched_cells.xlsx", rownames = T)

# (Plot all fibroblast CNs) ####

nonImmune <- data.frame()
p_values <- data.frame()
for (n in c(9, 12, 15, 18, 21, 24, 27, 30)){
  pct <- read.xlsx(paste0("GBM_spatial_analysis/outputs/LTS_analysis/neighbourhoods/CN_k10_K", n, "_pct.xlsx"), rowNames = T)
  colnames(pct) <- cellTypes
  rownames(pct) <- paste0(rownames(pct), " (K = ", n, ")")
  pct <- pct %>%
    filter(`Non-immune (other)` > 1.5 | Fibroblast > 1.5)
  nonImmune <- rbind(nonImmune, pct)
  
  inds <- rownames(pct)
  
  P <- read.xlsx(paste0("GBM_spatial_analysis/outputs/LTS_analysis/neighbourhoods/CN_k10_K", n, "_models_enriched_cells.xlsx"), rowNames = T)
  rownames(P) <- paste0(rownames(P), " (K = ", n, ")")
  P <- P[inds, "Pr(>|z|)", drop = FALSE]
  
  p_values <- rbind(p_values, P)
}

pval_column <- p_values %>%
  pull(`Pr(>|z|)`)

nonImmune <- nonImmune %>%
  mutate(pval = pval_column)

sig <- nonImmune %>%
  filter(pval <= 0.05)
sig$pval <- ifelse(sig$pval < 0.0001, paste("p < 0.0001"), paste0("p = ", signif(sig$pval, 2)))
sig <- sig %>%
  mutate(Axis_labels = paste0(rownames(.), ", ", pval))
rownames(sig) <- sig$Axis_labels
sig <- sig %>%
  select(-c(pval, Axis_labels))

nonsig <- nonImmune %>%
  filter(pval > 0.5)
nonsig$pval <- ifelse(nonsig$pval < 0.0001, paste("p < 0.0001"), paste0("p = ", signif(nonsig$pval, 2)))
nonsig <- nonsig %>%
  mutate(Axis_labels = paste0(rownames(.), ", ", pval))
rownames(nonsig) <- nonsig$Axis_labels
nonsig <- nonsig %>%
  select(-c(pval, Axis_labels))

palette <- colorRampPalette(
  c("blue",    
    "#E6E0EB",
    "firebrick3" 
  )
)(100)

sig_ht <- Heatmap(
  sig,
  name = "Enrichment score",
  col = palette,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 8.5),
  column_names_gp = gpar(fontsize = 8.5),
  column_names_rot = 45,
  border = FALSE,
  heatmap_legend_param = list(
    at = c(-6, -3, 0, +3, +6),
    labels = c("-6", "-3", "0", "+3", "+6"),
    title = "Enrichment score"
  ),
  column_title = "Fibroblast neighbourhoods enriched in LTS"
)

draw(sig_ht)

nonsig_ht <- Heatmap(
  nonsig,
  name = "Enrichment score",
  col = palette,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 8.5),
  column_names_gp = gpar(fontsize = 8.5),
  column_names_rot = 45,
  border = FALSE,
  heatmap_legend_param = list(
    at = c(-6, -3, 0, +3, +6),
    labels = c("-6", "-3", "0", "+3", "+6"),
    title = "Enrichment score"
  ),
  column_title = "Fibroblast neighbourhoods showing no cohort enrichment"
)

draw(nonsig_ht)

png("GBM_spatial_analysis/outputs/LTS_analysis/neighbourhoods/nonImmune_fibro_sig_heatmap.png",
    width = 2600, height = 1900, res = 300) 
draw(sig_ht, padding = unit(c(10, 15, 15, 15), "mm")) # add buffer
dev.off()

png("GBM_spatial_analysis/outputs/LTS_analysis/neighbourhoods/nonImmune_fibro_nonsig_heatmap.png",
    width = 2600, height = 1700, res = 300) 
draw(nonsig_ht, padding = unit(c(10, 15, 15, 15), "mm")) # add buffer
dev.off()

# Exploratory CN analysis #### 

# Re-annotate
m <- "MHC Class I"
cancer_cells <- c("GBM stem cell", "AstrocyteGBM cell", "Non-immune (other)", "Fibroblast")

temp_cellTypes <-
  c("MHC Class I+ GBM stem cell",
    "MHC Class I- GBM stem cell",
    "MHC Class I+ AstrocyteGBM cell",
    "MHC Class I- AstrocyteGBM cell",
    "Neuron", "Endothelial cell", "Vascular smooth muscle cell",
    "MHC Class I+ Fibroblast",
    "MHC Class I- Fibroblast",
    "MHC Class I+ Non-immune (other)",
    "MHC Class I- Non-immune (other)",
    "Macrophage", "Dendritic cell",
    "Microglia", "Neutrophil", "NK cell", "B cell","Plasma cell", "CD8+ T cell", "Th cell", "Treg", "T cell (other)", "Immune (other)")

keyExplore <- keyTumour %>%
  mutate(temp_celltype = case_when( 
    Cell.Type %in% cancer_cells &
      grepl(paste0(m, "($|:)"), Classifier) ~ paste0(m, "+ ", Cell.Type),
    Cell.Type %in% cancer_cells &
      !grepl(paste0(m, "($|:)"), Classifier) ~ paste0(m, "- ", Cell.Type),
    TRUE ~ Cell.Type  # leave other cell types unchanged
  )) %>%
  select(1:10, temp_celltype, everything())

comp_matrix <- list()
annots <- unique(keyExplore$Annotation_ID)

for (a in annots){
  
  cat("Assigning CNs for annotation", which(annots == a), "/", length(annots), "\n")
  
  annData <- filter(keyExplore, Annotation_ID == a)
  
  if (nrow(annData) < 50) {
    cat("Skipping annotation", a, "due to insufficient cell numbers \n")
    next }
  
  rownames(annData) <- annData$Cell_ID
  annData$Cell_ID <- NULL
  
  # create temp spe object
  
  im <- t(annData[, marker_start:marker_end])  # intensity matrix
  spe <- format_image_to_spe(format = "general", 
                             intensity_matrix = im,
                             coord_x = annData$X_coord,
                             coord_y = annData$Y_coord,
                             phenotype = NA)
  
  meta_data <- colnames(annData[, 1:which(colnames(annData) == "Y_coord")])
  remove <- c("Cell.ID", "Phenotype", "sample_id")  
  colData(spe) <- cbind(colData(spe), annData[, meta_data])
  colData(spe)[, remove] <- NULL
  rownames(spatialCoords(spe)) <- colnames(spe)
  
  # extract coordinates
  coords <- spatialCoords(spe) 
  cells <- rownames(coords)
  
  # Set number of NN
  
  k <- 10
  
  # Compute cell types of k-NN
  nn <- RANN::nn2(coords, k = 1 + k)
  nn_idx <- nn$nn.idx[, -1]
  rownames(nn_idx) <- cells
  neighb <- apply(nn_idx, 2, function(col){annData$temp_celltype[col]})
  rownames(neighb) <- cells
  
  # Adapt into composition matrix of all cell types
  
  comp_matrix[[a]] <- t(apply(neighb, 1, function(row) { # t() because apply() natively returns results as columns
    counts <- table(factor(row, levels = temp_cellTypes))  # counts all 19 types in right order
    pct <- counts / length(row) * 100               # as percentages
    as.numeric(pct)                           # return numeric vector
  }))
  colnames(comp_matrix[[a]]) <- temp_cellTypes
}

# Unlist final composition matrix across all annotations
comp_matrix <- do.call(rbind, comp_matrix)

# Set number of CNs for clustering
K <- 50

# Cluster
set.seed(2025)
mbk <- MiniBatchKmeans(data = comp_matrix, clusters = K, batch_size = 1024, num_init = 3)
cluster_assignments <- predict_MBatchKMeans(comp_matrix, mbk$centroids)

df <- data.frame(
  Cell_ID = rownames(comp_matrix),
  CN_MHC = cluster_assignments
)

# Add cluster assignments to metadata
keyExplore <- keyExplore %>%
  left_join(df, by = "Cell_ID") %>%
  select(1:10, CN_MHC, everything())

keyExplore$CN_MHC <- paste0("CN", keyExplore$CN_MHC)
keyExplore$CN_MHC <- factor(keyExplore$CN_MHC,
                            levels = sort(unique(keyExplore$CN_MHC)))
keyExplore$temp_celltype <- factor(keyExplore$temp_celltype,
                                   levels = temp_cellTypes)

# Plot heatmap of cell types per CN

hm_data <- keyExplore %>%
  filter(CN_MHC != "CNNA") %>%
  mutate(Overall_total = n()) %>%
  group_by(temp_celltype) %>%
  mutate(Overall_count = n()) %>%
  ungroup() %>%
  group_by(CN_MHC) %>%
  mutate(CN_total = n()) %>%
  group_by(CN_MHC, temp_celltype) %>%
  mutate(CN_count = n()) %>%
  ungroup() %>%
  mutate(Overall_pct = Overall_count / Overall_total * 100,
         CN_pct = CN_count / CN_total * 100) %>%
  select(CN_MHC, temp_celltype, CN_pct, Overall_pct) %>%
  distinct()

neighbourhoods <- unique(hm_data$CN_MHC)[order(as.numeric(gsub("\\D", "", unique(hm_data$CN_MHC))))]
pct_CN <- matrix(NA, nrow = length(neighbourhoods), ncol = length(temp_cellTypes),
                 dimnames = list(neighbourhoods, temp_cellTypes))

for (cn in neighbourhoods){
  
  cat("Calculating composition of", cn, "\n")
  
  cn_data <- hm_data %>%
    filter(CN_MHC == cn)
  
  for (ct in temp_cellTypes){
    
    if (ct %in% unique(cn_data$temp_celltype)){
      
      ct_data <- cn_data %>%
        filter(temp_celltype == ct)
      
      pct_CN[cn, ct] <-
        log2(ct_data$CN_pct / ct_data$Overall_pct)
    } else {
      next
    }
  }
}

pct_CN[is.na(pct_CN)] <- 0

write.xlsx(pct_CN, "GBM_spatial_analysis/outputs/LTS_analysis/neighbourhoods/CN_MHC_pct.xlsx", rowNames = TRUE)

palette <- colorRampPalette(
  c("blue",    
    "#E6E0EB",
    "firebrick3" 
  )
)(100)

ht <- Heatmap(
  pct_CN,
  name = "Enrichment score",
  col = palette,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 10),
  row_names_max_width = unit(200, "mm"),
  column_names_gp = gpar(fontsize = 10),
  column_names_rot = 45,
  border = FALSE,
  show_heatmap_legend = FALSE,
  heatmap_legend_param = list(
    at = c(-6, -3, 0, +3, +6),
    labels = c("-6", "-3", "0", "+3", "+6"),
    title = "Enrichment score"
  ),
  column_title = "MHC Class I-annotated cellular neighbourhoods"
)

draw(ht)

png("GBM_spatial_analysis/outputs/LTS_analysis/neighbourhoods/CN_MHC_heatmap.png",
    width = 2300, height = 1650, res = 300) 
draw(ht, padding = unit(c(20,35,5,5), "mm")) # add buffer
dev.off()


# Extract CN percentages, enrichment and significance

plotdata <- keyExplore %>%
  filter(CN_MHC != "CNNA") %>%
  group_by(Annotation_ID, Patient_ID, Cohort, CN_MHC) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Annotation_ID, Patient_ID, Cohort) %>%
  mutate(Total = sum(Count),
         pct = Count / Total * 100)

cn_model_output <- list()
for (cn in neighbourhoods){
  cn_data <- filter(plotdata, CN_MHC == cn)
  
  cn_model <- glmer(cbind(Count, Total - Count) ~ Cohort + (1 | Patient_ID),
                    family = binomial, data = cn_data)
  coefs <- summary(cn_model)$coefficients
  
  enriched_cells <- colnames(pct_CN)[pct_CN[cn, ] >= 1]
  enriched_cells_string <- paste(enriched_cells, collapse = "_")
  
  df <- as.data.frame(coefs["CohortSTS", c("Estimate", "Std. Error", "Pr(>|z|)")])
  cn_model_output[[cn]] <- t(df) %>%
    as.data.frame() %>%
    mutate(Enriched_cells = enriched_cells_string)
  rownames(cn_model_output[[cn]]) <- cn
}

cn_model_output <- do.call(rbind, cn_model_output)

write.xlsx(cn_model_output, "GBM_spatial_analysis/outputs/LTS_analysis/neighbourhoods/CN_MHC_models_enriched_cells.xlsx", rownames = T)

View(cn_model_output)

# Plot percentages

pval <- data.frame(
  CN_MHC = "CN5", 
  Pair = pbz_pairs,        # matches facet levels exactly
  group1 = "LTS",
  group2 = "STS",
  y.position = 62,
  p.signif = "****"
)

cn_bp <- ggplot(data = plotdata,
                mapping = aes(x = Cohort, y = pct, fill = Cohort)) + 
  geom_boxplot(outlier.color = "black", outlier.shape = 21, outlier.size = 1.3, color = "black") + 
  facet_wrap(~ CN_MHC, scales = "free_x", nrow = 1) +
  labs(y = "CN as percentage of total cells (%)") +
  theme_classic() +
  theme(
    plot.title = element_blank(),
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white"),
    panel.grid.major = element_line(color = "gray", size = 0.2),
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "black", size = 11),
    axis.text.x = element_blank(),
    axis.title = element_text(color = "black", size = 12),
    axis.title.x = element_blank(),
    legend.position = "right"
  ) +
  stat_pvalue_manual(data = pval, 
                     label = "p.signif",
                     tip.length = 0.01, inherit.aes = F)

ggsave("GBM_spatial_analysis/outputs/LTS_analysis/neighbourhoods/cn_mhc_props.png",
       cn_bp,
       height = 4, width = 8)

keyTumour <- keyExplore

# Analyse interactions driving LTS enrichment of MHC- CNs ####

# Set number of CNs for clustering
keyExplore <- keyTumour

K <- 50

# Cluster
set.seed(2025)
mbk <- MiniBatchKmeans(data = comp_matrix, clusters = K, batch_size = 1024, num_init = 3)
cluster_assignments <- predict_MBatchKMeans(comp_matrix, mbk$centroids)

df <- data.frame(
  Cell_ID = rownames(comp_matrix),
  CN_MHC_50 = cluster_assignments
)

# Add cluster assignments to metadata
keyExplore <- keyExplore %>%
  left_join(df, by = "Cell_ID") %>%
  select(1:10, CN_MHC_50, everything())

keyExplore$CN_MHC_50 <- paste0("CN", keyExplore$CN_MHC_50)
keyExplore$CN_MHC_50 <- factor(keyExplore$CN_MHC_50,
                               levels = sort(unique(keyExplore$CN_MHC_50)))
keyExplore$temp_celltype <- factor(keyExplore$temp_celltype,
                                   levels = temp_cellTypes)

# Plot heatmap of cell types per CN

hm_data <- keyExplore %>%
  filter(CN_MHC_50 != "CNNA") %>%
  mutate(Overall_total = n()) %>%
  group_by(temp_celltype) %>%
  mutate(Overall_count = n()) %>%
  ungroup() %>%
  group_by(CN_MHC_50) %>%
  mutate(CN_total = n()) %>%
  group_by(CN_MHC_50, temp_celltype) %>%
  mutate(CN_count = n()) %>%
  ungroup() %>%
  mutate(Overall_pct = Overall_count / Overall_total * 100,
         CN_pct = CN_count / CN_total * 100) %>%
  select(CN_MHC_50, temp_celltype, CN_pct, Overall_pct) %>%
  distinct()

neighbourhoods <- unique(hm_data$CN_MHC_50)[order(as.numeric(gsub("\\D", "", unique(hm_data$CN_MHC_50))))]
pct_CN <- matrix(NA, nrow = length(neighbourhoods), ncol = length(temp_cellTypes),
                 dimnames = list(neighbourhoods, temp_cellTypes))

for (cn in neighbourhoods){
  
  cat("Calculating composition of", cn, "\n")
  
  cn_data <- hm_data %>%
    filter(CN_MHC_50 == cn)
  
  for (ct in temp_cellTypes){
    
    if (ct %in% unique(cn_data$temp_celltype)){
      
      ct_data <- cn_data %>%
        filter(temp_celltype == ct)
      
      pct_CN[cn, ct] <-
        log2(ct_data$CN_pct / ct_data$Overall_pct)
    } else {
      next
    }
  }
}

pct_CN[is.na(pct_CN)] <- 0

write.xlsx(pct_CN, "GBM_spatial_analysis/outputs/LTS_analysis/neighbourhoods/CN_MHC_50_pct.xlsx", rowNames = TRUE)

palette <- colorRampPalette(
  c("blue",    
    "#E6E0EB",
    "firebrick3" 
  )
)(100)

ht <- Heatmap(
  pct_CN,
  name = "Enrichment score",
  col = palette,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 10),
  row_names_max_width = unit(200, "mm"),
  column_names_gp = gpar(fontsize = 10),
  column_names_rot = 45,
  border = FALSE,
  show_heatmap_legend = FALSE,
  heatmap_legend_param = list(
    at = c(-6, -3, 0, +3, +6),
    labels = c("-6", "-3", "0", "+3", "+6"),
    title = "Enrichment score"
  ),
  column_title = "Cellular neighbourhoods (K = 50)"
)

draw(ht)

png("GBM_spatial_analysis/outputs/LTS_analysis/neighbourhoods/CN_MHC_50_heatmap.png",
    width = 2000, height = 4000, res = 300) 
draw(ht, padding = unit(c(20,35,5,5), "mm")) # add buffer
dev.off()


# Extract CN percentages, enrichment and significance

plotdata <- keyExplore %>%
  filter(CN_MHC_50 != "CNNA") %>%
  group_by(Annotation_ID, Patient_ID, Cohort, CN_MHC_50) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Annotation_ID, Patient_ID, Cohort) %>%
  mutate(Total = sum(Count),
         pct = Count / Total * 100)

cn_model_output <- list()
for (cn in neighbourhoods){
  
  cn_data <- filter(plotdata, CN_MHC_50 == cn)
  
  # Skip if only one cohort present
  if (length(unique(cn_data$Cohort)) < 2) { next }
  
  cn_model <- glmer(cbind(Count, Total - Count) ~ Cohort + (1 | Patient_ID),
                    family = binomial, data = cn_data)
  coefs <- summary(cn_model)$coefficients
  
  enriched_cells <- names(sort(pct_CN[cn, pct_CN[cn, ] >= 1], decreasing = TRUE))
  enriched_cells_string <- paste(enriched_cells, collapse = "_")
  
  df <- as.data.frame(coefs["CohortSTS", c("Estimate", "Std. Error", "Pr(>|z|)")])
  cn_model_output[[cn]] <- t(df) %>%
    as.data.frame() %>%
    mutate(Enriched_cells = enriched_cells_string)
  rownames(cn_model_output[[cn]]) <- cn
}

cn_model_output <- do.call(rbind, cn_model_output)
table(keyExplore$Cohort, keyExplore$CN_MHC_50)
write.xlsx(cn_model_output, "GBM_spatial_analysis/outputs/LTS_analysis/neighbourhoods/CN_MHC_50_models_enriched_cells.xlsx", rownames = T)

temp_cellTypes_tidy <- c(
  "MHC_pos_GBM_stem_cell",
  "MHC_neg_GBM_stem_cell",
  "MHC_pos_AstrocyteGBM_cell",
  "MHC_neg_AstrocyteGBM_cell",
  "Neuron",
  "Endothelial_cell",
  "Vascular_smooth_muscle_cell",
  "MHC_pos_Fibroblast",
  "MHC_neg_Fibroblast",
  "MHC_pos_Non_immune_other",
  "MHC_neg_Non_immune_other",
  "Macrophage",
  "Dendritic_cell",
  "Microglia",
  "Neutrophil",
  "NK_cell",
  "B_cell",
  "Plasma_cell",
  "CD8_T_cell",
  "Th_cell",
  "Treg",
  "T_cell_other",
  "Immune_other"
)

colnames(pct_CN) <- temp_cellTypes_tidy


cn_regression <- data.frame(
  LTS_enriched = cn_model_output$Estimate,
  pct_CN
)

#Look for individual cell types predicting enrichment
results <- data.frame(
  CellType = character(),
  Estimate = numeric(),
  StdError = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for (ct in colnames(cn_regression)[colnames(cn_regression) != "LTS_enriched"]) {
  model <- lm(as.formula(paste("LTS_enriched ~", ct)), data = cn_regression)
  coef_summary <- summary(model)$coefficients
  
  if (nrow(coef_summary) >= 2) {  # Ensure there is a coefficient for the predictor
    results <- rbind(
      results,
      data.frame(
        CellType = ct,
        Estimate = coef_summary[2, "Estimate"],
        StdError = coef_summary[2, "Std. Error"],
        p_value = coef_summary[2, "Pr(>|t|)"]
      )
    )
  }
}

# Plot
overall_results <- results
overall_results$CellType <- temp_cellTypes

overall_cols <- c(my_cols[1],
                  my_cols[2],
                  my_cols[6],
                  "grey20",
                  my_cols[7],
                  my_cols[8],
                  my_cols[13],
                  my_cols[15]
)

overall_results <- overall_results %>%
  filter(p_value <= 0.05) %>%
  mutate(Upper = Estimate + (1.96 * StdError),
         Lower = Estimate - (1.96 * StdError),
         colours = overall_cols) %>%
  arrange(desc(Estimate))

overall_results$CellType <- factor(overall_results$CellType,
                                   levels = overall_results$CellType)
overall_results$colours <- factor(overall_results$colours,
                                  levels = overall_results$colours)

overall_results_plot <- ggplot(overall_results, aes(x = CellType, y = Estimate, color = colours)) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  geom_point(size = 5) +
  geom_linerange(aes(ymin = Lower, ymax = Upper),
                 size = 1) +
  geom_text(
    aes(
      label = ifelse(
        p_value <= 0.0001,
        "p < 0.0001",
        paste0("p = ", signif(p_value, 2))
      ),
      y = max(Upper) + 0.15  # adjust as needed for padding
    ),
    hjust = 0,
    color = "black",
    size = 5
  ) +
  coord_flip() +
  labs(y = "Effect size", x = NULL) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15),
    axis.title.x = element_text(size = 19, hjust = 0.43),
    plot.background = element_rect(fill = "white"),
    legend.position = "none") +
  scale_y_continuous(
    limits = c(min(mhc_results$Lower), max(mhc_results$Upper)),
    breaks = c(-0.5, -0.25, 0, 0.25, 0.5)
  ) +
  scale_color_manual(values = setNames(as.character(overall_results$colours), overall_results$colours))

ggsave("GBM_spatial_analysis/outputs/LTS_analysis/neighbourhoods/stats/celltype_contributions.png",
       overall_results_plot,
       height = 3, width = 7)

write.xlsx(overall_results, 
           "GBM_spatial_analysis/outputs/LTS_analysis/neighbourhoods/stats/data_for_celltype_contributions.xlsx",
           rowNames = F)

# look in CNs containing MHC-negative non-immune cells only
mhc_neg <- cn_model_output %>%
  filter(grepl("MHC Class I-", Enriched_cells)) %>%
  rownames()

cn_regression_mhc <- cn_regression[rownames(cn_regression) %in% mhc_neg, ]

results_mhc <- data.frame(
  CellType = character(),
  Estimate = numeric(),
  StdError = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for (ct in colnames(cn_regression_mhc)[colnames(cn_regression_mhc) != "LTS_enriched"]) {
  model <- lm(as.formula(paste("LTS_enriched ~", ct)), data = cn_regression_mhc)
  coef_summary <- summary(model)$coefficients
  
  if (nrow(coef_summary) >= 2) {  # Ensure there is a coefficient for the predictor
    results_mhc <- rbind(
      results_mhc,
      data.frame(
        CellType = ct,
        Estimate = coef_summary[2, "Estimate"],
        StdError = coef_summary[2, "Std. Error"],
        p_value = coef_summary[2, "Pr(>|t|)"]
      )
    )
  }
}

mhc_results <- results_mhc

mhc_results$CellType <- temp_cellTypes[-18]

mhc_results <- filter(mhc_results, p_value <= 0.05)

mhc_cols <- c(my_cols[4],
              my_cols[7],
              my_cols[8],
              my_cols[10],    
              my_cols[16],
              my_cols[18]
)

mhc_results <- mhc_results %>%
  mutate(Upper = Estimate + (1.96 * StdError),
         Lower = Estimate - (1.96 * StdError),
         colours = mhc_cols) %>%
  arrange(desc(Estimate))

mhc_results$CellType <- factor(mhc_results$CellType,
                               levels = mhc_results$CellType)
mhc_results$colours <- factor(mhc_results$colours,
                              levels = mhc_results$colours)

mhc_results_plot <- ggplot(mhc_results, aes(x = CellType, y = Estimate, color = colours)) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  geom_point(size = 5) +
  geom_linerange(aes(ymin = Lower, ymax = Upper),
                 size = 1) +
  geom_text(
    aes(
      label = ifelse(
        p_value <= 0.0001,
        "p < 0.0001",
        paste0("p = ", signif(p_value, 2))
      ),
      y = max(Upper) + 0.1  # adjust as needed for padding
    ),
    hjust = 0,
    color = "black",
    size = 5
  ) +
  coord_flip() +
  labs(y = "Effect size", x = NULL) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15),
    axis.title.x = element_text(size = 19, hjust = 0.36),
    plot.background = element_rect(fill = "white"),
    legend.position = "none") +
  scale_y_continuous(
    limits = c(min(mhc_results$Lower), max(mhc_results$Upper) + 0.3),
    breaks = c(-0.5, -0.25, 0, 0.25, 0.5, 0.75)
  ) +
  scale_color_manual(values = setNames(as.character(mhc_results$colours), mhc_results$colours))

ggsave("GBM_spatial_analysis/outputs/LTS_analysis/neighbourhoods/stats/mhc_celltype_contributions.png",
       mhc_results_plot,
       height = 3, width = 7)

write.xlsx(mhc_results, 
           "GBM_spatial_analysis/outputs/LTS_analysis/neighbourhoods/stats/data_for_mhc_celltype_contributions.xlsx")

legend_results <- rbind(overall_results, mhc_results) %>%
  select(CellType, colours) %>%
  distinct() %>%
  mutate(CellType = factor(CellType, levels = temp_cellTypes)) %>%
  arrange(CellType)


legend_results$CellType <- as.character(legend_results$CellType)
legend_results$colours <- as.character(legend_results$colours)

# Create named vector of colours: names must match CellType values
colour_map <- setNames(legend_results$colours, legend_results$CellType)

legend_results$CellType <- factor(legend_results$CellType, levels = legend_results$CellType)
legend_results$colours <- factor(legend_results$colours, levels = legend_results$colours)

# Dummy plot for generating the legend
legend_plot <- ggplot(legend_results, aes(x = CellType, y = CellType, fill = CellType)) +
  geom_bar(stat = "identity", show.legend = TRUE) +
  scale_fill_manual(values = colour_map) +
  guides(fill = guide_legend(title = "Cell Type")) +
  theme_void() +
  theme(legend.position = "right")
legend_plot
# Extract and display just the legend
legend_only <- get_legend(legend_plot)
plot_grid(legend_only)

ggsave("GBM_spatial_analysis/outputs/LTS_analysis/neighbourhoods/stats/legend_celltype_contributions.png",
       legend_only,
       height = 4, width = 3)

# Model effect of cell type combinations

temp_cellTypes_tidy

# Extract variable names as strings
var1 <- "MHC_pos_GBM_stem_cell"
var2 <- "Microglia"
ct1 <- temp_cellTypes[which(temp_cellTypes_tidy == var1)]
ct2 <- temp_cellTypes[which(temp_cellTypes_tidy == var2)]

# Build formula
fml <- as.formula(paste("LTS_enriched ~", var1, "*", var2))

# Fit model
model <- lm(fml, data = cn_regression_mhc)
summary(model)
coefs <- summary(model)$coefficients

p_raw <- interact_plot(model,
                       pred = !!sym(var1),
                       modx = !!sym(var2),
                       plot.points = TRUE,
                       interval = TRUE)

# Extract x and y ranges from plot data
gb <- ggplot_build(p_raw)
x_range <- gb$layout$panel_params[[1]]$x.range
y_range <- gb$layout$panel_params[[1]]$y.range

# Position text near top-right (or adjust as needed)
x_pos <- x_range[2] - 0.05 * diff(x_range)
y_neg <- y_range[2] + 0.05 * diff(y_range)

# Add annotation and clean styling
p <- p_raw +
  labs(title = paste("Effect of", ct2),
       y = "LTS ← → STS",
       x = paste0(ct1)) +
  annotate("text", x = x_pos - 1, y = y_neg -0.5, 
           label = paste0("p = ", signif(coefs[paste(var1, var2, sep = ":"), "Pr(>|t|)"], 2)), 
           size = 4) +
  theme(plot.title = element_text(size = 11, hjust = 0.5),
        legend.position = "none")

write.xlsx(coefs, paste0("GBM_spatial_analysis/outputs/LTS_analysis/neighbourhoods/stats/mhc/", var1, "_", var2, ".xlsx"),
           rowNames = T)

ggsave(paste0("GBM_spatial_analysis/outputs/LTS_analysis/neighbourhoods/stats/mhc/", var1, "_", var2, ".png"),
       p,
       height = 4, width = 4)

# Three-variable models

var1 <- "MHC_neg_Non_immune_other"
var2 <- "MHC_pos_GBM_stem_cell"
var3 <- "Endothelial_cell"

ct1 <- temp_cellTypes[which(temp_cellTypes_tidy == var1)]
ct2 <- temp_cellTypes[which(temp_cellTypes_tidy == var2)]
ct3 <- temp_cellTypes[which(temp_cellTypes_tidy == var3)]

fml <- as.formula(paste("LTS_enriched ~", var1, "*", var2, "*", var3))
model <- lm(fml, data = cn_regression_mhc)
summary(model)
coefs <- summary(model)$coefficients
p_raw <- interact_plot(model,
                       pred = !!sym(var1),
                       modx = !!sym(var2),
                       mod2 = !!sym(var3),
                       plot.points = TRUE,
                       interval = TRUE)
gb <- ggplot_build(p_raw)
x_range <- gb$layout$panel_params[[1]]$x.range
y_range <- gb$layout$panel_params[[1]]$y.range
x_pos <- x_range[2] - 0.05 * diff(x_range)
y_neg <- y_range[2] + 0.05 * diff(y_range)

p <- p_raw +
  labs(title = paste0(ct3, " modulation of MHC Class I- GBM stem cell tumour-promoting effect (p = ", signif(coefs[paste(var1, var2, var3, sep = ":"), "Pr(>|t|)"], 2), ")"),
       y = "LTS ← → STS",
       x = paste0(ct1),
       color = ct2) +
  theme(plot.title = element_text(size = 11, hjust = 0.5),
        legend.position = "none")
p
write.xlsx(coefs, paste0("GBM_spatial_analysis/outputs/LTS_analysis/neighbourhoods/stats/mhc/", var1, "_", var2, "_", var3, ".xlsx"),
           rowNames = T)

ggsave(paste0("GBM_spatial_analysis/outputs/LTS_analysis/neighbourhoods/stats/mhc/", var1, "_", var2, "_", var3, ".png"),
       p,
       height = 3.7, width = 8)

# Re-cluster CNs with different markers and analyse interactions ####

# Set terms
# (Note - example only; change separate model functions for two-part and three-part)

ct = "Neutrophil" 
tidy_ct = "Neutrophil" 

combinations <- list( 
   list(m = "CD45RO", tidy_m = "CD45RO")
)

for (combo in combinations){
  
  m <- combo$m
  tidy_m <- combo$tidy_m
  
  cat("Processing", m, "in", ct, "\n")
  temp_colname <- paste0("temp_celltype_", tolower(tidy_ct), "_", tidy_m)
  cn_colname <- paste0("CN_", tolower(tidy_ct), "_", tidy_m)
  
  cells <- c(
    "MHC Class I+ GBM stem cell", "MHC Class I- GBM stem cell",
    "MHC Class I+ AstrocyteGBM cell", "MHC Class I- AstrocyteGBM cell",
    "Neuron", "Endothelial cell", "Vascular smooth muscle cell",
    "MHC Class I+ Fibroblast", "MHC Class I- Fibroblast",
    "MHC Class I+ Non-immune (other)", "MHC Class I- Non-immune (other)",
    "Macrophage", "Dendritic cell", "Microglia", "Neutrophil", "NK cell",
    "B cell", "Plasma cell", "CD8+ T cell",
    "Th cell", "Treg", "T cell (other)", "Immune (other)"
  )
  
  cells_adj <- cells[cells != ct]
  
  temp_cellTypes <- c(paste0(m, "+ ", ct), paste0(m, "- ", ct), cells_adj)
  
  cellsTidy <- c(
    "MHC_pos_GBM_stem_cell",
    "MHC_neg_GBM_stem_cell",
    "MHC_pos_AstrocyteGBM_cell",
    "MHC_neg_AstrocyteGBM_cell",
    "Neuron",
    "Endothelial_cell",
    "Vascular_smooth_muscle_cell",
    "MHC_pos_Fibroblast",
    "MHC_neg_Fibroblast",
    "MHC_pos_Non_immune_other",
    "MHC_neg_Non_immune_other",
    "Macrophage",
    "Dendritic_cell",
    "Microglia",
    "Neutrophil",
    "NK_cell",
    "B_cell",
    "Plasma_cell",
    "CD8_T_cell",
    "Th_cell",
    "Treg",
    "T_cell_other",
    "Immune_other"
  )
  
  cellsTidy_adj <- cellsTidy[cellsTidy != tidy_ct]
  
  temp_cellTypes_tidy <- c(
    paste(tidy_m, "pos", tidy_ct, sep = "_"),
    paste(tidy_m, "neg", tidy_ct, sep = "_"),
    cellsTidy_adj
  )
  
  keyTumour[[temp_colname]] <- NULL
  keyTumour[[cn_colname]] <- NULL
  
  keyTumour[[temp_colname]] <- case_when(
    keyTumour$temp_celltype == ct & grepl(paste0(m, "($|:)"), keyTumour$Classifier) ~ paste0(m, "+ ", keyTumour$temp_celltype),
    keyTumour$temp_celltype == ct & !grepl(paste0(m, "($|:)"), keyTumour$Classifier) ~ paste0(m, "- ", keyTumour$temp_celltype),
    TRUE ~ keyTumour$temp_celltype
  )
  
  keyTumour <- keyTumour %>%
    select(1:10, .data[[temp_colname]], everything())
  
  if (anyNA(keyTumour[[temp_colname]])) {
    cat("NA introduced in", temp_colname, ":", sum(is.na(keyTumour[[temp_colname]])), "\n")
  }
  
  comp_matrix <- list()
  annots <- unique(keyTumour$Annotation_ID)
  
  for (a in annots){
    
    cat("Assigning nearest neighbour windows for annotation", which(annots == a), "/", length(annots), "\n")
    
    annData <- filter(keyTumour, Annotation_ID == a)
    
    if (nrow(annData) < 50) {
      cat("Skipping annotation", a, "due to insufficient cell numbers \n")
      next }
    
    rownames(annData) <- annData$Cell_ID
    annData$Cell_ID <- NULL
    
    # Create temp spe object
    im <- t(annData[, marker_start:marker_end])  # intensity matrix
    
    spe <- format_image_to_spe(
      format = "general", 
      intensity_matrix = im,
      coord_x = annData$X_coord,
      coord_y = annData$Y_coord,
      phenotype = NA
    )
    
    # Subset annData to match cells retained in spe
    retained_cells <- colnames(spe)
    annData <- annData[rownames(annData) %in% retained_cells, , drop = FALSE]
    
    # Continue as before
    meta_data <- colnames(annData[, 1:which(colnames(annData) == "Y_coord")])
    remove <- c("Cell.ID", "Phenotype", "sample_id")
    
    colData(spe) <- cbind(colData(spe), annData[, meta_data])
    colData(spe)[, remove] <- NULL
    rownames(spatialCoords(spe)) <- colnames(spe)
    
    # extract coordinates
    coords <- spatialCoords(spe) 
    cells <- rownames(coords)
    
    # Set number of NN
    
    k <- 10
    
    # Compute cell types of k-NN
    nn <- RANN::nn2(coords, k = 1 + k)
    nn_idx <- nn$nn.idx[, -1]
    rownames(nn_idx) <- cells
    neighb <- apply(nn_idx, 2, function(col){annData[[temp_colname]][col]})
    rownames(neighb) <- cells
    
    # Adapt into composition matrix of all cell types
    
    comp_matrix[[a]] <- t(apply(neighb, 1, function(row) {
      counts <- table(factor(row, levels = temp_cellTypes))
      pct <- counts / length(row) * 100
      as.numeric(pct)
    }))
    colnames(comp_matrix[[a]]) <- temp_cellTypes
  }
  
  if (any(sapply(neighb, function(x) any(is.na(x))))) {
    cat("NA in CN assignments for", temp_colname, "\n")
  }
  
  if (any(sapply(comp_matrix, function(x) any(is.na(x))))) {
    cat("NA in comp_matrix for", temp_colname, "\n")
  }
  
  # Unlist final composition matrix across all annotations
  comp_matrix_full <- do.call(rbind, comp_matrix)
  
  # Set number of CNs for clustering
  K <- 50
  
  cat("Clustering CNs... \n")
  
  # Cluster
  set.seed(2025)
  mbk <- MiniBatchKmeans(data = comp_matrix_full, clusters = K, batch_size = 1024, num_init = 3)
  cluster_assignments <- predict_MBatchKMeans(comp_matrix_full, mbk$centroids)
  
  df <- data.frame(Cell_ID = rownames(comp_matrix_full), cluster = paste0("CN", cluster_assignments))
  colnames(df)[2] <- cn_colname
  
  # Add cluster assignments to metadata
  
  keyTumour <- keyTumour %>%
    left_join(df, by = "Cell_ID") %>%
    select(1:10, .data[[cn_colname]], everything())
  
  keyTumour[[cn_colname]] <- factor(keyTumour[[cn_colname]],
                                    levels = sort(unique( keyTumour[[cn_colname]])))
  keyTumour[[temp_colname]] <- factor(keyTumour[[temp_colname]],
                                      levels = temp_cellTypes)
  
  # Plot heatmap of cell types per CN
  
  hm_data <- keyTumour %>%
    filter(.data[[cn_colname]] != "CNNA") %>%
    mutate(Overall_total = n()) %>%
    group_by(.data[[temp_colname]]) %>%
    mutate(Overall_count = n()) %>%
    ungroup() %>%
    group_by(.data[[cn_colname]]) %>%
    mutate(CN_total = n()) %>%
    group_by(.data[[cn_colname]], .data[[temp_colname]]) %>%
    mutate(CN_count = n()) %>%
    ungroup() %>%
    mutate(Overall_pct = Overall_count / Overall_total * 100,
           CN_pct = CN_count / CN_total * 100) %>%
    select(.data[[cn_colname]], .data[[temp_colname]], CN_pct, Overall_pct) %>%
    distinct()
  
  neighbourhoods <- unique(hm_data[[cn_colname]])[order(as.numeric(gsub("\\D", "", unique(hm_data[[cn_colname]]))))]
  pct_CN <- matrix(NA, nrow = length(neighbourhoods), ncol = length(temp_cellTypes),
                   dimnames = list(neighbourhoods, temp_cellTypes))
  
  
  
  for (cn in neighbourhoods){
    
    cat("Calculating composition of", cn, "\n")
    
    cn_data <- hm_data %>%
      filter(.data[[cn_colname]] == cn)
    
    for (celltype in temp_cellTypes){
      
      if (celltype %in% unique(cn_data[[temp_colname]])){
        
        ct_data <- cn_data %>%
          filter(.data[[temp_colname]] == celltype)
        
        ratio <- ct_data$CN_pct / ct_data$Overall_pct
        if (ratio > 0) {
          pct_CN[cn, celltype] <- log2(ratio)
        } else {
          pct_CN[cn, celltype] <- 0
        }
        
      } else {
        next
      }
    }
  }
  
  pct_CN[is.na(pct_CN)] <- 0
  
  write.xlsx(pct_CN, paste0("GBM_spatial_analysis/outputs/LTS_analysis/neighbourhoods/", cn_colname, "_pct.xlsx"), rowNames = TRUE)
  
  palette <- colorRampPalette(
    c("blue",    
      "#E6E0EB",
      "firebrick3" 
    )
  )(100)
  
  ht <- Heatmap(
    pct_CN,
    name = "Enrichment score",
    col = palette,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_gp = gpar(fontsize = 10),
    row_names_max_width = unit(200, "mm"),
    column_names_gp = gpar(fontsize = 10),
    column_names_rot = 45,
    border = FALSE,
    show_heatmap_legend = FALSE,
    heatmap_legend_param = list(
      at = c(-6, -3, 0, +3, +6),
      labels = c("-6", "-3", "0", "+3", "+6"),
      title = "Enrichment score"
    ),
    column_title = paste0(tidy_m, "+ ", ct)
  )
  
  draw(ht)
  
  png(paste0("GBM_spatial_analysis/outputs/LTS_analysis/neighbourhoods/", cn_colname, "_heatmap.png"),
      width = 2300, height = 1650, res = 300) 
  draw(ht, padding = unit(c(20,35,5,5), "mm")) # add buffer
  dev.off()
  
  # Extract CN percentages, enrichment and significance
  
  plotdata <- keyTumour %>%
    filter(.data[[cn_colname]] != "CNNA") %>%
    group_by(Annotation_ID, Patient_ID, Cohort, .data[[cn_colname]]) %>%
    summarise(Count = n(), .groups = "drop") %>%
    group_by(Annotation_ID, Patient_ID, Cohort) %>%
    mutate(Total = sum(Count),
           pct = Count / Total * 100)
  
  cn_model_output <- list()
  for (cn in neighbourhoods){
    cn_data <- filter(plotdata, .data[[cn_colname]] == cn)
    
    # Skip if only one cohort is present
    if (n_distinct(cn_data$Cohort) < 2) {
      message("Skipping ", cn, " due to single cohort present")
      next
    }
    
    cn_model <- glmer(cbind(Count, Total - Count) ~ Cohort + (1 | Patient_ID),
                      family = binomial, data = cn_data)
    coefs <- summary(cn_model)$coefficients
    
    enriched_cells <- colnames(pct_CN)[pct_CN[cn, ] >= 1]
    enriched_cells_string <- paste(enriched_cells, collapse = "_")
    
    df <- as.data.frame(coefs["CohortSTS", c("Estimate", "Std. Error", "Pr(>|z|)")])
    cn_model_output[[cn]] <- t(df) %>%
      as.data.frame() %>%
      mutate(Enriched_cells = enriched_cells_string)
    rownames(cn_model_output[[cn]]) <- cn
  }
  
  cn_model_output <- do.call(rbind, cn_model_output)
  
  write.xlsx(cn_model_output, 
             paste0("GBM_spatial_analysis/outputs/LTS_analysis/neighbourhoods", cn_colname, "_models_enriched_cells.xlsx"), rownames = T)
  
  # Analyse interactions driving LTS enrichment of MHC- CNs ####
  
  pct_CN <- pct_CN[rownames(pct_CN) %in% rownames(cn_model_output), ]
  colnames(pct_CN) <- temp_cellTypes_tidy
  
  cn_regression <- data.frame(
    LTS_enriched = cn_model_output$Estimate,
    pct_CN
  )
  
  # look in CNs containing MHC-negative non-immune cells only
  mhc_neg <- cn_model_output %>%
    filter(grepl("MHC Class I-", Enriched_cells)) %>%
    rownames()
  
  cn_regression_mhc <- cn_regression[rownames(cn_regression) %in% mhc_neg, ]
  
  results_mhc <- data.frame(
    CellType = character(),
    Estimate = numeric(),
    StdError = numeric(),
    p_value = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (column in colnames(cn_regression_mhc)[colnames(cn_regression_mhc) != "LTS_enriched"]) {
    model <- lm(as.formula(paste("LTS_enriched ~", column)), data = cn_regression_mhc)
    coef_summary <- summary(model)$coefficients
    
    if (nrow(coef_summary) >= 2) {  # Ensure there is a coefficient for the predictor
      results_mhc <- rbind(
        results_mhc,
        data.frame(
          CellType = column,
          Estimate = coef_summary[2, "Estimate"],
          StdError = coef_summary[2, "Std. Error"],
          p_value = coef_summary[2, "Pr(>|t|)"]
        )
      )
    }
  }
  
  # Model effect of cell type combinations 
  
  # Two-part model (sub in as needed)
  
  var1 <- "MHC_neg_Non_immune_other"
  var2 <- paste(tidy_m, "pos", tidy_ct, sep = "_")
  var3 <- paste(tidy_m, "neg", tidy_ct, sep = "_")
  
  ct1 <- "MHC Class I- Non-immune (other)"
  ct2 <- temp_cellTypes[temp_cellTypes == paste0(m, "+ ", ct)]
  ct3 <- temp_cellTypes[temp_cellTypes == paste0(m, "- ", ct)]
  
  output_dir <- paste0("GBM_spatial_analysis/outputs/LTS_analysis/neighbourhoods/stats/mhc/", tidy_ct, "/", tidy_m, "_", tidy_ct, "/")
  dir.create(output_dir, recursive=T, showWarnings=F)
  
  # Wrap function
  
  run_model2 <- function(modx, colour_label) {
    fml <- as.formula(paste("LTS_enriched ~", var1, "*", modx))
    model <- lm(fml, data = cn_regression_mhc)
    coefs <- summary(model)$coefficients
    
    p_raw <- interact_plot(model, pred = !!sym(var1), modx = !!sym(modx),
                           plot.points = TRUE, interval = TRUE)
    
    # Extract plot ranges for annotation placement
    gb <- ggplot_build(p_raw)
    x_range <- gb$layout$panel_params[[1]]$x.range
    y_range <- gb$layout$panel_params[[1]]$y.range
    x_pos <- x_range[2] - 0.05 * diff(x_range)
    y_pos <- y_range[2] + 0.05 * diff(y_range)
    
    # Extract p-value from interaction term
    interaction_term <- grep(paste0(var1, ".*", modx), rownames(coefs), value = TRUE)
    pval <- signif(coefs[interaction_term, "Pr(>|t|)"], 2)
    
    p <- p_raw +
      labs(title = paste0("Effect of ", colour_label),
           y = "LTS ← → STS",
           x = ct1,
           color = colour_label) +
      theme(plot.title = element_text(size = 11, hjust = 0.5),
            legend.position = "none") +
      annotate("text", x = x_pos - 1, y = y_pos -0.5, 
               label = paste0("p = ", pval),
               size = 4)
    
    write.xlsx(coefs, paste0(output_dir, modx, "_", var1, ".xlsx"), rowNames = TRUE)
    ggsave(paste0(output_dir, modx, "_", var1, "_", ".png"), p, height = 4, width = 4)
  }
  
  run_model2(var2, ct2) # model vars 1/2
  run_model2(var3, ct3) # model vars 1/3
  
}

# Three-part model (sub in as needed)

var1 <- "MHC_neg_Non_immune_other"
var2 <- "MHC_neg_GBM_stem_cell"
var3 <- "MHC_pos_GBM_stem_cell"
var4 <- paste(tidy_m, "pos", tidy_ct, sep = "_")
var5 <- paste(tidy_m, "neg", tidy_ct, sep = "_")

ct1 <- "MHC Class I- Non-immune (other)"
ct2 <- "MHC Class I- GBM stem cell"
ct3 <- "MHC Class I+ GBM stem cell"
ct4 <- temp_cellTypes[temp_cellTypes == paste0(m, "+ ", ct)]
ct5 <- temp_cellTypes[temp_cellTypes == paste0(m, "- ", ct)]

output_dir <- paste0("GBM_spatial_analysis/outputs/LTS_analysis/neighbourhoods/stats/mhc/", tidy_ct, "/", tidy_m, "_", tidy_ct, "/")
dir.create(output_dir, recursive=T, showWarnings=F)

# Wrap function

run_model <- function(modx, mod2, label, colour_label) {
  fml <- as.formula(paste("LTS_enriched ~", var1, "*", modx, "*", mod2))
  model <- lm(fml, data = cn_regression_mhc)
  coefs <- summary(model)$coefficients
  
  p_raw <- interact_plot(model, pred = !!sym(var1), modx = !!sym(modx), mod2 = !!sym(mod2),
                         plot.points = TRUE, interval = TRUE)
  
  gb <- ggplot_build(p_raw)
  x_range <- gb$layout$panel_params[[1]]$x.range
  y_range <- gb$layout$panel_params[[1]]$y.range
  
  p <- p_raw +
    labs(title = paste0(colour_label, " effect on ", label, " interaction (p = ", 
                        signif(coefs[grep(paste0(var1, ".*", modx, ".*", mod2), rownames(coefs)), "Pr(>|t|)"], 2), ")"),
         y = "LTS ← → STS",
         x = ct1,
         color = colour_label) +
    theme(plot.title = element_text(size = 11, hjust = 0.5),
          legend.position = "none")
  
  write.xlsx(coefs, paste0(output_dir, mod2, "_", var1, "_", modx, ".xlsx"), rowNames = TRUE)
  ggsave(paste0(output_dir, mod2, "_", var1, "_", modx, ".png"), p, height = 3.7, width = 8)
}

run_model(var2, var4, "MHC Class I- GSC", ct4) # model vars 1/2/4
run_model(var2, var5, "MHC Class I- GSC", ct5) # model vars 1/2/5 etc.
run_model(var3, var4, "MHC Class I+ GSC", ct4)
run_model(var3, var5, "MHC Class I+ GSC", ct5)


# Plot marker effects on interactions ####

# MHC- GSCs

mhc_neg <- read.xlsx("GBM_spatial_analysis/outputs/LTS_analysis/neighbourhoods/stats/mhc/MHC_neg_GBM_stem_cell/MHC_neg_GBM_stem_cell_markers.xlsx")

colnames(mhc_neg) <- c("Marker",
                       "Estimate",
                       "StdError",
                       "t_value",
                       "p_value")

mhc_neg <- mhc_neg %>%
  select(-t_value) %>%
  mutate(Upper = Estimate + (1.96 * StdError),
         Lower = Estimate - (1.96 * StdError),
         CellType = str_replace(Marker, ".*:\\s*", ""),
         Marker = str_extract(CellType, "^[^_]+_[^_]+"),   
         Marker = str_replace(Marker, "_pos", "+"),
         Marker = str_replace(Marker, "_neg", "-"),
         Marker = str_replace(Marker, "PDL", "PD-L"),
         CellType = str_replace(CellType, "^[^_]+_[^_]+_", ""),
         CellType = str_replace(CellType, "MHC_neg", "MHC Class I-"),
         CellType = str_replace_all(CellType, "_", " "),
         Colour = case_when(
           CellType == "MHC Class I- GBM stem cell" ~ "grey50"
         )) %>%
  arrange(Estimate)

mhc_neg$Marker <- factor(mhc_neg$Marker,
                                   levels = mhc_neg$Marker)

mhc_neg_plot <- ggplot(mhc_neg, aes(x = Marker, y = Estimate, color = Colour)) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  geom_point(size = 5) +
  geom_linerange(aes(ymin = Lower, ymax = Upper),
                 size = 1) +
  geom_text(
    aes(
      label = ifelse(
        p_value <= 0.0001,
        paste(Marker, "(p < 0.0001)"),
        paste0(Marker, " (p = ", signif(p_value, 2), ")")
      ),
      y = max(Upper + 0.05)  # adjust as needed for padding
    ),
    hjust = 0,
    color = "black",
    size = 5
  ) +
  coord_flip() +
  labs(y = "Effect size", x = NULL) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15),
    axis.title.x = element_text(size = 19, hjust = 0.2),
    plot.background = element_rect(fill = "white"),
    legend.position = "none") +
  scale_y_continuous(
    limits = c(min(mhc_neg$Lower) - 0.3, max(mhc_neg$Upper) + 0.35),
    breaks = c(-0.2, 0, 0.2, 0.4)
  ) +
  scale_color_manual(values = setNames(as.character(mhc_neg$Colour), mhc_neg$Colour))

mhc_neg_plot

ggsave("GBM_spatial_analysis/outputs/LTS_analysis/neighbourhoods/stats/mhc/MHC_neg_GBM_stem_cell/MHC_neg_GSC_markers.png",
       mhc_neg_plot,
       height = 3, width = 7)

# MHC+ GSCs

mhc_pos <- read.xlsx("GBM_spatial_analysis/outputs/LTS_analysis/neighbourhoods/stats/mhc/MHC_pos_GBM_stem_cell/MHC_pos_GBM_stem_cell_markers.xlsx")

colnames(mhc_pos) <- c("Marker",
                       "Estimate",
                       "StdError",
                       "t_value",
                       "p_value")

mhc_pos <- mhc_pos %>%
  select(-t_value) %>%
  mutate(Upper = Estimate + (1.96 * StdError),
         Lower = Estimate - (1.96 * StdError),
         CellType = str_replace(Marker, ".*:\\s*", ""),
         Marker = str_extract(CellType, "^[^_]+_[^_]+"),   
         Marker = str_replace(Marker, "_pos", "+"),
         Marker = str_replace(Marker, "_neg", "-"),
         Marker = str_replace(Marker, "PDL", "PD-L"),
         CellType = str_replace(CellType, "^[^_]+_[^_]+_", ""),
         CellType = str_replace(CellType, "MHC_pos", "MHC Class I+"),
         CellType = str_replace_all(CellType, "_", " "),
         Colour = case_when(
           CellType == "MHC Class I+ GBM stem cell" ~ "grey70"
         )) %>%
  arrange(Estimate)

mhc_pos$Marker <- factor(mhc_pos$Marker,
                         levels = mhc_pos$Marker)

mhc_pos_plot <- ggplot(mhc_pos, aes(x = Marker, y = Estimate, color = Colour)) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  geom_point(size = 5) +
  geom_linerange(aes(ymin = Lower, ymax = Upper),
                 size = 1) +
  geom_text(
    aes(
      label = ifelse(
        p_value <= 0.0001,
        paste(Marker, "(p < 0.0001)"),
        paste0(Marker, " (p = ", signif(p_value, 2), ")")
      ),
      y = max(Upper + 0.1)  # adjust as needed for padding
    ),
    hjust = 0,
    color = "black",
    size = 5
  ) +
  coord_flip() +
  labs(y = "Effect size", x = NULL) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15),
    axis.title.x = element_text(size = 19, hjust = 0.32),
    plot.background = element_rect(fill = "white"),
    legend.position = "none") +
  scale_y_continuous(
    limits = c(min(mhc_pos$Lower) - 0.3, max(mhc_pos$Upper) + 0.35),
    breaks = c(-0.2, 0, 0.2)
  ) +
  scale_color_manual(values = setNames(as.character(mhc_pos$Colour), mhc_pos$Colour))

mhc_pos_plot

ggsave("GBM_spatial_analysis/outputs/LTS_analysis/neighbourhoods/stats/mhc/MHC_pos_GBM_stem_cell/MHC_pos_GSC_markers.png",
       mhc_pos_plot,
       height = 3, width = 7)

# Immune on MHC- GSC

immune_mhc_neg <- read.xlsx("GBM_spatial_analysis/outputs/LTS_analysis/neighbourhoods/stats/mhc/MHC_neg_GBM_stem_cell/Immune_MHC_neg_GBM_stem_cell_markers.xlsx")

colnames(immune_mhc_neg) <- c("Marker",
                       "Estimate",
                       "StdError",
                       "t_value",
                       "p_value")

immune_mhc_neg <- immune_mhc_neg %>%
  select(-t_value) %>%
  mutate(Upper = Estimate + (1.96 * StdError),
         Lower = Estimate - (1.96 * StdError),
         CellType = str_replace(Marker,"^[^:]*:[^:]*:\\s*", ""),
         Marker = str_extract(CellType, "^[^_]+_[^_]+"),   
         Marker = str_replace(Marker, "_pos", "+"),
         Marker = str_replace(Marker, "_neg", "-"),
         Marker = str_replace(Marker, "PDL", "PD-L"),
         Marker = str_replace(Marker, "CTLA4", "CTLA-4"),
         Marker = str_replace(Marker, "LAG3", "LAG-3"),
         Marker = str_replace(Marker, "GZM", "Granzyme "),
         CellType = str_replace(CellType, "^[^_]+_[^_]+_", ""),
         CellType = str_replace(CellType, "MHC", "MHC Class I"),
         CellType = str_replace(CellType, "_neg", "-"),
         CellType = str_replace(CellType, "_pos", "+"),
         CellType = str_replace(CellType, "CD8", "CD8+"),
         CellType = str_replace_all(CellType, "_", " "),
         Colour = case_when(
           CellType == "MHC Class I- GBM stem cell" ~ "grey50",
           CellType == "MHC Class I+ GBM stem cell" ~ "grey70",
           CellType == "CD8+ T cell" ~ "#FFF200",
           CellType == "Th cell" ~ "#D35400",
           CellType == "Macrophage" ~ "#E60026",
           CellType == "Microglia" ~ "#FF6F61"
         ),
         CellType = factor(CellType, levels = c("CD8+ T cell",
                                                "Th cell",
                                                "Macrophage",
                                                "Microglia"))) %>%
  arrange(desc(CellType), desc(Estimate)) %>%
  mutate(MarkerLabel = Marker,
         Marker = paste0(Marker, "_", row_number()),
         Marker = factor(Marker, levels = Marker))

immune_mhc_neg$CellType <- factor(immune_mhc_neg$CellType)
immune_mhc_neg$Marker <- factor(immune_mhc_neg$Marker)
immune_mhc_neg$Colour <- factor(immune_mhc_neg$Colour)

immune_mhc_neg_plot <- ggplot(immune_mhc_neg, aes(x = Marker, y = Estimate, color = Colour)) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  geom_point(size = 5) +
  geom_linerange(aes(ymin = Lower, ymax = Upper),
                 size = 1) +
  geom_text(
    aes(
      label = ifelse(
        p_value <= 0.0001,
        paste(MarkerLabel, "(p < 0.0001)"),
        paste0(MarkerLabel, " (p = ", signif(p_value, 2), ")")
      ),
      y = max(Upper + 0.05)  # adjust as needed for padding
    ),
    hjust = 0,
    color = "black",
    size = 5
  ) +
  coord_flip() +
  labs(y = "Effect size", x = NULL) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15),
    axis.title.x = element_text(size = 19, hjust = 0.63),
    plot.background = element_rect(fill = "grey90"),
    legend.position = "none") +
  scale_y_continuous(
    limits = c(min(immune_mhc_neg$Lower), max(immune_mhc_neg$Upper) + 0.3),
    breaks = c(-0.4, -0.2, 0, 0.2)
  ) +
  scale_color_manual(values = setNames(as.character(immune_mhc_neg$Colour), immune_mhc_neg$Colour))

immune_mhc_neg_plot

ggsave("GBM_spatial_analysis/outputs/LTS_analysis/neighbourhoods/stats/mhc/MHC_neg_GBM_stem_cell/Immune_MHC_neg_GSC_markers.png",
       immune_mhc_neg_plot,
       height = 9, width = 6.6)

# Immune on MHC+ GSC

immune_mhc_pos <- read.xlsx("GBM_spatial_analysis/outputs/LTS_analysis/neighbourhoods/stats/mhc/MHC_pos_GBM_stem_cell/Immune_MHC_pos_GBM_stem_cell_markers.xlsx")

colnames(immune_mhc_pos) <- c("Marker",
                              "Estimate",
                              "StdError",
                              "t_value",
                              "p_value")

immune_mhc_pos <- immune_mhc_pos %>%
  select(-t_value) %>%
  mutate(Upper = Estimate + (1.96 * StdError),
         Lower = Estimate - (1.96 * StdError),
         Marker = str_replace(Marker, "MHC_II", "MHCII"),
         CellType = str_replace(Marker,"^[^:]*:[^:]*:\\s*", ""),
         Marker = str_extract(CellType, "^[^_]+_[^_]+"),   
         Marker = str_replace(Marker, "_pos", "+"),
         Marker = str_replace(Marker, "_neg", "-"),
         Marker = str_replace(Marker, "PDL", "PD-L"),
         Marker = str_replace(Marker, "PD1", "PD-1"),
         Marker = str_replace(Marker, "MHCII", "MHC Class II"),
         Marker = str_replace(Marker, "CTLA4", "CTLA-4"),
         Marker = str_replace(Marker, "LAG3", "LAG-3"),
         Marker = str_replace(Marker, "GZM", "Granzyme "),
         CellType = str_replace(CellType, "^[^_]+_[^_]+_", ""),
         CellType = str_replace(CellType, "MHC_neg", "MHC Class I-"),
         CellType = str_replace(CellType, "MHC_pos", "MHC Class I+"),
         CellType = str_replace(CellType, "_neg", "-"),
         CellType = str_replace(CellType, "_pos", "+"),
         CellType = str_replace(CellType, "CD8", "CD8+"),
         CellType = str_replace_all(CellType, "_", " "),
         Colour = case_when(
           CellType == "MHC Class I- GBM stem cell" ~ "grey50",
           CellType == "MHC Class I+ GBM stem cell" ~ "grey70",
           CellType == "CD8+ T cell" ~ "#FFF200",
           CellType == "Th cell" ~ "#D35400",
           CellType == "Macrophage" ~ "#E60026",
           CellType == "Microglia" ~ "#FF6F61"
         ),
         CellType = factor(CellType, levels = c("CD8+ T cell",
                                                "Th cell",
                                                "Macrophage",
                                                "Microglia"))) %>%
  arrange(desc(CellType), desc(Estimate)) %>%
  mutate(MarkerLabel = Marker,
         Marker = paste0(Marker, "_", row_number()),
         Marker = factor(Marker, levels = Marker))

immune_mhc_pos$CellType <- factor(immune_mhc_pos$CellType)
immune_mhc_pos$Marker <- factor(immune_mhc_pos$Marker)
immune_mhc_pos$Colour <- factor(immune_mhc_pos$Colour)

immune_mhc_pos_plot <- ggplot(immune_mhc_pos, aes(x = Marker, y = Estimate, color = Colour)) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  geom_point(size = 5) +
  geom_linerange(aes(ymin = Lower, ymax = Upper),
                 size = 1) +
  geom_text(
    aes(
      label = ifelse(
        p_value <= 0.0001,
        paste(MarkerLabel, "(p < 0.0001)"),
        paste0(MarkerLabel, " (p = ", signif(p_value, 2), ")")
      ),
      y = max(Upper + 0.05)  # adjust as needed for padding
    ),
    hjust = 0,
    color = "black",
    size = 5
  ) +
  coord_flip() +
  labs(y = "Effect size", x = NULL) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15),
    axis.title.x = element_text(size = 19, hjust = 0.45),
    plot.background = element_rect(fill = "grey90"),
    legend.position = "none") +
  scale_y_continuous(
    limits = c(min(immune_mhc_pos$Lower) - 0.02, max(immune_mhc_pos$Upper) + 0.32),
    breaks = c(-0.6, -0.4, -0.2, 0, 0.2, 0.4)
  ) +
  scale_color_manual(values = setNames(as.character(immune_mhc_pos$Colour), immune_mhc_pos$Colour))

immune_mhc_pos_plot

ggsave("GBM_spatial_analysis/outputs/LTS_analysis/neighbourhoods/stats/mhc/MHC_pos_GBM_stem_cell/Immune_MHC_pos_GSC_markers.png",
       immune_mhc_pos_plot,
       height = 6, width = 10)

# Legends

res1 <- mhc_pos %>%
  select(CellType, Colour)
res2 <- immune_mhc_pos %>%
  select(CellType, Colour)
legend_results <- rbind(res1, res2) %>%
  distinct() %>%
  mutate(CellType = factor(CellType, levels = c("MHC Class I+ GBM stem cell",
                                                "CD8+ T cell",
                                                "Th cell",
                                                "Macrophage",
                                                "Microglia"))) %>%
  arrange(CellType)
View(legend_results)

legend_results$CellType <- as.character(legend_results$CellType)
legend_results$Colour <- as.character(legend_results$Colour)

# Create named vector of colours: names must match CellType values
colour_map <- setNames(legend_results$Colour, legend_results$CellType)

legend_results$CellType <- factor(legend_results$CellType, levels = legend_results$CellType)
legend_results$Colour <- factor(legend_results$Colour, levels = legend_results$Colour)

# Dummy plot for generating the legend
legend_plot <- ggplot(legend_results, aes(x = CellType, y = CellType, fill = CellType)) +
  geom_bar(stat = "identity", show.legend = TRUE) +
  scale_fill_manual(values = colour_map) +
  guides(fill = guide_legend(title = "Cell Type")) +
  theme_void() +
  theme(legend.position = "right")
legend_plot

# Extract and display just the legend
legend_only <- get_legend(legend_plot)
plot_grid(legend_only)

ggsave("GBM_spatial_analysis/outputs/LTS_analysis/neighbourhoods/stats/mhc/MHC_pos_GBM_stem_cell/legend.png",
       legend_only,
       height = 4, width = 3)

# Example visualisation of CNs ####
marker_start <- which(colnames(annData) == "Cell_area") + 1
marker_end <- which(colnames(annData) == "CD39") - 1

View(cn4_props)
# 54 PT 1 shown to CHO
# 122 intravascular CN4 D12
# 004 F5 some haemorrhage and a vessel with perivascular macs
# 104 E11 haemorrhage but MHC+ non-immune


CN_colours <- c(
  "#e6ab02",
  "#1b9e77", 
  "#7570b3", 
  "#e7298a",
  "#66a61e",
  "#d95f02",
  "#1f78b4",
  "#a6761d",
  "#666666"
)

temp_celltype_colours <- c(
  # GBM stem cell
  "MHC Class I+ GBM stem cell" = "grey90", 
  "MHC Class I- GBM stem cell" = "grey60",   
  
  # AstrocyteGBM cell
  "MHC Class I+ AstrocyteGBM cell" = "#F0E1A0", 
  "MHC Class I- AstrocyteGBM cell" = "#C6A84D",  
  
  # Fibroblast
  "MHC Class I+ Fibroblast" = "#C7F1C4",     
  "MHC Class I- Fibroblast" = "#32CD32",     
  
  # Non-immune (other)
  "MHC Class I+ Non-immune (other)" = "grey30",
  "MHC Class I- Non-immune (other)" = "black", 
  
  # Others
  "Neuron" = "white",
  "Endothelial cell" = "#00FA9A",
  "Vascular smooth muscle cell" = "white",
  "Macrophage" = "#E60026",
  "Dendritic cell" = "white",
  "Microglia" = "#FF6F61",
  "Neutrophil" = "#D100D1",
  "NK cell" = "white",
  "B cell" = "white",
  "Plasma cell" = "white",
  "CD8+ T cell" = "#FFF200",
  "Th cell" = "#D35400",
  "Treg" = "white",
  "T cell (other)" = "white",
  "Immune (other)" = "white"
)

annots <- cn4_props %>%
  filter(Total > 3000) %>%
  group_by(Cohort) %>%
  slice_max(order_by = pct, n = 15) %>%
  ungroup() %>%
  pull(Annotation_ID)

annots <- c("132_Primary_Tumour_2",
            "122_Primary_Tumour_1",
            "004_Primary_Tumour_1",
            "104_Primary_Tumour_1",
            "054_Primary_Tumour_2")

for (a in annots){
  
  cat("Processing annotation", which(annots == a), "/", length(annots), "\n")
  
  annData <- keyTumour %>%
    filter(Annotation_ID == a)
  
  rownames(annData) <- annData$Cell_ID
  annData$Cell_ID <- NULL
  im <- t(annData[, marker_start:marker_end])  # intensity matrix
  spe <- format_image_to_spe(format = "general", 
                             intensity_matrix = im,
                             coord_x = annData$X_coord,
                             coord_y = annData$Y_coord,
                             phenotype = annData$temp_celltype)
  
  meta_data <- colnames(annData[, 1:which(colnames(annData) == "Y_coord")])
  remove <- c("Cell.ID", "Phenotype", "sample_id")  
  colData(spe) <- cbind(colData(spe), annData[, meta_data])
  colData(spe)[, remove] <- NULL
  rownames(spatialCoords(spe)) <- colnames(spe)
  
  # extract coordinates
  coords <- spatialCoords(spe) 
  cells <- rownames(coords)
  
  names(CN_colours) <- paste0("CN", 1:9)
  cns_in_core <- CN_colours[unique(annData$CN_k10_K9)]
  
  cn_plot <- plot_cell_categories(spe_object = spe, 
                                  feature_colname = "CN_k10_K9",
                                  categories_of_interest = sort(unique(annData$CN_k10_K9)),
                                  colour_vector = cns_in_core) +
    theme(plot.title = element_blank(),
          legend.position = "none")
  
  names(my_cols) <- levels(cellTypes)
  cols_for_core <- my_cols[unique(annData$Cell.Type)]
  
  temp_cols_for_core <- temp_celltype_colours[unique(annData$temp_celltype)]
  
  ct_plot <- plot_cell_categories(spe_object = spe, 
                                  feature_colname = "temp_celltype",
                                  categories_of_interest = unique(annData$temp_celltype),
                                  colour_vector = temp_cols_for_core) +
    theme(plot.title = element_blank(),
          legend.position = "none")
  
  output_dir <- paste0("GBM_spatial_analysis/outputs/LTS_analysis/neighbourhoods/pics/", a)
  dir.create(output_dir, recursive=T, showWarnings=F)
  
  ggsave(paste0(output_dir, "/cn_vis_", a, ".png"),
         cn_plot,
         height = 3, width = 3)
  
  ggsave(paste0(output_dir, "/ct_vis_", a, ".png"),
         ct_plot,
         height = 4.5, width = 4.5)
}

# Check for nuclei in MHC- non-immune ####

intensityMatrix <- read.csv("GBM_spatial_analysis/data/intensityMatrix_files_thesis/LTS_nucleus_intensities.csv")

colnames(intensityMatrix) <- c("Cell_ID",
                               "Mean_nuclear_DAPI",
                               "Median_nuclear_DAPI",
                               "Mean_cell_DAPI",
                               "Median_cell_DAPI")

keyTumour <- keyTumour %>%
  left_join(intensityMatrix, by = "Cell_ID")

nucleus_data <- keyTumour %>%
  filter(Cell.Type %in% c("Non-immune (other)", "GBM stem cell")) %>%
  select(Annotation_ID, Cell.Type, temp_celltype, Median_nuclear_DAPI.x, Nucleus_present, CN_k10_K9, Cell_area) %>%
  mutate(temp_CN = case_when(
    CN_k10_K9 == "CN4" ~ "CN4",
    CN_k10_K9 != "CN4" ~ "Other",
    TRUE ~ CN_k10_K9
  )) %>%
  group_by(Annotation_ID, 
           #temp_CN, 
           temp_celltype) %>%
  summarise(
    DAPI_level = mean(Median_nuclear_DAPI.x, na.rm = TRUE),
    CellArea = mean(Cell_area),
    N_total = n(),
    N_with_nucleus = sum(Nucleus_present == 1, na.rm = TRUE),
    Pct_with_nucleus = N_with_nucleus / N_total * 100,
    .groups = "drop"
  )

summ <- nucleus_data %>%
  group_by(temp_celltype) %>%
  summarise(
    pct_nuc_med = median(Pct_with_nucleus, na.rm = T),
    pct_nuc_upper = quantile(Pct_with_nucleus, 0.75, na.rm = T),
    pct_nuc_lower = quantile(Pct_with_nucleus, 0.25, na.rm = T),
    dapi_med = median(DAPI_level, na.rm = T),
    dapi_upper = quantile(DAPI_level, 0.75, na.rm = T),
    dapi_lower = quantile(DAPI_level, 0.25, na.rm = T),
    area_med = median(CellArea, na.rm = T),
    area_upper = quantile(CellArea, 0.75, na.rm = T),
    area_lower = quantile(CellArea, 0.25, na.rm = T),
  )

mhc_nucleus <- nucleus_data %>%
  filter(temp_celltype %in% c("MHC Class I+ GBM stem cell", "MHC Class I- Non-immune (other)")) %>%
  mutate(Donor = str_sub(Annotation_ID, 1, 3))

model <- glmer(cbind(N_with_nucleus, N_total - N_with_nucleus) ~ temp_celltype + (1 | Donor), family = binomial, data = mhc_nucleus)
summary(model)

model <- lmer(CellArea ~ temp_celltype + (1 | Donor), data = mhc_nucleus)
summary(model)

pval <- data.frame(
  group1 = "MHC Class I+ GBM stem cell",
  group2 = "MHC Class I- Non-immune (other)",
  y.position = 4900,
  p.signif = "****" 
)

p <- ggplot(mhc_nucleus, aes(y = DAPI_level, x = temp_celltype, fill = temp_celltype)) +
  geom_boxplot(outlier.color = "black", outlier.shape = 21, outlier.size = 1.3, color = "black") +
  labs(y = expression("Median nuclear DAPI intensity"),
       fill = "Cell type") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white"),
    panel.grid.major = element_line(color = "gray", size = 0.2),
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "black", size = 11),
    axis.text.x = element_blank(),
    axis.title = element_text(color = "black", size = 12),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    legend.position = "none",
    strip.text = element_text(size = 9)
  ) +
  scale_fill_manual(values = c("grey90", "grey10")) +
  stat_pvalue_manual(data = pval, 
                   label = "p.signif",
                   tip.length = 0.01,,
                   size = 5,
                   inherit.aes = F)


p
ggsave("GBM_spatial_analysis/outputs/LTS_analysis/neighbourhoods/pics/dapi_level.png",
       p,
       height = 4, width = 4)

# K-M analysis ####

lts_cohort <- lts_cohort %>%
  mutate(
    EOR_binary = case_when(
      EOR == "Complete" ~ "Complete",
      EOR == "No_imaging" ~ NA_character_,
      TRUE ~ "Subtotal"
    )
  )

lts_cohort <- lts_cohort %>%
  mutate(
    temp = if_else(
      Op_type1 == "R",
      "Complete",
      "Subtotal"),
    EOR_binary_2 = case_when(
      EOR_binary == "Complete" ~ "Complete",
      EOR_binary == "Subtotal" ~ "Subtotal",
      TRUE ~ temp)
  ) %>%
  select(-temp)

cn4_by_donor <- cn4_props %>%
  group_by(Patient_ID, Cohort) %>%
  summarise(
    SumCount = sum(Count),
    SumTotal = sum(Total),
    .groups = "drop"
  ) %>%
  mutate(SumPct = SumCount / SumTotal * 100,
         Over_1 = case_when(
           SumPct > 1 ~ 1,
           SumPct <= 1 ~ 0
         ),
         Over_3 = case_when(
           SumPct > 3 ~ 1,
           SumPct <= 3 ~ 0
         ),
         Over_5 = case_when(
           SumPct > 5 ~ 1,
           SumPct <= 5 ~ 0
         ),
         Over_10 = case_when(
           SumPct > 10 ~ 1,
           SumPct <= 10 ~ 0
         )
  )

survival <- lts_cohort %>%
  filter(Patient_ID %in% unique(cn4_by_donor$Patient_ID)) %>%
  mutate(OS_event = case_when(
    OS_censored == 0 ~ 1,
    OS_censored == 1 ~ 0),
    PFS_event = case_when(
      PFS_censored == 0 ~ 1,
      PFS_censored == 1 ~ 0)) %>%
  select(Patient_ID, Age_at_diagnosis, Gad_cortex, Gad_midline, Gad_vent,
         EOR_binary, EOR_binary_2, MGMT_final, PFS_months, PFS_event,
         OS_months, OS_event)

cn4_KM_data <- cn4_by_donor %>%
  left_join(survival, by = "Patient_ID")

# Fit survival object for baseline table

os_surv <- Surv(time = survival$OS_months, event = survival$OS_event)
pfs_surv <- Surv(time = survival$PFS_months, event = survival$PFS_event)
group_os <- survfit(os_surv ~ survival$Cohort)
group_pfs <- survfit(pfs_surv ~ survival$Cohort)
os_diff <- survdiff(os_surv ~ Cohort, data = survival)
pfs_diff <- survdiff(pfs_surv ~ Cohort, data = survival)
os_diff
pfs_diff

# Fit survival objects
os_surv_object <- Surv(time = cn4_KM_data$OS_months, event = cn4_KM_data$OS_event)
pfs_surv_object <- Surv(time = cn4_KM_data$PFS_months, event = cn4_KM_data$PFS_event)

# Fit Kaplan-Meier
os_fit <- survfit(os_surv_object ~ cn4_KM_data$Over_3)
pfs_fit <- survfit(pfs_surv_object ~ cn4_KM_data$Over_3)

# Plot
p <- ggsurvplot(os_fit, 
           data = cn4_KM_data,
           legend.labs = c("Below median (3%)", "Above median (3%)"),
           legend = "right",
           palette = c("grey70", "#e7298a"),
           ggtheme = theme_classic() +
             theme(legend.title = element_text(size = 7),
                   legend.text = element_text(size = 6),
                   axis.title = element_text(size = 13),
                   axis.text = element_text(size = 12))) +
  labs(x = "Time (months)",
       y = "Overall survival probability",
       color = "CN4 abundance")
p  
ggsave("GBM_spatial_analysis/outputs/LTS_analysis/neighbourhoods/stats/km_os3.png",
       p$plot,
       height = 4, width = 4)

cox_model <- coxph(Surv(OS_months, OS_event) ~ Over_1 + Age_at_diagnosis + MGMT_final,
                   data = cn4_KM_data)
summary(cox_model)

# Test other covariates ####

wilcox.test(PFS_days ~ Cohort, data = lts_cohort)

lts_cohort$ce_touching_ventricular_surface[lts_cohort$ce_touching_ventricular_surface == "y"] <- "Y"
lts_cohort$ce_touching_ventricular_surface[lts_cohort$ce_touching_ventricular_surface == "n"] <- "N"
ce_vent <- table(lts_cohort$Cohort, lts_cohort$ce_touching_ventricular_surface)
fisher.test(ce_vent)

# (Find marker expression on protective immune cells) ####

View(cn_regression_mhc)
CN_list <- paste0("CN",
                  c(12, 15, 19, 22, 25, 28, 30, 31, 35, 36, 38, 39, 40, 41, 42, 43, 46, 48))

target_CNs <- cn_regression_mhc %>%
  filter(MHC_neg_Non_immune_other > 1,
         MHC_neg_GBM_stem_cell > 1)

View(target_CNs)

macs <- keyTumour %>%
  filter(Cell.Type == "Macrophage",
         CN_k10_K50 %in% c("CN35", "CN36"))

for (m in markerList) {
  
  bpData <- macs %>%
    mutate(temp_celltype = case_when( 
      grepl(paste0(m, "($|:)"), Classifier) ~ paste0(m, "+ ", Cell.Type),
      !grepl(paste0(m, "($|:)"), Classifier) ~ paste0(m, "- ", Cell.Type),
      TRUE ~ Cell.Type  # leave other cell types unchanged
    )) %>%
    group_by(Annotation_ID, CN_k10_K50, Cell.Type) %>%
    mutate(Total = n()) %>%
    group_by(Annotation_ID, CN_k10_K50, temp_celltype) %>%
    reframe(Count = n(),
            Total = Total,
            CN_k10_K50 = CN_k10_K50,
            temp_celltype = temp_celltype) %>%
    mutate(pct = Count / Total * 100) %>%
    distinct() %>%
    filter(!grepl(paste0(m, "-"), temp_celltype),
           Total > 10)
  
  if (length(bpData$temp_celltype) == 0) {
    next } else {
      
      pvalues <- data.frame(
        group1 = "CN4",
        group2 = "Other",
        y.position = 73,
        p = "****",
        temp_celltype = c(
          paste0(m, "+ ", cancer_cells[1]),
          paste0(m, "+ ", cancer_cells[2])
        )
      )
      bp <- ggplot(data = bpData,
                   mapping = aes(x = CN_k10_K50, y = pct, fill = CN_k10_K50)) + 
        geom_boxplot(outlier.color = "black", outlier.shape = 21, outlier.size = 1.3, color = "black") +
        labs(y = "Percentage of cells expressing marker",
             title = paste0(m, "+ macrophages")) +
        theme_classic() +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold"),
          panel.background = element_rect(fill = "white", color = "white"),
          plot.background = element_rect(fill = "white", color = "white"),
          panel.grid.major = element_line(color = "gray", size = 0.2),
          panel.grid.minor = element_blank(),
          axis.text = element_text(color = "black", size = 11),
          axis.text.x = element_blank(),
          axis.title = element_text(color = "black", size = 12),
          axis.title.x = element_blank(),
          legend.position = "right",
          strip.text = element_text(size = 9)
        )
    bp  
    }
  
  ggsave(paste0("GBM_spatial_analysis/outputs/LTS_analysis/neighbourhoods/markers/CN_macs/", m, ".png"),
         bp,
         height = 4, width = 5.7)
}

# (Explore profile of interacting Th cells in LTS) ####

th_profiles <- c("CD38", "Granzyme A", "Granzyme B", "CD44", "CXADR", "PD-1", "CTLA-4", "MHC Class II", "CXCL13", "Ki67")
targets <- c("Endothelial cell", "Neuron", "AstrocyteGBM cell", "GBM stem cell", "CD8+ T cell", "Macrophage", "Microglia")

th_results <- data.frame()
for(th in th_profiles){
  
  cat("Assessing distances between", th, "Th cells and target cells \n")
  
  input <- keySpatial %>%
    filter(Cell.Type %in% c("Th cell", targets)) %>%
    mutate(temp_celltype = case_when( 
      Cell.Type == "Th cell" & grepl(paste0(th, "($|:)"), Classifier) ~ paste0(th, "+ Th cell"),
      Cell.Type == "Th cell" & !grepl(paste0(th, "($|:)"), Classifier) ~ paste0(th, "- Th cell"),
      TRUE ~ Cell.Type  # leave other cell types unchanged
    )) %>%
    select(temp_celltype, everything())
  
  output <- data.frame()
  annots <- unique(input$Annotation_ID)
  
  for (a in annots){
    
    cat("Processing annotation:", which(annots == a), "of", length(annots), "\n")
    
    annData <- filter(input,
                      Annotation_ID == a)
    rownames(annData) <- annData$Cell_ID
    annData$Cell_ID <- NULL
    im <- t(annData[, 15:ncol(annData)]) # may vary
    spe <- format_image_to_spe(format = "general", 
                               intensity_matrix = im,
                               coord_x = annData$X_coord,
                               coord_y = annData$Y_coord,
                               phenotype = NA)
    
    meta_data <- colnames(annData[, 1:11])
    remove <- c("Cell.ID", "Phenotype", "sample_id")  
    colData(spe) <- cbind(colData(spe), annData[, meta_data])
    colData(spe)[, remove] <- NULL
    rownames(spatialCoords(spe)) <- colnames(spe)
    
    #Pairwise distances:
    distances <- calculate_minimum_distances_between_celltypes(
      spe_object = spe,
      feature_colname = "temp_celltype",
      cell_types_of_interest = unique(input$temp_celltype)) %>%
      mutate(
        Annotation_ID = a,
        Patient_ID = unique(annData$Patient_ID),
        Annotation_Type = unique(annData$Annotation_Type),
        Cohort = unique(annData$Cohort)
      )
    
    if (nrow(distances) == 0) {
      cat("Skipping annotation:", a, "- no valid distance results.\n")
      next
    }
    
    summary_distances <- distances %>%
      calculate_summary_distances_between_celltypes() %>%
      mutate(
        Annotation_ID = a,
        Patient_ID = unique(annData$Patient_ID),
        Annotation_Type = unique(annData$Annotation_Type),
        Cohort = unique(annData$Cohort)
      )
    
    # Remove rows with missing Pair to avoid downstream errors
    cat("Unique Pair values in annotation", a, ":\n")
    print(unique(summary_distances$Pair))
    
    summary_distances <- summary_distances %>%
      filter(!is.na(Pair))
    
    if (nrow(summary_distances) == 0) {
      cat("Skipping annotation:", a, "- no valid cell pairs found after removing NAs.\n")
      next
    }
    
    summary_distances <- summary_distances %>%
      rowwise() %>%
      mutate(
        num_ref_cells = length(unique(distances$RefCell[distances$RefType == Reference]))
      ) %>%
      ungroup()
    
    cells_of_interest <- c(paste0(th, "+ Th cell"), paste0(th, "- Th cell"))
    output <- rbind(output, summary_distances) %>%
      filter(Target %in% cells_of_interest)
  }
  
  th_results <- rbind(th_results, output)
}

th_lts <- th_results %>%
  filter(num_ref_cells >= 5,
         Cohort == "LTS",
         !grepl("Th cell", Reference)) %>%
  group_by(Pair) %>%
  mutate(num_cores = n()) %>%
  ungroup()

th_lts_results <- th_lts %>%
  group_by(Pair) %>%
  summarise(Distance = mean(Mean),
            Pair = unique(Pair),
            Reference = unique(Reference),
            Target = unique(Target))

hm_data <- th_lts_results %>%
  select(Target, Reference, Distance) %>%
  pivot_wider(names_from = Reference, values_from = Distance) %>%
  column_to_rownames("Target")


# Within radius of Th

th_results <- data.frame()

output <- data.frame()
annots <- unique(input$Annotation_ID)

for (a in annots){
  
  cat("Processing annotation:", which(annots == a), "of", length(annots), "\n")
  
  annData <- filter(keySpatial,
                    Annotation_ID == a)
  rownames(annData) <- annData$Cell_ID
  annData$Cell_ID <- NULL
  im <- t(annData[, 13:ncol(annData)]) # may vary
  spe <- format_image_to_spe(format = "general", 
                             intensity_matrix = im,
                             coord_x = annData$X_coord,
                             coord_y = annData$Y_coord,
                             phenotype = NA)
  
  meta_data <- colnames(annData[, 1:11])
  remove <- c("Cell.ID", "Phenotype", "sample_id")  
  colData(spe) <- cbind(colData(spe), annData[, meta_data])
  colData(spe)[, remove] <- NULL
  rownames(spatialCoords(spe)) <- colnames(spe)
  
  for (t in targets){
    rad <- average_percentage_of_cells_within_radius(spe_object = spe, 
                                                     reference_celltype = t, 
                                                     target_celltype = "Th cell", 
                                                     radius=100,
                                                     feature_colname="Cell.Type")
  }
  
  rad  
  
  rad <- average_percentage_of_cells_within_radius(spe_object = spe, 
                                                   reference_celltype = "Th cell", 
                                                   target_celltype = c("Astrocyte/GBM cell"), 
                                                   radius=100,
                                                   feature_colname="Cell.Type") 
  
  
  
  
  distances <- calculate_minimum_distances_between_celltypes(
    spe_object = spe,
    feature_colname = "temp_celltype",
    cell_types_of_interest = unique(input$temp_celltype)) %>%
    mutate(
      Annotation_ID = a,
      Patient_ID = unique(annData$Patient_ID),
      Annotation_Type = unique(annData$Annotation_Type),
      Cohort = unique(annData$Cohort)
    )
  
  if (nrow(distances) == 0) {
    cat("Skipping annotation:", a, "- no valid distance results.\n")
    next
  }
  
  summary_distances <- distances %>%
    calculate_summary_distances_between_celltypes() %>%
    mutate(
      Annotation_ID = a,
      Patient_ID = unique(annData$Patient_ID),
      Annotation_Type = unique(annData$Annotation_Type),
      Cohort = unique(annData$Cohort)
    )
  
  # Remove rows with missing Pair to avoid downstream errors
  cat("Unique Pair values in annotation", a, ":\n")
  print(unique(summary_distances$Pair))
  
  summary_distances <- summary_distances %>%
    filter(!is.na(Pair))
  
  if (nrow(summary_distances) == 0) {
    cat("Skipping annotation:", a, "- no valid cell pairs found after removing NAs.\n")
    next
  }
  
  summary_distances <- summary_distances %>%
    rowwise() %>%
    mutate(
      num_ref_cells = length(unique(distances$RefCell[distances$RefType == Reference]))
    ) %>%
    ungroup()
  
  cells_of_interest <- c(paste0(th, "+ Th cell"), paste0(th, "- Th cell"))
  output <- rbind(output, summary_distances) %>%
    filter(Target %in% cells_of_interest)
}

th_results <- rbind(th_results, output)

average_percentage_of_cells_within_radius(spe_object = spe, 
                                          reference_celltype = "Th cell", 
                                          target_celltype = "Immune2", 
                                          radius=100, feature_colname="Cell.Type")



