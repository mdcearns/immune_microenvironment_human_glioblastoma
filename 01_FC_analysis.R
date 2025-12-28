library(tidyverse)
library(readxl)
library(lme4)
library(lmerTest)
library(compositions)
library(forcats)
library(cowplot)

# Read in and organise data 

data <- readxl::read_xls("GBM_single_cell_analysis/data/thesis_FC_data.xls")

flow <- data %>%
  mutate(
    Donor = str_extract(Sample, "^.{3}"),
    Timepoint = ifelse(grepl("Primary", Sample), "Primary", "Recurrence1"),
    Tissue = ifelse(grepl("Tumour", Sample), "Tumour", "PBZ"),
    T = CD4 + CD8 + DN,
    CTL = CD8 + DN,
    Myeloid = Macs + MG,
    Immune = T + Macs + MG + NK + Other,
    Ca = Live - Immune,
    p.total_Ca = ((Live - Immune) / Live) * 100,
    p.total_Immune = (Immune / Live) * 100,
    p.total_T = (T / Live) * 100,
    p.total_CD4 = (CD4 / Live) * 100,
    p.total_CTL = (CTL / Live) * 100,
    p.total_Macs = (Macs / Live) * 100,
    p.total_MG = (MG / Live) * 100,
    p.total_NK = (NK / Live) * 100,
    p.total_Other = (Other / Live) * 100,
    p.immune_T = (T / Immune) * 100,
    p.immune_CD4 = (CD4 / Immune) * 100,
    p.immune_CTL = (CTL / Immune) * 100,
    p.immune_Macs = (Macs / Immune) * 100,
    p.immune_MG = (MG / Immune) * 100,
    p.immune_NK = (NK / Immune) * 100,
    p.immune_Other = (Other / Immune) * 100,
    p.myeloid_Macs = (Macs / (Macs + MG) * 100),
    p.myeloid_MG = (MG / (Macs + MG) * 100),
    p.T_CD4 = (CD4 / T) * 100,
    p.T_CTL = (CTL / T) * 100
  ) %>%
  select(
    Sample, Donor, Timepoint, Tissue, Immune, Ca, Myeloid, T, CTL,
    everything()
  )

primary <- flow %>%
  filter(Timepoint == "Primary")

primary$Tissue <- factor(primary$Tissue, levels = c("PBZ", "Tumour"))

model <- glmer(cbind(Immune, Live - Immune) ~ Tissue + (1 | Donor),
               family = binomial,
               data = primary)
summary(model)

model <- glmer(cbind(MG, Myeloid - MG) ~ Tissue + (1 | Donor),
               family = binomial,
               data = primary)
summary(model)

model <- glmer(cbind(CD4, T - CD4) ~ Tissue + (1 | Donor),
               family = binomial,
               data = primary)
summary(model)

# Data for overall cell proportions by sample plot

immune <- flow %>%
  select(Sample, Donor, Timepoint, Tissue, p.total_Immune, p.total_Ca) %>%
  pivot_longer(cols = c(p.total_Immune, p.total_Ca), names_to = "Type", values_to = "Proportion")
  
donor_order <- c("N03", "N02", "N06", "N07", "N08", "N01", "N05")

immune <- immune %>%
  mutate(Donor = factor(Donor, levels = donor_order)) %>%
  arrange(Donor, desc(Tissue), desc(Type)) %>%
  mutate(Tissue = factor(Tissue, levels = unique(Tissue)),
         Sample = factor(Sample, levels = unique(Sample)),
         Type = factor(Type, levels = c("p.total_Immune", "p.total_Ca")))

# Plot

immunePlot <- ggplot(
  immune, aes(x = Sample, y = Proportion, fill = Type)
  ) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
  facet_wrap(~ Donor, scales = "free_x", nrow = 1) +  # to visually group by Donor
  scale_fill_manual(values = c("p.total_Immune" = "#33A02C", "p.total_Ca" = "gray70"),
                    labels = c("Immune (any)", "Non-immune")) +
  labs(x = "Sample", y = "Overall cell proportions", fill = "Cell Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, color = "black", vjust = 0.5, hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_text(color = "black", face = "bold"),
        strip.text = element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        legend.position = "none")

# Data for immune cell types plot

subtypes <- flow %>%
  select(Sample, Donor, Timepoint, Tissue, p.immune_CTL, p.immune_CD4, p.immune_Macs, p.immune_MG, p.immune_NK, p.immune_Other) %>%
  pivot_longer(cols = c(p.immune_CTL, p.immune_CD4, p.immune_Macs, p.immune_MG, p.immune_NK, p.immune_Other), names_to = "Type", values_to = "Proportion") %>%
  mutate(Donor = factor(Donor, levels = donor_order)) %>%
  arrange(Donor, desc(Tissue), desc(Type)) %>%
  mutate(Tissue = factor(Tissue, levels = unique(Tissue)),
         Sample = factor(Sample, levels = unique(Sample)),
         Type = factor(Type, levels = c("p.immune_Other", "p.immune_NK", "p.immune_MG", "p.immune_Macs", "p.immune_CTL","p.immune_CD4")))

# Plot
subtypesPlot <- ggplot(
  subtypes, aes(x = Sample, y = Proportion, fill = Type)
) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
  facet_wrap(~ Donor, scales = "free_x", nrow = 1) +  # optional, to visually group by Donor
  scale_fill_manual(values = c("lightgreen", "#FF7F00", "#CAB2D6", "#6A3D9A",  "#E31A1C", "#FB9A99"),
                    labels = c("Immune (other)", "NK cells", "Myeloid cells (CD45low)", "Myeloid cells (CD45high)", "CD4- T cells", "CD4+ T cells")) +
  labs(x = "Sample", y = "Immune cell proportions", fill = "Cell Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, color = "black", vjust = 0.5, hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_text(color = "black", face = "bold"),
        strip.text = element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        legend.position = "none")
subtypesPlot

output_dir <- "GBM_single_cell_analysis/outputs/flow/"
ggsave(paste0(output_dir, "overall_proportions.png"), immunePlot,
              height = 4, width = 4)
ggsave(paste0(output_dir, "subtypes.png"), subtypesPlot,
       height = 4, width = 4)

# Overall plot by region

model <- glmer(cbind(Immune, Live - Immune) ~ Tissue + (1 | Donor),
               family = binomial,
               data = primary)

overall_df <- data.frame(
  CellType = "Immune",
  Estimate = exp(fixef(model)["TissueTumour"]),
  SE = exp(summary(model)$coefficients["TissueTumour", "Std. Error"]),
  p_value = ifelse(summary(model)$coefficients["TissueTumour", "Pr(>|z|)"] < 0.0001, format(0.0001, scientific = F), signif(summary(model)$coefficients["TissueTumour", "Pr(>|z|)"], 2))
)

overall_df <- overall_df %>%
  mutate(
    Upper = Estimate + 1.96 * SE,
    Lower = Estimate - 1.96 * SE,
  )

overallDiff <- ggplot(overall_df, aes(x = CellType, y = Estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  geom_point(color = "#33A02C", size = 5) +
  geom_linerange(aes(ymin = Lower, ymax = Upper),
                 color = "#33A02C", size = 1) +
  geom_text(
    aes(label = ifelse(p_value == "0.0001", paste0("p < ", p_value), paste0("p = ", p_value)), y = Upper + 5),
    color = "black",
    size = 5
  ) +
  coord_flip() +
  labs(y = "Odds Ratio (relative abundance)", x = NULL) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15),
    axis.title.x = element_text(size = 19, hjust = 0.42),
    plot.background = element_rect(fill = "white")) +
  scale_y_continuous(limits = c(-20, 25),
                     breaks = seq(-20, 20, 10),
                     labels = seq(-20, 20, 10))

overallDiff

ggsave(paste0(output_dir, "overall_difference.png"), overallDiff,
              height = 4, width = 9)
  
# Subtype plot by region

celltypes <- c("p.total_CD4", "p.total_CTL", "p.total_Macs", "p.total_MG", "p.total_NK", "p.total_Other", "p.total_Ca")

# Get data formatted for CLR-transformation and LMM compositional analysis

primary_clr_input <- primary[, celltypes] /100 # switch back to proportions
min_nonzero <- min(primary_clr_input[primary_clr_input > 0])
pseudocount <- min_nonzero / 100
primary_clr_adj <- primary_clr_input
primary_clr_adj[primary_clr_adj == 0] <- pseudocount
primary_clr_adj <- primary_clr_adj / rowSums(primary_clr_adj)

clr_data <- clr(primary_clr_adj)
clr_df <- cbind(primary[, c("Donor", "Tissue")], as.data.frame(clr_data))

sub_df <- data.frame(
  CellType = celltypes,
  Estimate = NA,
  SE = NA,
  p_value = NA,
  Upper = NA,
  Lower = NA,
  Donor = NA,
  stringsAsFactors = FALSE
)

for (i in seq_along(celltypes)) {
  var <- celltypes[i]
  model_formula <- as.formula(paste(var, "~ Tissue + (1 | Donor)"))
  model <- lmer(model_formula, data = clr_df)
  
  est <- fixef(model)["TissueTumour"]
  se <- summary(model)$coefficients["TissueTumour", "Std. Error"]
  donor <- attr(summary(model)$varcor$Donor, "stddev")
  
  sub_df[i, "Estimate"] <- est
  sub_df[i, "SE"] <- se
  sub_df[i, "p_value"] <- ifelse(summary(model)$coefficients["TissueTumour", "Pr(>|t|)"] < 0.0001, "p < 0.0001", signif(summary(model)$coefficients["TissueTumour", "Pr(>|t|)"], 2))
  sub_df[i, "Upper"] <- est + 1.96 * se
  sub_df[i, "Lower"] <- est - 1.96 * se
  sub_df[i, "Donor"] <- donor
}

sub_df$FoldChange <- ifelse(sub_df$Estimate < 0, -exp(abs(sub_df$Estimate)), exp(sub_df$Estimate))
sub_df$Upper_FC <- ifelse(sub_df$Upper < 0, -exp(abs(sub_df$Upper)), exp(sub_df$Upper))
sub_df$Lower_FC <- ifelse(sub_df$Lower < 0, -exp(abs(sub_df$Lower)), exp(sub_df$Lower))
sub_df$Donor_FC <- exp(sub_df$Donor)

sub_df <- sub_df %>%
  arrange(FoldChange) %>%
  mutate(CellType = factor(CellType, levels = unique(CellType)))

subDiff <- ggplot(sub_df, aes(x = CellType, y = Estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  geom_point(color = rev(c("#FB9A99", "#E31A1C", "lightgreen", "#6A3D9A", "#CAB2D6", "#FF7F00","gray70")), size = 5) +
  geom_linerange(aes(ymin = Lower, ymax = Upper),
                 color = rev(c("#FB9A99", "#E31A1C", "lightgreen", "#6A3D9A", "#CAB2D6", "#FF7F00", "gray70")), size = 1) +
  geom_text(
    aes(label = ifelse(p_value == "0.0001", paste0("p < ", p_value), paste0("p = ", p_value)), y = max(Upper) + 1, hjust = 0),
    color = "black",
    size = 5
  ) +
  coord_flip() +
  labs(y = "Centered log-ratio change (compositional)", x = NULL) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15),
    axis.title.x = element_text(size = 19, hjust = 0.35),
    plot.background = element_rect(fill = "white")) +
  scale_y_continuous(limits = c(min(sub_df$Lower), max(sub_df$Upper) + 1.57),
                     breaks = seq(-4, 4, 1),
                     labels = seq(-4, 4, 1)
  )

ggsave(paste0(output_dir, "sub_difference.png"), subDiff,
       height = 4, width = 9)

write.xlsx(sub_df, paste0(output_dir, "sub_difference_models.xlsx"))

legend_df <- data.frame(
  Category = factor(c("Immune (all)", "CD4+ T cells", "CD4- T cells", "Macrophages", "Microglia", "NK cells", "Immune (other)", "Non-immune")),
  x = 1,
  y = 1
)

legend_cols <- c(
  "Immune (all)" = "#33A02C",
  "CD4+ T cells" = "#FB9A99",
  "CD4- T cells" = "#E31A1C",
  "Macrophages" = "#6A3D9A",
  "Microglia" = "#CAB2D6",
  "NK cells" =  "#FF7F00",
  "Immune (other)" = "lightgreen",
  "Non-immune" = "gray70"
)

legend_df$Category <- factor(legend_df$Category, levels = names(legend_cols))

legend_plot <- ggplot(legend_df, aes(x, y, fill = Category)) +
  scale_fill_manual(values = legend_cols) +
  geom_tile() +
  theme_void() +
  labs(fill = "Cell Type") +
  theme(legend.position = "right",
        plot.background = element_rect(fill = "white", color = NA))

# Extract the legend
legend_only <- get_legend(legend_plot)
grid.newpage()
grid.draw(legend_plot)

ggsave(paste0(output_dir, "legend.png"), legend_plot,
       height = 4, width = 7)
