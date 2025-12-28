# Load libraries ####

library(openxlsx)
library(dplyr)
library(stringr)
library(tidyr)
library(purrr)
library(ggplot2)
library(scales)
library(cowplot)
library(forcats)
library(janitor)

set.seed(2025)

# Read in TMA-level data ####

maps <- list.files("GBM_spatial_analysis/data/TMA_maps")

batchInfo <- data.frame(TMA = str_replace_all(maps, c("TMA_" = "", "_map.xlsx$" = "")))

data/TMA_maps_list <- list()
for (tma in seq_along(maps)){
  TMA_map <- read.xlsx(paste0("GBM_spatial_analysis/data/TMA_maps/", maps[[tma]]), colNames = FALSE)
  data/TMA_maps_list[[tma]] <- apply(TMA_map, 1, function(row) as.list(row))
  rm(TMA_map)
} # Read in TMA maps

cores_list <- unlist(data/TMA_maps_list, recursive = TRUE) # Convert to flat list of core names

coreData <- data.frame(Core_ID = cores_list,
                       Patient_ID = NA,
                       Timepoint = NA,
                       TMA = NA,
                       Num_timepoints = NA)

coreData <- coreData %>%
  mutate(
    Patient_ID = str_extract(Core_ID, "(?<=N)\\d{1}(\\d{2})|N\\d{1}(\\d{2})|\\d{3}"),
    Patient_ID = str_replace(Patient_ID, "^N\\d", "N"),
    Timepoint = case_when(
      grepl("aa|ab|ac", Core_ID) ~ "Primary",
      grepl("ba|bb|bc", Core_ID) ~ "Recurrence1",
      grepl("ca|cb|cc", Core_ID) ~ "Recurrence2",
      TRUE ~ "Other"
      ),
    extracted_char1 = str_extract(cores_list, "(?<=\\d{3}).*"), # extract Num_timepoints from core name
    extracted_char2 = substr(extracted_char1, 2, 2),
    TMA = str_extract(Core_ID, "^[^_]+"),
    Num_timepoints = case_when(
      extracted_char2 == "a" ~ 1,
      extracted_char2 == "b" ~ 2,
      extracted_char2 == "c" ~ 3,
      TRUE ~ NA_integer_
    )
  ) %>%
  select(-extracted_char1, -extracted_char2)

rm(cores_list, data/TMA_maps_list, maps)

control_coreData <- coreData %>% # Extract control cores
  filter(
    str_detect(Core_ID, "tonsil|placenta|liver|colon|small_bowel|cervix|spleen")
  )

control_coreData <- control_coreData %>%
  mutate(Core_ID = str_extract(Core_ID, "^[^_]+_[^_]+")) %>%
  select(Core_ID, Timepoint, TMA)

coreData <- coreData %>% # Remove blank/control cores
  mutate(Core_ID = paste(str_extract(TMA, ".*"),
                       str_extract(Core_ID, "(?<=_)[^_]+(?=_)"), 
                       sep = "_")) %>%
  filter(Timepoint != "Other") 

# Add number of cores and control cores to batchInfo:

core_counts <- coreData %>%
  group_by(TMA) %>%
  summarise(Num_cores = n()) %>%
  ungroup()

control_core_counts <- control_coreData %>%
  group_by(TMA) %>%
  summarise(Num_control_cores = n()) %>%
  ungroup()

batchInfo <- batchInfo %>%
  left_join(core_counts, by = "TMA") %>%
  left_join(control_core_counts, by = "TMA")

rm(core_counts, control_core_counts)

# Read in patient-level data ####

clinical_dataset <- read.xlsx("GBM_spatial_analysis/data/thesis_dataset.xlsx", colNames = TRUE)
names(clinical_dataset) <- make.names(names(clinical_dataset), unique = TRUE)

Num_timepoints_merge <- coreData %>%
  select(Patient_ID, Num_timepoints) %>%
  distinct(Patient_ID, .keep_all = TRUE) # Extract Num_timepoints to add to patientData

patientData <- clinical_dataset %>%
  filter(Patient_ID %in% coreData$Patient_ID) %>%
  left_join(Num_timepoints_merge, by = "Patient_ID") %>%
  select(Patient_ID, Cohort, Num_timepoints, Age_at_diagnosis, Sex, 
         PC_headache1, PC_vomiting1, PC_seizure1, PC_cognitive_change1, PC_focal_deficit1, Laterality1, Location1,
         WCC1, Neuts1, Lymphs1, Op_type1, X5.ALA1, RT_Gy1, RT_fractions1, Chemo1a, Concom1, Adj_cycles1,
         Diagnosis1, WHO_grade1, Postop_infection1, Diagnosis2, WHO_grade2, Postop_infection2,
         Diagnosis3, WHO_grade3, Postop_infection3, PFS_days, PFS_censored, OS_days, OS_censored, IDH_final, MGMT_final) %>%
  mutate(
    Age_at_diagnosis = round(Age_at_diagnosis, 1),
    WCC1 = round(WCC1, 1),
    Neuts1 = round(as.numeric(Neuts1), 2),
    Lymphs1 = round(as.numeric(Lymphs1), 2),
    OS_months = OS_days / 30.44,
    PFS_months = PFS_days / 30.44
  )

patientData$X5.ALA1[patientData$Patient_ID == "143"] <- 1

keyData <- patientData %>%
  filter(Cohort != "Other")

rm(clinical_dataset, Num_timepoints_merge)

coreData <- coreData %>%
  select(-Num_timepoints) %>% # Remove now passed to patientData
  select(Core_ID, everything(), TMA)

# Add in imaging data ####

imaging <- read.xlsx("GBM_spatial_analysis/data/LTS_imaging_data.xlsx") %>%
  clean_names()

names(imaging)

imaging <- imaging %>%
  dplyr::rename(
    Patient_ID = study_no,
    Complete = complete_resection_of_ce,
    Near_total = near_total_resection_of_ce_85_percent_99_9_percent_1_cm3_residual_tumour,
    Sub_total = subtotal_resection_of_ce_80_percent_94_9_percent_5cm3_residual_tumour,
    Partial = partial_resction_of_ce_tumour_80_percent_ce_5_cm3_residual_tumour,
    Biopsy = biopsy
  ) 

imaging[25, "Sub_total"] <- NA

imaging <- imaging %>%
  select(
    Patient_ID,
    pre_op_ce,
    pre_op_non_ce,
    Complete,
    Near_total,
    Sub_total,
    Partial,
    Biopsy,
    ce_contacting_cortex,
    ce_crossing_midline,
    ce_touching_ventricular_surface
  ) %>%
  mutate(
    Complete   = if_else(!is.na(Complete),   "Complete",   NA_character_),
    Near_total = if_else(!is.na(Near_total), "Near_total", NA_character_),
    Sub_total  = if_else(!is.na(Sub_total),  "Sub_total",  NA_character_),
    Partial    = if_else(!is.na(Partial),    "Partial",    NA_character_),
    Biopsy     = if_else(!is.na(Biopsy),     "Biopsy",     NA_character_),
    No_imaging = if_else(
      if_all(c(Complete, Near_total, Sub_total, Partial, Biopsy), is.na),
      "No_imaging", NA_character_),
    Gad_cortex = case_when(
      ce_contacting_cortex %in% c("Y", "y") ~ "Yes",
      ce_contacting_cortex %in% c("N", "n") ~ "No",
      TRUE ~ NA_character_),
    Gad_midline = case_when(
      ce_crossing_midline %in% c("Y", "y") ~ "Yes",
      ce_crossing_midline %in% c("N", "n") ~ "No",
      TRUE ~ NA_character_),
    Gad_vent = case_when(
      ce_touching_ventricular_surface %in% c("Y", "y") ~ "Yes",
      ce_touching_ventricular_surface %in% c("N", "n") ~ "No",
      TRUE ~ NA_character_)
  ) %>%
  pivot_longer(
    cols = c("Complete",
             "Near_total",
             "Sub_total",
             "Partial",
             "Biopsy",
             "No_imaging"),
    values_to = "EOR") %>%
  filter(!is.na(EOR)) %>%
  select(-name)

imaging[1:4, "Patient_ID"] <- c("004", "065", "026", "054")
imaging[8:9, "Patient_ID"] <- c("093", "016")
imaging[24, "pre_op_ce"] <- "11.2"
imaging[20, "pre_op_non_ce"] <- NA

imaging <- imaging %>%
  mutate(pre_op_ce = as.numeric(pre_op_ce))

imaging[26, "pre_op_non_ce"] <- "184.8"

imaging <- imaging %>%
  mutate(
    pre_op_ce = as.numeric(pre_op_ce),
    pre_op_non_ce = as.numeric(pre_op_non_ce),
    across(where(is.numeric), ~ round(.x, 1))
    )

keyData <- keyData %>%
  left_join(imaging, by = "Patient_ID"
)

keyData <- keyData %>%
  mutate(
    StuppRT = case_when(
      RT_Gy1 == "0" ~ "No",
      is.na(RT_Gy1) ~ "Unknown",
      TRUE ~ "Yes"),
    StuppConcomTMZ = case_when(
      Concom1 == "0" ~ "No",
      is.na(Concom1) ~ "Unknown",
      TRUE ~ "Yes"),
    StuppAdjCT = case_when(
      Adj_cycles1 == "0" ~ "None",
      Adj_cycles1 == "Y" ~ "CCNU",
      is.na(Adj_cycles1) ~ "Unknown",
      TRUE ~ "TMZ")
    )
  
lts_cohort <- filter(keyData, Cohort != "NAV")
nav_cohort <- filter(keyData, Cohort == "NAV")

count_summary <- lts_cohort %>%
  pivot_longer(cols = c(Sex, PC_headache1, PC_vomiting1, PC_seizure1, PC_cognitive_change1,
                        PC_focal_deficit1, Laterality1, Location1, Gad_cortex, Gad_midline, Gad_vent,
                        Op_type1, X5.ALA1, IDH_final, MGMT_final, EOR, StuppRT, StuppConcomTMZ, StuppAdjCT, PFS_censored, OS_censored),
               names_to = "Variable", values_to = "Category",
               values_transform = as.character) %>%
  count(Cohort, Variable, Category) %>%
  arrange(Cohort, Variable, desc(n))
  
contin_summary <- lts_cohort %>%
  group_by(Cohort) %>%
  summarise(
    across(c(WCC1, Neuts1, Lymphs1, pre_op_ce, pre_op_non_ce, OS_months, PFS_months),
           list(median = ~median(.x, na.rm = T),
                lower = ~quantile(.x, 0.25, na.rm = T),
                upper = ~quantile(.x, 0.75, na.rm = T)),
           .names = "{.col}_{.fn}"),
    across(Age_at_diagnosis,
           list(mean = ~mean(.x, na.rm = T),
                sd = ~sd(.x, na.rm = T)),
           .names = "{.col}_{.fn}")
  )
  
write.xlsx(count_summary,
           "GBM_spatial_analysis/outputs/cohort_count_data.xlsx")

write.xlsx(contin_summary,
           "GBM_spatial_analysis/outputs/cohort_continuous_data.xlsx")


table <- read.xlsx("GBM_spatial_analysis/outputs/cohort_table.xlsx")
View(table)

table$X4[!grepl("Median", table$X2)] <- round(as.numeric(table$X4), 1)
table$X6[!grepl("Median", table$X2)] <- round(as.numeric(table$X6), 1)
table$`Long-term.survivors.(13)` <- paste0(table$`Long-term.survivors.(13)`, " (", table$X4, ")")
table$`Short-term.survivors.(14)` <- paste0(table$`Short-term.survivors.(14)`, " (", table$X6, ")")

write.xlsx(table, "GBM_spatial_analysis/outputs/cohort_table.xlsx")

rm(contin_summary, count_summary, imaging, table)

# Statistical tests:

median(lts_cohort$Age_at_diagnosis[lts_cohort$Cohort == "STS"])
quantile(lts_cohort$Age_at_diagnosis[lts_cohort$Cohort == "STS"], 0.75)
quantile(lts_cohort$Age_at_diagnosis[lts_cohort$Cohort == "STS"], 0.25)
wilcox.test(Age_at_diagnosis ~ Cohort, data = lts_cohort)
wilcox.test(pre_op_non_ce ~ Cohort, data = lts_cohort)

lat_cohort <- lts_cohort %>%
  filter(Laterality1 %in% c("R", "L"))
fisher.test(table(lat_cohort$Cohort, lat_cohort$Laterality1))

op_cohort <- lts_cohort %>%
  mutate(optype = if_else(Op_type1 == "B", "B", "R"))
fisher.test(table(op_cohort$Cohort, op_cohort$optype))
tbl <- table(op_cohort$Cohort, op_cohort$optype)
tbl_adj <- tbl + 0.5
OR_adj <- (tbl_adj[1,1]*tbl_adj[2,2]) / (tbl_adj[1,2]*tbl_adj[2,1])
OR_adj

radio_cohort <- lts_cohort %>%
  mutate(radio = case_when(
    RT_Gy1 == "0" ~ "0",
    is.na(RT_Gy1) ~ "NA",
    TRUE ~ "1"
  )
  ) %>%
  filter(radio != "NA")

fisher.test(table(radio_cohort$Cohort, radio_cohort$radio), conf.int = TRUE, workspace = 2e8)
tbl <- table(radio_cohort$Cohort, radio_cohort$radio)
tbl_adj <- tbl + 0.5
OR_adj <- (tbl_adj[1,1]*tbl_adj[2,2]) / (tbl_adj[1,2]*tbl_adj[2,1])
OR_adj

chemo_cohort <- lts_cohort %>%
  mutate(chemo = case_when(
    Chemo1a == "TMZ" ~ "1",
    is.na(Chemo1a) ~ "NA",
    TRUE ~ "0"
  )
  ) %>%
  filter(chemo != "NA")

fisher.test(table(chemo_cohort$Cohort, chemo_cohort$chemo), conf.int = TRUE, workspace = 2e8)
tbl <- table(chemo_cohort$Cohort, chemo_cohort$chemo)
tbl_adj <- tbl + 0.5
OR_adj <- (tbl_adj[1,1]*tbl_adj[2,2]) / (tbl_adj[1,2]*tbl_adj[2,1])
OR_adj

wilcox.test(OS_days ~ Cohort, data = lts_cohort)
wilcox.test(PFS_days ~ Cohort, data = lts_cohort)

# For separate batch import:

intensityMatrix_list <- list()

#intensityMatrix_files <- list.files("GBM_spatial_analysis/data/intensityMatrix_files",
#                                    full.names = TRUE)

#for (tma in seq_along(intensityMatrix_files)){
 # intensityMatrix <- read.csv(paste(intensityMatrix_files[[tma]]))
 # intensityMatrix <- intensityMatrix %>%
 # intensityMatrix_list[[tma]] <- intensityMatrix
 # rm(intensityMatrix)
#}

# For global import:

intensityMatrix <- read.csv("GBM_spatial_analysis/data/intensityMatrix_files_thesis/thesis_measurements.csv")

intensityMatrix_list[[1]] <- intensityMatrix %>%
  select(-3, -9) %>%
  filter(Image == "LTS")

intensityMatrix_list[[2]] <- intensityMatrix %>%
  select(-3, -9) %>%
  filter(Image == "NAV")

cellData_list <- list() # Create dataframe lists
annotationData_list <- list()
control_cellData_list <- list()
na.cells_list <- list()
batchInfo$Cells_segmented <- NA
batchInfo$After_controls_removed <- NA
batchInfo$After_annotation <- NA

for (tma in seq_along(batchInfo$TMA)){

  print(paste("Processing file", tma))
  
  cellData <- intensityMatrix_list[[tma]] %>% 
    dplyr::rename(TMA = Image,
                  Cell_ID = Object.ID,
                  Classifier = Classification,
                  Core_ID = TMA.core,
                  Annotation_Type = Parent,
                  X_coord = Centroid.X.µm,
                  Y_coord = Centroid.Y.µm,
                  Cell_area = Cell..Area.µm.2) %>%
    dplyr::rename_with(~ gsub(".*\\.\\.(.*?)\\.\\..*", "\\1", .), contains("..")) %>%
    dplyr::rename("MHC Class I" = MHC.Class.I,
                  "MHC Class II" = MHC.Class.II,
                  "PD-L2" = PD.L2,
                  "CTLA-4" = CTLA.4,
                  "Granzyme K" = Granzyme.K,
                  "Granzyme B" = Granzyme.B,
                  "PD-1" = PD.1,
                  "LAG-3" = LAG.3,
                  "PD-L1" = PD.L1,
                  "Granzyme A" = Granzyme.A) %>%
    mutate(Core_ID = paste0(str_extract(TMA, ".*"), "_", gsub("-", "", Core_ID)),
           Annotation_Type = str_extract(Annotation_Type, "(?<=\\().*(?=\\))"))
  
  annotationData <- cellData %>% # Extract annotation level details
    select(TMA, Core_ID, Annotation_Type) %>%
    distinct() # Remove duplicate rows
  
  Annotation_ID_merge <- coreData %>%
    select(-TMA)
  
  annotationData <- annotationData %>%
    left_join(Annotation_ID_merge, by = "Core_ID") %>% # Get details for Annotation_ID
    mutate(Original_Order = row_number(),
           Annotation_ID = paste0(Patient_ID, "_", Timepoint, "_", Annotation_Type)) %>%
    group_by(Patient_ID, Timepoint, Annotation_Type) %>% # Identify duplicate Annotation_Type per patient
    mutate(n = n(),
           Annotation_ID = paste0(Annotation_ID, "_", row_number())) %>% # Dynamically add suffix
    ungroup() %>%
    arrange(Original_Order) %>%
    select(Annotation_ID, Core_ID, TMA, Annotation_Type) # Re-order and remove unnecessary columns
  
  cellData <- cellData %>% # Add Annotation_ID to cells
    left_join(annotationData, by = c("Core_ID", "Annotation_Type", "TMA")) %>%
    select(-Annotation_Type) %>%
    select(Cell_ID, Annotation_ID, Core_ID, TMA, everything()) 
  
  batchInfo$Cells_segmented[tma] <- nrow(cellData)
  
  control_cellData <- cellData %>% # Extract out and then remove cells from control cores
    filter(grepl("Control", Annotation_ID)) 
  cellData <- cellData %>%
    filter(!grepl("Control", Annotation_ID)) %>% 
    mutate_at(vars(8:ncol(cellData)), ~ replace(., is.na(.), 0)) # assuming standardised data export
  
  batchInfo$After_controls_removed[tma] <- nrow(cellData)
  
  na.cells <- cellData %>% # Extract and remove unannotated cells
    filter(grepl("NA", Annotation_ID))
  cellData <- cellData %>%
    filter(!grepl("NA", Annotation_ID))
  annotationData <- annotationData %>% # Remove NA and Control annotations
    filter(!is.na(Annotation_Type),
           !grepl("Control", Annotation_Type))
 
  batchInfo$After_annotation[tma] <- nrow(cellData)

  cellData_list[[tma]] <- cellData # Re-embed in lists
  annotationData_list[[tma]] <- annotationData
  control_cellData_list[[tma]] <- control_cellData
  na.cells_list[[tma]] <- na.cells
} # Organise data

# NAV annotations re-named to match single cell dataset (assuming [[2]] in list)
# NOTE: in navigate patients (N01 - N09), Annotation_ID matches tumour/PBZ allocations and
# Annotation_Type indicates the actual histology of the core (inc. necrosis)
# in some cases this results in two Tumour annotations from one core - must use
# Annotation_Type if interested in histology
# hyphenated Annotation_ID = image-guided sample, otherwise underscored
# Specific note: "N05_Primary_" annotations - spatially separated regions chosen from
# historical blocks to mimic image-guided sampling method in Recurrence1 sample

nav_renamed_annotations <- read.xlsx("GBM_spatial_analysis/data/NAV_renamed_annotations.xlsx")
nav_renamed_annotations <- nav_renamed_annotations %>%
  rename(Assigned = Assigned.Annotation_ID, Final = Final.Annotation_ID)

annotationData_list[[2]] <- annotationData_list[[2]] %>%
  left_join(nav_renamed_annotations, by = c("Annotation_ID" = "Assigned")) %>% # left_join by Annotation_ID which matches "Assigned"
  mutate(Annotation_ID = ifelse(!is.na(Final), Final, Annotation_ID)) %>% # replaces Annotation_ID with contents of "Final" if populated
  select(-Final)

cellData_list[[2]] <- cellData_list[[2]] %>%
  left_join(nav_renamed_annotations, by = c("Annotation_ID" = "Assigned")) %>%
  mutate(Annotation_ID = ifelse(!is.na(Final), Final, Annotation_ID)) %>%
  select(-Final)

rm(cellData, control_cellData, na.cells, intensityMatrix_list,
   intensityMatrix, Annotation_ID_merge, nav_renamed_annotations)

# Integrate annotation and na.cells data (leave cell level data until after QC)
annotationData <- dplyr::bind_rows(annotationData_list)
na.cells <- dplyr::bind_rows(na.cells_list)
rm(annotationData_list, na.cells_list)

gc()

# Batch quality control ####

markerList <- colnames(cellData_list[[1]])[9:ncol(cellData_list[[1]])] # NB check if re-export

long_cellData_list <- lapply(cellData_list, function(cellData){
  cellData %>%
    pivot_longer(cols = 9:ncol(cellData),
                 names_to = "Marker",
                 values_to = "Intensity")
})

# Plot cell area by TMA and filter cells

batchInfo$After_cell_area_filtering <- NA

for (tma in seq_along(batchInfo$TMA)){
  
  print(paste("Plotting cell area for TMA", batchInfo$TMA[[tma]]))
  
  cellData <- cellData_list[[tma]]
  
  cell_area_data <- long_cellData_list[[tma]] %>%
    filter(TMA == batchInfo$TMA[tma]) %>%
    arrange(Cell_area)
  
  cell_area_violin <- ggplot(cell_area_data, aes(x = batchInfo$TMA[[tma]], y = Cell_area)) +
    geom_violin(trim = FALSE, fill = "steelblue", alpha = 0.7) + 
    theme_minimal() +
    labs(title = "Cell area by TMA",
         x = paste("TMA",
                   batchInfo$TMA[[tma]]),
         y = "Cell area (µm2)") +
    theme(plot.title = element_text(hjust = 0.5, size = 26, face = "bold"),
          axis.title.x = element_text(vjust = 0.5, size = 22, face = "bold"),
          axis.title.y = element_text(vjust = 0.5, size = 22, face = "bold"),
          axis.text.y = element_text(size = 18),
          axis.text.x = element_blank(),
          panel.background = element_rect(fill = "white", color = NA),
          plot.background = element_rect(fill = "white", color = NA),
          strip.text = element_blank()) +
    facet_wrap(~ TMA) +
    geom_hline(yintercept = 3.5, color = "red", linetype = "dashed", linewidth = 0.5)
  
  print("Saving...")
  
  output_dir <- paste0("GBM_spatial_analysis/outputs/QC/", batchInfo$TMA[[tma]])
  dir.create(output_dir, recursive=T, showWarnings=F)
  ggsave(paste0(output_dir, "/", batchInfo$TMA[[tma]], "_cell_area_violin.png"),
         cell_area_violin,
         height = 12, width = 8)
  
  rm(cell_area_violin)
  
  print("Filtering...")
  
  cellData <- cellData %>%
    filter(Cell_area > 3.5)
  
  batchInfo$After_cell_area_filtering[tma] <- nrow(cellData)
  cellData_list[[tma]] <- cellData
  rm(cellData, cell_area_data)
  
  print("Complete!")
}

rm(cell_area_data, cell_area_violin)

# Marker intensity plots:

Marker_Violin_PlotList <- list() # structure: Marker_Violin_Plotlist[[tma]][[marker]]

top0.1_cells <- list() # for Elbow plot

# Note, adjusted to Elbow Plot top 0.1% of each marker, original QC
# deleted top 0.1% and Elbow Plot of top 0.5%. Plots adjusted accordingly

for (tma in seq_along(batchInfo$TMA)){

  Marker_Violin_PlotList[[tma]] <- list()
  top0.1_cells[[tma]] <- list()
  
  for (marker in seq_along(markerList)){

    print(paste0("Processing TMA ", batchInfo$TMA[[tma]], ", marker ", markerList[[marker]]))
    
    marker_data <- long_cellData_list[[tma]] %>%
      filter(TMA == batchInfo$TMA[tma], # filter by TMA and by marker
             Marker == markerList[marker]) %>%
      arrange(Intensity)
    
    # Marker Violin Plots:
    
    print("Generating Marker Violin Plot...")
    
    # Find 99.9th percentile of cells by intensity
    
    max_intensity <- round(max(marker_data$Intensity, na.rm = TRUE))
    
    p99.9_row <- round(0.999 * nrow(marker_data))
    p99.9_intensity <- marker_data$Intensity[p99.9_row]
    p99.9_cells <- marker_data %>%
      filter(Intensity >= p99.9_intensity)
    
    top0.1_cells[[tma]][[marker]] <- p99.9_cells$Cell_ID
    
    # Set the ylim and y_annotation placement for consistency
    
    ylim_lower <- min(marker_data$Intensity, na.rm = TRUE)
    intensity_diff <- p99.9_intensity - ylim_lower
    ylim_upper <- intensity_diff * 1.4
    y_annotation <- p99.9_intensity + 0.05 * (ylim_upper - p99.9_intensity)
    
    # Plot over distribution
    
    p <- ggplot(
      marker_data,
      aes(x = Marker,
          y = Intensity)
    ) +
      geom_violin(trim = FALSE, fill = "steelblue", alpha = 0.7) + 
      theme_minimal() +
      coord_flip() +  # Flip axes for better readability
      labs(x = markerList[marker]) +
      theme(axis.title.x = element_text(hjust = 0.5, size = 24, face = "bold"),
            axis.title.y = element_text(vjust = 0.5, size = 24, face = "bold", angle = 0),
            axis.text.x = element_text(size = 15),
            axis.text.y = element_blank(),
            panel.background = element_rect(fill = "white", color = NA),
            plot.background = element_rect(fill = "white", color = NA)) +
      scale_y_continuous(expand = c(0, 0), limits = c(ylim_lower, ylim_upper)) +  # Limit the y-axis
      scale_x_discrete(expand = c(0, 0)) +
      geom_hline(yintercept = p99.9_intensity, color = "darkorange", linetype = "dashed", linewidth = 1) + # Add line for 99.9th percentile
      annotate("text", x = 1.4, y = y_annotation,
               label = "top 0.1% of cells",
               color = "darkorange", hjust = 0, size = 6) +
      annotate("text", x = 1.3, y = y_annotation,
               label = paste0("(intensity = ", round(p99.9_intensity, 0), ")"),
               color = "darkorange", hjust = 0, size = 6, fontface = "italic") +
      # Annotate the maximum intensity value
      annotate("text", x = 0.9, y = y_annotation,
               label = paste("to max intensity:", round(max_intensity, 0)),
               color = "black", hjust = 0, size = 6, fontface = "italic")
    p
    
    output_dir <- paste0("GBM_spatial_analysis/outputs/QC/",
                         batchInfo$TMA[[tma]],
                         "/marker_plots/")
    dir.create(output_dir, recursive=T, showWarnings=F)
    
    ggsave(paste0(output_dir,
                  batchInfo$TMA[[tma]], "_",
                  markerList[[marker]], "_",
                  "marker_violin.png"),
           p,
           height = 4, width = 12)
    
    Marker_Violin_PlotList[[tma]][[marker]] <- p
    
  } # Plot marker violin plots with 99.9th percentile intensity
  
  rm(positive_cells, negative_cells)
  
  top0.1_cells[[tma]] <- unlist(top0.1_cells[[tma]])
  
  print("Complete!")
  
  }

rm(intensity_diff, p99.9_intensity, p99.9_row,
   y_annotation_lower, y_annotation_mid, y_annotation_upper,
   ylim_lower, ylim_upper, Marker_Violin_PlotList, output_dir,
   marker_data, p)

# Elbow Plot of number of cells in top 0.1% by number of markers

top0.1_df_list <- list()
top0.1_counts_list <- list()
top0.1_ElbowPlot_List <- list()

for (tma in seq_along(batchInfo$TMA)){
  
  top0.1_df_list[[tma]] <- as.data.frame(top0.1_cells[[tma]]) %>%
    rename(Cell_ID = 1) %>%  # Rename the single column to a known name
    group_by(Cell_ID) %>%   # Group by actual column name, not an external variable
    summarise(marker_number = n(), .groups = "drop") %>%
    ungroup() %>%
    arrange(marker_number)
  
  top0.1_counts_list[[tma]] <- top0.1_df_list[[tma]] %>%
    count(marker_number) %>%
    rename(cell_number = n) %>%
    arrange(marker_number)
  
  p <- ggplot(
    top0.1_counts_list[[tma]],
    aes(x = marker_number,
        y = cell_number)) +
    labs(title = paste0("Highest-expressing cells by number of markers (TMA ",
                        batchInfo$TMA[[tma]],
                        ")"),
         x = "Number of markers",
         y = "Number of cells in top 0.1%") +
    geom_point() +
    geom_line() +
    scale_y_continuous(labels = label_number(),
                       breaks = seq(0,
                                    max(top0.1_counts_list[[tma]]$cell_number),
                                    by = 5000)) +
    scale_x_continuous(labels = seq(0,
                                    max(top0.1_counts_list[[tma]]$marker_number),
                                    by = 5),
                       breaks = seq(0,
                                    max(top0.1_counts_list[[tma]]$marker_number),
                                    by = 5)) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 14,
                                face = "bold",
                                hjust = 0.5),
      axis.title = element_text(size = 10,
                                face = "bold"),
      axis.text = element_text(size = 8)
    )
  
  p
  
  top0.1_ElbowPlot_List[[tma]] <- p
}

# Inspect Elbow Plot and adjust as necessary

# Plot with line and save, add cells to cells_to_exclude:

cells_to_exclude <- list()

for (tma in seq_along(batchInfo$TMA)){
  top0.1_ElbowPlot_List[[tma]] <- top0.1_ElbowPlot_List[[tma]] +
    geom_vline(xintercept = 4.5,
               color = "red",
               linetype = "dashed",
               linewidth = 0.5)
  
  output_dir <- paste0("GBM_spatial_analysis/outputs/QC/", batchInfo$TMA[[tma]])
  dir.create(output_dir, recursive=T, showWarnings=F)
  ggsave(paste0(output_dir, "/", batchInfo$TMA[[tma]], "_top0.1_elbow.png"),
         top0.1_ElbowPlot_List[[tma]],
         height = 5, width = 8)
  
  cells_to_exclude[[tma]] <- top0.1_df_list[[tma]] %>%
    filter(marker_number > 4) %>%
    distinct(Cell_ID) %>%
    pull(Cell_ID)
  
}

# Elbow Plot of number of cells by number of positive classifications

multi_classifier <- list()
multi_classifier_ElbowPlot_List <- list()

for (tma in seq_along(batchInfo$TMA)){

  multi_classifier[[tma]] <- cellData_list[[tma]] %>%
    mutate(Classifier_length = ifelse(Classifier == "", 0, str_count(Classifier, ":") + 1)) %>%
    arrange(Classifier_length)
  
  multi_classifier_table <- as.data.frame(table(multi_classifier[[tma]]$Classifier_length))
  colnames(multi_classifier_table) <- c("marker_number", "cell_number")
  multi_classifier_table$marker_number <- as.numeric(multi_classifier_table$marker_number)
  
  p <- ggplot(
    multi_classifier_table,
    aes(x = marker_number,
        y = cell_number,
        group = 1)) +
    labs(title = paste0("Positively-classified cells by number of markers (TMA ",
                        batchInfo$TMA[[tma]],
                        ")"),
         x = "Number of markers",
         y = "Number of cells classified positive") +
    geom_point() +
    geom_line() +
    scale_y_continuous(labels = label_number(),
                       breaks = seq(0,
                                    max(multi_classifier_table$cell_number),
                                    by = 20000)) +
    scale_x_continuous(labels = seq(0,
                                    max(as.numeric(multi_classifier_table$marker_number)),
                                    by = 5),
                       breaks = seq(0,
                                    max(as.numeric(multi_classifier_table$marker_number)),
                                    by = 5)) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 14,
                                face = "bold",
                                hjust = 0.5),
      axis.title = element_text(size = 10,
                                face = "bold"),
      axis.text = element_text(size = 8)
    )
  
  p
  
  multi_classifier_ElbowPlot_List[[tma]] <- p
}

# Inspect Elbow Plot and adjust as necessary

# Plot with line and save, add cells to cells_to_exclude:

for (tma in seq_along(batchInfo$TMA)){
  
  multi_classifier_ElbowPlot_List[[tma]] <- multi_classifier_ElbowPlot_List[[tma]] +
    geom_vline(xintercept = 20.5,
               color = "red",
               linetype = "dashed",
               linewidth = 0.5)
  
  output_dir <- paste0("GBM_spatial_analysis/outputs/QC/", batchInfo$TMA[[tma]])
  dir.create(output_dir, recursive=T, showWarnings=F)
  
  ggsave(paste0(output_dir, "/", batchInfo$TMA[[tma]], "_multi_classifier_elbow.png"),
         multi_classifier_ElbowPlot_List[[tma]],
         height = 5, width = 8)
  
  multi_classifier[[tma]] <- multi_classifier[[tma]] %>%
    filter(Classifier_length > 20) %>%
    distinct(Cell_ID) %>%
    pull(Cell_ID)
}

# Filter cells by cells_to_exclude

filtered_cells <- list()
batchInfo$After_marker_filtering <- NA

for (tma in seq_along(batchInfo$TMA)){

  cells_to_exclude[[tma]] <- unique(c(cells_to_exclude[[tma]], multi_classifier[[tma]]))
  
  print(paste0("Number of cells to exclude for TMA ",
               batchInfo$TMA[[tma]], ": ",
               length(cells_to_exclude[[tma]])))
  
  filtered_cells[[tma]] <- cellData_list[[tma]] %>%
    filter(!Cell_ID %in% cells_to_exclude[[tma]])
  
  batchInfo$After_marker_filtering[tma] <- nrow(filtered_cells[[tma]])
  
}

batchInfo

rm(p99.9_cells, top0.1_cells, top0.1_counts_list,
   top0.1_df_list, top0.1_ElbowPlot_List, cells_to_exclude,
   multi_classifier, multi_classifier_table, multi_classifier_ElbowPlot_List, p,
   max_intensity, y_annotation, output_dir, cellData_list)

# Classifier box plots:
# Note - plotted across entire TMA, original QC in LTS was across sample cores

long_filtered_cells <- list()
Classifier_Box_PlotList <- list()

for (tma in seq_along(batchInfo$TMA)){
  long_filtered_cells[[tma]] <- filtered_cells[[tma]] %>%
    pivot_longer(cols = 9:ncol(filtered_cells[[tma]]), # assuming standardised data export
                 names_to = "Marker",
                 values_to = "Intensity")
  
  Classifier_Box_PlotList[[tma]] <- list()
  
  for (marker in seq_along(markerList)){
    
    print(paste0("Processing TMA ", batchInfo$TMA[[tma]], ", marker ", markerList[[marker]]))

    classifier_data <- long_filtered_cells[[tma]] %>%
      filter(TMA == batchInfo$TMA[tma], # filter by TMA and by marker
             Marker == markerList[marker]) %>%
      arrange(Intensity) %>%
      mutate(marker_status = ifelse(grepl(paste0(markerList[[marker]], "($|:)"), Classifier), "Positive", "Negative"))
    
    # Classifier Box Plots:
    
    print("Generating Classifier Box Plot...")
    
    Classifier_Box_PlotList[[tma]][[marker]] <- list()

    # Create the plot
    q <- ggplot(classifier_data, aes(
      x = marker_status, 
      y = Intensity, 
      color = marker_status  # Keep TRUE/FALSE
    )) +
      geom_boxplot(outlier.shape = NA,
                   alpha = 0.5,
                   width = 0.8) +  # Boxplot without outliers
      geom_jitter(aes(
        color = marker_status),
        width = 0.15,
        alpha = 0.4,
        size = 0.5) +  # Jittered points
      stat_summary(fun.min = min, fun.max = max, geom = "errorbar", width = 0.6, linewidth = 0.7) +  # Whiskers
      labs(
        x = paste(markerList[marker], "Classifier"),
        y = paste(markerList[marker], "Intensity"),
        title = paste0("Classifier intensity (TMA ", batchInfo$TMA, ")"),
        color = "Classifier Status"  # Change legend title
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
        axis.title = element_text(face = "bold", size = 10),
        axis.text.x = element_blank(),
        legend.title = element_text("Classifier status"),
        legend.position = "bottom",  # Moves legend below
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)
      ) +
      scale_fill_manual(values = c("Positive" = "steelblue", "Negative" = "#FF6F61")) +  # Adjust fill colors
      scale_color_manual(values = c("Positive" = "steelblue", "Negative" = "#FF6F61")) +  # Adjust jitter colors
      guides(color = guide_legend(nrow = 1))  # Forces legend to be one row

    Classifier_Box_PlotList[[tma]][[marker]] <- q
    
    output_dir <- paste0("GBM_spatial_analysis/outputs/QC/",
                         batchInfo$TMA[[tma]],
                         "/classifier_plots/")
    dir.create(output_dir, recursive=T, showWarnings=F)
    
    ggsave(paste0(output_dir,
                  batchInfo$TMA[[tma]], "_",
                  markerList[[marker]], "_classifier_plot.png"),
           q,
           height = 5, width = 7)
  }

  print("Complete!")
  
}

rm(Classifier_Box_PlotList)

dyn_range_values <- data.frame(matrix(ncol = length(batchInfo$TMA),
                                      nrow = length(markerList)))
rownames(dyn_range_values) <- markerList

for (tma in seq_along(batchInfo$TMA)){

  colnames(dyn_range_values)[tma] <- batchInfo$TMA[[tma]]
  
  for (marker in seq_along(markerList)){
    
    intensities <- filtered_cells[[tma]][[markerList[marker]]]
    signal <- mean(intensities[intensities >= quantile(intensities, 0.8, na.rm = TRUE)], na.rm = TRUE)
    noise <- mean(intensities[intensities <= quantile(intensities, 0.1, na.rm = TRUE)], na.rm = TRUE)
    dyn_range <- signal / noise
    dyn_range_values[marker, tma] <- dyn_range
    rm(dyn_range, intensities, signal, noise)
  }
}

dyn_range_values[] <- apply(dyn_range_values, 2, function(x){
  max_finite <- max(x[is.finite(x)])
  x[!is.finite(x)] <- max_finite
  return(x)
}) # Replace infinite values with max finite

# Plot dynamic range by marker

dyn_range_values$Marker <- markerList

# Set marker order (categorical) as a factor
dyn_range_values$Marker <- factor(dyn_range_values$Marker, 
                                  levels = dyn_range_values$Marker)

for (tma in seq_along(batchInfo$TMA)){
  dyn_range_values <- dyn_range_values %>%
    arrange(desc(!!sym(batchInfo$TMA[[tma]]))) %>%
    mutate(Value = case_when(
      !!sym(batchInfo$TMA[[tma]]) > 30 ~ "Capped",         
      !!sym(batchInfo$TMA[[tma]]) > 10 ~ "Ideal",          
      !!sym(batchInfo$TMA[[tma]]) > 2  ~ "Acceptable",     
      !!sym(batchInfo$TMA[[tma]]) <= 2 ~ "Unacceptable",   
      TRUE ~ NA_character_  # In case of missing or non-numeric values, return NA
    )) %>%
    mutate(!!sym(batchInfo$TMA[[tma]]) := pmin(!!sym(batchInfo$TMA[[tma]]), 30))
  
  p <- ggplot(
    dyn_range_values %>%
      mutate(Marker = fct_reorder(Marker, !!sym(batchInfo$TMA[[tma]]), .desc = TRUE)),  # Order markers
    aes(x = Marker,
        y = !!sym(batchInfo$TMA[[tma]]),
        fill = Value)
  ) +
    geom_bar(stat = "identity") +
    labs(title = paste0("Dynamic Range by Marker (TMA ", batchInfo$TMA[[tma]], ")"),
         y = "Dynamic Range") +
    scale_fill_manual(values = c("Capped" = "darkgreen",
                                 "Ideal" = "#66BB6A",
                                 "Acceptable" = "orange",
                                 "Unacceptable" = "red"),
                      breaks = c("Capped", "Ideal", "Acceptable", "Unacceptable")) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.title = element_blank(),
      legend.text = element_text(size = 8)) +
    geom_hline(yintercept = 2, color = "red", linetype = "longdash", linewidth = 0.5)
  
  output_dir <- paste0("GBM_spatial_analysis/outputs/QC/",
                       batchInfo$TMA[[tma]], "/dynamic_ranges/")
  dir.create(output_dir, recursive=T, showWarnings=F)
  ggsave(paste0(output_dir, batchInfo$TMA[[tma]], "_dynamic_ranges.png"),
         p,
         height = 4, width = 8)
}

# Calculate classifier ratio per marker per TMA:


classifier_ratio_values <- data.frame(matrix(ncol = length(batchInfo$TMA),
                                      nrow = length(markerList)))
rownames(classifier_ratio_values) <- markerList

marker_abundance_values <- classifier_ratio_values

for (tma in seq_along(batchInfo$TMA)){

  print(paste("Processing TMA", batchInfo$TMA[[tma]]))
  
  for (marker in seq_along(markerList)){
    
    print(paste("Calculating for", markerList[[marker]]))

    pos_neg <- long_filtered_cells[[tma]] %>%
      filter(TMA == batchInfo$TMA[[tma]], # filter by TMA and by marker
             Marker == markerList[[marker]]) %>%
      mutate(marker_status = ifelse(grepl(paste0(markerList[[marker]], "($|:)"), Classifier), "Positive", "Negative"))
    
    pos <- pos_neg %>%
      filter(marker_status == "Positive") %>%
      arrange(Intensity)
    
    neg <- pos_neg %>%
      filter(marker_status == "Negative") %>%
      arrange(Intensity)
    
    colnames(classifier_ratio_values)[tma] <- batchInfo$TMA[[tma]]
    intensities <- filtered_cells[[tma]][[markerList[marker]]]
    positive <- mean(pos$Intensity, na.rm = TRUE) 
    negative <- mean(neg$Intensity, na.rm = TRUE)
    classifier_ratio <- positive / negative
    classifier_ratio_values[marker, tma] <- classifier_ratio
    
    colnames(marker_abundance_values)[tma] <- batchInfo$TMA[[tma]]
    
    marker_abundance <- (nrow(pos) / nrow(pos_neg)) * 100
    marker_abundance_values[marker, tma] <- marker_abundance
  
    }
}

rm(p, q, pos, neg, pos_neg, classifier_data, intensities,
   marker_abundance, positive, negative, output_dir, 
   dyn_range_values, long_filtered_cells)


# Lineage assignment ####

# Average the CCR and abundance values across TMAs (doesn't use raw intensity so no normalisation required):

CR_mean <- classifier_ratio_values %>%
  mutate(Mean = rowMeans(., na.rm = TRUE)) %>%
  select(Mean)

A_mean <- marker_abundance_values %>%
  mutate(Mean = rowMeans(., na.rm = TRUE)) %>%
  select(Mean)

# Read in expression profiles

exp <- read.xlsx("GBM_spatial_analysis/outputs/lineage_assignment/thesis_expression_profiles.xlsx", colNames=T)

exp <- exp %>%
  # Keep track of expression profiles
  mutate(Exp = apply(., 1, function(x) paste(na.omit(x[2:ncol(exp)]), collapse = ","))) %>%
  select(Cell.Type, Exp, everything())

CCR <- exp %>%
  # Replace NAs with values from CR_mean or with 1 (so they don't affect the product)
  mutate(across(3:ncol(.), 
                ~ ifelse(!is.na(.), CR_mean$Mean[match(., rownames(CR_mean))], 1))) %>%
  # Calculate the product, replacing NAs with 1
  mutate(Product = reduce(across(3:ncol(.), as.numeric), 
                          `*`, 
                          .init = 1, 
                          .dir = "forward"))

A_min <- exp %>%
  mutate(across(3:ncol(.),
                ~ ifelse(!is.na(.), A_mean$Mean[match(., rownames(CR_mean))], NA))) %>%
  mutate(Min = apply(across(3:ncol(.)), 1, function(x) min(x, na.rm = TRUE)))

options(scipen = 999)
CS <- exp %>%
  mutate(CS = CCR$Product / A_min$Min) %>%
  select(Cell.Type, Exp, CS) %>%
  arrange(-CS)
  
# Assess and CS document to create delay_call document:
# Identify expression profiles to delay
# Identify marker to delay
# Input Delay_marker and its Marker_number
# Input Remaining_markers such that identical sets are inputted identically

write.xlsx(CS, "GBM_spatial_analysis/outputs/lineage_assignment/thesis_CS_initial.xlsx")

delay_call <- read.xlsx("GBM_spatial_analysis/outputs/lineage_assignment/thesis_delay_call.xlsx", colNames=T)

# Delay calls by excluding markers from calculations:

delay_call <- delay_call %>%
  group_by(Remaining_markers) %>%  # Group by Remaining_markers
  mutate(Adj = dense_rank(CS)) %>%  # Assign ranks based on CS values (1 to x)
  mutate(Adj = 1 + Adj / 1000) %>%
  ungroup()  # Remove the grouping after the mutate operation

CCR_delay <- CCR
for (i in seq_len(nrow(delay_call))){
  # For each delay_call, find the right row in CCR
  row_idx <- which(CCR_delay$Cell.Type == delay_call$Cell.Type[i] & 
                     CCR_delay$Exp == delay_call$Exp[i])
  # Find the right marker to adjust
  marker_col <- paste0("Marker", delay_call$Marker_number[i])
  
  if (marker_col %in% colnames(CCR_delay)){
    CCR_delay[row_idx, marker_col] <- delay_call$Adj[i]
  }
}

CCR_delay <- CCR_delay %>%
  mutate(Product = reduce(across(3:(ncol(.)-1), as.numeric), 
                          `*`, 
                          .init = 1, 
                          .dir = "forward"))

A_min_delay <- A_min
  for (i in seq_len(nrow(delay_call))){
    # For each delay_call, find the right row in CCR
    row_idx <- which(A_min_delay$Cell.Type == delay_call$Cell.Type[i] & 
                       CCR_delay$Exp == delay_call$Exp[i])
    # Find the right marker to adjust
    marker_col <- paste0("Marker", delay_call$Marker_number[i])
    
    if (marker_col %in% colnames(A_min_delay)){
      A_min_delay[row_idx, marker_col] <- Inf
    }
  }

A_min_delay <- A_min_delay %>%
  rowwise() %>%
  mutate(Min = min(c_across(3:(ncol(.) - 1)), na.rm = TRUE)) %>%
  mutate(Min = ifelse(is.infinite(Min), Inf, Min)) %>%
  ungroup()

CS_delay <- exp %>%
  mutate(CS = CCR_delay$Product / A_min_delay$Min) %>%
  select(Cell.Type, Exp, CS) %>%
  arrange(-CS)

write.xlsx(CS_delay, "GBM_spatial_analysis/outputs/lineage_assignment/thesis_CS_final.xlsx")

# Assign lineages

cellData <- rbind(filtered_cells[[1]], filtered_cells[[2]])

cellData$Cell.Type <- NA

cellData <- cellData %>%
  select(Cell_ID, Annotation_ID, Core_ID, TMA, Classifier, Cell.Type, everything())

for (i in seq_len(nrow(CS_delay))) {
  exp <- unlist(strsplit(CS_delay$Exp[i], ","))  # Split markers into a vector
  
  # Identify unassigned rows where all markers are present
  match_rows <- is.na(cellData$Cell.Type) & 
    Reduce(`&`, lapply(exp, function(m) grepl(paste0(m, "($|:)"), cellData$Classifier)))
  
  # Assign the corresponding Cell.Type from CS_delay
  cellData$Cell.Type[match_rows] <- CS_delay$Cell.Type[i]
  
}

cellData$Cell.Type[is.na(cellData$Cell.Type)] <- "Non-immune (other)"

rm(A_mean, A_min, A_min_delay, CCR, CCR_delay,
   CR_mean, CS, CS_delay, delay_call, exp, match_rows,
   row_idx, sorted_idx, marker_col)

# Need metadata attributed at cell level for stats:

lts_donors <- patientData %>%
  filter(Cohort == "LTS") %>%
  pull(Patient_ID)

sts_donors <- patientData %>%
  filter(Cohort == "STS") %>%
  pull(Patient_ID)

nav_donors <- patientData %>%
  filter(Cohort == "NAV") %>%
  pull(Patient_ID)

other_donors <- patientData %>%
  filter(Cohort == "Other") %>%
  pull(Patient_ID)

annIDs <- as.data.frame(unique(cellData$Annotation_ID))
colnames(annIDs) <- "Annotation_ID"

annIDs <- cellData %>%
  distinct(Annotation_ID) %>%
  mutate(Cohort = case_when(
    str_detect(Annotation_ID, str_c(lts_donors, collapse = "|")) ~ "LTS",
    str_detect(Annotation_ID, str_c(sts_donors, collapse = "|")) ~ "STS",
    str_detect(Annotation_ID, str_c(nav_donors, collapse = "|")) ~ "NAV",
    str_detect(Annotation_ID, str_c(other_donors, collapse = "|")) ~ "Other",
    TRUE ~ NA_character_
  ))

cellData <- cellData %>%
  mutate(
    Timepoint = case_when(
      grepl("Primary", Annotation_ID) ~ "Primary",
      grepl("Recurrence1", Annotation_ID) ~ "Recurrence1",
      grepl("Recurrence2", Annotation_ID) ~ "Recurrence2",
      TRUE ~ NA
    ),
    Annotation_Type = case_when(
      grepl("Tumour", Annotation_ID) ~ "Tumour",
      grepl("PBZ", Annotation_ID) ~ "PBZ",
      grepl("PPN", Annotation_ID) ~ "PPN",
      grepl("Necrosis", Annotation_ID) ~ "Necrosis",
      TRUE ~ NA
    )
  ) %>%
  left_join(annIDs, by = "Annotation_ID") %>%
  select(Cell_ID, Annotation_ID, Cohort, Timepoint, Annotation_Type,
         Core_ID, TMA, Classifier, Cell.Type,
         X_coord, Y_coord, everything())

# save ####
save.image("GBM_single_cell_analysis/outputs/thesis_spatial_environment.RData")

write.xlsx(batchInfo, "GBM_spatial_analysis/outputs/batchInfo.xlsx")
saveRDS(coreData, "GBM_spatial_analysis/outputs/coreData.RDS")
saveRDS(annotationData, "GBM_spatial_analysis/outputs/annotationData.RDS")
saveRDS(patientData, "GBM_spatial_analysis/outputs/patientData.RDS")
saveRDS(cellData, "GBM_spatial_analysis/outputs/cellData.RDS")
