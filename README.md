# The immune microenvironment of human glioblastoma

## Purpose

This serves as a public repository of code used for all analyses in PhD thesis 'The immune microenvironment of human glioblastoma', submitted to University of Liverpool January 2026.

## Raw data availability

Due to patient confidentiality and institutional governance restrictions, raw and processed data used in this thesis are not publicly available and are therefore not included in this repository.

## Requirements

This project was carried out using the following programs, softwares and algorithms:

FlowJo v10.10.0, BD Biosciences	(https://www.flowjo.com)

Cell Ranger	v8.0.1,	Zheng et al., 2017	(https://10xgenomics.com/software)

Ensembl	v7.1.0, Dyer et al., 2025	(https://www.ensembl.org)

R	v4.4.1, R Core Team	(https://www.R-project.org)

Seurat v5.1.0, Hao et al., 2024	(https://satijalab.org/seurat/)

ab_capture, Meckiff et al., 2020	(https://github.com/vijaybioinfo/ab_capture)

harmony	v1.2.3, Korsunsky et al., 2019	(https://CRAN.R-project.org/package=harmony)

clustree v0.5.1, Zappia & Oshlack, 2018	(https://github.com/lazappi/clustree)

ScType, Ianevski et al., 2022	(https://github.com/IanevskiAleksandr/sc-type)

celldex	v1.14.0, Aran et al., 2019	(10.18129/B9.bioc.celldex)

SingleR	v2.6.0,	Aran et al., 2019	(10.18129/B9.bioc.SingleR)

CelliD	v1.12.0,	Akira et al., 2020	(10.18129/B9.bioc.CelliD)

DESeq2	v1.44.0,	Love et al., 2014	(https://github.com/thelovelab/DESeq2)

clusterProfiler	v4.12.6,	Xu et al., 2024	(https://github.com/YuLab-SMU/clusterProfiler)

slingshot	v2.12.0,	Street et al., 2018	(https://github.com/kstreet13/slingshot)

scRepertoire	v2.0.8,	Borcherding et al., 2020	(https://github.com/BorchLab/scRepertoire)

QuPath	v0.6.0,	Bankhead et al., 2017	(https://qupath.github.io)

InstanSeg,	Goldsborough et al., 2024	(https://github.com/instanseg)

StarDist	v0.7.3,	Schmidt et al., 2018	(https://github.com/stardist)

SpatialExperiment	v1.14.0,	Righelli et al., 2022	(10.18129/B9.bioc.SpatialExperiment)

SPIAT	v1.7.2,	Feng et al., 2023	(https://trigosteam.github.io/SPIAT/)

ClusterR	v1.3.3,	Sanderson & Curtin, 2016	(https://CRAN.R-project.org/package=ClusterR)

RANN	v2.6.2,	Jefferis et al., 2024	(https://CRAN.R-project.org/package=RANN)

lme4	v1.1.37,	Bates et al., 2015	(https://CRAN.R-project.org/package=lme4)

lmerTest	v3.1.3,	Kuznetsova et al., 2017	(https://CRAN.R-project.org/package=lmerTest)

emmeans	v1.11.1,	Lenth, 2025	(https://CRAN.R-project.org/package=emmeans)

survival	v3.8.3,	Therneau & Grambsch, 2000	(https://CRAN.R-project.org/package=survival)

ggplot2 v.3.5.3, Wickham, 2016 (https://ggplot2.tidyverse.org)


## Content

The code is provided to document the analytical approach, processing steps, and figure generation logic, including during exploratory analyses. Analysis scripts are therefore not intended as a fully automated pipeline or to run in the absence of the raw data, and may require manual adaptation to run in a new environment.

Scripts reference local data directories used during the PhD, including the following examples, as well as additional subdirectories specifically referenced within the code:

GBM_single_cell_analysis/

├── data/

└── outputs/

GBM_spatial_analysis/

├── data/

└── outputs/

These directories are expected to exist outside the repository root and contain processed data objects generated during the analysis pipeline.

## Contact

Please e-mail Michael Cearns (michael.cearns@liverpool.ac.uk).
