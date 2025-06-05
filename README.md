# Artemisinin Partial Resistance Mutations in Zanzibar and Tanzania Suggest Regional Spread and African Origins (2023)

**DOI:** [10.1101/2025.03.22.25323829](https://doi.org/10.1101/2025.03.22.25323829)  

## Overview

This repository contains the **processed data** and **analysis scripts** for the study titled  
*“Artemisinin Partial Resistance Mutations in Zanzibar and Tanzania Suggest Regional Spread and African Origins (2023)”*

The study investigates the genetic signatures of **artemisinin partial resistance** in *Plasmodium falciparum* across **Zanzibar** and **mainland Tanzania**, providing evidence of regional connectivity and the emergence of resistance mutations in parasite populations.

## Repository Structure

```plaintext
.
├── data/
│ └── processed/ # Final processed data used in analyses
│ ├── metadata/ # Sample metadata
│
├── analysis/
│ ├── TableS1.R # Script to generate Table S1
│ ├── TableS2.R # Scripts to generate Table S2
│ ├── FigureS2_TableS4.R # Script to generate Figure S2 and Table S4
│ ├── TableS5.R # Script to generate Table S5
├── README.md
└── LICENSE
```

## Clone the Repository

To clone the repository, run:

```bash
git clone https://github.com/sconnelly007/ZIM_ART_R.git
cd ZIM_ART_R
```

## Set Up the Environment with renv

```R
# Install renv if not already installed
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}
# Initialize renv
renv::init()
# Restore the environment
renv::restore()
```

## Citation
Connelly S.V., et al. (2025) Artemisinin Partial Resistance Mutations in Zanzibar and Tanzania Suggest Regional Spread and African Origins (2023). medRxiv. https://doi.org/10.1101/2025.03.22.25323829