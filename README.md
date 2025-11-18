# 🪸 AR-Producer 🐠

Using shipwreck fish assemblage data to test whether newly deployed artificial reefs act as **producers** of fish from nearby natural reefs in the Gulf of Thailand.  
This repository covers only the *producer* half of the broader **Attractor vs Producer** project by [Global Reef](https://global-reef.com).

---

## 🐠 Overview  
Artificial reefs can increase local fish abundance either by producing new biomass or by attracting existing fish from surrounding reefs.
This repository focuses on the production pathway using repeated fish surveys from wreck deployments near Koh Tao, Thailand.
Analysis uses size-structured data to test whether artificial reefs support growth, recruitment, or long-term increases in biomass over time.


The complementary “Attractor” half of this project is hosted separately in [`AR-attractor`](https://github.com/global-reef/AR-attractor) and focuses on long-term productivity and recruitment processes.

---

## 🎯 Objectives  
- Quantify fish assemblage structure and abundance on newly deployed wrecks vs nearby natural reefs.  
- Test whether observed assemblages indicate attraction of existing reef fish.  
- Explore responses by **functional group** (e.g., grazers, invertivores, mesopredators, higher-trophic predators) and **species-specific**.  
- Provide management-relevant guidance for future AR deployments under Thailand’s DMCR program.

---

## 📁 Repository structure   
```
DATA/                   # cleaned and raw data (fish counts, site metadata)
00_RUN.R                # master script to reproduce full analysis
01_CLEAN.R              # data cleaning and preprocessing
02_EXPLORE.R            # exploratory summaries and visualizations
03_TOTALMODEL.R         # m1 total density models (Type × Size_Class_f)
04_JUVPROPS.R           # juvenile proportion by reef type and pair
05_FUNGROUPS.R          # m_fg_stage and related contrasts
06_SPP_MODELS.R         # species-level GLMMs and diagnostics (model_diag)
07_FINAL_PRODUCTION.R   # unified summary of AR − NR contrasts
999_scratch.R           # testing / notes


outputs/
  figs/                         # saved figures
  tables/                       # csv tables for manuscript
  eda/                          # exploratory analysis following zuur et al (2010)
  rds/                          # model fits, predictions, posterior draws
scratch/                        # temporary work
AR_Producer.Rproj               # RStudio project
```

Outputs (figures, tables, RDS objects) are generated within each script and can be directed to an `outputs/` folder.

---

## ⚙️ Getting started  
1. Clone the repository:
   ```bash
   git clone https://github.com/global-reef/AR-Producer.git
   ```
2. Open `AR_Producer.Rproj` in RStudio.  
3. Run `00_RUN.R` to assemble cleaned data.  
4. Execute subsequent scripts (`01_*` → `07_*`) in order.  
5. All figures and summary tables will be created automatically.

---

## 📦 Dependencies  
Required R packages:
- tidyverse
- glmmTMB
- emmeans
- ggeffects
- broom.mixed
- performance
- vegan (optional for multivariate summaries)

Each script lists additional dependencies at the top.

---

## 📊 Outputs  
- Total abundance contrasts (AR − NR)
- Juvenile vs adult production signatures
- Functional group contrasts across life stages
- Species-level models and diagnostics
- Unified numerical summary from 07_FINAL_PRODUCTION_SUMM.R for manuscript text
- Figures suitable for publication

---

## 🔍 Relation to “Attractor” work  
This repository isolates the *producer* question:  
> Do new wrecks  wrecks generate *new* fish biomass through recruitment and growth?

The **attractor** repository will assess whether new wrecks pull fish from existing reef systems, altering local assemblage balance. 

 [`AR-attractor`](https://github.com/global-reef/AR-attractor)


## Notes
-----

- Survey data are collected by Global Reef researchers based in Koh Tao, Thailand.
- Fieldwork and data processing are ongoing; results may be updated as more surveys are completed.
- This project contributes to understanding the ecological role of mid-shelf reef features in regional conservation planning.

## License
-------

This project is private and not licensed for redistribution. For collaboration inquiries, please contact scarlett@global-reef.com.
---

**Affiliation:** [Global Reef](https://global-reef.com), Koh Tao, Thailand  
