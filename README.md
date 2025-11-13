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
- Explore responses by **functional group** (e.g., grazers, invertivores, mesopredators, higher-trophic predators).  
- Provide management-relevant guidance for future AR deployments under Thailand’s DMCR program.

---

## 📁 Repository structure   ### TBA
```
DATA/                   # cleaned and raw data (fish counts, site metadata)
00.RUN.R                # master script to reproduce full analysis
01.CLEAN.R              # data cleaning and preprocessing
02.EXPLORE.R            # exploratory summaries and visualizations
03.DIVERSITY.R          # diversity and assemblage metrics
04.FUNGROUPS.R          # functional group models
05.PAIRED.R             # paired comparisons: wreck vs reef
06.NMDS_BRAYCURT.R      # multivariate ordination (NMDS)
scratch.R               # testing / notes
AR_Producer.Rproj       # R project file
```

Outputs (figures, tables, RDS objects) are generated within each script and can be directed to an `outputs/` folder.

---

## ⚙️ Getting started  
1. Clone the repository:
   ```bash
   git clone https://github.com/global-reef/AR-Producer.git
   ```
2. Open `AR_Producer.Rproj` in RStudio.  
3. Run `01.CLEAN.R` to assemble cleaned data.  
4. Execute subsequent scripts (`02.*` → `06.*`) in order.  
5. All figures and summary tables will be created automatically.

---

## 📦 Dependencies  
Required R packages:
- **tidyverse** (dplyr, ggplot2, tidyr)
- **vegan** (for NMDS and Bray-Curtis dissimilarities)
- **emmeans**, **glmmTMB**, **broom.mixed**, **performance** (for model evaluation)

Each script lists additional dependencies at the top.

---

## 📊 Outputs  
- Fish abundance summaries by site, date, and functional group  
- Diversity and NMDS ordination plots  
- Paired comparisons between wreck and reef sites  
- Figures for reporting to DMCR and Publication

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
