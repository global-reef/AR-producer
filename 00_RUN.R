
#########  Set Analysis Date & Create Output Folder ########################
library(ggplot2)

# Enter the date for this analysis
analysis_date <- "2025.11.13"  # update each run !!
# File path (automatically join folder + date + filename)
file_path <- paste0(
  "~/Documents/1_GLOBAL REEF/0_PROJECTS/AR_Producer_Attractor/AR_Producer/DATA/",
  analysis_date,
  "_FishSize_MASTER.csv"
)
raw_fish <- read.csv(file_path, stringsAsFactors=TRUE, strip.white=TRUE) 


# Create a folder named with the date inside the working directory
output_dir <- file.path(getwd(), paste0("Analysis_", analysis_date))
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

fits_dir  <- file.path(output_dir, "fits")
plots_dir <- file.path(output_dir, "plots")
stats_dir <- file.path(output_dir, "stats")
summ_dir  <- file.path(output_dir, "summaries")
dir.create(fits_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(stats_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(summ_dir,  showWarnings = FALSE, recursive = TRUE)

### custom theme and colour palettes ###############################################################
theme_clean <- theme_minimal(base_family = "Arial") +
  theme(
    legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.grid = element_blank()
  )

# colour palettes
fg_cols <- c(
  "Grazer" = "#66c2a4",
  "Invertivore" = "#41b6c4",
  "Mesopredator" = "#2c7fb8",
  "HTLP" = "#253494"
) 
reef_cols <- c("Natural" = "#66BFA6", "Artificial" = "#007A87")

