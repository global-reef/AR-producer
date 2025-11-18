######### 02_EXPLORE.R #######################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(lubridate)
  library(glmmTMB)
  library(lme4)
  library(DHARMa)
  library(readr)
  library(RColorBrewer)
})

## Inputs ###################################################################

stopifnot(exists("output_dir"))

eda_dir <- file.path(output_dir, "EDA")
dir.create(eda_dir, showWarnings = FALSE, recursive = TRUE)

fish_size <- readRDS(file.path(output_dir, "fish_size_cleaned.rds"))

## Helper: base plotting theme ##############################################

base_theme <- if (exists("theme_clean")) theme_clean else theme_minimal()

## Zuur-style EDA function ##################################################

zuur_eda <- function(dat,
                     output_dir = NULL,
                     save_plots = TRUE,
                     save_csv = TRUE,
                     show_plots = FALSE,
                     species_min_n = 50) {
  
  d <- dat %>%
    mutate(
      Date   = if (inherits(Date, "Date")) Date else as.Date(as.character(Date)),
      day_id = interaction(site, Date, drop = TRUE)
    )
  
  ## Effort: unique surveys per site ########################################
  surv <- d %>%
    dplyr::distinct(site, type, Date, survey_id)
  
  effort_by_site <- surv %>%
    count(site, name = "n_surveys")
  
  p_effort_hist <- ggplot(effort_by_site, aes(x = n_surveys)) +
    geom_histogram(bins = 30, boundary = 0) +
    labs(
      title = "Surveys per site",
      x = "Number of surveys",
      y = "Frequency"
    ) +
    base_theme
  
  ## Species-level zero inflation ###########################################
  species_zero <- d %>%
    group_by(Species) %>%
    summarise(
      n              = n(),
      detections     = sum(Count > 0, na.rm = TRUE),
      zero_rate      = mean(Count == 0, na.rm = TRUE),
      mean_count     = mean(Count, na.rm = TRUE),
      mean_count_pos = ifelse(detections > 0,
                              mean(Count[Count > 0], na.rm = TRUE),
                              NA_real_),
      .groups        = "drop"
    ) %>%
    arrange(desc(zero_rate), desc(n))
  
  species_flags <- species_zero %>%
    filter(n >= species_min_n, zero_rate >= 0.9)
  
  p_species_detect_hist <- ggplot(species_zero, aes(x = detections)) +
    geom_histogram(bins = 40, boundary = 0) +
    labs(
      title = "Detections per species",
      x = "Number of nonzero observations",
      y = "Frequency"
    ) +
    base_theme
  
  ## Per survey richness ####################################################
  richness <- d %>%
    group_by(survey_id) %>%
    summarise(
      richness = n_distinct(Species[Count > 0]),
      .groups  = "drop"
    )
  
  p_richness <- ggplot(richness, aes(x = richness)) +
    geom_histogram(bins = 40, boundary = 0) +
    labs(
      title = "Species richness per survey",
      x = "Number of species (Count > 0)",
      y = "Frequency"
    ) +
    base_theme
  
  ## Core count diagnostics #################################################
  zero_rate_all <- mean(d$Count == 0, na.rm = TRUE)
  
  m_pois <- glm(Count ~ 1, data = d, family = poisson())
  disp   <- sum(residuals(m_pois, type = "pearson")^2) / df.residual(m_pois)
  
  d$log1pC <- log1p(d$Count)
  
  m_day  <- lme4::lmer(log1pC ~ 1 + (1 | day_id), data = d)
  vc     <- lme4::VarCorr(m_day)
  sd_day <- as.data.frame(vc)$sdcor[1]
  sd_res <- attr(vc, "sc")
  icc_day <- sd_day^2 / (sd_day^2 + sd_res^2)
  
  # simple NB model for DHARMa residuals
  m_nb <- glmmTMB::glmmTMB(
    Count ~ 1 + (1 | site) + (1 | day_id),
    data   = d,
    family = glmmTMB::nbinom2()
  )
  sim  <- DHARMa::simulateResiduals(m_nb, n = 500)
  
  p_hist <- ggplot(d, aes(x = Count)) +
    geom_histogram(bins = 50, boundary = 0) +
    coord_cartesian(xlim = c(0, stats::quantile(d$Count, 0.99, na.rm = TRUE))) +
    labs(
      title = "Count distribution",
      x = "Count",
      y = "Frequency"
    ) +
    base_theme
  
  by_day <- d %>%
    group_by(day_id) %>%
    summarise(
      mean_count = mean(Count, na.rm = TRUE),
      var_count  = var(Count, na.rm = TRUE),
      .groups    = "drop"
    )
  
  p_mv <- ggplot(by_day, aes(x = mean_count, y = var_count)) +
    geom_point(alpha = 0.6) +
    geom_smooth(se = FALSE, method = "lm") +
    labs(
      title = "Mean vs variance by day",
      x = "Mean count per day",
      y = "Variance in counts per day"
    ) +
    base_theme
  
  ## Extra covariate EDA ####################################################
  d <- d %>%
    mutate(hour = lubridate::hour(hms(Time)))
  
  # Depth
  p_depth <- ggplot(d, aes(Depth_m, Count)) +
    geom_boxplot(outlier.size = 0.8) +
    labs(title = "Count by depth", x = "Depth (m)", y = "Count") +
    base_theme
  
  # Visibility
  p_vis <- ggplot(d, aes(Visibility_m, Count)) +
    geom_boxplot(outlier.size = 0.8) +
    labs(title = "Count by visibility", x = "Visibility (m)", y = "Count") +
    base_theme
  
  # Current
  p_current <- ggplot(d, aes(Current, Count)) +
    geom_boxplot(outlier.size = 0.8) +
    labs(title = "Count by current", x = "Current", y = "Count") +
    base_theme
  
  # Boats
  p_boats <- d %>%
    mutate(Boats = factor(Boats)) %>%
    ggplot(aes(Boats, Count)) +
    geom_boxplot(outlier.size = 0.8) +
    labs(title = "Count by number of boats", x = "Boats", y = "Count") +
    base_theme
  
  # Hour of day
  p_hour <- ggplot(d, aes(hour, Count)) +
    geom_point(alpha = 0.15) +
    geom_smooth(se = FALSE) +
    labs(title = "Count vs hour", x = "Hour", y = "Count") +
    base_theme
  
  # Month-Year
  p_month <- ggplot(d, aes(Month_Year, Count)) +
    geom_boxplot(outlier.size = 0.8) +
    labs(title = "Count by month", x = "Month-Year", y = "Count") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    base_theme
  
  ## Save plots and CSVs ####################################################
  if (!is.null(output_dir) && (save_plots || save_csv)) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  }
  
  if (!is.null(output_dir) && save_plots) {
    ggsave(file.path(output_dir, "EDA_hist_counts.png"),
           p_hist, width = 6, height = 4, dpi = 300)
    ggsave(file.path(output_dir, "EDA_mean_variance_by_day.png"),
           p_mv, width = 6, height = 4, dpi = 300)
    ggsave(file.path(output_dir, "EDA_surveys_per_site.png"),
           p_effort_hist, width = 6, height = 4, dpi = 300)
    ggsave(file.path(output_dir, "EDA_species_detections_hist.png"),
           p_species_detect_hist, width = 6, height = 4, dpi = 300)
    ggsave(file.path(output_dir, "EDA_richness_per_survey.png"),
           p_richness, width = 6, height = 4, dpi = 300)
    ggsave(file.path(output_dir, "EDA_depth.png"),      p_depth,  width = 6, height = 4, dpi = 300)
    ggsave(file.path(output_dir, "EDA_visibility.png"), p_vis,    width = 6, height = 4, dpi = 300)
    ggsave(file.path(output_dir, "EDA_current.png"),    p_current,width = 6, height = 4, dpi = 300)
    ggsave(file.path(output_dir, "EDA_boats.png"),      p_boats,  width = 6, height = 4, dpi = 300)
    ggsave(file.path(output_dir, "EDA_hour.png"),       p_hour,   width = 6, height = 4, dpi = 300)
    ggsave(file.path(output_dir, "EDA_month.png"),      p_month,  width = 8, height = 5, dpi = 300)
    
    png(file.path(output_dir, "EDA_DHARMa_nb.png"),
        width = 1200, height = 900, res = 150)
    plot(sim)
    dev.off()
  }
  
  if (!is.null(output_dir) && save_csv) {
    readr::write_csv(effort_by_site,
                     file.path(output_dir, "effort_by_site.csv"))
    readr::write_csv(species_zero,
                     file.path(output_dir, "species_zero_summary.csv"))
    readr::write_csv(species_flags,
                     file.path(output_dir, "species_zero_flags.csv"))
    readr::write_csv(richness,
                     file.path(output_dir, "survey_richness.csv"))
  }
  
  if (show_plots) {
    print(p_effort_hist)
    print(p_species_detect_hist)
    print(p_richness)
    print(p_hist)
    print(p_mv)
  }
  
  list(
    metrics = tibble::tibble(
      n_rows             = nrow(d),
      n_surveys          = nrow(surv),
      n_sites            = dplyr::n_distinct(surv$site),
      n_days             = dplyr::n_distinct(surv$Date),
      zero_rate          = zero_rate_all,
      poisson_dispersion = disp,
      icc_day            = icc_day
    ),
    effort  = list(by_site = effort_by_site),
    species = list(summary = species_zero, flags = species_flags),
    survey  = list(richness = richness)
  )
}

## Run EDA ###################################################################

# fish_long should already exist from 01_CLEAN.R
eda <- zuur_eda(
  dat        = fish_long,
  output_dir = eda_dir,
  save_plots = TRUE,
  save_csv   = TRUE,
  show_plots = FALSE,
  species_min_n = 50
)

print(eda$metrics)
eda$species$flags

## Basic size structure summaries ###########################################

# Surveys per site
survey_counts <- fish_long %>%
  dplyr::distinct(site, Date, survey_id) %>%
  count(site, name = "n_surveys")

readr::write_csv(
  survey_counts,
  file.path(eda_dir, "survey_counts_by_site.csv")
)

# Total abundance by reef type
total_by_type <- fish_long %>%
  group_by(type) %>%
  summarise(total_abundance = sum(Count, na.rm = TRUE), .groups = "drop")

readr::write_csv(
  total_by_type,
  file.path(eda_dir, "total_abundance_by_type.csv")
)

# Total abundance by pair and site
total_by_pair_site <- fish_long %>%
  group_by(pair, site) %>%
  summarise(total_abundance = sum(Count, na.rm = TRUE), .groups = "drop")

readr::write_csv(
  total_by_pair_site,
  file.path(eda_dir, "total_abundance_by_pair_site.csv")
)

## Size class distribution plots ############################################

# Normalised by number of surveys per site
fish_long_norm <- fish_long %>%
  left_join(survey_counts, by = "site") %>%
  mutate(Count_perSurvey = Count / n_surveys)

size_site_summary <- fish_long_norm %>%
  group_by(site, Size_Class) %>%
  summarise(
    Mean_Count = sum(Count_perSurvey, na.rm = TRUE),
    .groups    = "drop"
  )

p_size_by_site <- ggplot(size_site_summary,
                         aes(x = Size_Class, y = Mean_Count, fill = Size_Class)) +
  geom_col(position = "dodge", na.rm = TRUE) +
  facet_wrap(~ site, nrow = 3, ncol = 2, scales = "free_y") +
  scale_fill_brewer(palette = "GnBu") +
  labs(
    title = "Mean fish count per size class (normalised by survey)",
    x     = "Size class (cm)",
    y     = "Mean count per survey",
    fill  = "Size class"
  ) +
  base_theme +
  theme(
    axis.text.x    = element_text(angle = 45, hjust = 1),
    strip.text     = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(
  file.path(eda_dir, "size_class_mean_by_site.png"),
  p_size_by_site,
  width = 9, height = 7, dpi = 300
)

# Mean size classes (not by site)
survey_counts2 <- fish_long %>%
  distinct(Site, Date) %>%
  count(Site, name = "n_surveys")

size_class_summary <- fish_long %>%
  left_join(survey_counts2, by = c("Site" = "Site")) %>%
  mutate(Count_perSurvey = Count / n_surveys) %>%
  group_by(Size_Class) %>%
  summarise(
    Mean_Count = sum(Count_perSurvey, na.rm = TRUE),
    .groups    = "drop"
  )

p_size_classes <- ggplot(size_class_summary,
                         aes(x = Size_Class, y = Mean_Count, fill = Size_Class)) +
  geom_col(position = "dodge") +
  scale_fill_brewer(palette = "GnBu") +
  labs(
    title = "Mean fish count per size class (normalised across all sites)",
    x     = "Size class (cm)",
    y     = "Mean count per survey",
    fill  = "Size class"
  ) +
  base_theme +
  theme(
    axis.text.x    = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

ggsave(
  file.path(eda_dir, "size_class_mean.png"),
  p_size_classes,
  width = 9, height = 7, dpi = 300
)

## Histograms of counts #####################################################

ggplot(fish_long, aes(x = Count)) +
  geom_histogram(binwidth = 1, color = "white", fill = "grey40") +
  labs(
    title = "Histogram of counts (including zeros)",
    x     = "Count",
    y     = "Frequency"
  ) +
  coord_cartesian(xlim = c(0, 10)) +
  theme_minimal()

ggplot(filter(fish_long, Count > 0), aes(x = Count)) +
  geom_histogram(binwidth = 1, color = "white", fill = "steelblue") +
  labs(
    title = "Histogram of counts (positive values only)",
    x     = "Count",
    y     = "Frequency"
  ) +
  coord_cartesian(xlim = c(0, 50)) +
  theme_minimal()

## Proportion of fish in each size class by reef type #######################

size_prop_type <- fish_long %>%
  group_by(type, Size_Class) %>%
  summarise(Total = sum(Count, na.rm = TRUE), .groups = "drop") %>%
  group_by(type) %>%
  mutate(Proportion = Total / sum(Total))

p_size_prop_type <- ggplot(size_prop_type,
                           aes(x = Size_Class, y = Proportion, fill = type)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("Natural" = "#66BFA6",
                               "Artificial" = "#007A87")) +
  labs(
    title = "Proportion of fish in each size class by reef type",
    x     = "Size class (cm)",
    y     = "Proportion of total count",
    fill  = "Reef type"
  ) +
  base_theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  file.path(eda_dir, "size_class_proportion_by_type.png"),
  p_size_prop_type,
  width = 8, height = 5, dpi = 300
)

## Extra diagnostics (as comments in your original) left out here for brevity ##
