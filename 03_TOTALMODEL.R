## Total density, size structure, life stage #############################################


## Libraries #############################################
library(glmmTMB)
library(dplyr)
library(stringr)
library(tibble)
library(ggeffects)
library(ggplot2)

## 1. Add time variables and clean size factor ################

origin_date <- min(fish_long$Date, na.rm = TRUE)

fish_long <- fish_long %>%
  mutate(
    Size_Class_f = factor(Size_Class, ordered = FALSE),
    date_num     = as.numeric(Date - origin_date),
    date_s       = as.numeric(scale(date_num, center = TRUE, scale = TRUE))
  )

## 2. Main Model: size structure by reef type ################

# “Are size distributions different on artificial vs natural reefs,
# within pairs, across all months pooled?”

m1 <- glmmTMB(
  Count ~ Type * Size_Class_f +
    Pair +
    Count.Type +
    offset(log(Inclusion_m)) +
    (1 | Site) +
    (1 | survey_id),
  ziformula = ~ Type + Size_Class_f,   # zero part by reef type and size class
  family    = nbinom2,
  data      = fish_long
)

# Save model and summary
saveRDS(m1, file.path(fits_dir, paste0("m1_size_structure_", analysis_date, ".rds")))
capture.output(
  summary(m1),
  file = file.path(stats_dir, paste0("m1_size_structure_summary_", analysis_date, ".txt"))
)

## 3. Time model: do size distributions change over time? ################

m_time <- glmmTMB(
  Count ~ Type * Size_Class_f +           # size structure by reef type
    Type * date_s +                       # time trend differs by reef type
    Size_Class_f * date_s +               # time trend differs by size class
    Pair +
    Count.Type +
    offset(log(Inclusion_m)) +
    (1 | Site) +
    (1 | survey_id),
  ziformula = ~ Type + Size_Class_f,      # simpler zero model, no time in zero part
  family    = nbinom2,
  data      = fish_long
)

saveRDS(m_time, file.path(fits_dir, paste0("m_time_size_time_", analysis_date, ".rds")))
capture.output(
  summary(m_time),
  file = file.path(stats_dir, paste0("m_time_size_time_summary_", analysis_date, ".txt"))
)

# Optional: site random effects check (interactive, not saved)
# ranef(m_time)$cond$Site


## 4. SUPPLEMENTARY: deterministic life-stage model (with mixed) ################

# This section is kept for robustness checks / supplement only

# Species table (used by maturity lookup script)
spp_tbl <- fish_size %>%
  distinct(Sci_Name) %>%
  mutate(
    Sci_Name = as.character(Sci_Name),
    genus    = word(Sci_Name, 1),
    species  = word(Sci_Name, 2)
  )

# Maturity lookup
source("~/Documents/1_GLOBAL REEF/0_PROJECTS/AR_Producer_Attractor/AR_Producer/02.1_matlookup.R")

size_bins <- tibble::tribble(
  ~Size_Class, ~lower, ~upper,
  "0-1",       0,       1,
  "1-2",       1,       2,
  "2-5",       2,       5,
  "5-10",      5,      10,
  "10-15",    10,      15,
  "15-20",    15,      20,
  "20-50",    20,      50,
  "50-100",   50,     100,
  "100+",    100,    Inf
)

fish_long_life <- fish_long %>%
  left_join(maturity_lookup, by = "Sci_Name") %>%
  left_join(size_bins,       by = "Size_Class") %>%
  mutate(
    life_stage = case_when(
      is.na(Lmat_cm) ~ NA_character_,
      upper   < Lmat_cm ~ "juvenile",     # whole bin below Lmat
      lower   >= Lmat_cm ~ "adult",       # whole bin above Lmat
      TRUE             ~ "mixed"          # bin crosses Lmat
    ),
    life_stage = factor(
      life_stage,
      levels = c("juvenile", "mixed", "adult")
    )
  )

# Deterministic life stage model (main effect of life_stage)
m_stage <- glmmTMB(
  Count ~ Type * life_stage +
    date_num + Pair +
    Count.Type +
    offset(log(Inclusion_m)) +
    (1 | Site) +
    (1 | survey_id),
  family = nbinom2,
  data   = fish_long_life
)

saveRDS(m_stage, file.path(fits_dir, paste0("m_stage_deterministic_", analysis_date, ".rds")))
capture.output(
  summary(m_stage),
  file = file.path(stats_dir, paste0("m_stage_deterministic_summary_", analysis_date, ".txt"))
)

# With Pair interaction (supplementary final deterministic model)
m_stage_fx <- glmmTMB(
  Count ~ Type * life_stage * Pair +
    date_num +
    Count.Type +
    offset(log(Inclusion_m)) +
    (1 | Site) +
    (1 | survey_id),
  family = nbinom2,
  data   = fish_long_life
)

saveRDS(m_stage_fx, file.path(fits_dir, paste0("m_stage_deterministic_pair_", analysis_date, ".rds")))
capture.output(
  summary(m_stage_fx),
  file = file.path(stats_dir, paste0("m_stage_deterministic_pair_summary_", analysis_date, ".txt"))
)


## 5. MAIN: probabilistic life stage model (juvenile vs adult) ################

# Re-use size_bins from above and fish_long with date_num already added

fish_long_prob <- fish_long %>%
  left_join(maturity_lookup, by = "Sci_Name") %>%   # has Lmat_cm
  left_join(size_bins,       by = "Size_Class") %>%
  mutate(
    # probability that a random fish in the bin is juvenile
    p_juv = case_when(
      is.na(Lmat_cm)   ~ NA_real_,
      upper <= Lmat_cm ~ 1,  # whole bin below Lmat
      lower >= Lmat_cm ~ 0,  # whole bin above Lmat
      TRUE             ~ (Lmat_cm - lower) / (upper - lower)
    ),
    p_adult = if_else(is.na(p_juv), NA_real_, 1 - p_juv)
  )

# Turn probabilities into juvenile/adult counts
fish_long_life_prob <- fish_long_prob %>%
  filter(!is.na(p_juv), !is.na(Count)) %>%
  mutate(
    n_juv   = round(Count * p_juv),
    n_adult = Count - n_juv
  ) %>%
  tidyr::pivot_longer(
    cols      = c(n_juv, n_adult),
    names_to  = "life_stage",
    values_to = "stage_Count"
  ) %>%
  filter(stage_Count > 0) %>%   # fix: filter using stage_Count, not Count
  mutate(
    life_stage = recode(
      life_stage,
      "n_juv"   = "juvenile",
      "n_adult" = "adult"
    ),
    life_stage = factor(life_stage, levels = c("juvenile", "adult"))
  )

# Probabilistic life stage model (MAIN inference model)
m_stage_prob <- glmmTMB(
  stage_Count ~ Type * life_stage * Pair +
    date_num +
    Count.Type +
    offset(log(Inclusion_m)) +
    (1 | Site) +
    (1 | survey_id),
  family = nbinom2,
  data   = fish_long_life_prob
)

saveRDS(m_stage_prob, file.path(fits_dir, paste0("m_stage_prob_", analysis_date, ".rds")))
capture.output(
  summary(m_stage_prob),
  file = file.path(stats_dir, paste0("m_stage_prob_summary_", analysis_date, ".txt"))
)


## 6. PREDICTIONS AND PLOTS ###############################

### Figure 1. Total fish density by reef type and pair ################

pred_tot <- ggpredict(
  m1,
  terms     = c("Type", "Pair"),
  condition = c(Inclusion_m = 10000)    # per hectare
) %>% as.data.frame()

p_tot <- ggplot(pred_tot, aes(x = x, y = predicted, color = group)) +
  geom_point(size = 3, position = position_dodge(width = 0.4)) +
  geom_errorbar(
    aes(ymin = conf.low, ymax = conf.high),
    width    = 0.15,
    position = position_dodge(width = 0.4)
  ) +
  scale_color_manual(values = reef_cols) +
  theme_clean +
  labs(
    x     = "Reef type",
    y     = "Total fish density (ind. ha⁻¹)",
    color = "Reef type"
  )

ggsave(
  filename = file.path(plots_dir, paste0("Fig1_total_density_", analysis_date, ".png")),
  plot     = p_tot,
  width    = 6.5,
  height   = 4.5,
  dpi      = 300
)


### Figure 2. Size structure by reef type ################

pred_size <- ggpredict(
  m1,
  terms     = c("Size_Class_f", "Type"),
  condition = c(Inclusion_m = 10000)
) %>% as.data.frame()

p_size <- ggplot(pred_size, aes(x = x, y = predicted, color = group)) +
  geom_line(aes(group = group), linewidth = 0.9) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high, fill = group),
    alpha = 0.15,
    color = NA
  ) +
  geom_point(size = 2.4) +
  scale_color_manual(values = reef_cols) +
  scale_fill_manual(values = reef_cols) +
  theme_clean +
  labs(
    x     = "Size class (cm)",
    y     = "Predicted density (ind. ha⁻¹)",
    color = "Reef type",
    fill  = "Reef type"
  )

ggsave(
  filename = file.path(plots_dir, paste0("Fig2_size_structure_", analysis_date, ".png")),
  plot     = p_size,
  width    = 6.5,
  height   = 4.5,
  dpi      = 300
)


### Figure 3. Time trend in total density by reef type ################

pred_time <- ggpredict(
  m_time,
  terms     = c("date_s [minmax]", "Type"),
  condition = c(Inclusion_m = 10000)
) %>% as.data.frame()

pred_time$x_date <- origin_date + pred_time$x

p_time <- ggplot(pred_time, aes(x = x_date, y = predicted, color = group)) +
  geom_line(linewidth = 0.9) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high, fill = group),
    alpha = 0.15,
    color = NA
  ) +
  scale_color_manual(values = reef_cols) +
  scale_fill_manual(values = reef_cols) +
  theme_clean +
  labs(
    x     = "Date",
    y     = "Predicted density (ind. ha⁻¹)",
    color = "Reef type",
    fill  = "Reef type"
  )

ggsave(
  filename = file.path(plots_dir, paste0("Fig3_time_trend_", analysis_date, ".png")),
  plot     = p_time,
  width    = 6.5,
  height   = 4.5,
  dpi      = 300
)


### SUPP Figure: deterministic life stage (juvenile / mixed / adult) ################

pred_stage_det <- ggpredict(
  m_stage_fx,
  terms     = c("life_stage", "Type", "Pair"),
  condition = c(Inclusion_m = 100)
) %>% as.data.frame()

pred_stage_det$x <- factor(
  pred_stage_det$x,
  levels = c("juvenile", "mixed", "adult"),
  ordered = TRUE
)

p_pair_det <- ggplot(
  pred_stage_det,
  aes(x = x, y = predicted, color = group, fill = group, group = group)
) +
  geom_errorbar(
    aes(ymin = conf.low, ymax = conf.high),
    width    = 0.15,
    position = position_dodge(width = 0.5)
  ) +
  geom_point(
    size     = 2.4,
    position = position_dodge(width = 0.5)
  ) +
  facet_wrap(~ facet, nrow = 1) +
  scale_color_manual(values = reef_cols) +
  scale_fill_manual(values = reef_cols) +
  theme_clean +
  labs(
    x     = "Life stage",
    y     = "Predicted count per 100 m²",
    color = "Reef type",
    fill  = "Reef type"
  )

ggsave(
  filename = file.path(plots_dir, paste0("Supp_Fig_lifestage_deterministic_", analysis_date, ".png")),
  plot     = p_pair_det,
  width    = 7,
  height   = 4.5,
  dpi      = 300
)


### MAIN Figure: probabilistic life stage (juvenile vs adult) ################

pred_stage_prob <- ggpredict(
  m_stage_prob,
  terms     = c("life_stage", "Type", "Pair"),
  condition = c(Inclusion_m = 100)
) %>% as.data.frame()

p_pair_prob <- ggplot(
  pred_stage_prob,
  aes(x = x, y = predicted, color = group, fill = group, group = group)
) +
  geom_errorbar(
    aes(ymin = conf.low, ymax = conf.high),
    width    = 0.15,
    position = position_dodge(width = 0.5)
  ) +
  geom_point(
    size     = 2.4,
    position = position_dodge(width = 0.5)
  ) +
  facet_wrap(~ facet, nrow = 1) +
  scale_color_manual(values = reef_cols) +
  scale_fill_manual(values = reef_cols) +
  theme_clean +
  labs(
    x     = "Life stage",
    y     = "Predicted count per 100 m²",
    color = "Reef type",
    fill  = "Reef type"
  )

ggsave(
  filename = file.path(plots_dir, paste0("Fig4_lifestage_prob_", analysis_date, ".png")),
  plot     = p_pair_prob,
  width    = 7,
  height   = 4.5,
  dpi      = 300
)
