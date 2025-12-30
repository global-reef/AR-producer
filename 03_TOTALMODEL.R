### Total density, size structure, life stage #################################

### Libraries ##################################################################

suppressPackageStartupMessages({
  library(glmmTMB)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(tidyr)
  library(ggeffects)
  library(ggplot2)
  library(emmeans)
})

### 1) Add time variables and clean size factor ################################

origin_date <- min(fish_long$Date, na.rm = TRUE)

fish_long <- fish_long %>%
  mutate(
    Size_Class_f = factor(Size_Class, ordered = FALSE),
    date_num     = as.numeric(Date - origin_date),
    date_s       = as.numeric(scale(date_num, center = TRUE, scale = TRUE))
  )

### 2) Life-stage preprocessing (deterministic + probabilistic) #################

spp_tbl <- fish_size %>%
  distinct(Sci_Name) %>%
  mutate(
    Sci_Name = as.character(Sci_Name),
    genus    = stringr::word(Sci_Name, 1),
    species  = stringr::word(Sci_Name, 2)
  )

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
      is.na(Lmat_cm)       ~ NA_character_,
      upper   <  Lmat_cm   ~ "juvenile",
      lower   >= Lmat_cm   ~ "adult",
      TRUE                 ~ "mixed"
    ),
    life_stage = factor(life_stage, levels = c("juvenile", "mixed", "adult"))
  )

fish_long_prob <- fish_long %>%
  left_join(maturity_lookup, by = "Sci_Name") %>%
  left_join(size_bins,       by = "Size_Class") %>%
  mutate(
    p_juv = case_when(
      is.na(Lmat_cm)    ~ NA_real_,
      upper <= Lmat_cm  ~ 1,
      lower >= Lmat_cm  ~ 0,
      TRUE              ~ (Lmat_cm - lower) / (upper - lower)
    ),
    p_adult = if_else(is.na(p_juv), NA_real_, 1 - p_juv)
  )

fish_long_life_prob <- fish_long_prob %>%
  filter(!is.na(p_juv), !is.na(Count)) %>%
  mutate(
    n_juv   = round(Count * p_juv),
    n_adult = Count - n_juv
  ) %>%
  pivot_longer(
    cols      = c(n_juv, n_adult),
    names_to  = "life_stage",
    values_to = "stage_Count"
  ) %>%
  filter(stage_Count > 0) %>%
  mutate(
    life_stage = recode(life_stage, "n_juv" = "juvenile", "n_adult" = "adult"),
    life_stage = factor(life_stage, levels = c("juvenile", "adult"))
  )

### 3) MAIN inference: probabilistic life-stage model ###########################

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

model_export(
  model       = m_stage_prob,
  model_name  = "T2_PrimaryInference_P",
  output_dir  = results_dir
)

### 4) SUPPLEMENTARY: deterministic life-stage models ###########################

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

model_export(
  model       = m_stage,
  model_name  = "S99_DeterministicLS_P",
  output_dir  = results_dir
)

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

### 5) Size structure model #####################################################

m1 <- glmmTMB(
  Count ~ Type * Size_Class_f +
    Pair +
    Count.Type +
    offset(log(Inclusion_m)) +
    (1 | Site) +
    (1 | survey_id),
  ziformula = ~ Type + Size_Class_f,
  family    = nbinom2,
  data      = fish_long
)

saveRDS(m1, file.path(fits_dir, paste0("m1_size_structure_", analysis_date, ".rds")))
capture.output(
  summary(m1),
  file = file.path(stats_dir, paste0("m1_size_structure_summary_", analysis_date, ".txt"))
)

model_export(
  model       = m1,
  model_name  = "S5_sizestruct_P",
  output_dir  = results_dir
)

emm_m1 <- emmeans(m1, ~ Type | Size_Class_f, offset = log(100))
capture.output(
  pairs(emm_m1),
  file = file.path(stats_dir, paste0("m1_emm_summary_", analysis_date, ".txt"))
)

### 6) Predictions and plots ####################################################

make_totalfish_plots <- function(m1,
                                 m_stage_fx,
                                 m_stage_prob,
                                 origin_date,
                                 plots_dir,
                                 analysis_date,
                                 reef_cols,
                                 theme_clean) {
  
  pred_tot <- ggpredict(
    m1,
    terms     = c("Type", "Pair"),
    condition = c(Inclusion_m = 100)
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
      y     = "Total fish density per 100 m²",
      color = "Reef type"
    )
  
  ggsave(
    filename = file.path(plots_dir, paste0("Fig1_total_density_", analysis_date, ".png")),
    plot     = p_tot,
    width    = 6.5,
    height   = 4.5,
    dpi      = 300
  )
  
  pred_size <- ggpredict(
    m1,
    terms     = c("Size_Class_f", "Type"),
    condition = c(Inclusion_m = 100)
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
      y     = "Predicted density per 100 m²",
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
  
  pred_size_pub <- pred_size %>%
    mutate(
      x = factor(
        x,
        levels = c("0-1","1-2","2-5","5-10","10-15","15-20","20-50","50-100","100+"),
        ordered = TRUE
      ),
      group = factor(group, levels = c("Natural","Artificial"))
    )
  
  p_size_pub <- ggplot(
    pred_size_pub,
    aes(x = x, y = predicted, color = group, fill = group, group = group)
  ) +
    geom_ribbon(
      aes(ymin = conf.low, ymax = conf.high),
      alpha = 0.12,
      color = NA
    ) +
    geom_line(linewidth = 1.1) +
    geom_point(size = 2.6) +
    scale_color_manual(values = reef_cols) +
    scale_fill_manual(values = reef_cols) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
    theme_clean +
    theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
    labs(
      x     = "Size class (cm)",
      y     = "Predicted density per 100 m²",
      color = "Reef type",
      fill  = "Reef type"
    )
  
  ggsave(
    filename = file.path(plots_dir, paste0("Fig2_size_structure_pub_", analysis_date, ".png")),
    plot     = p_size_pub,
    width    = 7,
    height   = 4.5,
    dpi      = 300
  )
  
  pred_stage_det <- ggpredict(
    m_stage_fx,
    terms     = c("life_stage", "Type", "Pair"),
    condition = c(Inclusion_m = 100)
  ) %>% as.data.frame() %>%
    mutate(
      x = factor(x, levels = c("juvenile", "mixed", "adult"), ordered = TRUE)
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
  
  invisible(list(
    p_tot       = p_tot,
    p_size      = p_size,
    p_size_pub  = p_size_pub,
    p_pair_det  = p_pair_det,
    p_pair_prob = p_pair_prob
  ))
}

plots_totalfish <- make_totalfish_plots(
  m1            = m1,
  m_stage_fx    = m_stage_fx,
  m_stage_prob  = m_stage_prob,
  origin_date   = origin_date,
  plots_dir     = plots_dir,
  analysis_date = analysis_date,
  reef_cols     = reef_cols,
  theme_clean   = theme_clean
)

plots_totalfish$p_pair_prob

### 7) Contrasts and exports (primary inference) ################################

em_ls_pair <- emmeans(
  m_stage_prob,
  ~ Type | life_stage * Pair,
  at = list(
    Inclusion_m = 100,
    date_num    = 0,
    Count.Type  = "Belt"
  ),
  type = "response"
)

ct_ls_pair <- as.data.frame(contrast(em_ls_pair, method = "revpairwise"))

out_ct <- file.path(results_dir, paste0("T3_PrimaryInference_emmeans_contrasts_", analysis_date, ".csv"))
write.csv(ct_ls_pair, out_ct, row.names = FALSE)

out_em <- file.path(results_dir, paste0("T3_PrimaryInference_emmeans_grid_", analysis_date, ".csv"))
write.csv(as.data.frame(em_ls_pair), out_em, row.names = FALSE)

em_juv_AR <- emmeans(
  m_stage_prob,
  ~ Pair,
  at = list(
    Type        = "Artificial",
    life_stage  = "juvenile",
    Inclusion_m = 100,
    date_num    = 0,
    Count.Type  = "Belt"
  ),
  type = "response"
)

ct_juv_AR <- as.data.frame(contrast(em_juv_AR, method = "pairwise"))

out_ct2 <- file.path(results_dir, paste0("T3b_PrimaryInference_emmeans_contrasts_", analysis_date, ".csv"))
write.csv(ct_juv_AR, out_ct2, row.names = FALSE)

message("Main models done. Plots and results saved to: ", output_dir)
