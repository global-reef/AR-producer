### 05_FUNCTIONAL_GROUP_MODELS.R ###

library(glmmTMB)
library(dplyr)
library(ggeffects)
library(ggplot2)

### Prepare functional group x life-stage data ###

# Aggregate stage-level counts per survey x functional group
fish_fg_ls <- fish_long_life_prob %>%
  group_by(
    survey_id, Site, Pair, Type, Date,
    Count.Type, Inclusion_m, fgroup,
    life_stage, date_num, date_s
  ) %>%
  summarise(
    stage_Count = sum(stage_Count, na.rm = TRUE),
    .groups     = "drop"
  ) %>%
  # keep only rows with valid group, life stage, and area
  filter(
    !is.na(fgroup),
    !is.na(life_stage),
    !is.na(Inclusion_m),
    stage_Count > 0
  ) %>%
  mutate(
    life_stage = factor(life_stage, levels = c("juvenile", "adult"))
  )

### Fit multigroup functional group x life-stage model (main text) ###

m_fg_stage <- glmmTMB(
  stage_Count ~ Type * life_stage * fgroup +
    Pair +
    date_s +
    Count.Type +
    offset(log(Inclusion_m)) +
    (1 | Site) +
    (1 | survey_id),
  family = nbinom2,
  data   = fish_fg_ls
)

summary(m_fg_stage)


saveRDS(
  m_fg_stage,
  file.path(fits_dir, paste0("m_fg_stage_", analysis_date, ".rds"))
)

capture.output(
  summary(m_fg_stage),
  file = file.path(stats_dir, paste0("m_fg_stage_summary_", analysis_date, ".txt"))
)

### Plot function for functional group x life-stage model ###

make_fg_stage_plots <- function(m_fg_stage,
                                plots_dir,
                                analysis_date,
                                reef_cols,
                                theme_clean) {
  
  ## Main figure: Type x life_stage x fgroup (per 100 m²) ##
  
  pred_fg_ls <- ggpredict(
    m_fg_stage,
    terms     = c("life_stage", "Type", "fgroup"),
    condition = c(Inclusion_m = 100)   # 100 m²
  ) %>% as.data.frame()
  
  p_fg_stage <- ggplot(
    pred_fg_ls,
    aes(x = x, y = predicted, color = group, fill = group)
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
      y     = "Predicted density (per 100 m²)",
      color = "Reef type",
      fill  = "Reef type"
    )
  
  ggsave(
    filename = file.path(plots_dir, paste0("Fig_FG1_LifeStage_", analysis_date, ".png")),
    plot     = p_fg_stage,
    width    = 7,
    height   = 4.5,
    dpi      = 300
  )
  
  ## Supplementary: add Pair to see within-pair patterns (optional) ##
  
  pred_fg_ls_pair <- ggpredict(
    m_fg_stage,
    terms     = c("life_stage", "Type", "fgroup", "Pair"),
    condition = c(Inclusion_m = 100)
  ) %>% as.data.frame()
  
  p_fg_stage_pair <- ggplot(
    pred_fg_ls_pair,
    aes(x = x, y = predicted, color = group, fill = group)
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
    facet_grid(facet ~ panel) +  # facet = fgroup, panel = Pair
    scale_color_manual(values = reef_cols) +
    scale_fill_manual(values = reef_cols) +
    theme_clean +
    labs(
      x     = "Life stage",
      y     = "Predicted density (per 100 m²)",
      color = "Reef type",
      fill  = "Reef type"
    )
  
  ggsave(
    filename = file.path(plots_dir, paste0("Fig_FG2_LifeStage_Pair_", analysis_date, ".png")),
    plot     = p_fg_stage_pair,
    width    = 8,
    height   = 6,
    dpi      = 300
  )
  
  
  pred_fg_ls_pair$x <- factor(
    pred_fg_ls_pair$x,
    levels = c("juvenile", "adult")
  )
  
  dodge <- position_dodge(width = 0.5)
  
  p_fg_stage_pair2 <- ggplot(
    pred_fg_ls_pair,
    aes(x = x, y = predicted, color = group)
  ) +
    geom_errorbar(
      aes(ymin = conf.low, ymax = conf.high),
      width    = 0.15,
      position = dodge
    ) +
    geom_point(
      size     = 2.4,
      position = dodge
    ) +
    facet_grid(facet ~ panel, scales = "free_y") +   # facet = fgroup (rows), panel = Pair (cols)
    scale_color_manual(values = reef_cols) +
    theme_clean +
    labs(
      x     = "Life stage",
      y     = "Predicted density (per 100 m²)",
      color = "Reef type"
    )
  
  ggsave(
    filename = file.path(plots_dir, paste0("Fig_FG3_LifeStage_Pair_clean_", analysis_date, ".png")),
    plot     = p_fg_stage_pair2,
    width    = 8,
    height   = 6,
    dpi      = 300
  )
  
  
  ### Δ (AR − NR) per functional group, life stage, and pair ###
  
  delta_df <- pred_fg_ls_pair %>%
    select(facet, panel, x, group, predicted) %>%
    tidyr::pivot_wider(
      names_from  = group,
      values_from = predicted
    ) %>%
    mutate(
      delta  = Artificial - Natural,
      fgroup = factor(facet, levels = names(fg_cols)),
      x      = factor(x, levels = c("juvenile", "adult"))
    )
  # old 
  p_delta <- ggplot(
    delta_df,
    aes(x = x, y = delta, color = fgroup)
  ) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_point(size = 2.4) +
    facet_wrap(~ panel, nrow = 1) +  # one panel per pair
    scale_color_manual(values = fg_cols) +
    theme_clean +  
    labs(
      x     = "Life stage",
      y     = "Δ density (AR − NR, per 100 m²)",
      color = "Functional group"
    )
  
  p_delta 
  
  
  ggsave(
    file.path(plots_dir, paste0("Fig_FG_LifeStage_delta_", analysis_date, ".png")),
    p_delta,
    width = 7,
    height = 4.5,
    dpi = 300
  )
  
  # trying % change 
  
  delta_df2 <- pred_fg_ls_pair %>%
    select(facet, panel, x, group, predicted) %>%
    tidyr::pivot_wider(
      names_from  = group,
      values_from = predicted
    ) %>%
    mutate(
      pct_change = 100 * ((Artificial / Natural) - 1),
      fgroup     = factor(facet, levels = names(fg_cols)),
      x          = factor(x, levels = c("juvenile", "adult"))
    )
  p_delta2 <- ggplot(
    delta_df2,
    aes(x = x, y = pct_change, color = fgroup)
  ) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_point(size = 2.4, alpha = 0.9, position = position_dodge(width = 0.25)) +
    scale_color_manual(values = fg_cols) +
    theme_clean +
    labs(
      x     = "Life stage",
      y     = "Percent change (AR relative to NR)",
      color = "Functional group"
    )
  
  ggsave(
    file.path(plots_dir, paste0("Fig_FG_LifeStage_delta_percent_", analysis_date, ".png")),
    p_delta2,
    width = 7,
    height = 4.5,
    dpi = 300
  )
  
  invisible(list(
    p_fg_stage      = p_fg_stage,
    p_fg_stage_pair = p_fg_stage_pair,
    p_fg_stage_pair2 = p_fg_stage_pair2,
    p_delta = p_delta2
  ))
}

### Call the functional group stage-level plot function ###

plots_fg_stage <- make_fg_stage_plots(
  m_fg_stage   = m_fg_stage,
  plots_dir    = plots_dir,
  analysis_date = analysis_date,
  reef_cols    = reef_cols,
  theme_clean  = theme_clean
)

message("✅ Functional group stage-level modelling done! Plots and results saved to: ", output_dir)








