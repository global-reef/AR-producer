### 07_FINAL_PRODUCTION_SUMM.R ###
### Final producer vs attractor inference ###
suppressPackageStartupMessages({
library(dplyr)
library(glmmTMB)
library(ggeffects)
library(ggplot2)
library(emmeans)
library(tidyr)
})
message("### Starting final production summary ###")

## Inputs already created earlier
## fish_long_life_prob
## m_stage_prob
## origin_date
## theme_clean, reef_cols
## output_dir, plots_dir, stats_dir, summ_dir

### 1. Build survey level juvenile and adult totals ###

juv_adult_total <- fish_long_life_prob %>%
  group_by(
    survey_id, Site, Pair, Type, Date,
    Count.Type, Inclusion_m, date_num
  ) %>%
  summarise(
    juv   = sum(stage_Count[life_stage == "juvenile"]),
    adult = sum(stage_Count[life_stage == "adult"]),
    .groups = "drop"
  ) %>%
  mutate(
    juv_per100   = juv   / (Inclusion_m / 100),
    adult_per100 = adult / (Inclusion_m / 100)
  )

readr::write_csv(
  juv_adult_total,
  file.path(summ_dir, paste0("survey_juv_adult_totals_", analysis_date, ".csv"))
)


### 2. Fit final J and A models (clean and simple NB2) ###

m_juv <- glmmTMB(
  juv ~ Type + Pair + date_num +
    offset(log(Inclusion_m)) +
    (1 | Site) + (1 | survey_id),
  family = nbinom2,
  data   = juv_adult_total
)

m_adult <- glmmTMB(
  adult ~ Type + Pair + date_num +
    offset(log(Inclusion_m)) +
    (1 | Site) + (1 | survey_id),
  family = nbinom2,
  data   = juv_adult_total
)

saveRDS(m_juv,  file.path(fits_dir,  paste0("m_final_juv_",  analysis_date, ".rds")))
saveRDS(m_adult,file.path(fits_dir,  paste0("m_final_adult_",analysis_date, ".rds")))

capture.output(summary(m_juv),
               file = file.path(stats_dir, paste0("m_final_juv_summary_", analysis_date, ".txt")))
capture.output(summary(m_adult),
               file = file.path(stats_dir,paste0("m_final_adult_summary_",analysis_date, ".txt")))


### 3. Contrasts: AR vs NR for juveniles and adults ###

em_juv   <- emmeans(m_juv,   ~ Type)
em_adult <- emmeans(m_adult, ~ Type)

ct_juv   <- contrast(em_juv,   method = "revpairwise")
ct_adult <- contrast(em_adult, method = "revpairwise")

readr::write_csv(
  as.data.frame(ct_juv),
  file.path(summ_dir, paste0("contrast_juv_AR_NR_", analysis_date, ".csv"))
)

readr::write_csv(
  as.data.frame(ct_adult),
  file.path(summ_dir, paste0("contrast_adult_AR_NR_", analysis_date, ".csv"))
)


### 4. Prediction plots for clean inference ###

# Juveniles
pred_juv <- ggpredict(
  m_juv,
  terms = c("Type"),
  condition = list(Inclusion_m = 100)
) %>% as.data.frame()

p_juv <- ggplot(pred_juv, aes(x = x, y = predicted, color = x)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.15) +
  scale_color_manual(values = reef_cols) +
  theme_clean +
  labs(
    x = "Reef type",
    y = "Juvenile density (per 100 m2)",
    color = "Type"
  )

# Adults
pred_adult <- ggpredict(
  m_adult,
  terms = c("Type"),
  condition = list(Inclusion_m = 100)
) %>% as.data.frame()

p_adult <- ggplot(pred_adult, aes(x = x, y = predicted, color = x)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.15) +
  scale_color_manual(values = reef_cols) +
  theme_clean +
  labs(
    x = "Reef type",
    y = "Adult density (per 100 m2)",
    color = "Type"
  )

ggsave(
  file.path(plots_dir, paste0("Final_Juvenile_density_", analysis_date, ".png")),
  p_juv, width = 6, height = 4.2, dpi = 300
)

ggsave(
  file.path(plots_dir, paste0("Final_Adult_density_", analysis_date, ".png")),
  p_adult, width = 6, height = 4.2, dpi = 300
)

library(patchwork)

p_combine <- p_juv + p_adult
ggsave(
  file.path(plots_dir, paste0("FIG_99_Final_combined_density_", analysis_date, ".png")),
  p_combine, width = 6, height = 4.2, dpi = 300
)


### 5. Δ (AR minus NR) summaries ###

ct_juv   <- contrast(em_juv,   method = "revpairwise", infer = c(TRUE, TRUE))
ct_adult <- contrast(em_adult, method = "revpairwise", infer = c(TRUE, TRUE))

juv_df <- as.data.frame(ct_juv)   %>% mutate(stage = "juvenile")
adult_df <- as.data.frame(ct_adult) %>% mutate(stage = "adult")

delta_final <- bind_rows(juv_df, adult_df) %>%
  dplyr::select(
    stage, contrast, estimate, SE, df,
    asymp.LCL, asymp.UCL, p.value
  )

readr::write_csv(
  delta_final,
  file.path(summ_dir, paste0("Final_AR_NR_delta_summary_", analysis_date, ".csv"))
)


print_production_summary <- function(delta_final, total_df = NULL) {
message("### Final production summary complete ###")

cat("\n===== PRODUCTION SUMMARY =====\n")

# 1. Juvenile effect
juv_row <- delta_final %>% filter(stage == "juvenile")
cat("\nJuvenile AR minus NR:\n")
cat("  Estimate:", round(juv_row$estimate, 3), "\n")
cat("  CI:", round(juv_row$asymp.LCL, 3), "to", round(juv_row$asymp.UCL, 3), "\n")
cat("  p:", signif(juv_row$p.value, 3), "\n")

# 2. Adult effect
adult_row <- delta_final %>% filter(stage == "adult")
cat("\nAdult AR minus NR:\n")
cat("  Estimate:", round(adult_row$estimate, 3), "\n")
cat("  CI:", round(adult_row$asymp.LCL, 3), "to", round(adult_row$asymp.UCL, 3), "\n")
cat("  p:", signif(adult_row$p.value, 3), "\n")


cat("\n===== END SUMMARY =====\n")
} 

print_production_summary(delta_final, total_df)

summary(m_juv)
summary(m_adult)

## Extra summaries for main models ###########################################

suppressPackageStartupMessages({
  library(emmeans)
  library(ggeffects)
  library(dplyr)
  library(tidyr)
})






