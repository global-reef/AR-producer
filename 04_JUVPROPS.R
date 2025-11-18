### 04_JUVENILE_PROPORTIONS.R ###
suppressPackageStartupMessages({
library(dplyr)
library(ggplot2)
})
### Build survey-level juvenile proportions ###

juv_prop_survey <- fish_long_life_prob %>%
  group_by(
    survey_id, Site, Pair, Type, Date,
    Count.Type, Inclusion_m
  ) %>%
  summarise(
    juv_count   = sum(stage_Count[life_stage == "juvenile"], na.rm = TRUE),
    adult_count = sum(stage_Count[life_stage == "adult"],    na.rm = TRUE),
    total_la    = juv_count + adult_count,
    .groups     = "drop"
  ) %>%
  filter(total_la > 0) %>%              # keep only surveys with some juv+adult
  mutate(
    juv_prop = juv_count / total_la,
    date_num = as.numeric(Date - min(Date, na.rm = TRUE))
  )

# Save table for later use
readr::write_csv(
  juv_prop_survey,
  file.path(summ_dir, paste0("juvenile_proportion_survey_", analysis_date, ".csv"))
)

### Plot 1 - juvenile proportion by reef type and pair over time ###

p_juv_time <- ggplot(
  juv_prop_survey,
  aes(x = Date, y = juv_prop, color = Type)
) +
  geom_point(alpha = 0.4, position = position_jitter(width = 0, height = 0.02)) +
  geom_smooth(se = TRUE, method = "loess", span = 0.8) +
  facet_wrap(~ Pair, nrow = 1) +
  scale_color_manual(values = reef_cols) +
  scale_x_date(date_breaks = "4 month",  date_labels = "%b\n%Y") +
  theme_clean +
  labs(
    x     = "Date (months)",
    y     = "Juvenile proportion",
    color = "Reef type"
  )
p_juv_time
ggsave(
  filename = file.path(plots_dir, paste0("Supp_JuvProp_time_", analysis_date, ".png")),
  plot     = p_juv_time,
  width    = 7,
  height   = 4.5,
  dpi      = 300
)

### Plot 2 - juvenile proportion by reef type and pair (boxplot) ###

p_juv_box <- ggplot(
  juv_prop_survey,
  aes(x = Type, y = juv_prop, fill = Type)
) +
  geom_boxplot(outlier.alpha = 0.4) +
  facet_wrap(~ Pair, nrow = 1) +
  scale_fill_manual(values = reef_cols) +
  theme_clean +
  labs(
    x    = "Reef type",
    y    = "Juvenile proportion",
    fill = "Reef type"
  )

ggsave(
  filename = file.path(plots_dir, paste0("Supp_JuvProp_box_", analysis_date, ".png")),
  plot     = p_juv_box,
  width    = 7,
  height   = 4.5,
  dpi      = 300
)

message("✅ 3A: Juvenile Proportions done! Plots and resules saved to to: ", output_dir)
