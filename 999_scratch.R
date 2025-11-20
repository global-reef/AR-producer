### trying stuff out page #### 
# checking size strcture of S rivulatus 
library(dplyr)

riv_tab <- fish_long_norm %>%
  filter(Sci_Name == "Scarus rivulatus") %>%
  group_by(Type, Size_Class) %>%
  summarise(
    total_n = sum(Count, na.rm = TRUE),   # raw counts
    .groups = "drop"
  ) %>%
  arrange(Type, Size_Class)

riv_tab



p_riv_hist <- ggplot(riv_tab, aes(x = Size_Class, y = total_n, fill = Type)) +
  geom_col(position = "dodge") +
  labs(
    title = "Size structure of Scarus rivulatus by reef type",
    x = "Size class (cm)",
    y = "Total count"
  ) + 
  scale_color_manual(values = reef_cols) +
  scale_fill_manual(values = reef_cols) +
  theme_clean

p_riv_hist


triggers <- 
  riv_tab <- fish_long_norm %>%
  filter(Sci_Name == "Ballistoides viridescens") %>%
  group_by(Type, Size_Class) %>%
  summarise(
    total_n = sum(Count, na.rm = TRUE),   # raw counts
    .groups = "drop"
  ) %>%
  arrange(Type, Size_Class)

triggers

triggers_under10_simple <- fish_long_norm %>%
  filter(
    Sci_Name == "Ballistoides viridescens",
    Size_Class %in% c("0-1", "1-2", "2-5", "5-10"),
    Count > 0
  ) %>%
  distinct(Date, Site, Pair, Type, survey_id, Researchers) %>%
  arrange(Date, Site, Type, Researchers)

triggers_under10_simple

