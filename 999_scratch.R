### trying stuff out page #### 
# checking size strcture of S rivulatus 


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

