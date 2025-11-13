## exploratory anlaysis 
library(dplyr)
library(ggplot2)
library(RColorBrewer)



# mutate to long for exploratory analysis 
size_class_cols <- new_size_cols
fish_long <- fish_size %>%
  pivot_longer(
    cols = all_of(size_class_cols),
    names_to = "Size_Class",
    values_to = "Count"
  )
# Define ordered levels
ordered_levels <- c("0-1", "1-2", "2-5", "5-10", "10-20", "20-50", "50-100", "100+")
# Reorder the factor
fish_long <- fish_long %>%
  mutate(Size_Class = factor(Size_Class, levels = ordered_levels, ordered = TRUE))




##### exploratory analysis 
library(ggplot2)

# hist of size classes 
ggplot(fish_long, aes(x = Size_Class, y = Count, fill = Site_Type)) +
  geom_col(position = "dodge", na.rm = TRUE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "Fish Size Class Distribution by Habitat Type",
    x = "Size Class (cm)",
    y = "Fish Count",
    fill = "Site Type"
  ) +
  theme_minimal(base_size = 14) +  scale_fill_brewer("BuGn")

# standardize by site and survey 
fish_long <- fish_long %>%
  mutate(Site_Pair = case_when(
    Site %in% c("Aow Mao", "Aow Mao Wreck") ~ "Aow Mao Pair",
    Site %in% c("No Name Pinnacle", "No Name Wreck") ~ "No Name Pair",
    Site %in% c("Hin Pee Wee", "Sattakut") ~ "Sattakut Pair"
  ))
ggplot(fish_long, aes(x = Size_Class, y = Count, fill = Site_Type)) +
  geom_col(position = "dodge", na.rm = TRUE) +
  facet_wrap(~ Site_Pair, ncol = 1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "Size Class Distribution by Site Pair",
    x = "Size Class (cm)",
    y = "Fish Count",
    fill = "Site Type"
  ) +
  theme_minimal(base_size = 14) + scale_fill_brewer("BuGn") +
  theme(
    strip.text = element_text(size = 13, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# normalize by survey effort 
library(dplyr)

# Count number of unique surveys per site
survey_counts <- fish_long %>%
  distinct(Site, Date_mm.dd.yy) %>%
  count(Site, name = "n_surveys")
fish_long_norm <- fish_long %>%
  left_join(survey_counts, by = "Site") %>%
  mutate(Count_perSurvey = Count / n_surveys)
# size class hist
ggplot(fish_long_norm, aes(x = Size_Class, y = Count_perSurvey, fill = Site_Type)) +
  geom_col(position = "dodge", na.rm = TRUE) +
  facet_wrap(~ Site_Pair, ncol = 1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "Normalized Fish Size Class Distribution by Site Pair",
    x = "Size Class (cm)",
    y = "Mean Count per Survey",
    fill = "Site Type"
  ) +
  theme_minimal(base_size = 14) + scale_fill_brewer("BuGn") +
  theme(
    strip.text = element_text(size = 13, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )



fish_long %>%
  distinct(Site, Site_Type, Date_mm.dd.yy) %>%
  count(Site, Site_Type, name = "n_surveys")
fish_long %>%
  group_by(Site_Type, Sci_Name) %>%
  summarise(Total = sum(Count, na.rm = TRUE), .groups = "drop") %>%
  arrange(Site_Type, desc(Total))
# proportion of size class per site type 
fish_long %>%
  group_by(Site_Type, Size_Class) %>%
  summarise(Total = sum(Count, na.rm = TRUE), .groups = "drop") %>%
  group_by(Site_Type) %>%
  mutate(Proportion = Total / sum(Total)) %>%
  ggplot(aes(x = Size_Class, y = Proportion, fill = Site_Type)) +
  geom_col(position = "dodge") +
  scale_fill_brewer(palette = "BuGn") + theme_minimal() +
  labs(title = "Proportion of Fish in Each Size Class by Habitat",
       x = "Size Class", y = "Proportion")

# total fish per size class by sites to check for outliers 

# Set site order to control facet layout
site_order <- c("Aow Mao", "Aow Mao Wreck",
                "No Name Pinnacle", "No Name Wreck",
                "Hin Pee Wee", "Sattakut")

fish_long_norm <- fish_long %>%
  mutate(Site = factor(Site, levels = site_order))

# Count surveys per site
survey_counts <- fish_long_norm %>%
  distinct(Site, Date_mm.dd.yy) %>%
  count(Site, name = "n_surveys")

gnbu_colors <- brewer.pal(n = 9, name = "GnBu")
gnbu_shifted <- gnbu_colors[2:9]  # keeps 8 colors, drops lightest
size_order <- c("0-1", "1-2", "2-5", "5-10", "10-20", "20-50", "50-100", "100+")
fish_long_norm <- fish_long_norm %>%
  mutate(Size_Class = factor(Size_Class, levels = size_order))
# Normalize and plot
fish_long_norm %>%
  left_join(survey_counts, by = "Site") %>%
  mutate(Count_perSurvey = Count / n_surveys) %>%
  group_by(Site, Size_Class) %>%
  summarise(Mean_Count = sum(Count_perSurvey, na.rm = TRUE), .groups = "drop") %>%
  
  ggplot(aes(x = Size_Class, y = Mean_Count, fill = Size_Class)) +
  geom_col(position = "dodge", na.rm = TRUE) +
  facet_wrap(~ Site, nrow = 3, ncol = 2, scales = "free_y") + # free_y allows y axis to vary 
  scale_fill_manual(values = setNames(gnbu_shifted, size_order)) + 
  labs(
    title = "Mean Fish Count per Size Class (Normalized by Survey)",
    x = "Size Class (cm)",
    y = "Mean Count per Survey",
    fill = "Size Class"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )


# Prepare the data
fish_agg <- fish_long_norm %>%
  left_join(survey_counts, by = "Site") %>%
  mutate(Count_perSurvey = Count / n_surveys) %>%
  group_by(Site_Pair, Site, Size_Class) %>%
  summarise(Mean_Count = sum(Count_perSurvey, na.rm = TRUE), .groups = "drop")

# Filter and plot each site pair individually

# 1. Aow Mao Pair
aow_mao_plot <- fish_agg %>%
  filter(Site_Pair == "Aow Mao Pair") %>%
  ggplot(aes(x = Site, y = Mean_Count, fill = Size_Class)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = setNames(gnbu_shifted, size_order)) +
  labs(
    title = "Aow Mao Pair",
    x = "Site", y = "Mean Count per Survey", fill = "Size Class"
  ) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 2. No Name Pair
no_name_plot <- fish_agg %>%
  filter(Site_Pair == "No Name Pair") %>%
  ggplot(aes(x = Site, y = Mean_Count, fill = Size_Class)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = setNames(gnbu_shifted, size_order)) +
  labs(
    title = "No Name Pair",
    x = "Site", y = "Mean Count per Survey", fill = "Size Class"
  ) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 3. State Pair
satta_plot <- fish_agg %>%
  filter(Site_Pair == "Sattakut Pair") %>%
  ggplot(aes(x = Site, y = Mean_Count, fill = Size_Class)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = setNames(gnbu_shifted, size_order)) +
  labs(
    title = "Sattakut Pair",
    x = "Site", y = "Mean Count per Survey", fill = "Size Class"
  ) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Print them one by one
print(aow_mao_plot)
print(no_name_plot)
print(satta_plot)





# check 0s and NAs 
fish_long %>%
  group_by(Size_Class, Site_Type) %>%
  summarise(
    n_zero = sum(Count == 0, na.rm = TRUE),
    n_na = sum(is.na(Count)),
    n_total = n()
  )


fish_long %>%
  filter(str_detect(Sci_Name, regex("Gnathanodon speciosus", ignore_case = TRUE))) %>%    
  filter(!is.na(Count) & Count > 0)
