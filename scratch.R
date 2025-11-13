# Install if you haven't already
# install.packages("rfishbase")

library(rfishbase)
library(dplyr)
library(brms)
library(tidyverse) 

# Step 1: Your list of species
unique_fish <- tibble::tribble(
  ~Original_Name,                 ~FishBase_Name,
  "Neopomacentrus cyanos",        "Neopomacentrus cyanomos",
  "Pomacentrus alexanderae",      "Pomacentrus alexanderae",
  "Scarus rivulatus",             "Scarus rivulatus",
  "Siganus javus",                "Siganus javus",
  "Chaetodon wiebeli",            "Chaetodon wiebeli",
  "Heniochus acuminatus",         "Heniochus acuminatus",
  "Chaetodon octofasciatus",      "Chaetodon octofasciatus",
  "Pomacanthus annularis",        "Pomacanthus annularis",
  "Ephippidae spp.",              "Platax teira",
  "Caesio xanthonota",            "Caesio xanthonota",
  "Labroides dimidiatus",         "Labroides dimidiatus",
  "Cheilinus fasciatus",          "Cheilinus fasciatus",
  "Thalassoma lunare",            "Thalassoma lunare",
  "Ballistoides viridescens",     "Balistoides viridescens",
  "Lutjanus vitta",               "Lutjanus vitta",
  "Lutjanus russellii",           "Lutjanus russellii",
  "Lutjanus griseus",             "Lutjanus griseus",
  "Diagramma pictum",             "Diagramma pictum",
  "Plectorhinchus chaetodonoides","Plectorhinchus chaetodonoides",
  "Plectorhinchus gibbosus",      "Plectorhinchus gibbosus",
  "Ephinphelus fasciatus",        "Epinephelus fasciatus",
  "Plectropomus spp.",            NA,
  "Epinephelus fuscoguttatus",    "Epinephelus fuscoguttatus",
  "Carangoides bajad",            "Carangoides bajad",
  "Carangoides fulvoguttatus",    "Turrum fulvoguttatus",
  "Gnanthanodon speciosus",       "Gnathanodon speciosus",
  "Sphyraena flavicauda",         "Sphyraena flavicauda",
  "Sphyraena qenie",              "Sphyraena qenie",
  "Sphyraena jello",              "Sphyraena jello",
  "Sphyraena barracuda",          "Sphyraena barracuda"
)






# Step 2: Get Lmax values from popgrowth
species_data <- species(unique_fish$FishBase_Name) %>%
  select(Species, Length, LTypeMaxM)



estimate_Lmat_direct <- function(Lmax) {
  10^(-0.282 + 0.8979 * log10(Lmax))
}


species_life_history <- species_data %>%
  filter(!is.na(Length)) %>%
  mutate(
    Lmax_cm = Length,
    Lmat_cm = estimate_Lmat(Lmax_cm)
  ) %>%
  select(Species, Lmax_cm, Lmat_cm)


# First ensure both data frames have the same join key and then join them 
lookup <- unique_fish %>%
  left_join(species_life_history, by = c("FishBase_Name" = "Species"))







##### use these to apply juvenile and adult labels to fish size 
long_fish <- fish_model %>%
  pivot_longer(cols = starts_with("bin_"), names_to = "Size_Class", values_to = "Count")

long_fish <- long_fish %>%
  mutate(
    Size_Min = as.numeric(str_extract(Size_Class, "(?<=bin_)[0-9]+")),
    Size_Max = as.numeric(str_extract(Size_Class, "(?<=_)[0-9]+|(?<=plus).*")) %>% 
      replace_na(120),  # assume upper bound of 120 cm for 100+ bin
    Midpoint_cm = (Size_Min + Size_Max) / 2
  )
long_fish <- long_fish %>%
  left_join(lookup, by = c("Sci_Name" = "Original_Name")) %>%
  mutate(Maturity = case_when(
    Midpoint_cm < Lmat_cm ~ "juvenile",
    Midpoint_cm >= Lmat_cm ~ "adult",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Maturity))  # drop unmatched or NA cases


survey_summary <- long_fish %>%
  group_by(Site, Date_mm.dd.yy, Researchers, Sci_Name, Maturity) %>%
  summarise(n = sum(Count, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = Maturity,
    values_from = n,
    values_fill = 0
  ) %>%
  mutate(total = juvenile + adult)

brm_juv_adult <- brm(
  juvenile | trials(total) ~ Site_Type + (1 | Site) + (1 | Sci_Name),
  data = survey_summary %>%
    left_join(fish_model %>% select(Site, Date_mm.dd.yy, Site_Type) %>% distinct(), 
              by = c("Site", "Date_mm.dd.yy")),
  family = binomial(), 
  backend = "cmdstanr"
)
summary(brm_juv_adult)


bbfish <- brm(
  juvenile | trials(total) ~ Site_Type + (Site_Type | Sci_Name) + (1 | Site),
  data = survey_summary %>%
    left_join(fish_model %>% select(Site, Date_mm.dd.yy, Site_Type) %>% distinct(), 
              by = c("Site", "Date_mm.dd.yy")),
  family = binomial(), 
  backend = "cmdstanr"
)

summary(bbfish)

#### extract effects 

library(brms)
library(dplyr)
library(ggplot2)

# Extract random slopes for Site_TypeWreck grouped by species
species_slopes <- ranef(bbfish)$Sci_Name[, , "Site_TypeWreck"] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Species") %>%
  rename(
    estimate = Estimate,
    lower = Q2.5,
    upper = Q97.5
  )
species_slopes <- species_slopes %>%
  arrange(estimate) %>%
  mutate(Species = factor(Species, levels = Species))
ggplot(species_slopes, aes(x = estimate, y = Species)) +
  geom_point() +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    x = "Effect of Wreck on Juvenile Proportion (log-odds)",
    y = "Species",
    title = "Species-specific juvenile/adult slope by Site Type"
  ) +
  theme_minimal()




# Summarise total counts per species
library(dplyr)

species_counts <- fish_model %>%
  group_by(Sci_Name) %>%
  summarise(total_individuals = sum(total_count, na.rm = TRUE),
            n_surveys = n()) %>%
  arrange(desc(total_individuals))

# View top and bottom
print(species_counts, n=30)

# Filter for species with > X total individuals (e.g., 50 as a working threshold)
well_sampled_species <- species_counts %>%
  filter(total_individuals >= 50) %>%
  pull(Sci_Name)

# Subset the main data for only those species
fish_model_filtered <- fish_model %>%
  filter(Sci_Name %in% well_sampled_species)



# Assign group based on Inclusion_m and Count.Type
fish_model_filtered <- fish_model_filtered %>%
  mutate(
    Ecological_Group = case_when(
      Inclusion_m == 1 ~ "Small Reef Fish",
      Inclusion_m == 5 & Count.Type == "Stationary" ~ "Mid-size Invertivores",
      Inclusion_m == 5 & Count.Type == "Belt" ~ "Large Predators",
      TRUE ~ NA_character_
    )
  )
table(fish_model_filtered$Ecological_Group, useNA = "ifany")


long_fish <- fish_model_filtered %>%
  pivot_longer(cols = starts_with("bin_"), names_to = "Size_Class", values_to = "Count") %>%
  mutate(
    Size_Min = as.numeric(str_extract(Size_Class, "(?<=bin_)[0-9]+")),
    Size_Max = as.numeric(str_extract(Size_Class, "(?<=_)[0-9]+|(?<=plus).*")) %>% 
      replace_na(120),
    Midpoint_cm = (Size_Min + Size_Max) / 2
  ) %>%
  left_join(lookup, by = c("Sci_Name" = "Original_Name")) %>%
  mutate(Maturity = case_when(
    Midpoint_cm < Lmat_cm ~ "juvenile",
    Midpoint_cm >= Lmat_cm ~ "adult",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Maturity))


survey_summary <- long_fish %>%
  group_by(Site, Site_Type, Date_mm.dd.yy, Researchers, Sci_Name, Ecological_Group, Maturity) %>%
  summarise(n = sum(Count, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = Maturity,
    values_from = n,
    values_fill = 0
  ) %>%
  mutate(total = juvenile + adult)

# fit with new interaction 
brm_ecogroup <- brm(
  juvenile | trials(total) ~ Site_Type  * Ecological_Group + (1 | Sci_Name) + (1 | Site),
  data = survey_summary,
  family = binomial(),
  chains = 4, cores = 4, iter = 2000
)

# or 
split_data <- survey_summary %>%
  group_split(Ecological_Group)



summary(brm_ecogroup)
# Extract random slopes for Site_TypeWreck grouped by species
species_slopes <- ranef(brm_ecogroup)$Sci_Name[, , "Site_TypeWreck"] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Species") %>%
  rename(
    estimate = Estimate,
    lower = Q2.5,
    upper = Q97.5
  )
species_slopes <- species_slopes %>%
  arrange(estimate) %>%
  mutate(Species = factor(Species, levels = Species))
ggplot(species_slopes, aes(x = estimate, y = Species)) +
  geom_point() +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    x = "Effect of Wreck on Juvenile Proportion (log-odds)",
    y = "Species",
    title = "Species-specific juvenile/adult slope by Site Type"
  ) +
  theme_minimal()


split_data <- survey_summary %>%
  group_split(Ecological_Group)



# Fit a separate model per group (you can loop this if you like) 
model_small <- brm(
  juvenile | trials(total) ~ Site_Type,
  data = split_data[[1]], 
  iter = 4000, chains =4, cores = 4, backend = "cmdstanr", control = list(adapt_delta = 0.95), 
  prior = prior,
  family = binomial()
)

summary(model_small)


data_mid <- split_data[[2]] %>% filter(total > 0)

model_mid <- brm(
  juvenile | trials(total) ~ Site_Type,
  data = data_mid,
  family = binomial(),
  control = list(adapt_delta = 0.95)
)
summary(model_mid)



data_big <- split_data[[3]] %>% filter(total > 0)

prior <- c(
  prior(normal(0, 5), class = "b"),
  prior(normal(0, 5), class = "Intercept")
)

model_big <- brm(
  juvenile | trials(total) ~ Site_Type,
  data = data_big,
  family = binomial(),
  prior = prior,
  control = list(adapt_delta = 0.99)
)
summary(model_big)
plot(model_big)



library(dplyr)
library(ggplot2)

# Filter out rows with total == 0 to avoid NaNs
survey_summary_filtered <- survey_summary %>%
  filter(total > 0) %>%
  mutate(prop_juv = juvenile / total)

# Summarize mean proportion juvenile per species
species_props <- survey_summary_filtered %>%
  group_by(Sci_Name) %>%
  summarise(
    mean_prop_juv = mean(prop_juv, na.rm = TRUE),
    n = n()
  ) %>%
  arrange(desc(mean_prop_juv))

# Plot
ggplot(species_props, aes(x = reorder(Sci_Name, mean_prop_juv), y = mean_prop_juv)) +
  geom_col(fill = "#66BFA6") +
  coord_flip() +
  labs(
    x = "Species",
    y = "Mean Proportion Juvenile",
    title = "Juvenile Proportion by Species"
  ) +
  theme_minimal(base_size = 14)

