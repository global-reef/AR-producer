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


##### 
library(dplyr)

# Distinct size-structure surveys
size_surveys <- fish_long_norm %>%
  distinct(survey_id, Site, Pair, Type, Date)

size_surveys %>% count()       # <N_size_surveys>

# Number of species with size-structured data
fish_long_norm %>%
  distinct(Sci_Name) %>%
  nrow()                       # <N_species_size>


# Distinct abundance surveys
abund_surveys <- fish_long %>%
  distinct(survey_id, site, pair, type, Date)

# Total surveys and by type
abund_surveys %>% count()                      # <N_abund_surveys>
abund_surveys %>% count(type)                 # <N_abund_AR>, <N_abund_NR>
abund_surveys %>% count(pair)                 # <N_AowMao>, <N_NoName>, <N_Sattakut>

# Date range
abund_surveys %>%
  summarise(
    start_date = min(Date, na.rm = TRUE),
    end_date   = max(Date, na.rm = TRUE)
  )


# Location: file.path(summ_dir, paste0("summ_species_overall_", analysis_date, ".csv"))
# If you are not sure about the date, inspect the directory:
list.files(summ_dir, pattern = "summ_species_overall_")

summ_spp_overall <- readr::read_csv(
  file.path(summ_dir, "summ_species_overall_2025.11.18.csv")
)

summ_spp_overall

# finding \Beta values for m_juv and m_adult 
coefs_juv <- summary(m_final_juv)$coefficients$cond
coefs_adult <- summary(m_final_adult)$coefficients$cond

beta_type_juv   <- coefs_juv["TypeArtificial", c("Estimate","Std. Error","Pr(>|z|)")]
beta_type_adult <- coefs_adult["TypeArtificial", c("Estimate","Std. Error","Pr(>|z|)")]

print(beta_type_juv)
print(beta_type_adult)



