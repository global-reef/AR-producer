######### 01 CLEANING DATA   ############################################

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tibble)
})

## Main cleaning function ####################################################

clean_fish_size <- function(file_path, output_dir) {
  
  ## 1. Read raw data --------------------------------------------------------
  raw_fish <- read.csv(file_path, stringsAsFactors = TRUE, strip.white = TRUE)
  
  ## 2. Standardise site names and attach site metadata ----------------------
  fish_size <- raw_fish %>%
    mutate(
      Site_raw = Site,
      Site = dplyr::recode(
        Site,
        "Aow Mao Wall" = "Aow Mao",
        "No Name"      = "No Name Pinnacle"
      ),
      Site = factor(Site)
    )
  
  site_lookup <- tibble::tribble(
    ~Site,              ~Type,        ~Pair,
    "Aow Mao",          "Natural",    "Aow Mao",
    "Aow Mao Wreck",    "Artificial", "Aow Mao",
    "Hin Pee Wee",      "Natural",    "Sattakut",
    "Sattakut",         "Artificial", "Sattakut",
    "No Name Pinnacle", "Natural",    "No Name",
    "No Name Wreck",    "Artificial", "No Name"
  )
  
  fish_size <- fish_size %>%
    left_join(site_lookup, by = "Site") %>%
    mutate(
      Type = factor(Type, levels = c("Natural", "Artificial")),
      Pair = factor(Pair),
      Species = factor(Species),
      Site = factor(Site)
    )
  
  ## 3. Dates and Month_Year -------------------------------------------------
  fish_size <- fish_size %>%
    mutate(
      Date = as.Date(as.character(Date_mm.dd.yy), format = "%m/%d/%Y"),
      # Fix mangled years (< 2024 assumed to be 2025, same month/day)
      Date = if_else(
        Date < as.Date("2024-01-01") & !is.na(Date),
        as.Date(paste0("2025-", format(Date, "%m-%d"))),
        Date
      ),
      Month_Year = format(Date, "%Y-%m")
    ) %>%
    dplyr::select(-Date_mm.dd.yy, -Site_raw) %>% 
    mutate(Month_Year = factor(Month_Year))
  
  ## 4. Rename size bins to bin_*_* -----------------------------------------
  size_cols <- grep("^X", names(fish_size), value = TRUE)
  
  new_size_names <- size_cols |>
    sub("^X\\.?", "bin_", x = _) |>
    gsub("\\.", "_", x = _)
  
  names(fish_size)[match(size_cols, names(fish_size))] <- new_size_names
  
  ## 5. Structural size-bin NAs ----------------------------------------------
  ### Species that should never exceed 20 cm
  spp_max20 <- c(
    "Damsels - Regal Demoiselle",
    "Damsels - Alexanders"
  )
  
  ### Species that can be up to 50 cm, but never > 50
  spp_max50 <- c(
    "Parrotfish - Surf",
    "Rabbit - Java",
    "Butterfly - Weibels",
    "Butterfly - Longfin bannerfish",
    "Butterfly - Eight banded",
    "Angel - Blue-ringed"
  )
  
  fish_size <- fish_size %>%
    mutate(
      bin_20_50  = ifelse(Species %in% spp_max20, NA_integer_, bin_20_50),
      bin_50_100 = ifelse(Species %in% c(spp_max20, spp_max50),
                          NA_integer_, bin_50_100),
      bin_100    = ifelse(Species %in% c(spp_max20, spp_max50),
                          NA_integer_, bin_100)
    )
  
  ## 6. Functional group lookup ----------------------------------------------
  spp_lookup <- tibble::tribble(
    ~Species,                          ~fgroup,
    "Angel - Blue-ringed",            "Invertivore",
    "Barracuda - Chevron",            "Mesopredator",
    "Barracuda - Great",              "HTLP",
    "Barracuda - Pickhandle",         "HTLP",
    "Barracuda - Yellowtail",         "Mesopredator",
    "Batfish - ALL",                  "Invertivore",
    "Butterfly - Eight banded",       "Grazer",
    "Butterfly - Longfin bannerfish", "Grazer",
    "Butterfly - Weibels",            "Grazer",
    "Cleaner - Blue-streaked",        "Invertivore",
    "Damsels - Alexanders",           "Grazer",
    "Damsels - Regal Demoiselle",     "Grazer",
    "Fusiliers - Yellowback",         "Grazer",
    "Grouper - Blacktip",             "HTLP",
    "Grouper - Brown marbled",        "HTLP",
    "Grouper - Coral groupers (all)", "HTLP",
    "Parrotfish - Surf",              "Grazer",
    "Rabbit - Java",                  "Grazer",
    "Rabbit - Virgate",               "Grazer",
    "Snapper - Brownstripe",          "Mesopredator",
    "Snapper - Mangrove",             "HTLP",
    "Snapper - Russells",             "Mesopredator",
    "Sweetlips - Harlequin",          "Mesopredator",
    "Sweetlips - Harry hotlips",      "HTLP",
    "Sweetlips - Painted",            "Mesopredator",
    "Tigger - Titan",                 "Invertivore",
    "Trevally - Gold spotted",        "HTLP",
    "Trevally - Golden",              "HTLP",
    "Trevally - Orange spotted",      "HTLP",
    "Wrasse - Moon",                  "Invertivore",
    "Wrasse - Redbreasted",           "Invertivore"
  )
  
  fish_size <- fish_size %>%
    left_join(spp_lookup, by = "Species") %>%
    mutate(
      fgroup = factor(
        fgroup,
        levels = c("Grazer", "Invertivore", "Mesopredator", "HTLP")
      )
    )
  
  # fix some names 
  fish_size <- fish_size %>%
    mutate(
      Sci_Name = case_when(
        Sci_Name == "Carangoides bajad"         ~ "Flavocaranx bajad",
        Sci_Name == "Carangoides fulvoguttatus" ~ "Turrum fulvoguttatum",
        TRUE                                    ~ Sci_Name
      )
    )
  
  fix_dict <- c(
    "Neopomacentrus cyanos"     = "Neopomacentrus cyanomos",
    "Chaetodon wiebeli"         = "Chaetodon weibeli",
    "Ephinphelus fasciatus"     = "Epinephelus fasciatus",
    "Gnanthanodon speciosus"    = "Gnathanodon speciosus",
    "Siganus Virgatus"          = "Siganus virgatus",
    "Epiphelus spp."            = "Ephippidae spp."
  )
  
  # Apply recodes, then drop the invalid Epiphelus spp.
  fish_size <- fish_size %>%
    mutate(Sci_Name = recode(Sci_Name, !!!fix_dict)) 
  
  fish_size <- fish_size %>% filter(Sci_Name != "Epinephelus spp.")
  
  fish_size <- fish_size %>% 
    mutate(
      Species = factor(Species),
      Sci_Name = factor(Sci_Name))
  
  ## 7. Time cleaning (Time (hh:mm)) -----------------------------------------
  fish_size <- fish_size %>%
    mutate(
      Time = as.character(Time_start),
      Time = if_else(
        str_detect(Time, "^[0-9]:"),
        str_pad(Time, width = 5, side = "left", pad = "0"),
        Time
      )
    )
  
  ## 8. Reorder columns: survey info, covariates, counts ---------------------
  survey_cols <- c(
    "Site", "Pair", "Type",
    "Date", "Month_Year",
    "Time", "Count.Type",
    "Inclusion_m", "Researchers"
  )
  
  covar_cols <- c(
    "Depth_m", "Visibility_m", "Weather",
    "Current", "Boats",
    "Species", "Sci_Name", "fgroup"
  )
  
  count_cols <- grep("^bin_", names(fish_size), value = TRUE)
  
  fish_size <- fish_size %>%
    dplyr::select(
      all_of(survey_cols),
      all_of(covar_cols),
      all_of(count_cols),
      -any_of("Time_start")
    )
  

  
  ## 9. Save cleaned file ----------------------------------------------------
  out_path <- file.path(output_dir, "fish_size_cleaned.rds")
  saveRDS(fish_size, out_path)
  message("✅ Saved 'fish_size_cleaned.rds' to: ", out_path)
  
  return(fish_size)
}

fish_size <- clean_fish_size(file_path, output_dir)

## Quick check
str(fish_size)


fish_size <- readRDS(file.path(output_dir, "fish_size_cleaned.rds"))

## Build long format for size classes #######################################

# Identify size class columns
size_cols <- grep("^bin_", names(fish_size), value = TRUE)

# Map bin_* names to human readable size classes
size_map <- c(
  bin_0_1    = "0-1",
  bin_1_2    = "1-2",
  bin_2_5    = "2-5",
  bin_5_10   = "5-10",
  bin_10_15  = "10-15",
  bin_15_20  = "15-20",
  bin_20_50  = "20-50",
  bin_50_100 = "50-100",
  bin_100    = "100+"
)

size_levels <- c("0-1", "1-2", "2-5", "5-10", "10-15", "15-20",
                 "20-50", "50-100", "100+")

fish_long <- fish_size %>%
  tidyr::pivot_longer(
    cols      = all_of(size_cols),
    names_to  = "Size_Class_raw",
    values_to = "Count"
  ) %>%
  mutate(
    Size_Class = factor(size_map[Size_Class_raw],
                        levels = size_levels,
                        ordered = TRUE),
    site       = Site,
    type       = Type,
    pair       = Pair,
    # survey_id = unique survey per site x date x count type
    survey_id  = interaction(Site, Date, Count.Type, drop = TRUE)
  ) %>%
  # keep rows where size class exists
  filter(!is.na(Size_Class))

fish_long <- fish_long %>%
  mutate(
    Size_Class_f = factor(Size_Class, ordered = FALSE)  # take away factor orders
  )

