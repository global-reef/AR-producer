## data clean up 


library(dplyr)
library(stringr)
library(tidyr)


fish_size <- read.csv("~/Documents/1_GLOBAL REEF/0_PROJECTS/AR_Producer_Attractor/AR_Producer/2025.06.18_FishSize_MASTER.csv", stringsAsFactors=TRUE)



# clean up size columns and replace blanks with 0s 
old_size_cols <- c("X0.1", "X.1.2", "X.2.5", "X.5.10", 
                   "X.10.20", "X.20.50", "X.50.100", "X.100")

new_size_cols <- c("0-1", "1-2", "2-5", "5-10", 
                   "10-20", "20-50", "50-100", "100+")

# Only convert "" to 0, keep true NAs intact
fish_size[old_size_cols] <- lapply(fish_size[old_size_cols], function(x) {
  x <- as.character(x)
  x[x == ""] <- "0"         # replace only blank strings
  as.numeric(x)             # NA stays NA
})

# Rename columns
names(fish_size)[names(fish_size) %in% old_size_cols] <- new_size_cols
# drop unneeded columns for this analysis 
fish_size <- fish_size %>%
  select(-Time_start, -Depth_m, -Weather, -Current, -Boats)


# Clean and relabel the Site names
fish_size <- fish_size %>%
  mutate(Site = as.character(Site)) %>%
  mutate(Site = case_when(
    Site %in% c("313 (No Name Wreck)", "No Name Wreck") ~ "No Name Wreck",
    Site %in% c("No Name Pinnacle", "No Name") ~ "No Name Pinnacle",
    Site %in% c("Sattakut") ~ "Sattakut",
    Site %in% c("Hin Pee Wee") ~ "Hin Pee Wee",
    Site %in% c("Aow Mao Wreck") ~ "Aow Mao Wreck",
    Site %in% c("Aow Mao", "Aow Mao Natural") ~ "Aow Mao",
    TRUE ~ NA_character_  # Remove any unmatched sites
  ))
# Label each site as Wreck or Natural
fish_size <- fish_size %>%
  mutate(Site_Type = case_when(
    Site %in% c("No Name Wreck", "Sattakut", "Aow Mao Wreck") ~ "Wreck",
    Site %in% c("No Name Pinnacle", "Hin Pee Wee", "Aow Mao") ~ "Natural"
  )) %>%
  mutate(Site_Type = factor(Site_Type, levels = c("Natural", "Wreck"))) %>% 
  filter(
    !(is.na(Site) | Site == "") &              # Drop if Site is NA or blank
      !(is.na(Date_mm.dd.yy) | Date_mm.dd.yy == "")) %>%
  mutate(Site_Pair = case_when(
    Site %in% c("Aow Mao", "Aow Mao Wreck") ~ "Aow Mao Pair",
    Site %in% c("No Name Pinnacle", "No Name Wreck") ~ "No Name Pair",
    Site %in% c("Hin Pee Wee", "Sattakut") ~ "Sattakut Pair",
    TRUE ~ NA_character_
  )) %>%
  mutate(Site_Pair = factor(Site_Pair))

fish_size <- fish_size %>%
  mutate(Date_mm.dd.yy = as.Date(as.character(Date_mm.dd.yy), format = "%m/%d/%Y")) %>%
  filter(!(Site == "Aow Mao Wreck" & Date_mm.dd.yy == as.Date("2025-03-27") & Researchers == "Georgia"))


