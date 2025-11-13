library(brms)
library(tidyverse)

# define the ordered size class columns
size_bins <- c("0-1", "1-2", "2-5", "5-10", "10-20", "20-50", "50-100", "100+")


# Make a copy of the dataset for modeling
fish_model <- fish_size %>%
  rename(
    bin_0_1 = `0-1`,
    bin_1_2 = `1-2`,
    bin_2_5 = `2-5`,
    bin_5_10 = `5-10`,
    bin_10_20 = `10-20`,
    bin_20_50 = `20-50`,
    bin_50_100 = `50-100`,
    bin_100_plus = `100+`
  ) %>%
  mutate(total_count = rowSums(across(starts_with("bin_"))))


# fulllllll mulitvariate model 
fit_multinom <- brm(
  cbind(bin_0_1, bin_1_2, bin_2_5, bin_5_10, bin_10_20, bin_20_50, bin_50_100, bin_100_plus) |
    trials(total_count) ~ Site_Type + (1 | Site) + (1 | Species),
  data = fish_model,
  family = multinomial(),
  chains = 4, cores = 4, iter = 2000, 
  backend = "cmdstanr"
)


summary(fit_multinom)



# snappers only but all spp together 
library(stringr)
snappers <- fish_model %>%
  filter(str_detect(Sci_Name, "Lutjanus")) %>%
  filter(rowSums(across(starts_with("bin_"))) > 0)  # drop rows with 0 total count if needed


estimate_Lmat <- function(Linf_cm) {
  log_Lmat <- -0.1189 + 0.9157 * log10(Linf_cm)
  Lmat_cm <- 10^log_Lmat
  return(Lmat_cm)
}


# Add Lmat estimates using the function
snappers <- snappers %>%
  mutate(Lmat_cm = case_when(
    Sci_Name == "Lutjanus vitta" ~ estimate_Lmat(40),
    Sci_Name == "Lutjanus russellii" ~ 29,  # observed from FishBase
    Sci_Name == "Lutjanus griseus" ~ estimate_Lmat(37),
    TRUE ~ NA_real_
  ))

# change size bins into midpoints, and then categorize if bigger or smaller than mat 
size_bins <- tibble(
  bin_name = c("bin_0_1", "bin_1_2", "bin_2_5", "bin_5_10", "bin_10_20", "bin_20_50", "bin_50_100", "bin_100_plus"),
  midpoint = c(0.5, 1.5, 3.5, 7.5, 15, 35, 75, 125)
)

snappers <- snappers %>%
  rowwise() %>%
  mutate(
    juvenile = sum(
      c_across(all_of(size_bins$bin_name))[size_bins$midpoint < Lmat_cm],
      na.rm = TRUE
    ),
    adult = sum(
      c_across(all_of(size_bins$bin_name))[size_bins$midpoint >= Lmat_cm],
      na.rm = TRUE
    )
  ) %>%
  ungroup() %>%
  mutate(
    juvenile = as.integer(juvenile),
    adult = as.integer(adult),
    total = juvenile + adult )



juv_snap <- brm(
  juvenile | trials(total) ~ Site_Type + (1 | Site) + (1 | Sci_Name),
  family = binomial(),
  data = snappers
)


summary(juv_snap)
plot(juv_snap)
# Do snappers, accounting for differences between species and sites, tend to shift size structure across habitat types?

