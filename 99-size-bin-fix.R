
###############################################################################
# Size-bin reconstruction for early surveys
#
# Early versions of the size-estimation protocol used a single 10–20 cm bin
# (X.10.20). Later surveys refined this into two bins: 10–15 cm (X.10.15) and
# 15–20 cm (X.15.20). To make the entire dataset comparable, we reconstruct the
# early size structure by splitting the old 10–20 cm counts into 10–15 and
# 15–20 cm portions.
#
# The split is done *species specifically*:
#   - We use only the later surveys (where both 10–15 and 15–20 exist) to
#     calculate the observed proportion of individuals falling into each sub-bin
#     for each species.
#   - For species without enough later data, we fall back to global proportions
#     calculated across all species.
#
#   For each species s:
#       p10_s = proportion in 10–15 cm
#       p15_s = proportion in 15–20 cm
#
#   For early rows with only a 10–20 cm count (n_10_20):
#       X.10.15  = round(n_10_20 * p10_s)
#       X.15.20  = n_10_20 - X.10.15
#
# After reconstruction, the original X.10.20 column is removed, leaving all
# surveys with a consistent and comparable set of size bins.
###############################################################################

########  4. Size bins to numeric ###########################################
########  4. Size bins to numeric ###########################################

size_cols <- c("X0.1", "X.1.2", "X.2.5", "X.5.10",
               "X.10.15", "X.15.20", "X.20.50",
               "X.50.100", "X.100", "X.10.20")

clean_int <- function(x) {
  x_chr <- as.character(x)
  x_chr[x_chr %in% c("", " ", "NA", "N/A")] <- NA
  suppressWarnings(as.integer(x_chr))
}

fish_size <- fish_size %>%
  mutate(across(all_of(size_cols), clean_int),
         X.10.20_orig = X.10.20)   # keep a copy for QC

#########  5. Species specific split of 10–20 into 10–15 and 15–20 ###########

sp_max <- tibble::tribble(
  ~Species,                        ~max_size_cm,
  "Damsels - Regal Demoiselle",     20,
  "Damsels - Alexanders",           20,
  "Parrotfish - Surf",              50,
  "Rabbit - Java",                  50,
  "Butterfly - Weibels",            50,
  "Butterfly - Longfin bannerfish", 50,
  "Butterfly - Eight banded",       50,
  "Angel - Blue-ringed",            50
  # add more if needed, others default to NA = no restriction
)
fish_size <- fish_size %>%
  left_join(sp_max, by = "Species")

# flag rows where any size data exists 
fish_size <- fish_size %>%
  mutate(
    size_data_present = if_any(
      c(X0.1, X.1.2, X.2.5, X.5.10, X.20.50, X.50.100, X.100),
      ~ !is.na(.x)
    )
  )


# compute species-specific proportions for the 10-15 and 15-20 size classes from later data 
# Only rows with true 10–15 and 15–20 values
split_base <- fish_size %>%
  filter(!is.na(X.10.15), !is.na(X.15.20)) %>%
  mutate(total_10_20 = X.10.15 + X.15.20) %>%
  filter(total_10_20 > 0)

# Global fallback proportions
global_props <- split_base %>%
  summarise(
    sum_10_15 = sum(X.10.15, na.rm = TRUE),
    sum_15_20 = sum(X.15.20, na.rm = TRUE)
  ) %>%
  mutate(
    p10_global = sum_10_15 / (sum_10_15 + sum_15_20),
    p15_global = 1 - p10_global
  )

p10_global <- global_props$p10_global
p15_global <- global_props$p15_global

# Species specific proportions
spp_props <- split_base %>%
  group_by(Species) %>%
  summarise(
    sum_10_15 = sum(X.10.15, na.rm = TRUE),
    sum_15_20 = sum(X.15.20, na.rm = TRUE),
    total_10_20 = sum_10_15 + sum_15_20,
    .groups = "drop"
  ) %>%
  mutate(
    p10 = ifelse(total_10_20 > 0, sum_10_15 / total_10_20, NA_real_),
    p15 = 1 - p10
  )

# Combine with global fallback
spp_props_full <- tibble(Species = unique(fish_size$Species)) %>%
  left_join(spp_props, by = "Species") %>%
  mutate(
    p10 = ifelse(is.na(p10), p10_global, p10),
    p15 = ifelse(is.na(p15), p15_global, p15)
  )

# apply proportional split ONLY when valid 
fish_size <- fish_size %>%
  left_join(spp_props_full, by = "Species") %>%
  mutate(
    
    # CASE A: valid old 10–20 bin that needs splitting
    split_allowed =
      size_data_present &
      !is.na(X.10.20_orig) &
      X.10.20_orig > 0 &                       # must be >0, include 1
      is.na(X.10.15) & is.na(X.15.20) &        # no sub-bin data yet
      (is.na(max_size_cm) | max_size_cm >= 15),
    
    # CASE B: X.10.20_orig == 0 and diver DID estimate size
    # (a true zero; must produce 0/0 in sub bins)
    true_zero =
      size_data_present &
      !is.na(X.10.20_orig) &
      X.10.20_orig == 0 &
      is.na(X.10.15) & is.na(X.15.20),
    
    # APPLY CASE A
    X.10.15 = ifelse(split_allowed,
                     as.integer(round(X.10.20_orig * p10)),
                     X.10.15),
    
    X.15.20 = ifelse(split_allowed,
                     as.integer(X.10.20_orig - X.10.15),
                     X.15.20),
    
    # APPLY CASE B
    X.10.15 = ifelse(true_zero, 0L, X.10.15),
    X.15.20 = ifelse(true_zero, 0L, X.15.20)
  ) %>%
  select(-p10, -p15, -split_allowed, -true_zero)


# quick console check
summary(fish_size$X.10.15)
summary(fish_size$X.15.20)

# remove old 10-20 bin 
fish_size <- fish_size %>%
  select(-X.10.20, 
         # -X.10.20_orig, 
         -max_size_cm)
na_10_15 <- fish_size %>%
  filter(is.na(X.10.15)) %>%
  select(Site, Date, Species, X.10.20_orig, X0.1:X.100) %>%
  head(30)

na_15_20 <- fish_size %>%
  filter(is.na(X.15.20)) %>%
  select(Site, Date, Species, X.10.20_orig, X0.1:X.100) %>%
  head(30)

print(na_15_20) 