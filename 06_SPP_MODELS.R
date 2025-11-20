### 06_SPECIES_MODELS.R ###
suppressPackageStartupMessages({
library(glmmTMB)
library(dplyr)
library(purrr)
library(tidyr)
library(tibble)
library(ggeffects)
library(ggplot2)
  library(patchwork)
})
### Prepare species-level dataset ####

# We use the probabilistic life stage data to keep the juvenile vs adult story
fish_species_ls <- fish_long_life_prob %>%
  # keep only species with non missing life stage counts
  filter(stage_Count > 0) %>%
  # make sure Sci_Name is character
  mutate(Sci_Name = as.character(Sci_Name))

### Species diagnostics helper ####

diagnose_species <- function(data) {
  
  all_species <- sort(unique(data$Sci_Name))
  
  map_dfr(all_species, function(sp) {
    
    df_sp <- data %>%
      filter(Sci_Name == sp) %>%
      droplevels()
    
    if (nrow(df_sp) == 0) {
      return(tibble(
        Sci_Name         = sp,
        n_rows           = 0L,
        total_count      = 0,
        n_Type           = 0L,
        n_stage          = 0L,
        n_Pair           = 0L,
        n_cells_nonzero  = 0L
      ))
    }
    
    tab_3d <- xtabs(stage_Count ~ Type + life_stage + Pair, df_sp)
    
    tibble(
      Sci_Name        = sp,
      n_rows          = nrow(df_sp),
      total_count     = sum(df_sp$stage_Count, na.rm = TRUE),
      n_Type          = nlevels(df_sp$Type),
      n_stage         = nlevels(df_sp$life_stage),
      n_Pair          = nlevels(df_sp$Pair),
      n_cells_nonzero = sum(tab_3d > 0)
    )
  })
}

species_diag <- diagnose_species(fish_species_ls)

# Set thresholds for candidate species

candidate_species <- species_diag %>%
  mutate(
    ok_candidate = n_rows >= 10 &
      total_count >= 40 &
      n_Type  >= 2 &
      n_stage >= 1 &
      n_Pair  >= 2 &
      n_cells_nonzero >= 4
  )

candidate_species %>%
  arrange(desc(ok_candidate), desc(total_count)) %>%
  print(n = 40)

key_species <- candidate_species %>%
  filter(ok_candidate) %>%
  pull(Sci_Name)

# ones that didn't make the cut 
candidate_species %>% filter(!ok_candidate)


### Data driven choice of model structure per species ####

choose_species_formula <- function(df_sp, min_cell = 5, min_total = 50) {
  tab_ts  <- xtabs(stage_Count ~ Type + life_stage, df_sp)
  total   <- sum(tab_ts)
  min_ts  <- if (length(tab_ts) > 0) min(tab_ts) else 0
  
  has_T   <- nlevels(df_sp$Type)       > 1
  has_S   <- nlevels(df_sp$life_stage) > 1
  has_P   <- nlevels(df_sp$Pair)       > 1
  has_CT  <- nlevels(df_sp$Count.Type) > 1
  
  rhs <- c("offset(log(Inclusion_m))",
           "(1 | Site)",
           "(1 | survey_id)")
  
  if (total < min_total) {
    # rare species, no time, simple structure
    if (has_T)  rhs <- c("Type", rhs)
    if (has_S)  rhs <- c("life_stage", rhs)
    if (has_P)  rhs <- c("Pair", rhs)
  } else {
    # allow time and Count.Type when they vary
    rhs <- c("date_num", rhs)
    if (has_CT) rhs <- c("Count.Type", rhs)
    if (has_T)  rhs <- c("Type", rhs)
    if (has_S)  rhs <- c("life_stage", rhs)
    if (has_P)  rhs <- c("Pair", rhs)
    if (has_T && has_S && min_ts >= min_cell) {
      rhs <- c("Type:life_stage", rhs)
    }
  }
  
  f_str <- paste("stage_Count ~", paste(rhs, collapse = " + "))
  as.formula(f_str)
}


### Helper to fit one species model ####


fit_species_model <- function(sp_name, data) {
  message("Fitting species model for: ", sp_name)
  
  df_sp <- data %>%
    filter(Sci_Name == sp_name) %>%
    droplevels() %>%
    mutate(
      Type       = factor(Type),
      life_stage = factor(life_stage),
      Pair       = factor(Pair),
      Count.Type = factor(Count.Type),
      Site       = factor(Site),
      survey_id  = factor(survey_id)
    )
  
  if (nrow(df_sp) == 0) {
    warning("No rows for species: ", sp_name)
    return(NULL)
  }
  
  f_final <- choose_species_formula(df_sp)
  
  m_sp <- glmmTMB(f_final, family = nbinom2, data = df_sp)
  
  list(
    species = sp_name,
    model   = m_sp,
    formula = f_final
  )
}


### Fit models for all key species ####

species_models <- key_species %>%
  map(~ fit_species_model(.x, fish_species_ls)) %>%
  compact()

### Model quality diagnostics ####
summarise_model_structure <- function(formula_obj) {
  f <- paste(deparse(formula_obj), collapse = " ")
  
  rhs <- strsplit(f, "~")[[1]][2]
  rhs <- gsub("\\s+", " ", rhs)
  
  # classify based on presence of interaction and life_stage
  has_int   <- grepl("Type:life_stage", rhs)
  has_stage <- grepl("life_stage", rhs)
  has_time  <- grepl("date_num", rhs)
  
  if (has_int) {
    desc <- "Type x life_stage interaction"
  } else if (has_stage) {
    desc <- "Additive Type + life_stage"
  } else {
    desc <- "Type only (life stage removed)"
  }
  
  if (!has_time) desc <- paste0(desc, " (rare species model)")
  
  desc
}
check_model_quality <- function(sp_obj) {
  sp   <- sp_obj$species
  mod  <- sp_obj$model
  fml  <- sp_obj$formula
  
  coefs <- summary(mod)$coefficients$cond
  
  tibble(
    species        = sp,
    model_type     = summarise_model_structure(fml),
    AIC            = AIC(mod),
    any_na_est     = any(!is.finite(coefs[, "Estimate"])),
    any_na_se      = any(!is.finite(coefs[, "Std. Error"])),
    ok_for_plots   = !any_na_est & !any_na_se & is.finite(AIC)
  )
}

model_diag <- species_models %>%
  map_dfr(check_model_quality)

print(model_diag, n=Inf)

### Save species models and summaries ####

walk(species_models, function(x) {
  sp_name <- x$species
  m_sp    <- x$model
  
  safe_name <- gsub(" ", "_", sp_name)
  
  saveRDS(
    m_sp,
    file.path(fits_dir, paste0("m_species_", safe_name, "_", analysis_date, ".rds"))
  )
  
  capture.output(
    summary(m_sp),
    file = file.path(stats_dir, paste0("m_species_", safe_name, "_summary_", analysis_date, ".txt"))
  )
})

### Plotting function for species models ####

make_species_plots <- function(species_models,
                               plots_dir,
                               analysis_date,
                               reef_cols,
                               theme_clean,
                               fish_long_norm) {
  
  # canonical bin order from size_bins
  bin_levels <- size_bins$Size_Class
  
  plot_list <- list()
  
  for (x in species_models) {
    
    sp_name   <- x$species
    m_sp      <- x$model
    safe_name <- gsub(" ", "_", sp_name)
    
    message("Plotting species: ", sp_name)
    
    ## 1) PREDICTED JUV / ADULT PANEL ########################################
    
    pred_sp <- try(
      ggpredict(
        m_sp,
        terms     = c("life_stage", "Type", "Pair"),
        condition = list(Inclusion_m = 100)
      ),
      silent = TRUE
    )
    
    if (inherits(pred_sp, "try-error")) {
      message("  ggpredict failed for ", sp_name, ": ",
              conditionMessage(attr(pred_sp, "condition")))
      next
    }
    
    pred_sp <- as.data.frame(pred_sp)
    
    if (nrow(pred_sp) == 0) {
      message("  Empty prediction grid for ", sp_name)
      next
    }
    
    if (!all(c("conf.low", "conf.high") %in% names(pred_sp))) {
      if ("std.error" %in% names(pred_sp)) {
        pred_sp <- pred_sp %>%
          mutate(
            conf.low  = predicted - 1.96 * std.error,
            conf.high = predicted + 1.96 * std.error
          )
        message("  CI rebuilt from std.error for ", sp_name)
      } else {
        pred_sp <- pred_sp %>%
          mutate(
            conf.low  = predicted,
            conf.high = predicted
          )
      }
    }
    
    pred_sp <- pred_sp %>%
      mutate(
        x     = factor(x,     levels = c("juvenile", "adult")),
        group = factor(group, levels = c("Natural", "Artificial"))
      )
    
    p_sp <- ggplot(
      pred_sp,
      aes(x = x, y = predicted, color = group, fill = group, group = group)
    ) +
      geom_errorbar(
        aes(ymin = conf.low, ymax = conf.high),
        width    = 0.15,
        position = position_dodge(width = 0.5)
      ) +
      geom_point(
        size     = 2.4,
        position = position_dodge(width = 0.5)
      ) +
      facet_wrap(~ facet, nrow = 1) +
      scale_color_manual(values = reef_cols) +
      scale_fill_manual(values = reef_cols) +
      theme_clean +
      labs(
        x     = "Life stage",
        y     = "Predicted count per 100 m²",
        color = "Reef type",
        fill  = "Reef type"
      ) +
      theme(
        plot.title       = element_text(face = "italic"),
        strip.background = element_blank()
      )
    
    ## 2) RAW SIZE STRUCTURE PANEL ###########################################
    
    size_tab <- fish_long_norm %>%
      filter(Sci_Name == sp_name) %>%
      group_by(Type, Size_Class) %>%
      summarise(
        total_n = sum(Count, na.rm = TRUE),
        .groups = "drop"
      )
    
    if (nrow(size_tab) == 0) {
      message("  No size data for ", sp_name, " skipping size panel")
      p_combined <- p_sp
    } else {
      
      # canonical bin order from size_bins
      bin_levels <- size_bins$Size_Class
      
      # join bin bounds, THEN enforce factor order
      size_tab <- size_tab %>%
        left_join(size_bins, by = "Size_Class") %>%
        mutate(
          Size_Class = factor(Size_Class, levels = bin_levels)
        ) %>%
        arrange(Type, Size_Class)
      
      # base histogram
      p_hist <- ggplot(size_tab,
                       aes(x = Size_Class, y = total_n, fill = Type)) +
        geom_col(position = "dodge") +
        scale_fill_manual(values = reef_cols) +
        labs(
          x = "Size class (cm)",
          y = "Total count"
        ) +
        theme_clean +
        theme(
          axis.text.x     = element_text(angle = 45, hjust = 1),
          legend.position = "none"
        )
      
      # look up Lmat
      Lmat_val <- maturity_lookup %>%
        filter(Sci_Name == sp_name) %>%
        pull(Lmat_cm)
      
      if (length(Lmat_val) == 1 && is.finite(Lmat_val)) {
        
        # bins object for Lmat logic, using same factor levels
        bins <- size_bins %>%
          mutate(
            Size_Class   = factor(Size_Class, levels = bin_levels),
            crosses_Lmat = lower < Lmat_val & upper > Lmat_val
          )
        
        idx <- which(bins$crosses_Lmat)
        if (length(idx) == 0) {
          if (Lmat_val <= min(bins$lower)) {
            idx <- 1L
          } else {
            idx <- nrow(bins)
          }
        }
        
        Lmat_bin <- bins$Size_Class[idx[1]]
        x_pos    <- which(bin_levels == as.character(Lmat_bin))
        y_max    <- max(size_tab$total_n, na.rm = TRUE)
        
        # label with actual size at maturity
        Lmat_lab <- paste0("Lmat = ", round(Lmat_val, 1), " cm")
        
        p_hist <- p_hist +
          geom_vline(
            xintercept = x_pos,
            linetype   = 2,
            linewidth  = 0.6,
            color      = "red"
          ) +
          annotate(
            "text",
            x     = x_pos + 1.0,
            y     = y_max * 1.05,
            label = Lmat_lab,
            color = "red",
            size  = 3,
            vjust = 0
          ) +
          expand_limits(y = y_max * 1.1)
      }
      
      ## 3) COMBINE PANELS (vertical) #######################################
      p_combined <- (p_sp / p_hist) + # use / for above/below or + for side to side but much change 
        plot_layout(heights = c(1.3, 1)) + # use for /
        # plot_layout(widths = c(2, 1)) + # use for + 
        plot_annotation(
          title = sp_name,
          theme = theme(plot.title = element_text(face = "italic"))
        )
    }
    
    
    ggsave(
      filename = file.path(
        plots_dir,
        paste0("Fig_species_", safe_name, "_", analysis_date, ".png")
      ),
      plot   = p_combined,
      width  = 7,
      height = 6,
      dpi    = 300
    )
    
    plot_list[[safe_name]] <- p_combined
  }
  
  invisible(plot_list)
}




### Restrict plots to good models and make plots ####

good_species_for_plots <- model_diag %>%
  filter(ok_for_plots) %>%
  pull(species)

species_models_plot <- species_models %>%
  keep(~ .x$species %in% good_species_for_plots)

plots_species <- make_species_plots(
  species_models = species_models_plot,
  plots_dir      = plots_dir,
  analysis_date  = analysis_date,
  reef_cols      = reef_cols,
  theme_clean    = theme_clean,
  fish_long_norm = fish_long_norm
)

## Optional check for a specific species
fish_species_ls %>%
   filter(Sci_Name == "Caesio xanthonota") %>%
   group_by(Type, life_stage, Pair) %>%
   summarise(n = sum(stage_Count), .groups = "drop") %>%
   pivot_wider(names_from = Pair, values_from = n, values_fill = 0)




##### Movement of preds by time ########## 
library(dplyr)
library(ggplot2)
library(patchwork)

## 1) Define the predator subset from spp_lookup ####

pred_spp <- spp_lookup %>%
  filter(grepl("Snapper|Grouper|Sweetlips", Species))

# Optional: check which species are included
pred_spp
# Snapper - Brownstripe / Mangrove / Russells
# Grouper - Blacktip / Brown marbled / Coral groupers (all)
# Sweetlips - Harlequin / Harry hotlips / Painted

# 2) filter fish_long_norm by those species  #### 
fish_pred <- fish_long_norm %>%
  filter(
    Species %in% pred_spp$Species,
    !is.na(Date)
  )


## 3) Time series: predator abundance through time on wrecks #### 

pred_time_df <- fish_pred %>%
  filter(
    Type == "Artificial",
    !is.na(Date)
  ) %>%
  group_by(Pair, Site, Date) %>%
  summarise(
    total_pred = sum(Count_perSurvey, na.rm = TRUE),
    .groups    = "drop"
  )

p_pred_time <- ggplot(pred_time_df,
                      aes(x = Date, y = total_pred, color = Site, group = Site)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ Pair, scales = "free_y") +
  labs(
    x = "Date",
    y = "Mean predator count per survey",
    color = "Site"
  ) +
  theme_clean

## 4) First-arrival dates for each predator species on each wreck ----------

pred_arrivals <- fish_pred %>%
  filter(
    Type == "Artificial",
    Count > 0,
    !is.na(Date)
  ) %>%
  group_by(Pair, Site, Sci_Name, Species) %>%
  summarise(
    first_date = min(Date, na.rm = TRUE),
    .groups    = "drop"
  )

p_arrivals <- ggplot(pred_arrivals,
                     aes(x = first_date, y = Species, color = Site)) +
  geom_point(size = 2) +
  facet_wrap(~ Pair) +
  labs(
    x = "First detection date",
    y = "Predator species",
    color = "Site"
  ) +
  theme_clean

## 5) Optional combined figure ---------------------------------------------

p_pred_time / p_arrivals + plot_layout(heights = c(2, 1))

