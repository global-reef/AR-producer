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
  library(tidyverse)
})
### Prepare species-level dataset ####

# We use the probabilistic life stage data to keep the juvenile vs adult counts
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
      total_count >= 30 &
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
      
      
      ## 3) COMBINE PANELS (hortizontal) #######################################
      p_hor <- (p_sp + p_hist) + # use / for above/below or + for side to side but much change 
        # plot_layout(heights = c(1.3, 1)) + # use for /
       plot_layout(widths = c(2, 1)) + # use for + 
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
    
    ggsave(
      filename = file.path(
        plots_dir,
        paste0("Fig_species_hor_", safe_name, "_", analysis_date, ".png")
      ),
      plot   = p_hor,
      width  = 7,
      height = 4,
      dpi    = 300
    )
    
    plot_list[[safe_name]] <- p_hor
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



#### new plots and summary tables ##### 

### Table S7: Species-level reef-type contrasts (AR:NR) #########################

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tibble)
  library(stringr)
  library(emmeans)
})

### Focal species list (only those that are ok_for_plots) #######################

focal_spp <- c(
  "Diagramma pictum",
  "Lutjanus russellii",
  "Lutjanus argentimaculatus",
  "Neopomacentrus cyanomos",
  "Epinephelus fasciatus",
  "Siganus javus"
)

### Helper: detect model structure #############################################

detect_model_structure <- function(formula_obj) {
  ftxt <- paste(deparse(formula_obj), collapse = " ")
  has_int   <- grepl("Type:life_stage", ftxt, fixed = TRUE)
  has_stage <- grepl("life_stage", ftxt)
  if (has_int) "Type x life_stage"
  else if (has_stage) "Type + life_stage"
  else "Type only"
}

### Helper: compute AR:NR ratio with CI and p ##################################

### Robust extractor: always returns ratio + CI + p ################################

extract_type_ratio <- function(mod, fml, prefer_stage = "juvenile", inclusion_m = 100) {
  
  ftxt      <- paste(deparse(fml), collapse = " ")
  has_int   <- grepl("Type:life_stage", ftxt, fixed = TRUE)
  has_stage <- grepl("life_stage", ftxt)
  
  off_val <- log(inclusion_m)
  
  # helper to compute AR - NR on link scale
  link_contrast <- function(emm_df) {
    
    # ensure ordering
    emm_df <- emm_df %>% dplyr::arrange(Type)
    
    # Artificial - Natural on log scale
    diff_est <- emm_df$emmean[emm_df$Type == "Artificial"] -
      emm_df$emmean[emm_df$Type == "Natural"]
    
    # SE of difference
    se_diff <- sqrt(
      emm_df$SE[emm_df$Type == "Artificial"]^2 +
        emm_df$SE[emm_df$Type == "Natural"]^2
    )
    
    z <- qnorm(0.975)
    
    tibble::tibble(
      ratio    = exp(diff_est),
      conf.low = exp(diff_est - z * se_diff),
      conf.high= exp(diff_est + z * se_diff),
      p.value  = 2 * pnorm(-abs(diff_est / se_diff))
    )
  }
  
  if (has_int && has_stage) {
    
    emm <- emmeans::emmeans(mod, ~ Type | life_stage, offset = off_val) %>%
      as.data.frame() %>%
      dplyr::mutate(life_stage = as.character(life_stage))
    
    if (nrow(emm) == 0) {
      return(tibble::tibble(
        contrast_reported = NA_character_,
        ratio = NA_real_, conf.low = NA_real_, conf.high = NA_real_, p.value = NA_real_
      ))
    }
    
    # choose juvenile if present
    if (prefer_stage %in% emm$life_stage) {
      emm_sub <- emm %>% dplyr::filter(life_stage == prefer_stage)
      stage   <- prefer_stage
    } else {
      emm_sub <- emm %>% dplyr::filter(life_stage == emm$life_stage[1])
      stage   <- emm_sub$life_stage[1]
    }
    
    out <- link_contrast(emm_sub)
    out$contrast_reported <- stage
    out %>% dplyr::select(contrast_reported, ratio, conf.low, conf.high, p.value)
    
  } else {
    
    emm <- emmeans::emmeans(mod, ~ Type, offset = off_val) %>%
      as.data.frame()
    
    out <- link_contrast(emm)
    out$contrast_reported <- "Marginal"
    out %>% dplyr::select(contrast_reported, ratio, conf.low, conf.high, p.value)
  }
}



### Build Table S7 #############################################################

# species_models is your list of lists:
# list(species = sp_name, model = m_sp, formula = f_final)
table_s7 <- species_models %>%
  # keep(~ .x$species %in% focal_spp) %>%
  map_dfr(function(x) {
    
    sp  <- x$species
    mod <- x$model
    fml <- x$formula
    
    out <- try(
      extract_type_ratio(mod, fml, prefer_stage = "juvenile", inclusion_m = 100),
      silent = TRUE
    )
    
    if (inherits(out, "try-error") || is.null(out) || nrow(out) == 0) {
      return(tibble(
        Species = sp,
        Model_structure = detect_model_structure(fml),
        Contrast = NA_character_,
        ar_nr_ratio = NA_real_,
        ci_low = NA_real_,
        ci_high = NA_real_,
        p_value = NA_real_,
        note = "contrast failed"
      ))
    }
    
    tibble(
      Species = sp,
      Model_structure = detect_model_structure(fml),
      Contrast = out$contrast_reported,
      ar_nr_ratio = as.numeric(out$ratio),
      ci_low  = as.numeric(out$conf.low),
      ci_high = as.numeric(out$conf.high),
      p_value = as.numeric(out$p.value),
      note = NA_character_
    )
  }) %>%
  mutate(
    ar_nr_ratio = round(ar_nr_ratio, 2),
    ci_low      = round(ci_low, 2),
    ci_high     = round(ci_high, 2),
    p_value_fmt = case_when(
      is.na(p_value)       ~ NA_character_,
      p_value < 0.001      ~ "<0.001",
      TRUE                 ~ formatC(p_value, format = "f", digits = 3)
    )
  ) %>%
  select(Species, Model_structure, Contrast,
         ar_nr_ratio, ci_low, ci_high, p_value_fmt, -note) %>%
  rename(
    `AR:NR ratio` = ar_nr_ratio,
    `95% CI low`  = ci_low,
    `95% CI high` = ci_high,
    `p-value`     = p_value_fmt
  )
table_s7_clean <- table_s7 %>%
  mutate(
    `95% CI low`  = ifelse(is.na(`95% CI low`) | is.nan(`95% CI low`), NA, `95% CI low`),
    `95% CI high` = ifelse(is.na(`95% CI high`) | is.nan(`95% CI high`), NA, `95% CI high`),
    `p-value`     = ifelse(is.na(`p-value`), NA, `p-value`)
  ) %>%
  mutate(
    `95% CI` = ifelse(
      is.na(`95% CI low`),
      "n.e.",
      paste0("[", `95% CI low`, ", ", `95% CI high`, "]")
    ),
    `p-value` = ifelse(is.na(`p-value`), "n.e.", `p-value`)
  ) %>%
  select(Species, Model_structure, Contrast, `AR:NR ratio`, `95% CI`, `p-value`)

print(table_s7_clean, n = Inf)

# Save for manuscript assembly
write.csv(table_s7, file.path(results_dir, paste0("S7_species_contrasts_", analysis_date, ".csv")), row.names = FALSE)




### Figure S3: Species case studies (6 spp, two-panel each, ncol = 2) ###########

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(ggeffects)
  library(ggplot2)
  library(patchwork)
})

theme_clean2 <- theme_minimal(base_family = "Arial") +
  theme(
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    plot.title        = element_text(face = "italic", size = 10, margin = margin(0, 0, 2, 0)),
    panel.grid = element_blank(),
    plot.margin       = margin(3, 10, 3, 3)
  )



# focal species (order matters)
focal_spp <- c(
  "Neopomacentrus cyanomos",
  "Siganus javus",
  "Epinephelus fasciatus",
  "Diagramma pictum",
  "Lutjanus russellii",
  "Lutjanus argentimaculatus"
)

# subset model objects (only these 6)
species_models_plot6 <- species_models %>%
  keep(~ .x$species %in% focal_spp)


# --- One combined panel per species (pred left, size right) -------------------
plot_species_panel <- function(sp_obj,
                               reef_cols,
                               theme_clean2,
                               fish_long_norm,
                               size_bins,
                               maturity_lookup,
                               inclusion_m = 100) {
  
  sp_name <- sp_obj$species
  m_sp    <- sp_obj$model
  
  # ---- predictions ----
  pred_sp <- try(
    ggeffects::ggpredict(
      m_sp,
      terms     = c("life_stage", "Type", "Pair"),
      condition = list(Inclusion_m = inclusion_m)
    ),
    silent = TRUE
  )
  if (inherits(pred_sp, "try-error")) return(NULL)
  
  pred_sp <- as.data.frame(pred_sp)
  if (nrow(pred_sp) == 0) return(NULL)
  
  if (!all(c("conf.low", "conf.high") %in% names(pred_sp))) {
    if ("std.error" %in% names(pred_sp)) {
      pred_sp <- pred_sp %>%
        mutate(conf.low = predicted - 1.96 * std.error,
               conf.high = predicted + 1.96 * std.error)
    } else {
      pred_sp <- pred_sp %>% mutate(conf.low = predicted, conf.high = predicted)
    }
  }
  
  pred_sp <- pred_sp %>%
    mutate(
      life_stage = factor(x, levels = c("juvenile", "adult")),
      Type       = factor(group, levels = c("Natural", "Artificial"))
    )
  
  p_sp <- ggplot(pred_sp, aes(x = life_stage, y = predicted, fill = Type, group = Type)) +
    geom_errorbar(
      aes(ymin = conf.low, ymax = conf.high, colour = Type),
      width = 0.15,
      position = position_dodge(width = 0.5),
      show.legend = FALSE
    ) +
    geom_point(
      shape = 21,
      size = 2.6,
      stroke = 0.25,
      position = position_dodge(width = 0.5),
      show.legend = FALSE
    ) +
    facet_wrap(~ facet, nrow = 1) +
    scale_fill_manual(values = reef_cols, name = "Reef type") +
    scale_colour_manual(values = reef_cols, guide = "none") +
    theme_clean2 +
    labs(x = NULL, y = expression("Predicted count per 100 m"^2)) +
    theme(
      plot.margin = margin(6, 6, 4, 4)  # give the species title room
    )
  p_sp <- p_sp +
    labs(title = sp_name) +
    theme(
      plot.title = element_text(face = "italic", size = 11, hjust = 0.5, margin = margin(0, 0, 6, 0)),
      plot.margin = margin(14, 6, 4, 4)
    )
  
  # ---- size structure ----
  size_tab <- fish_long_norm %>%
    filter(Sci_Name == sp_name) %>%
    group_by(Type, Size_Class) %>%
    summarise(total_n = sum(Count, na.rm = TRUE), .groups = "drop")
  
  if (nrow(size_tab) == 0) {
    return(
      p_sp +
        patchwork::plot_annotation(
          title = sp_name,
          theme = theme(plot.title = element_text(face = "italic", size = 10),
                        plot.margin = margin(10, 6, 4, 4))
        )
    )
  }
  
  bin_levels <- size_bins$Size_Class
  keep_labs  <- c("0-1", "2-5", "10-15", "20-50", "100+")
  
  size_tab <- size_tab %>%
    left_join(size_bins, by = "Size_Class") %>%
    mutate(
      Size_Class = factor(Size_Class, levels = bin_levels),
      Type       = factor(Type, levels = c("Natural", "Artificial"))
    )
  
  p_hist <- ggplot(size_tab, aes(x = Size_Class, y = total_n, fill = Type)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = reef_cols, name = "Reef type") +
    scale_x_discrete(labels = function(x) ifelse(x %in% keep_labs, x, "")) +
    theme_clean2 +
    labs(x = NULL, y = "Total count") +
    theme(
      legend.position = "none",
      axis.text.x     = element_text(angle = 45, hjust = 1, size = 7)
    )
  
  # Lmat
  Lmat_val <- maturity_lookup %>% filter(Sci_Name == sp_name) %>% pull(Lmat_cm)
  
  if (length(Lmat_val) == 1 && is.finite(Lmat_val)) {
    
    bins <- size_bins %>%
      mutate(
        Size_Class   = factor(Size_Class, levels = bin_levels),
        crosses_Lmat = lower < Lmat_val & upper > Lmat_val
      )
    
    idx <- which(bins$crosses_Lmat)
    if (length(idx) == 0) idx <- if (Lmat_val <= min(bins$lower)) 1L else nrow(bins)
    
    x_pos <- which(bin_levels == as.character(bins$Size_Class[idx[1]]))
    y_max <- max(size_tab$total_n, na.rm = TRUE)
    
    p_hist <- p_hist +
      geom_vline(xintercept = x_pos, linetype = 2, linewidth = 0.6) +
      expand_limits(y = y_max * 1.10) +
      labs(title = paste0("Lmat = ", round(Lmat_val, 1), " cm")) +
      theme(
        plot.title = element_text(size = 8, hjust = 1, margin = margin(0, 0, 2, 0)),
        plot.margin = margin(6, 6, 4, 4)  # you can usually drop the huge right margin now
      )
  } else {
    p_hist <- p_hist + labs(title = NULL)
  }
  
  
  # ---- combine (NO guide collection here) ----
  (p_sp + p_hist) +
    plot_layout(widths = c(2.2, 1)) +
    plot_annotation(
      title = sp_name,
      tag_levels = NULL,   # <- CRITICAL: stop inner tagging
      theme = theme(
        plot.title  = element_text(
          face = "italic",
          size = 11,
          hjust = 0.5,
          margin = margin(0, 0, 6, 0)
        ),
        plot.margin = margin(14, 6, 6, 6)
      )
    )
}


# Build panels in the focal order (no Ballistoides etc)
plots_species6 <- map(focal_spp, function(sp) {
  sp_obj <- species_models_plot6 %>% keep(~ .x$species == sp) %>% pluck(1)
  plot_species_panel(
    sp_obj          = sp_obj,
    reef_cols       = reef_cols,
    theme_clean2    = theme_clean2,
    fish_long_norm  = fish_long_norm,
    size_bins       = size_bins,
    maturity_lookup = maturity_lookup,
    inclusion_m     = 100
  )
})
names(plots_species6) <- focal_spp
plots_species6 <- compact(plots_species6)

# Assemble S3 with 6 tags only (a–f)
fig_s3 <- wrap_plots(plots_species6, ncol = 2) +
  plot_layout(guides = "collect")  &
  theme(
    legend.position = "bottom",
    plot.tag = element_text(face = "bold", size = 12),
    plot.tag.position = c(0, 1)
  )

fig_s3

ggsave(
  filename = file.path(plots_dir, paste0("Fig_9_species_case_studies_6spp_", analysis_date, ".png")),
  plot     = fig_s3,
  width = 12, height = 8.45,
  dpi      = 300
)


### size structure summary for species #### 
make_table_S8_size_maturity <- function(fish_long_norm,
                                        maturity_lookup) {
  
  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(tibble)
  })
  
  # size-class midpoints (cm)
  size_midpoints <- tibble(
    Size_Class_f = factor(
      c("0-1","1-2","2-5","5-10","10-15","15-20","20-50","50-100",">100"),
      levels = c("0-1","1-2","2-5","5-10","10-15","15-20","20-50","50-100",">100")
    ),
    Size_mid_cm = c(0.5, 1.5, 3.5, 7.5, 12.5, 17.5, 35, 75, 125)
  )
  
  # expand to individual fish sizes
  fish_sizes <- fish_long_norm %>%
    filter(!is.na(Size_Class_f), !is.na(Count), Count > 0) %>%
    left_join(size_midpoints, by = "Size_Class_f") %>%
    tidyr::uncount(weights = Count)
  
  # join maturity info
  fish_sizes <- fish_sizes %>%
    left_join(
      maturity_lookup %>% select(Sci_Name, Lmat_cm),
      by = "Sci_Name"
    ) %>%
    mutate(mature = Size_mid_cm >= Lmat_cm)
  
  # summarise by species × reef type
  tab <- fish_sizes %>%
    group_by(Sci_Name, Type) %>%
    summarise(
      median_size = median(Size_mid_cm, na.rm = TRUE),
      pct_mature  = mean(mature, na.rm = TRUE) * 100,
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from  = Type,
      values_from = c(median_size, pct_mature),
      names_sep   = "_"
    ) %>%
    left_join(
      maturity_lookup %>% select(Sci_Name, Lmat_cm),
      by = "Sci_Name"
    ) %>%
    relocate(Lmat_cm, .after = Sci_Name) %>%
    arrange(Sci_Name)
  
  # round for presentation
  tab %>%
    mutate(
      Lmat_cm = round(Lmat_cm, 1),
      across(starts_with("median_size"), ~ round(.x, 1)),
      across(starts_with("pct_mature"),  ~ round(.x, 1))
    )
}
table_s8 <- make_table_S8_size_maturity(
  fish_long_norm   = fish_long_norm,
  maturity_lookup  = maturity_lookup
) %>%
dplyr::rename(
  `L-mat (cm)`          = Lmat_cm,
  `AR median size (cm)` = median_size_Artificial,
  `NR median size (cm)` = median_size_Natural,
  `AR ≥ L-mat (%)`      = pct_mature_Artificial,
  `NR ≥ L-mat (%)`      = pct_mature_Natural
)

print(table_s8, n = Inf)


# Save for manuscript assembly
write.csv(table_s8, file.path(results_dir, paste0("S8_species_mat_", analysis_date, ".csv")), row.names = FALSE)
