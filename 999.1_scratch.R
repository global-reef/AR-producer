# ---- MAKE SPECIES PANELS (updated per your notes) ----
make_species_plots <- function(species_models,
                               plots_dir,
                               analysis_date,
                               reef_cols,
                               theme_clean,
                               fish_long_norm,
                               bottom_row_species = NULL) {
  
  bin_levels <- as.character(size_bins$Size_Class)
  
  size_bins_join <- size_bins %>%
    mutate(Size_Class = ordered(as.character(Size_Class), levels = bin_levels))
  
  plot_list <- list()
  
  for (x in species_models) {
    
    sp_name   <- x$species
    m_sp      <- x$model
    safe_name <- gsub(" ", "_", sp_name)
    
    ## 1) PREDICTED JUV / ADULT PANEL ########################################
    
    pred_sp <- try(
      ggpredict(
        m_sp,
        terms     = c("life_stage", "Type", "Pair"),
        condition = list(Inclusion_m = 100)
      ),
      silent = TRUE
    )
    if (inherits(pred_sp, "try-error")) next
    
    pred_sp <- as.data.frame(pred_sp)
    if (nrow(pred_sp) == 0) next
    
    if (!all(c("conf.low", "conf.high") %in% names(pred_sp))) {
      if ("std.error" %in% names(pred_sp)) {
        pred_sp <- pred_sp %>%
          mutate(
            conf.low  = predicted - 1.96 * std.error,
            conf.high = predicted + 1.96 * std.error
          )
      } else {
        pred_sp <- pred_sp %>%
          mutate(conf.low = predicted, conf.high = predicted)
      }
    }
    
    pred_sp <- pred_sp %>%
      mutate(
        x     = factor(x, levels = c("juvenile", "adult")),
        group = factor(group, levels = c("Natural", "Artificial"))
      )
    
    p_sp <- ggplot(pred_sp, aes(x = x, y = predicted, color = group, group = group)) +
      geom_errorbar(
        aes(ymin = conf.low, ymax = conf.high),
        width = 0.14,
        position = position_dodge(width = 0.5)
      ) +
      geom_point(
        size = 2.8,
        position = position_dodge(width = 0.5)
      ) +
      facet_wrap(~ facet, nrow = 1) +
      scale_color_manual(values = reef_cols, guide = "none") +
      theme_clean +
      labs(
        x = NULL,
        y = expression("Predicted count per 100 m"^2)
      ) +
      theme(
        strip.background = element_blank(),
        strip.text       = element_text(size = 10),
        axis.title.y     = element_text(size = 10),
        axis.text        = element_text(size = 9),
        plot.margin      = margin(4, 4, 4, 4)
      )
    
    if (!is.null(bottom_row_species) && !(sp_name %in% bottom_row_species)) {
      p_sp <- p_sp + theme(axis.text.x = element_blank())
    }
    
    ## 2) RAW SIZE STRUCTURE PANEL ###########################################
    
    size_tab <- fish_long_norm %>%
      filter(Sci_Name == sp_name) %>%
      group_by(Type, Size_Class) %>%
      summarise(total_n = sum(Count, na.rm = TRUE), .groups = "drop") %>%
      mutate(Size_Class = ordered(as.character(Size_Class), levels = bin_levels)) %>%
      tidyr::complete(
        Type,
        Size_Class = ordered(bin_levels, levels = bin_levels),
        fill = list(total_n = 0)
      ) %>%
      left_join(size_bins_join, by = "Size_Class") %>%
      arrange(Type, Size_Class)
    
    p_hist <- ggplot(size_tab, aes(x = Size_Class, y = total_n, fill = Type)) +
      geom_col(position = position_dodge(width = 0.85), width = 0.8) +
      scale_fill_manual(values = reef_cols, name = "Reef type") +
      labs(x = NULL, y = "Total count") +
      theme_clean +
      theme(
        axis.text.x     = element_text(angle = 35, hjust = 1, size = 8),
        axis.text.y     = element_text(size = 9),
        axis.title.y    = element_text(size = 10),
        legend.position = "bottom",
        plot.margin     = margin(4, 4, 4, 4)
      )
    
    if (!is.null(bottom_row_species) && !(sp_name %in% bottom_row_species)) {
      p_hist <- p_hist + theme(axis.text.x = element_blank())
    }
    
    Lmat_val <- maturity_lookup %>%
      filter(Sci_Name == sp_name) %>%
      pull(Lmat_cm)
    
    if (length(Lmat_val) == 1 && is.finite(Lmat_val)) {
      
      bins <- size_bins_join %>%
        mutate(crosses_Lmat = lower < Lmat_val & upper > Lmat_val)
      
      idx <- which(bins$crosses_Lmat)
      if (length(idx) == 0) idx <- if (Lmat_val <= min(bins$lower)) 1L else nrow(bins)
      
      Lmat_bin <- bins$Size_Class[idx[1]]
      x_pos    <- which(bin_levels == as.character(Lmat_bin))
      Lmat_lab <- paste0("Lmat = ", round(Lmat_val, 1), " cm")
      
      y_max <- max(size_tab$total_n, na.rm = TRUE)
      if (!is.finite(y_max) || y_max <= 0) y_max <- 1
      
      if (x_pos >= (length(bin_levels) - 1)) {
        x_lab <- x_pos - 0.35
        hjust_lab <- 1
      } else {
        x_lab <- x_pos + 0.35
        hjust_lab <- 0
      }
      
      p_hist <- p_hist +
        geom_vline(xintercept = x_pos, linetype = 2, linewidth = 0.6, color = "red") +
        annotate(
          "text",
          x = x_lab,
          y = y_max * 1.05,
          label = Lmat_lab,
          color = "red",
          size = 3.1,
          hjust = hjust_lab
        ) +
        expand_limits(y = y_max * 1.12) +
        coord_cartesian(clip = "off")
    }
    
    ## 3) COMBINE (title strip smaller + tighter) #############################
    
    p_title <- ggplot() +
      theme_void() +
      labs(title = sp_name) +
      theme(
        plot.title  = element_text(face = "italic", size = 10.5, hjust = 0, margin = margin(b = 0)),
        plot.margin = margin(0, 3, 0, 3)
      )
    
    p_body <- (p_sp + p_hist) + plot_layout(widths = c(2.2, 1))
    
    p_hor <- (p_title / p_body) +
      plot_layout(heights = c(0.07, 1))  # less whitespace above panels
    p_hor <- patchwork::wrap_elements(full = p_hor)
    
    ggsave(
      filename = file.path(plots_dir, paste0("Fig_species_hor_", safe_name, "_", analysis_date, ".png")),
      plot   = p_hor,
      width  = 7.4,
      height = 3.6,
      dpi    = 300
    )
    
    plot_list[[safe_name]] <- p_hor
  }
  
  invisible(plot_list)
}

# ---- 6-species composite ----

bottom_row_spp <- focal_spp

plots_species6 <- make_species_plots(
  species_models     = species_models_plot6,
  plots_dir          = plots_dir,
  analysis_date      = analysis_date,
  reef_cols          = reef_cols,
  theme_clean        = theme_clean,
  fish_long_norm     = fish_long_norm,
  bottom_row_species = bottom_row_spp
)

# Only 6 tags (a–f), no overall title/subtitle
fig_spp6 <- wrap_plots(plots_species6, ncol = 2, guides = "collect", tag_level = "new") +
  plot_annotation(tag_levels = "a") &
  theme(
    legend.position   = "bottom",
    legend.title      = element_text(size = 10),
    legend.text       = element_text(size = 9),
    plot.tag          = element_text(face = "bold", size = 13),
    plot.tag.position = c(0.01, 0.99),
    plot.margin       = margin(6, 6, 6, 6)
  )

fig_spp6
