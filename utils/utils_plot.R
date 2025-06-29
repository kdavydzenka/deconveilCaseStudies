performance_plot <- function(plot_df,
                             method_colors,
                             x_label = "sample size",
                             y_label = "Performance metric",
                             x_labels = c("10", "20", "40", "100"),
                             y_limits = c(0.7, 1),
                             y_breaks = seq(0, 1, by = 0.1),
                             legend_position = "bottom") {
  
  ggplot(plot_df, aes(x = SampleSize, y = Mean, color = Method, group = Method)) +
    geom_line(size = 1.4) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.15, size = 0.7) +
    ggh4x::facet_nested(factor(GeneSize) ~ factor(Metric)) +
    labs(title = "", x = x_label, y = y_label) +
    scale_x_discrete(labels = x_labels) +   
    scale_y_continuous(limits = y_limits, breaks = y_breaks) + 
    scale_color_manual(values = method_colors) +
    theme_bw() +
    theme(
      strip.text = element_text(size = 12, face = "plain", color = "black"),       
      axis.title.x = element_text(size = 12, color = "black"),                     
      axis.title.y = element_text(size = 12, color = "black"),                     
      axis.text.x = element_text(size = 12),                      
      axis.text.y = element_text(size = 12),                      
      legend.text = element_text(size = 12, color = "black"),                      
      legend.title = element_text(size = 12, color = "black"),
      legend.position = legend_position
    )
}


plot_metrics_robustness <- function(df, 
                                    metric_col, 
                                    y_label, 
                                    y_lim = c(0, 1), 
                                    noise_labels = c("noise_10" = "10", "noise_15" = "15", "noise_20" = "20", "noise_25" = "25"),
                                    sample_labels = c("10" = "# samples 10", "20" = "# samples 20", "40" = "# samples 40", "60" = "# samples 60"),
                                    gene_group_colors = c("DSG" = "#D55E00", "DIG" = "#771C19", "DCG" = "#FFC300", "non_DEG" = "#888888")) {
  
  mean_col <- paste0("mean_", metric_col)
  sd_col <- paste0("sd_", metric_col)
  
  # Ensure factor levels are correct
  df <- df %>%
    mutate(
      noise_level = factor(noise_level, levels = names(noise_labels)),
      sample_size = factor(sample_size, levels = names(sample_labels)),
      group_clean = factor(group_clean, levels = names(gene_group_colors))
    )
  
  p <- ggplot(df, aes(x = noise_level, y = .data[[mean_col]], color = group_clean, group = group_clean)) +
    geom_line(size = 1.4) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = .data[[mean_col]] - .data[[sd_col]],
                      ymax = .data[[mean_col]] + .data[[sd_col]]),
                  width = 0.2, size = 0.6) +
    facet_wrap(~ sample_size, nrow = 1, labeller = as_labeller(sample_labels)) +
    labs(
      title = "",
      x = "Noise level, %",
      y = y_label,
      color = "Gene Group"
    ) +
    scale_x_discrete(labels = noise_labels) +
    scale_color_manual(values = gene_group_colors) +
    ylim(y_lim) +
    theme_bw() +
    theme(
      strip.text = element_text(size = 12, face = "plain"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12),
      legend.position = "right"
    )
  
  return(p)
}


forest_plot <- function(sel_genes_data, title = NULL) {
  plot_data <- sel_genes_data %>%
    mutate(
      gene = factor(Gene, levels = rev(unique(Gene))),
      is_summary = FALSE,
      box_color = ifelse(HR > 1, "#D60C00FF", "blue3")
    )
  
  plot_data$is_summary[1] <- TRUE  
  
  ggplot(plot_data, aes(x = HR, y = gene)) +
    geom_pointrange(aes(xmin = CI_lower, xmax = CI_upper, color = box_color),
                    size = 0.9, show.legend = FALSE) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "black") +
    labs(
      x = "Hazard Ratio",
      y = NULL,
      title = title
    ) +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 15, color = "black"),
      axis.text.x = element_text(size = 14, color = "black"),
      axis.title.x = element_text(size = 15),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank()
    ) +
    scale_color_identity() +
    coord_cartesian(xlim = c(min(plot_data$CI_lower), max(plot_data$CI_upper))) +
    scale_x_continuous(trans = "identity")
}


venn_plot <- function(data, group_name) {
  title_map <- c(
    d_sensitive = "DSGs",      
    d_insensitive = "DIGs",    
    d_compensated = "DCGs"     
  )
  
  # Prepare the wide-format presence matrix
  venn_data <- data %>%
    dplyr::mutate(present = 1) %>%
    tidyr::pivot_wider(
      names_from = tumor_type,
      values_from = present,
      values_fill = list(present = 0)
    ) %>%
    distinct(geneID, .keep_all = TRUE) %>%
    dplyr::select(geneID, any_of(c("LUAD", "LUSC", "BRCA"))) 
  
  # Generate list of gene IDs per tumor type
  venn_list <- lapply(c("LUAD", "LUSC", "BRCA"), function(tt) {
    venn_data %>% filter(.data[[tt]] == 1) %>% pull(geneID)
  })
  names(venn_list) <- c("LUAD", "LUSC", "BRCA")
  
  # Plot
  ggVennDiagram(venn_list) +
    scale_fill_gradient(low = "lightgray", high = "#B2474599") +
    theme_minimal() +
    labs(title = title_map[[group_name]]) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.background = element_rect(fill = "white", color = NA),
      legend.position = "right"
    )
}


enrichment_dotplot <- function(dotplot_data) {
  # Ensure proper factor levels if not already set
  dotplot_data <- dotplot_data %>%
    mutate(
      Tumor = factor(Tumor, levels = c("BRCA", "LUAD", "LUSC")),
      Category = factor(Category, levels = c("DSGs", "DIGs", "DCGs")),
      Description = factor(Description, levels = unique(Description))
    )
  
  ggplot(dotplot_data, aes(x = Tumor, y = Description)) +
    geom_point(aes(size = Count, color = log_padjust)) +
    facet_wrap(~Category, scales = "free", ncol = 3) +
    scale_color_gradient(
      low = "mediumpurple",
      high = "coral1",
      name = expression(-log[10](p.adjust))
    ) +
    scale_size_continuous(name = "Gene Count", range = c(2, 6)) +
    scale_x_discrete(name = "Tumor type") +
    scale_y_discrete(name = "Biological Process GO term") +
    theme_bw(base_size = 12) +
    theme(
      panel.spacing = unit(1, "lines"),
      axis.text.x = element_text(angle = 30, hjust = 1),
      legend.key.size = unit(0.6, "cm"),
      legend.spacing.y = unit(0.2, 'cm'),
      strip.text = element_text(face = "plain")
    )
}


gene_barplot <- function(distribution_summary, gene_type_colors) {
  # Ensure gene_category is ordered correctly
  distribution_summary <- distribution_summary %>%
    mutate(gene_category = factor(gene_category, levels = c("DSGs", "DIGs", "DCGs")))
  
  ggplot(distribution_summary, aes(x = gene_category, y = gene_count, fill = gene_type)) +
    geom_bar(stat = "identity") +
    geom_text(
      aes(label = gene_count),
      position = position_stack(vjust = 0.5),
      size = 4,
      color = "white"
    ) +
    scale_fill_manual(values = gene_type_colors) +
    labs(
      x = "Gene category",
      y = "Gene count",
      fill = "Gene type"
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text = element_text(color = "black"),
      axis.title = element_text(size = 12),
      legend.position = "right",
      legend.key.size = unit(0.6, "cm"),
      legend.spacing.y = unit(0.2, "cm"),
      legend.title = element_text(color = "black"),
      legend.text = element_text(color = "black")
    )
}

sankey_plot <- function(data,
                        x_var = "x",
                        id_var = "id",
                        split_var = "y",
                        value_var = "freq",
                        fill_var = "CN_naive",
                        group_colors = NULL,
                        x_labels = c("CN-naive", "CN-aware"),
                        title = NULL) {
  
  ggplot(data, aes_string(x = x_var, id = id_var, split = split_var, value = value_var)) +
    geom_parallel_sets(aes_string(fill = fill_var), alpha = 0.9, axis.width = 0.2,
                       n = 4415, strength = 0.5, color = "black", linewidth = 0.3) +
    geom_parallel_sets_axes(axis.width = 0.25, fill = "gray93",
                            color = "gray", linewidth = 0.5) +
    
    # geom_parallel_sets_labels(colour = 'gray35', size = 2.0, angle = 0, fontface = "plain") +
    scale_fill_manual(values = group_colors, name = "Gene group") +
    scale_color_manual(values = group_colors) +
    scale_x_continuous(breaks = seq_along(x_labels), labels = x_labels) +
    theme_classic() +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 15, face = "plain"),
      legend.text = element_text(size = 13),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank()
    ) +
    ggtitle(title)
}
