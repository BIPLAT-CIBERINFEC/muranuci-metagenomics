# ===================================================
# FUNCIONES COMUNES PARA GRÁFICOS TAXONÓMICOS
# MURANUCI
# ===================================================

# ---------------------------------------------------
# Paleta fija usada en los plots originales
# 15 colores para top15 + gris para Others
# ---------------------------------------------------
fixed_taxonomy_palette <- c(
  "#E6194B", "#3CB44B", "#FFE119", "#4363D8",
  "#F58231", "#911EB4", "#42D4F4", "#FF69B4",
  "#BFEF45", "#FABEBE", "#808000", "#E6BEFF",
  "#9A6324", "#FFFAC8", "#800000"
)

others_colour <- "#A9A9A9"

# Formato científico para etiquetas taxonómicas de género/especie.
# Las categorías agregadas permanecen en redonda y, en el complejo de
# E. cloacae, solo el nombre taxonómico aparece en cursiva.
taxon_plotmath_labels <- function(x) {
  x <- as.character(x)
  aggregate_labels <- c(
    "Others", "Other NFGNB", "Gram-positive bacteria",
    "Other Enterobacterales"
  )

  label_code <- vapply(x, function(label) {
    escaped <- gsub("'", "\\\\'", label, fixed = TRUE)

    if (label %in% aggregate_labels) {
      return(paste0("'", escaped, "'"))
    }
    if (label == "Enterobacter cloacae complex") {
      return("italic('Enterobacter cloacae')~'complex'")
    }
    paste0("italic('", escaped, "')")
  }, character(1))

  parse(text = label_code)
}

taxon_plotmath_strings <- function(x) {
  vapply(as.character(x), function(label) {
    escaped <- gsub("'", "\\\\'", label, fixed = TRUE)
    if (label == "Enterobacter cloacae complex") {
      return("italic('Enterobacter cloacae')~'complex'")
    }
    paste0("italic('", escaped, "')")
  }, character(1))
}


# ---------------------------------------------------
# Crear paleta fija con nombres de taxones
# Mantiene Others siempre en gris y al final
# ---------------------------------------------------
make_fixed_cols <- function(taxa, others_label = "Others") {
  
  taxa <- unique(as.character(taxa))
  taxa_no_others <- setdiff(taxa, others_label)
  
  if (length(taxa_no_others) <= length(fixed_taxonomy_palette)) {
    cols_taxa <- fixed_taxonomy_palette[seq_along(taxa_no_others)]
  } else {
    cols_taxa <- grDevices::colorRampPalette(fixed_taxonomy_palette)(
      length(taxa_no_others)
    )
  }
  
  cols <- setNames(cols_taxa, taxa_no_others)
  cols <- c(cols, setNames(others_colour, others_label))
  
  return(cols)
}


# ---------------------------------------------------
# Normalizar IDs de paciente y semana para figuras
# Usa Patient = PAC1, PAC2...
# Usa follow_up_week = W0, W1, W2
# ---------------------------------------------------
prepare_patient_week_labels <- function(tabla) {
  
  if (!"follow_up_week" %in% names(tabla)) {
    tabla <- tabla %>%
      mutate(follow_up_week = paste0("W", Week))
  }
  
  tabla <- tabla %>%
    mutate(
      Patient = as.character(Patient),
      Patient = factor(Patient, levels = paste0("PAC", 1:8)),
      
      follow_up_week = as.character(follow_up_week),
      follow_up_week = if_else(
        stringr::str_detect(follow_up_week, "^W"),
        follow_up_week,
        paste0("W", follow_up_week)
      ),
      follow_up_week = factor(follow_up_week, levels = c("W0", "W1", "W2"))
    )
  
  return(tabla)
}


# ---------------------------------------------------
# Tema común para todos los gráficos taxonómicos
# ---------------------------------------------------
theme_taxonomy <- function(base_size = 16) {
  theme_minimal(base_size = base_size) +
    theme(
      axis.title = element_text(face = "bold"),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      axis.text.x = element_text(size = 12),
      legend.position = "top",
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12),
      strip.text = element_text(size = 14, face = "bold", color = "black"),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5)
    )
}


# ---------------------------------------------------
# Gráfico rank-abundance con umbrales
# Ejemplo: Family 10% y 0.5%
# ---------------------------------------------------
plot_rank_abundance_thresholds <- function(
  tabla,
  rank_code,
  taxon_label,
  thresholds = c(0.10, 0.005)
) {
  
  rank_global <- tabla %>%
    filter(grepl("Bacteria", Fullname)) %>%
    filter(TaxRank == rank_code) %>%
    rename(Taxon = Name2) %>%
    group_by(Taxon) %>%
    summarise(
      TotalReads = sum(Reads, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      GlobalRel = TotalReads / sum(TotalReads, na.rm = TRUE)
    ) %>%
    arrange(desc(GlobalRel)) %>%
    mutate(Rank = row_number())
  
  p <- ggplot(rank_global, aes(x = Rank, y = GlobalRel)) +
    geom_point() +
    geom_line() +
    scale_y_continuous(labels = scales::percent) +
    labs(
      title = paste0(
        "Rank-abundance of ",
        ifelse(taxon_label == "Family", "families", paste0(tolower(taxon_label), "s"))
      ),
      subtitle = paste0(
        "Global read-based abundance; thresholds: ",
        paste(scales::percent(thresholds, accuracy = 0.1), collapse = " and "),
        ""
      ),
      x = paste0(taxon_label, " rank (most to least abundant)"),
      y = "Global relative abundance (%)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.title = element_text(face = "bold", size = 12),
      axis.text = element_text(color = "black", size = 11),
      plot.title = element_text(face = "bold", color = "black", size = 15),
      plot.subtitle = element_text(size = 11, margin = margin(b = 10)),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(12, 16, 12, 12)
    )
  
  for (thr in thresholds) {
    p <- p +
      geom_hline(yintercept = thr, linetype = "dashed", color = "red") +
      annotate(
        "text",
        x = Inf,
        y = thr,
        label = scales::percent(thr, accuracy = 0.1),
        color = "red",
        hjust = 1.1,
        vjust = -0.5,
        size = 4
      )
  }
  
  return(p)
}


# ---------------------------------------------------
# Top 15 global por Reads + Others
# Orden alfabético en leyenda y stack
# ---------------------------------------------------
plot_top15_global <- function(
  tabla,
  rank_code,
  taxon_label,
  taxon_label_plural
) {
  
  tabla <- prepare_patient_week_labels(tabla)
  
  # ---------------------------------------------------
  # 0) Filtrar a Bacteria + rango taxonómico
  # ---------------------------------------------------
  tabla_rank_raw <- tabla %>%
    filter(grepl("Bacteria", Fullname)) %>%
    filter(TaxRank == rank_code) %>%
    rename(Taxon = Name2)
  
  # ---------------------------------------------------
  # 1) Top 15 taxones globales por reads
  # ---------------------------------------------------
  global_reads <- tabla_rank_raw %>%
    group_by(Taxon) %>%
    summarise(
      TotalReads = sum(Reads, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(TotalReads))
  
  top15_taxa <- global_reads %>%
    slice_head(n = 15) %>%
    pull(Taxon)
  
  # ---------------------------------------------------
  # 2) Orden alfabético de los top15 + Others al final
  # ---------------------------------------------------
  taxon_levels <- c(sort(top15_taxa), "Others")
  
  # ---------------------------------------------------
  # 3) Preparar datos para plot: colapsar a Others
  # ---------------------------------------------------
  tabla_plot <- tabla_rank_raw %>%
    mutate(
      Taxon = if_else(Taxon %in% top15_taxa, Taxon, "Others")
    ) %>%
    group_by(Patient, follow_up_week, Taxon) %>%
    summarise(
      Reads = sum(Reads, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      Patient = factor(Patient, levels = paste0("PAC", 1:8)),
      follow_up_week = factor(follow_up_week, levels = c("W0", "W1", "W2")),
      Taxon = factor(Taxon, levels = taxon_levels)
    )
  
  # ---------------------------------------------------
  # 4) Paleta fija siguiendo ese orden alfabético
  # ---------------------------------------------------
  cols <- make_fixed_cols(taxon_levels)
  
  # ---------------------------------------------------
  # 5) Plot
  # ---------------------------------------------------
  p <- ggplot(
    tabla_plot,
    aes(x = follow_up_week, y = Reads, fill = Taxon)
  ) +
    geom_col(position = "fill", color = "black", width = 0.8) +
    scale_fill_manual(
      values = cols,
      drop = FALSE,
      labels = if (rank_code %in% c("G", "G1", "S")) taxon_plotmath_labels else waiver()
    ) +
    scale_y_continuous(labels = scales::percent) +
    labs(
      title = paste0(
        "Taxonomic Composition by Sample at ",
        taxon_label,
        " level"
      ),
      subtitle = paste0(
        "Top 15 ",
        taxon_label_plural,
        " by total reads across the cohort (others grouped as 'Others')"
      ),
      x = "Follow-up week",
      y = "Relative abundance (%)",
      fill = taxon_label
    ) +
    theme_taxonomy() +
    facet_wrap(~ Patient)
  
  return(p)
}


# ---------------------------------------------------
# Taxones por umbral global de Reads + Others
# Orden alfabético en leyenda y stack
# Ejemplo: Family >= 10% o Family >= 0.5%
# ---------------------------------------------------
plot_global_threshold <- function(
  tabla,
  rank_code,
  taxon_label,
  taxon_label_plural,
  threshold
) {
  
  tabla <- prepare_patient_week_labels(tabla)
  
  # ---------------------------------------------------
  # 0) Filtrar a Bacteria + rango taxonómico
  # ---------------------------------------------------
  tabla_rank_raw <- tabla %>%
    filter(grepl("Bacteria", Fullname)) %>%
    filter(TaxRank == rank_code) %>%
    rename(Taxon = Name2)
  
  # ---------------------------------------------------
  # 1) Abundancia global por reads
  # ---------------------------------------------------
  global_reads <- tabla_rank_raw %>%
    group_by(Taxon) %>%
    summarise(
      TotalReads = sum(Reads, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      GlobalRel = TotalReads / sum(TotalReads, na.rm = TRUE)
    ) %>%
    arrange(desc(GlobalRel))
  
  taxa_keep <- global_reads %>%
    filter(GlobalRel >= threshold) %>%
    pull(Taxon)
  
  # ---------------------------------------------------
  # 2) Orden alfabético de los taxones seleccionados + Others al final
  # ---------------------------------------------------
  taxon_levels <- c(sort(taxa_keep), "Others")
  
  # ---------------------------------------------------
  # 3) Preparar datos para plot
  # ---------------------------------------------------
  tabla_plot <- tabla_rank_raw %>%
    mutate(
      Taxon = if_else(Taxon %in% taxa_keep, Taxon, "Others")
    ) %>%
    group_by(Patient, follow_up_week, Taxon) %>%
    summarise(
      Reads = sum(Reads, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      Patient = factor(Patient, levels = paste0("PAC", 1:8)),
      follow_up_week = factor(follow_up_week, levels = c("W0", "W1", "W2")),
      Taxon = factor(Taxon, levels = taxon_levels)
    )
  
  # ---------------------------------------------------
  # 4) Paleta fija siguiendo ese orden alfabético
  # ---------------------------------------------------
  cols <- make_fixed_cols(taxon_levels)
  
  threshold_label <- scales::percent(threshold, accuracy = 0.1)
  
  # ---------------------------------------------------
  # 5) Plot
  # ---------------------------------------------------
  p <- ggplot(
    tabla_plot,
    aes(x = follow_up_week, y = Reads, fill = Taxon)
  ) +
    geom_col(position = "fill", color = "black", width = 0.8) +
    scale_fill_manual(
      values = cols,
      drop = FALSE,
      labels = if (rank_code %in% c("G", "G1", "S")) taxon_plotmath_labels else waiver()
    ) +
    scale_y_continuous(labels = scales::percent) +
    labs(
      title = paste0(
        "Taxonomic Composition by Sample at ",
        taxon_label,
        " level"
      ),
      subtitle = paste0(
        taxon_label_plural,
        " with ≥",
        threshold_label,
        " of total reads across the cohort (others grouped as 'Others')"
      ),
      x = "Follow-up week",
      y = "Relative abundance (%)",
      fill = taxon_label
    ) +
    theme_taxonomy() +
    facet_wrap(~ Patient)
  
  return(p)
}
