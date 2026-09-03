suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
})

# -----------------------------------------------------------------------------
# Plot di una singola condizione lungo la posizione proteica di NF1.
# Tutte le guide/varianti sono colorate in base allo score renovo
# e annotate con il cambio aminoacidico/HGVSp.
# -----------------------------------------------------------------------------

project_dir <- Sys.getenv("ANALYSIS_DIR", unset = file.path(getwd(), "analysis"))

# Parametri modificabili da RStudio.
# Da terminale si possono sovrascrivere con 2 argomenti, per esempio:
#   Rscript --vanilla scripts/plot_nf1_probability_class_single_condition.R UN_vs_T0 enriched
condition_to_plot <- "starv_vs_UN"
direction_to_plot <- "depleted"

input_dir <- file.path(project_dir, "after_musa", "output")
out_dir <- file.path(input_dir, "plots_nf1_probability_class_single_condition")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 2) {
  condition_to_plot <- args[[1]]
  direction_to_plot <- args[[2]]
} else if (length(args) != 0) {
  stop(
    "Uso: Rscript --vanilla plot_nf1_probability_class_single_condition.R ",
    "<condition> <direction>\n",
    "Esempio: Rscript --vanilla plot_nf1_probability_class_single_condition.R ",
    "UN_vs_T0 enriched"
  )
}

if (!direction_to_plot %in% c("depleted", "enriched")) {
  stop("direction must be one of: depleted, enriched")
}

input_file <- file.path(
  input_dir,
  paste0(condition_to_plot, "_", direction_to_plot, "_all_selected_all_missense.csv")
)

if (!file.exists(input_file)) {
  stop("File condizione non trovato: ", input_file)
}

standardize_comparison <- function(x, comparison) {
  lfc_col <- paste0(comparison, "_LFC")
  fdr_col <- paste0(comparison, "_FDR")
  sgrna_col <- paste0(comparison, "_sgrna")
  
  missing_cols <- setdiff(c(lfc_col, fdr_col, sgrna_col), names(x))
  if (length(missing_cols) > 0) {
    stop(
      "Nel file ", comparison, " mancano colonne: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  x %>%
    mutate(
      LFC = suppressWarnings(as.numeric(.data[[lfc_col]])),
      FDR = suppressWarnings(as.numeric(.data[[fdr_col]])),
      sgrna = .data[[sgrna_col]]
    )
}

domains <- data.frame(
  domain = c("Ras-GAP", "CRAL-TRIO", "Lipid binding", "Disordered"),
  start  = c(1251, 1580, 1580, 2787),
  end    = c(1482, 1738, 1837, 2839)
)

domain_colors <- c(
  "Ras-GAP" = "#E41A1C",
  "CRAL-TRIO" = "#984EA3",
  "Lipid binding" = "#377EB8",
  "Disordered" = "#4DAF4A",
  "No domain" = "grey70"
)

probability_score_colors <- c("#377EB8", "grey30", "#E41A1C")

prepare_plot_df <- function(x) {
  if (!"curated_Domain" %in% names(x)) {
    x$curated_Domain <- "No domain"
  }
  if (!"HGVSp" %in% names(x)) {
    x$HGVSp <- NA_character_
  }
  if (!"Amino_acids" %in% names(x)) {
    x$Amino_acids <- NA_character_
  }
  if (!"Probability_Class_1" %in% names(x)) {
    x$Probability_Class_1 <- NA_character_
  }
  
  x %>%
    mutate(
      Protein_position_num = suppressWarnings(as.numeric(Protein_position)),
      Probability_Class_1_num = suppressWarnings(as.numeric(Probability_Class_1)),
      curated_Domain = ifelse(
        is.na(curated_Domain) | curated_Domain == "",
        "No domain",
        curated_Domain
      ),
      variant_label = case_when(
        !is.na(HGVSp) & HGVSp != "" ~ HGVSp,
        !is.na(Amino_acids) & Amino_acids != "" ~ Amino_acids,
        TRUE ~ as.character(Protein_position)
      )
    ) %>%
    filter(!is.na(Protein_position_num))
}

plot_df <- read.csv(input_file, check.names = FALSE, stringsAsFactors = FALSE) %>%
  standardize_comparison(condition_to_plot) %>%
  mutate(confronto = paste(condition_to_plot, direction_to_plot)) %>%
  prepare_plot_df()

if (nrow(plot_df) == 0) {
  stop("Nessuna variante con Protein_position numerica da plottare.")
}

plot_label <- paste(condition_to_plot, direction_to_plot)
file_label <- str_replace_all(plot_label, "[^A-Za-z0-9]+", "_")
file_label <- str_replace_all(file_label, "^_|_$", "")

p <- ggplot(plot_df, aes(x = Protein_position_num, y = LFC)) +
  geom_rect(
    data = domains,
    aes(xmin = start, xmax = end, ymin = -0.25, ymax = 0.25, fill = domain),
    inherit.aes = FALSE,
    alpha = 0.25
  ) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_segment(
    aes(xend = Protein_position_num, y = 0, yend = LFC, color = Probability_Class_1_num),
    linewidth = 0.55,
    alpha = 0.9
  ) +
  geom_point(
    aes(size = -log10(FDR), color = Probability_Class_1_num),
    alpha = 0.85,
    stroke = 0.7
  ) +
  geom_text(
    aes(label = variant_label),
    size = 2.6,
    vjust = -0.8,
    check_overlap = TRUE,
    show.legend = FALSE
  ) +
  scale_color_gradientn(
    colors = probability_score_colors,
    limits = c(0, 1),
    na.value = "grey75"
  ) +
  scale_fill_manual(values = domain_colors, drop = FALSE) +
  scale_size_continuous(range = c(1.8, 4.5)) +
  scale_x_continuous(limits = c(1, 2839), expand = c(0.01, 0.01)) +
  labs(
    title = paste0("NF1 protein positions colored by renovo (", plot_label, ")"),
    x = "NF1 protein position",
    y = "LFC",
    fill = "Protein domains",
    color = "renovo",
    size = expression(-log[10](FDR))
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.direction = "horizontal",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9),
    legend.key.size = unit(0.35, "cm"),
    legend.spacing.x = unit(0.2, "cm"),
    legend.margin = margin(t = 2, r = 5, b = 2, l = 5),
    plot.margin = margin(t = 10, r = 25, b = 15, l = 25)
  ) +
  guides(
    fill = guide_legend(nrow = 1, byrow = TRUE, title.position = "top"),
    color = guide_colorbar(title.position = "top", barwidth = unit(4, "cm")),
    size = guide_legend(
      nrow = 1,
      byrow = TRUE,
      title.position = "top",
      override.aes = list(color = "black", alpha = 0.8)
    )
  )

png_file <- paste0(file_label, "_nf1_renovo_labeled.png")
pdf_file <- paste0(file_label, "_nf1_renovo_labeled.pdf")

ggsave(
  file.path(out_dir, png_file),
  p,
  device = "png",
  width = 42,
  height = 12,
  units = "cm",
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(out_dir, pdf_file),
  p,
  device = "pdf",
  width = 42,
  height = 12,
  units = "cm"
)

message("Condition: ", plot_label)
message("Plot salvati in: ", out_dir)
message("PNG: ", png_file)
message("PDF: ", pdf_file)


##############################################################################################################
##############################################################################################################


# Reuse the public configurable analysis directory for the second plotting block.
base_dir <- Sys.getenv("ANALYSIS_DIR", unset = file.path(getwd(), "analysis"))

args <- commandArgs(trailingOnly = TRUE)
comparison <- if (length(args) >= 1) args[[1]] else "starv_vs_UN"

input_dir <- file.path(base_dir, "after_musa", "output")
plot_dir <- file.path(input_dir, paste0("plots_", comparison, "_all_missense"))
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

input_files <- c(
  enriched = file.path(input_dir, paste0(comparison, "_enriched_all_missense.csv")),
  depleted = file.path(input_dir, paste0(comparison, "_depleted_all_missense.csv"))
)

missing_input_files <- input_files[!file.exists(input_files)]
if (length(missing_input_files) > 0) {
  stop(
    paste0(
      "Missing input file(s) for comparison '", comparison, "':\n",
      paste(missing_input_files, collapse = "\n"),
      "\nRun 10_filter_musa_variants.R first for this comparison to create enriched/depleted all_missense tables."
    )
  )
}

required_columns <- c("AlphaMissense_score", "Probability_Class_1")
plot_columns <- c("AlphaMissense_score_max", "Probability_Class_1")
score_labels <- c(
  AlphaMissense_score_max = "AlphaMissense",
  Probability_Class_1 = "Renovo"
)

max_semicolon_numeric <- function(x) {
  vapply(
    strsplit(as.character(x), split = ";", fixed = TRUE),
    function(values) {
      values <- trimws(values)
      values <- values[!is.na(values) & values != "" & values != "." & values != "NA"]
      numeric_values <- suppressWarnings(as.numeric(values))
      numeric_values <- numeric_values[!is.na(numeric_values)]
      
      if (length(numeric_values) == 0) {
        return(NA_real_)
      }
      
      max(numeric_values)
    },
    numeric(1)
  )
}

read_all_missense <- function(file_path, condition) {
  data <- read.csv(
    file_path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  missing_columns <- setdiff(required_columns, names(data))
  if (length(missing_columns) > 0) {
    stop(
      paste0(
        "Missing columns in ", basename(file_path), ": ",
        paste(missing_columns, collapse = ", ")
      )
    )
  }
  
  data %>%
    mutate(
      condition = condition,
      AlphaMissense_score_max = max_semicolon_numeric(AlphaMissense_score),
      Probability_Class_1 = as.numeric(Probability_Class_1)
    ) %>%
    select(condition, all_of(plot_columns))
}

all_missense <- bind_rows(
  read_all_missense(input_files[["enriched"]], "enriched"),
  read_all_missense(input_files[["depleted"]], "depleted")
) %>%
  mutate(condition = factor(condition, levels = c("enriched", "depleted")))

all_missense_long <- all_missense %>%
  pivot_longer(
    cols = all_of(plot_columns),
    names_to = "score",
    values_to = "value"
  ) %>%
  filter(!is.na(value)) %>%
  mutate(score = factor(score, levels = plot_columns, labels = score_labels[plot_columns]))

summary_table <- all_missense_long %>%
  group_by(condition, score) %>%
  summarise(
    n = n(),
    mean = mean(value),
    median = median(value),
    sd = sd(value),
    min = min(value),
    q25 = quantile(value, 0.25),
    q75 = quantile(value, 0.75),
    max = max(value),
    .groups = "drop"
  )

write.csv(
  summary_table,
  file = file.path(plot_dir, paste0(comparison, "_all_missense_score_distribution_summary.csv")),
  row.names = FALSE,
  quote = TRUE
)

density_plot <- ggplot(
  all_missense_long,
  aes(x = value, fill = condition, colour = condition)
) +
  geom_density(alpha = 0.35, linewidth = 0.8, na.rm = TRUE) +
  facet_wrap(~ score) +
  theme_classic(base_size = 12) +
  labs(
    title = paste0(comparison, " all missense variants"),
    x = "Score value",
    y = "Density",
    fill = "LFC group",
    colour = "LFC group"
  )

ggsave(
  filename = file.path(plot_dir, paste0(comparison, "_all_missense_density_distributions.png")),
  plot = density_plot,
  width = 9,
  height = 5,
  dpi = 300
)
