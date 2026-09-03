library(dplyr)
library(stringr)

# Filter variant-prediction tables after joining MAGeCK/BEEstimate statistics.
# Override ANALYSIS_DIR instead of editing local/HPC paths in the script.
base_dir <- Sys.getenv("ANALYSIS_DIR", unset = file.path(getwd(), "analysis"))

input_dir <- file.path(base_dir, "after_musa")
output_dir <- file.path(input_dir, "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

save_variant_tables <- function(tables, comparison, direction, output_dir) {
  for (table_name in names(tables)) {
    write.csv(
      tables[[table_name]],
      file = file.path(output_dir, paste0(comparison, "_", direction, "_", table_name, ".csv")),
      row.names = FALSE,
      quote = TRUE
    )
  }
}


filter_variant_tables <- function(
    file_path,
    comparison = NULL,
    gene = "NF1",
    mane_label = "MANE_Select",
    variant_class = "SNV",
    lfc_min = 1,
    lfc_direction = c("enriched", "depleted"),
    fdr_max = 0.05,
    pathogenic_label = "Pathogenic",
    probability_min = 0.95,
    consequence_separator = ","
) {
  
  # Importazione
  data <- read.csv(
    file_path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  # Nomi dinamici delle colonne
  lfc_column <- paste0(comparison, "_LFC")
  fdr_column <- paste0(comparison, "_FDR")
  lfc_direction <- match.arg(lfc_direction)
  
  # Controllo delle colonne necessarie
  required_columns <- c(
    "MANE",
    "Consequence",
    "Prediction_Label",
    "Probability_Class_1",
    "Hugo_Symbol",
    "VARIANT_CLASS",
    lfc_column,
    fdr_column
  )
  
  missing_columns <- setdiff(required_columns, names(data))
  
  if (length(missing_columns) > 0) {
    stop(
      paste0(
        "Nel file mancano queste colonne: ",
        paste(missing_columns, collapse = ", ")
      )
    )
  }
  
  # Funzione per cercare termini esatti nella colonna Consequence
  has_consequence <- function(x, wanted_terms) {
    
    x[is.na(x)] <- ""
    
    vapply(
      strsplit(
        x,
        split = consequence_separator,
        fixed = TRUE
      ),
      function(terms) {
        any(trimws(terms) %in% wanted_terms)
      },
      logical(1)
    )
  }
  
  # Conseguenze mantenute indipendentemente dalla predizione
  priority_non_missense <- c(
    "start_lost",
    "stop_lost",
    "splice_donor_variant",
    "splice_acceptor_variant",
    "NMD_transcript_variant"
  )
  
  missense_consequence <- "missense_variant"
  
  # Filtri comuni
  common_filtered <- data %>%
    filter(
      .data[["MANE"]] == mane_label,
      if (lfc_direction == "enriched") {
        .data[[lfc_column]] > lfc_min
      } else {
        .data[[lfc_column]] < -lfc_min
      },
      .data[[fdr_column]] < fdr_max,
      .data[["Hugo_Symbol"]] == gene,
      .data[["VARIANT_CLASS"]] == variant_class
    )
  
  # Tabella 1:
  # contiene almeno una conseguenza prioritaria non-missense.
  # Può contenere anche missense se associata, per esempio,
  # a NMD_transcript_variant.
  non_missense_priority <- common_filtered %>%
    filter(
      has_consequence(
        .data[["Consequence"]],
        priority_non_missense
      )
    )
  
  # Tabella 2:
  # contiene missense e supera i filtri di patogenicità
  pathogenic_missense <- common_filtered %>%
    filter(
      has_consequence(
        .data[["Consequence"]],
        missense_consequence
      ),
      .data[["Prediction_Label"]] == pathogenic_label,
      .data[["Probability_Class_1"]] > probability_min
    )
  
  # Tabella 3:
  # contiene tutti i missense, senza filtri su patogenicità
  # o score/probabilità della predizione.
  all_missense <- common_filtered %>%
    filter(
      has_consequence(
        .data[["Consequence"]],
        missense_consequence
      )
    )
  
  # Tabella 4:
  # unione finale, rimuovendo le righe eventualmente presenti
  # in entrambe le tabelle
  all_selected <- bind_rows(
    non_missense_priority,
    pathogenic_missense
  ) %>%
    distinct()
  
  # Tabella 5:
  # nuova unione finale che mantiene tutti i missense,
  # senza applicare i filtri su patogenicità o score.
  all_selected_all_missense <- bind_rows(
    non_missense_priority,
    all_missense
  ) %>%
    distinct()
  
  return(
    list(
      non_missense_priority = non_missense_priority,
      pathogenic_missense = pathogenic_missense,
      all_missense = all_missense,
      all_selected = all_selected,
      all_selected_all_missense = all_selected_all_missense
    )
  )
}

run_filtering_for_comparison <- function(file_name, comparison) {
  enriched_tables <- filter_variant_tables(
    file_path = file.path(input_dir, file_name),
    comparison = comparison,
    lfc_direction = "enriched"
  )
  
  depleted_tables <- filter_variant_tables(
    file_path = file.path(input_dir, file_name),
    comparison = comparison,
    lfc_direction = "depleted"
  )
  
  save_variant_tables(
    tables = enriched_tables,
    comparison = comparison,
    direction = "enriched",
    output_dir = output_dir
  )
  
  save_variant_tables(
    tables = depleted_tables,
    comparison = comparison,
    direction = "depleted",
    output_dir = output_dir
  )
  
  list(
    enriched = enriched_tables,
    depleted = depleted_tables
  )
}

starv_vs_UN_tables <- run_filtering_for_comparison(
  file_name = "musa_predictions_with_starv_vs_UN_screen_stats_expanded.csv",
  comparison = "starv_vs_UN"
)

message("Saved filtering outputs to: ", output_dir)