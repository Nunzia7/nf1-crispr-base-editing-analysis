suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
  library(scales)
})

# ---------------------------
# paths
# ---------------------------
# Override PROJECT_DIR, ANALYSIS_DIR, or BESTIMATE_DIR with environment variables
# instead of editing hard-coded local/HPC paths in the script.
project_dir <- Sys.getenv("PROJECT_DIR", unset = getwd())
analysis_dir <- Sys.getenv("ANALYSIS_DIR", unset = file.path(project_dir, "analysis"))
in_dir  <- Sys.getenv("BESTIMATE_DIR", unset = file.path(analysis_dir, "bestimate_out_plus_vep_full"))
out_dir <- file.path(in_dir, "plots_bestimate_nf1")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

crispr_file <- file.path(in_dir, "NF1_NGN_A2G_act4-8_crispr_df.csv")
edit_file   <- file.path(in_dir, "NF1_NGN_A2G_act4-8_edit_df.csv")

# ---------------------------
# read data
# ---------------------------
crispr_df <- fread(crispr_file)
edit_df   <- fread(edit_file)

# ---------------------------
# helpers
# ---------------------------
extract_start <- function(x) {
  # Converts strings like "17:31377624-31377646" -> 31377624
  as.numeric(sub(".*:(\\d+)-.*", "\\1", x))
}

extract_end <- function(x) {
  # Converts strings like "17:31377624-31377646" -> 31377646
  as.numeric(sub(".*-(\\d+)$", "\\1", x))
}

safe_bool <- function(x) {
  # Handles logical / character / integer robustly
  if (is.logical(x)) return(replace(x, is.na(x), FALSE))
  x <- as.character(x)
  x <- trimws(tolower(x))
  x %in% c("true", "t", "1", "yes")
}

# ---------------------------
# clean guide-level table
# ---------------------------
crispr_clean <- crispr_df %>%
  mutate(
    guide_id = CRISPR_PAM_Sequence,
    guide_start = extract_start(Location),
    guide_end   = extract_end(Location),
    guide_mid   = (guide_start + guide_end) / 2,
    guide_in_CDS = safe_bool(guide_in_CDS),
    Poly_T       = safe_bool(Poly_T),
    GC_percent   = suppressWarnings(as.numeric(`GC%`))
  ) %>%
  distinct(guide_id, .keep_all = TRUE)

# ---------------------------
# summarize edit-level info to one row per guide
# ---------------------------
edit_summary <- edit_df %>%
  mutate(
    guide_id = CRISPR_PAM_Sequence,
    Edit_in_Exon = safe_bool(Edit_in_Exon),
    Edit_in_CDS  = safe_bool(Edit_in_CDS),
    mutation_on_guide  = safe_bool(mutation_on_guide),
    guide_change_mutation = safe_bool(guide_change_mutation),
    mutation_on_window = safe_bool(mutation_on_window),
    mutation_on_PAM    = safe_bool(mutation_on_PAM),
    Edit_Location = suppressWarnings(as.numeric(Edit_Location)),
    n_edits_from_row = suppressWarnings(as.numeric(`# Edits/guide`))
  ) %>%
  group_by(guide_id) %>%
  summarise(
    n_edit_rows        = n(),  # actual number of editable positions observed
    n_edits_reported   = suppressWarnings(max(n_edits_from_row, na.rm = TRUE)),
    has_edit_in_exon   = any(Edit_in_Exon, na.rm = TRUE),
    has_edit_in_CDS    = any(Edit_in_CDS, na.rm = TRUE),
    any_mut_on_guide   = any(mutation_on_guide, na.rm = TRUE),
    any_mut_on_window  = any(mutation_on_window, na.rm = TRUE),
    min_edit_pos       = suppressWarnings(min(Edit_Location, na.rm = TRUE)),
    max_edit_pos       = suppressWarnings(max(Edit_Location, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    n_edits_reported = ifelse(is.infinite(n_edits_reported), NA, n_edits_reported),
    n_edit_positions = ifelse(!is.na(n_edits_reported), n_edits_reported, n_edit_rows)
  )

# ---------------------------
# merge final guide-level table
# ---------------------------
guide_tbl <- crispr_clean %>%
  left_join(edit_summary, by = "guide_id") %>%
  mutate(
    has_any_edit      = !is.na(n_edit_positions) & n_edit_positions > 0,
    has_edit_in_exon  = replace_na(has_edit_in_exon, FALSE),
    has_edit_in_CDS   = replace_na(has_edit_in_CDS, FALSE),
    category = case_when(
      has_edit_in_CDS  ~ "Edit in CDS",
      has_edit_in_exon ~ "Edit in exon only",
      has_any_edit     ~ "Editable, non-exonic",
      TRUE             ~ "No editable base"
    ),
    category = factor(
      category,
      levels = c("Edit in CDS", "Edit in exon only", "Editable, non-exonic", "No editable base")
    )
  )

# Save summary table
fwrite(guide_tbl, file.path(out_dir, "NF1_guides_summary_one_row_per_guide.csv"))

# ---------------------------
# plotting theme
# ---------------------------
theme_natureish <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(size = base_size + 1, face = "bold"),
      plot.subtitle = element_text(size = base_size),
      axis.title = element_text(size = base_size, face = "bold"),
      axis.text = element_text(size = base_size - 1, colour = "black"),
      legend.title = element_text(size = base_size, face = "bold"),
      legend.text = element_text(size = base_size - 1),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.line = element_line(linewidth = 0.4),
      axis.ticks = element_line(linewidth = 0.4)
    )
}

save_plot <- function(filename, plot, width = 16, height = 10) {
  ggsave(
    filename = file.path(out_dir, filename),
    plot = plot,
    width = width,
    height = height,
    units = "cm",
    dpi = 600,
    device = "png",
    bg = "white"
  )
}

# ---------------------------
# Plot 1: distribution of editable positions per guide
# ---------------------------
p1_df <- guide_tbl %>%
  filter(has_any_edit)

p1 <- ggplot(p1_df, aes(x = n_edit_positions)) +
  geom_histogram(binwidth = 1, boundary = 0, closed = "left",
                 linewidth = 0.3, fill = "grey70", colour = "black") +
  scale_x_continuous(breaks = pretty_breaks()) +
  labs(
    title = "Editable positions per NF1 guide",
    x = "Number of editable positions per guide",
    y = "Number of guides"
  ) +
  theme_natureish()

save_plot("plot1_editable_positions_histogram.png", p1)

# ---------------------------
# Plot 2: percentage of guides by category
# ---------------------------
p2_df <- guide_tbl %>%
  count(category, name = "n") %>%
  mutate(percent = 100 * n / sum(n))

p2 <- ggplot(p2_df, aes(x = category, y = percent)) +
  geom_col(width = 0.7, fill = "grey70", colour = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f%%", percent)), vjust = -0.35, size = 3.5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(
    title = "Functional annotation of NF1 guides",
    x = NULL,
    y = "Guides (%)"
  ) +
  theme_natureish() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

save_plot("plot2_guide_categories_barplot.png", p2)

# ---------------------------
# Plot 3: guide distribution along NF1
# ---------------------------
# ordered rank on y just to spread points
p3_df <- guide_tbl %>%
  arrange(guide_mid) %>%
  mutate(guide_rank = row_number())

p3 <- ggplot(p3_df, aes(x = guide_mid, y = 0, shape = category)) +
  geom_point(size = 2.1, stroke = 0.3, colour = "black") +
  labs(
    title = "Distribution of NF1 guides along genomic coordinates",
    x = "Genomic position (midpoint of guide + PAM)",
    y = NULL
  ) +
  theme_natureish()

save_plot("plot3_guides_along_nf1.png", p3, width = 18, height = 11)


# ---------------------------
# optional: top candidate table
# ---------------------------
top_candidates <- guide_tbl %>%
  filter(
    has_edit_in_CDS,
    !Poly_T,
    !is.na(GC_percent),
    !is.na(n_edit_positions),
    n_edit_positions <= 2,
    GC_percent >= 30,
    GC_percent <= 75
  ) %>%
  arrange(desc(has_edit_in_CDS), n_edit_positions, desc(GC_percent)) %>%
  select(
    Hugo_Symbol, guide_id, gRNA_Target_Sequence, Location, Direction,
    guide_in_CDS, has_any_edit, has_edit_in_exon, has_edit_in_CDS,
    n_edit_positions, GC_percent, Poly_T
  )

fwrite(top_candidates, file.path(out_dir, "NF1_top_candidates_CDS_guides.csv"))

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#conditions = c("DM1_1_vs_DM1_5","DM1_1_vs_UN","DM1_5_vs_UN","RMC_vs_UN","UN_vs_T0")
conditions = c("starv_vs_UN")

for (condition in conditions){
merged <- read.delim(file.path(analysis_dir, "mageck_test_x_bestimate_nf1", condition, "merged_mageck_test_bestimate.tsv"), sep = "\t")
enriched <- merged %>% filter(LFC>1,FDR<0.05,Gene=="NF1",!is.na(n_edit_positions))

crispr_enriched <- enriched %>%
  mutate(
    guide_id = CRISPR_PAM_Sequence,
    guide_start = extract_start(Location),
    guide_end   = extract_end(Location),
    guide_mid   = (guide_start + guide_end) / 2,
    guide_in_CDS = safe_bool(guide_in_CDS),
  ) %>%
  distinct(guide_id, .keep_all = TRUE)

enriched_positions <- ggplot(crispr_enriched, aes(x = guide_mid, y = 0)) +
  geom_point(size = 2.1, stroke = 0.3, colour = "black") +
  labs(
    title = "Enriched NF1 guides along genomic coordinates",
    x = "Genomic position (midpoint of guide + PAM)",
    y = NULL
  ) +
  theme_classic()


#subset for BEstimate
ngn_library_tot <- read.csv(file.path(project_dir, "input", "NGN_library.csv"))
ngn_library_enriched <- ngn_library_tot %>% filter(ID %in% enriched$sgrna)
#write.csv(ngn_library_enriched,"${PROJECT_DIR}/NGN_library_enriched.csv", row.names = FALSE)
#write.csv(ngn_library_enriched,file.path(analysis_dir, "mageck_test_x_bestimate_nf1", condition, "NGN_library_enriched.csv"), row.names = FALSE)

#protein_df <- read.csv("${PROJECT_DIR}/bestimate_out_plus_vep_enriched/NF1_NGN_A2G_act4-8_protein_df.csv")

protein_df <- read.csv(
  file.path(in_dir, "NF1_NGN_A2G_act4-8_protein_df.csv")
)

table(protein_df$Domain == "" | is.na(protein_df$Domain))
table(protein_df$curated_Domain == "" | is.na(protein_df$curated_Domain))
unique(na.omit(protein_df$Domain))
unique(na.omit(protein_df$curated_Domain))

plot_df <- merge(
  enriched,
  protein_df,
  by.x = c("CRISPR_PAM_Sequence", "Location", "Direction"),
  by.y = c("CRISPR_PAM_Sequence", "CRISPR_PAM_Location", "Direction"),
  all.x = TRUE
)

# Tieni solo varianti con posizione proteica
plot_df <- plot_df %>%
  mutate(
    Protein_Position = suppressWarnings(as.numeric(Protein_Position)),
    curated_Domain = ifelse(is.na(curated_Domain) | curated_Domain == "", "No domain", curated_Domain),
    label_me = curated_Domain != "No domain"
  ) %>%
  filter(!is.na(Protein_Position))


#write.csv(plot_df,"${PROJECT_DIR}/bestimate_out/mageck_test_x_bestimate_nf1/mageck_enriched_bestimate_plus_vep.csv")
write.csv(plot_df,file.path(analysis_dir, "mageck_test_x_bestimate_nf1", condition, "mageck_enriched_bestimate_plus_vep.csv"))


# tabella domini minimale da disegnare sotto
# puoi cambiarla se hai intervalli migliori
domains <- data.frame(
  domain = c("Ras-GAP", "CRAL-TRIO", "Lipid binding", "Disordered"),
  start  = c(1251, 1580, 1580, 2787),
  end    = c(1482, 1738, 1837, 2839)
)

domain_colors <- c(
  "Ras-GAP" = "#E41A1C",
  "Lipid binding" = "#377EB8",
  "Disordered" = "#4DAF4A",
  "No domain" = "grey70"
)

p <- ggplot(plot_df, aes(x = Protein_Position, y = LFC)) +
  
  geom_rect(
    data = domains,
    aes(xmin = start, xmax = end, ymin = -0.25, ymax = 0.25, fill = domain),
    inherit.aes = FALSE,
    alpha = 0.25
  ) +
  
  geom_hline(yintercept = 0, linewidth = 0.4) +
  
  geom_segment(aes(xend = Protein_Position, y = 0, yend = LFC),
               linewidth = 0.35, alpha = 0.8) +
  
  geom_point(aes(color = curated_Domain, size = -log10(FDR)), alpha = 0.9) +
  
  scale_color_manual(values = domain_colors) +
  scale_fill_manual(values = domain_colors) +
  
  scale_x_continuous(limits = c(1, 2839), expand = c(0.01, 0.01)) +
  
  labs(
    x = "NF1 protein position",
    y = "LFC",
    color = "Annotated domain",
    fill = "Protein domains",
    size = expression(-log[10](FDR))
  ) +
  
  theme_classic(base_size = 15) +
  theme(
    legend.position = "bottom"
  )

p_amm <- ggplot(plot_df, aes(x = Protein_Position, y = LFC)) +
  
  geom_rect(
    data = domains,
    aes(xmin = start, xmax = end, ymin = -0.25, ymax = 0.25, fill = domain),
    inherit.aes = FALSE,
    alpha = 0.25
  ) +
  
  geom_hline(yintercept = 0, linewidth = 0.4) +
  
  geom_segment(aes(xend = Protein_Position, y = 0, yend = LFC),
               linewidth = 0.35, alpha = 0.8) +
  
  geom_point(aes(color = curated_Domain, size = -log10(FDR)), alpha = 0.9) +
  geom_text(
    data = subset(plot_df, label_me),
    aes(label = Protein_Change),
    size = 5,
    vjust = -0.8,
    check_overlap = TRUE
  ) +
  
  scale_color_manual(values = domain_colors) +
  scale_fill_manual(values = domain_colors) +
  
  scale_x_continuous(limits = c(1, 2839), expand = c(0.01, 0.01)) +
  
  labs(
    x = "NF1 protein position",
    y = "LFC",
    color = "Annotated domain",
    fill = "Protein domains",
    size = expression(-log[10](FDR))
  ) +
  
  theme_classic(base_size = 15) +
  theme(
    legend.position = "bottom"
  )

#ggsave("${PROJECT_DIR}/bestimate_out/mageck_test_x_bestimate_nf1/domains.png", p,device = "png", width = 20, height = 5)
#ggsave("${PROJECT_DIR}/bestimate_out/mageck_test_x_bestimate_nf1/domains.pdf", p,device = "pdf", width = 20, height = 5)
#ggsave("${PROJECT_DIR}/bestimate_out/mageck_test_x_bestimate_nf1/domains_amm.png", p_amm,device = "png", width = 20, height = 5)

ggsave(file.path(analysis_dir, "mageck_test_x_bestimate_nf1", condition, "domains.png"), p,device = "png", width = 20, height = 5)
ggsave(file.path(analysis_dir, "mageck_test_x_bestimate_nf1", condition, "domains.pdf"), p,device = "pdf", width = 20, height = 5)
ggsave(file.path(analysis_dir, "mageck_test_x_bestimate_nf1", condition, "domains_amm.png"), p_amm,device = "png", width = 20, height = 5)
}


conditions = c("starv_vs_UN")

ngn_library_tot <- read.csv(file.path(project_dir, "input", "NGN_library.csv"))

protein_df <- read.csv(
  file.path(in_dir, "NF1_NGN_A2G_act4-8_protein_df.csv")
)

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

for (condition in conditions) {
  
  out_dir <- file.path(analysis_dir, "mageck_test_x_bestimate_nf1", condition)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  merged <- read.delim(
    paste0(out_dir, "merged_mageck_test_bestimate.tsv"),
    sep = "\t"
  )
  
  merged_clean <- merged %>%
    mutate(
      CRISPR_PAM_Sequence = trimws(CRISPR_PAM_Sequence),
      Location = trimws(Location),
      Direction = trimws(Direction)
    )
  
  protein_df_clean <- protein_df %>%
    mutate(
      CRISPR_PAM_Sequence = trimws(CRISPR_PAM_Sequence),
      CRISPR_PAM_Location = trimws(CRISPR_PAM_Location),
      Direction = trimws(Direction)
    )
  
  merged_plus_vep <- merged_clean %>%
    left_join(
      protein_df_clean,
      by = c(
        "CRISPR_PAM_Sequence" = "CRISPR_PAM_Sequence",
        "Location" = "CRISPR_PAM_Location",
        "Direction" = "Direction"
      )
    )
  
  write.csv(
    merged_plus_vep,
    file.path(out_dir, "merged_mageck_test_bestimate_plus_vep.csv"),
    row.names = FALSE
  )
  
  for (direction_type in c("enriched", "depleted")) {
    
    if (direction_type == "enriched") {
      selected <- merged %>%
        filter(LFC > 1, FDR < 0.05, Gene == "NF1", !is.na(n_edit_positions))
    } else {
      selected <- merged %>%
        filter(LFC < -1, FDR < 0.05, Gene == "NF1", !is.na(n_edit_positions))
    }
    
    if (nrow(selected) == 0) {
      message("No ", direction_type, " guides found for ", condition)
      next
    }
    
    crispr_selected <- selected %>%
      mutate(
        guide_id = CRISPR_PAM_Sequence,
        guide_start = extract_start(Location),
        guide_end   = extract_end(Location),
        guide_mid   = (guide_start + guide_end) / 2,
        guide_in_CDS = safe_bool(guide_in_CDS)
      ) %>%
      distinct(guide_id, .keep_all = TRUE)
    
    positions_plot <- ggplot(crispr_selected, aes(x = guide_mid, y = 0)) +
      geom_point(size = 2.1, stroke = 0.3, colour = "black") +
      labs(
        title = paste0(direction_type, " NF1 guides along genomic coordinates"),
        x = "Genomic position (midpoint of guide + PAM)",
        y = NULL
      ) +
      theme_classic()
    
    ggsave(
      file.path(out_dir, paste0("positions_", direction_type, ".png")),
      positions_plot,
      device = "png",
      width = 12,
      height = 3
    )
    
    ngn_library_selected <- ngn_library_tot %>%
      filter(ID %in% selected$sgrna)
    
    write.csv(
      ngn_library_selected,
      file.path(out_dir, paste0("NGN_library_", direction_type, ".csv")),
      row.names = FALSE
    )
    
    plot_df <- merge(
      selected,
      protein_df,
      by.x = c("CRISPR_PAM_Sequence", "Location", "Direction"),
      by.y = c("CRISPR_PAM_Sequence", "CRISPR_PAM_Location", "Direction"),
      all.x = TRUE
    )
    
    plot_df <- plot_df %>%
      mutate(
        Protein_Position = suppressWarnings(as.numeric(Protein_Position)),
        curated_Domain = ifelse(
          is.na(curated_Domain) | curated_Domain == "",
          "No domain",
          curated_Domain
        ),
        label_me = curated_Domain != "No domain"
      ) %>%
      filter(!is.na(Protein_Position))
    
    write.csv(
      plot_df,
      file.path(out_dir, paste0("mageck_", direction_type, "_bestimate_plus_vep.csv")),
      row.names = FALSE
    )
    
    p <- ggplot(plot_df, aes(x = Protein_Position, y = LFC)) +
      geom_rect(
        data = domains,
        aes(xmin = start, xmax = end, ymin = -0.25, ymax = 0.25, fill = domain),
        inherit.aes = FALSE,
        alpha = 0.25
      ) +
      geom_hline(yintercept = 0, linewidth = 0.4) +
      geom_segment(
        aes(xend = Protein_Position, y = 0, yend = LFC),
        linewidth = 0.35,
        alpha = 0.8
      ) +
      geom_point(aes(color = curated_Domain, size = -log10(FDR)), alpha = 0.9) +
      scale_color_manual(values = domain_colors) +
      scale_fill_manual(values = domain_colors) +
      scale_x_continuous(limits = c(1, 2839), expand = c(0.01, 0.01)) +
      labs(
        x = "NF1 protein position",
        y = "LFC",
        color = "Annotated domain",
        fill = "Protein domains",
        size = expression(-log[10](FDR))
      ) +
      theme_classic(base_size = 15) +
      theme(legend.position = "bottom")
    
    p_amm <- p +
      geom_text(
        data = subset(plot_df, label_me),
        aes(label = Protein_Change),
        size = 5,
        vjust = ifelse(direction_type == "enriched", -0.8, 1.3),
        check_overlap = TRUE
      )
    
    ggsave(
      file.path(out_dir, paste0("domains_", direction_type, ".png")),
      p,
      device = "png",
      width = 20,
      height = 5
    )
    
    ggsave(
      file.path(out_dir, paste0("domains_", direction_type, ".pdf")),
      p,
      device = "pdf",
      width = 20,
      height = 5
    )
    
    ggsave(
      file.path(out_dir, paste0("domains_amm_", direction_type, ".png")),
      p_amm,
      device = "png",
      width = 20,
      height = 5
    )
  }
}

