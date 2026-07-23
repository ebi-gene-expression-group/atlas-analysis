#!/usr/bin/env Rscript

# Identification of marker genes for baseline RNA-seq experiments in Expression Atlas.

# -----------------------
# Setup & dependencies
# -----------------------
required_packages <- c("xml2", "dplyr", "data.table", "this.path")

require_or_fail <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required but not installed.", pkg), call. = FALSE)
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

invisible(lapply(required_packages, require_or_fail))
message("All required packages are loaded successfully.")

# -----------------------
# Argument parsing
# -----------------------
args <- commandArgs(trailingOnly = TRUE)

usage <- paste(
  "Usage:",
  "  Rscript script.R <config_xml> <expression_data.undecorated> <expression_data> <output_file> [<specificity_score_cutoff> <expression_cutoff>]",
  "",
  "Notes:",
  "  - specificity_score_cutoff in [0,1], default 0.25",
  "  - expression_cutoff in TPMS/FPKMs (or ppb), default 0.5",
  sep = "\n"
)

if (length(args) < 4 || length(args) > 6) {
  stop(usage, call. = FALSE)
}

# Default values
default_specificity_cutoff <- 0.25   # between 0 and 1
default_expression_cutoff  <- 0.5    # in tpms/fpkms/ppb

# Normalize args to length 6 (fill optional cutoffs when missing)
if (length(args) == 4) {
  args <- c(args, default_specificity_cutoff, default_expression_cutoff)
} else if (length(args) == 5) {
  args <- c(args, default_expression_cutoff)
}

arg_names <- c(
  "config_xml", "expr_undecorated", "expr_decorated",
  "output_file", "specificity_cutoff", "expression_cutoff"
)
names(args) <- arg_names

message("Arguments to be used:")
for (nm in names(args)) message("  ", nm, ": ", args[[nm]])

# Coerce cutoffs to numeric (they arrive as character)
SPECIFICITY_SCORE_CUTOFF <- suppressWarnings(as.numeric(args[["specificity_cutoff"]]))
EXPRESSION_CUTOFF        <- suppressWarnings(as.numeric(args[["expression_cutoff"]]))

if (!is.finite(SPECIFICITY_SCORE_CUTOFF) || SPECIFICITY_SCORE_CUTOFF < 0 || SPECIFICITY_SCORE_CUTOFF > 1) {
  stop("specificity_cutoff must be a number in [0,1].", call. = FALSE)
}
if (!is.finite(EXPRESSION_CUTOFF) || EXPRESSION_CUTOFF < 0) {
  stop("expression_cutoff must be a non-negative number.", call. = FALSE)
}

# Determine data type and accession
metric   <- sub(".*-(.*)\\.tsv\\.undecorated$", "\\1", args[["expr_undecorated"]])
is_rnaseq <- tolower(metric) %in% c("tpms", "fpkms")
message("Detected data type: ", ifelse(is_rnaseq, "RNA-seq (MGFR)", "Proteomics (specificity)"))

if (is_rnaseq) {
  accession <- sub("-(tpms|fpkms)\\.tsv\\.undecorated$", "", args[["expr_undecorated"]])
} else {
  accession   <- sub("\\.tsv\\.undecorated\\.aggregated$", "", args[["expr_undecorated"]])
  script_path <- normalizePath(this.path::this.path())
  script_dir  <- dirname(script_path)
  metric      <- "ppb"
  utils_file  <- file.path(script_dir, "marker_gene_utils.R")
  if (file.exists(utils_file)) {
    source(utils_file, local = TRUE)
  } else {
    message("Note: 'marker_gene_utils.R' not found at: ", utils_file, " (continuing without it)")
  }
}

# -----------------------
# Read config XML
# -----------------------
doc <- xml2::read_xml(args[["config_xml"]])

groups <- xml2::xml_find_all(doc, ".//assay_group")

assay_df <- dplyr::bind_rows(lapply(groups, function(group) {
  group_id <- xml2::xml_attr(group, "id")
  group_label <- xml2::xml_attr(group, "label")
  assays <- xml2::xml_find_all(group, ".//assay")
  data.frame(
    group    = group_id,
    label    = group_label,
    assay    = xml2::xml_text(assays),
    tech_rep = xml2::xml_attr(assays, "technical_replicate_id"),
    stringsAsFactors = FALSE
  )
}))

print(assay_df)

summary_df <- assay_df |>
  dplyr::group_by(group, label) |>
  dplyr::summarise(n_assays = dplyr::n(), .groups = "drop")
print(summary_df)

# -----------------------
# Read undecorated expression matrix
# -----------------------
dt <- data.table::fread(args[["expr_undecorated"]])
data.table::setDT(dt)

if (!is_rnaseq) {
  # Keep only Gene ID and *WithInSampleAbundance* columns (case-insensitive)
  keep_cols <- c("Gene ID", grep("WithInSampleAbundance", names(dt), ignore.case = TRUE, value = TRUE))
  dt <- dt[, ..keep_cols]

  # Keep rows with valid Gene ID (not NA/empty/purely numeric)
  dt <- dt[!is.na(`Gene ID`) & `Gene ID` != "" & !grepl("^\\d+$", `Gene ID`)]

  # Rename abundance columns from group codes to assay IDs found in XML
  # (strip everything after first dot in undecorated column names)
  old_cols <- names(dt)[-1]
  new_cols <- vapply(
    old_cols,
    function(nm) {
      grp <- sub("\\..*$", "", nm)
      # map group -> first assay id under that group
      cand <- assay_df$assay[assay_df$group == grp]
      if (length(cand) && !is.na(cand[1])) cand[1] else nm
    },
    character(1)
  )
  data.table::setnames(dt, old = old_cols, new = new_cols)
}

# Remove duplicated columns (keep first occurrence)
dt <- dt[, !duplicated(names(dt)), with = FALSE]

# Map assay -> label; rename matching columns to labels
assay_to_label <- setNames(assay_df$label, assay_df$assay)
columns_to_rename <- intersect(names(assay_to_label), colnames(dt))
if (length(columns_to_rename)) {
  data.table::setnames(dt, columns_to_rename, assay_to_label[columns_to_rename])
}

# Keep only columns present in config XML (plus Gene ID first column)
cfg_labels <- unique(assay_df$label)
keep_idx <- c(1L, which(names(dt)[-1] %in% cfg_labels) + 1L)
keep_idx <- unique(keep_idx[keep_idx <= ncol(dt)])
dt <- dt[, ..keep_idx]

# -----------------------
# Collapse technical reps by label (mean)
# -----------------------
dt_numeric <- dt[, -1, with = FALSE]
gene_ids   <- dt[[1]]

col_names    <- colnames(dt_numeric)
unique_names <- unique(col_names)

averaged_list <- lapply(unique_names, function(nm) {
  cols <- which(col_names == nm)
  rowMeans(dt_numeric[, ..cols], na.rm = TRUE)
})

result_dt <- data.table::data.table(`Gene ID` = gene_ids, data.table::as.data.table(averaged_list))
data.table::setnames(result_dt, old = names(result_dt)[-1], new = unique_names)

# Replace NA with 0 (skip Gene ID)
if (ncol(result_dt) > 1) {
  result_dt[, (2:ncol(result_dt)) := lapply(.SD, function(x) data.table::fifelse(is.na(x), 0, x)), .SDcols = 2:ncol(result_dt)]
}

# Drop genes with zero expression across all groups
if (ncol(result_dt) > 1) {
  result_dt <- result_dt[rowSums(result_dt[, -1, with = FALSE] != 0) > 0]
}

# -----------------------
# Build numeric matrix
# -----------------------
mat <- as.matrix(result_dt[, !"Gene ID", with = FALSE])
rownames(mat) <- result_dt[["Gene ID"]]

# -----------------------
# Marker detection
# -----------------------
if (is_rnaseq) {
  require_or_fail("MGFR")
  markers.list <- getMarkerGenes.rnaseq(
    mat,
    class.vec       = colnames(mat),
    samples2compare = "all",
    annotate        = FALSE,
    gene.ids.type   = "ensembl",
    score.cutoff    = 1
  )
} else {
  message("Running MGFR-style specificity scoring for proteomics…")

  require_or_fail("tidyr")
  group_names <- colnames(mat)
  SPECIFICITY_MARGIN <- 0.02   # min gap to second-best to call exclusivity
  MIN_LOG2FC         <- 0.5
  eps <- 1e-8

  compute_specificity <- function(x) {
    sapply(seq_along(x), function(j) {
      target <- x[j]
      others <- if (length(x) > 1) mean(x[-j]) else 0
      # lower = more specific
      others / (target + others + eps)
    })
  }

  spec_mat <- t(apply(mat, 1, compute_specificity))
  colnames(spec_mat) <- paste0(group_names, "_markers")

  markers.list <- setNames(vector("list", length(group_names)), paste0(group_names, "_markers"))
  for (g in seq_along(group_names)) {
    gname  <- group_names[g]
    scores <- spec_mat[, g]
    exprs  <- mat[, g]
    keep   <- exprs > EXPRESSION_CUTOFF
    scores_f <- scores[keep]
    if (length(scores_f) == 0) {
      markers.list[[paste0(gname, "_markers")]] <- character(0)
      next
    }
    ord   <- order(scores_f, decreasing = FALSE, na.last = NA)
    ids   <- rownames(spec_mat)[keep][ord]
    scvec <- scores_f[ord]
    markers.list[[paste0(gname, "_markers")]] <- paste0(ids, " : ", format(scvec, scientific = FALSE, trim = TRUE))
  }

  parse_markers_list <- function(mk) {
    if (!length(mk)) return(data.frame())
    pieces <- lapply(names(mk), function(nm) {
      v <- mk[[nm]]
      if (!length(v)) return(NULL)
      spl <- strsplit(v, " : ")
      data.frame(
        FEATURE_ID = vapply(spl, `[`, "", 1),
        SPEC_SCORE = suppressWarnings(as.numeric(vapply(spl, `[`, "", 2))),
        GROUP_NAME = sub("_markers$", "", nm),
        stringsAsFactors = FALSE
      )
    })
    dplyr::bind_rows(pieces)
  }

  cand_df <- parse_markers_list(markers.list)
  if (nrow(cand_df) == 0) stop("No candidates in markers.list.", call. = FALSE)

  get_mean_others <- function(fid, gname) {
    others <- setdiff(colnames(mat), gname)
    if (!length(others)) 0 else mean(mat[fid, others], na.rm = TRUE)
  }

  cand_df <- cand_df |>
    dplyr::rowwise() |>
    dplyr::mutate(
      TARGET_EXPR = as.numeric(mat[FEATURE_ID, GROUP_NAME]),
      MEAN_OTHERS = get_mean_others(FEATURE_ID, GROUP_NAME),
      LOG2FC      = log2((TARGET_EXPR + eps) / (MEAN_OTHERS + eps))
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(
      !is.na(SPEC_SCORE),
      TARGET_EXPR > EXPRESSION_CUTOFF,
      SPEC_SCORE <= SPECIFICITY_SCORE_CUTOFF,
      LOG2FC >= MIN_LOG2FC
    )

  if (nrow(cand_df) == 0) stop("No candidates after filters; relax thresholds.", call. = FALSE)

  cand_df <- cand_df |>
    dplyr::arrange(FEATURE_ID, SPEC_SCORE, dplyr::desc(TARGET_EXPR), dplyr::desc(LOG2FC), GROUP_NAME) |>
    dplyr::group_by(FEATURE_ID) |>
    dplyr::mutate(RANK_WITHIN_FEAT = dplyr::row_number()) |>
    dplyr::ungroup()

  second_best_tbl <- cand_df |>
    dplyr::group_by(FEATURE_ID) |>
    dplyr::summarise(SECOND_BEST = dplyr::nth(SPEC_SCORE, 2, order_by = SPEC_SCORE), .groups = "drop")

  best_df <- cand_df |>
    dplyr::filter(RANK_WITHIN_FEAT == 1) |>
    dplyr::left_join(second_best_tbl, by = "FEATURE_ID") |>
    dplyr::mutate(
      SECOND_BEST = ifelse(is.na(SECOND_BEST), SPEC_SCORE + SPECIFICITY_MARGIN + 1, SECOND_BEST),
      MARGIN_OK   = (SECOND_BEST - SPEC_SCORE) >= SPECIFICITY_MARGIN
    ) |>
    dplyr::filter(MARGIN_OK) |>
    dplyr::select(-MARGIN_OK, -RANK_WITHIN_FEAT)

  markers.list <- setNames(vector("list", length(group_names)), paste0(group_names, "_markers"))
  ranked_df <- best_df |>
    dplyr::group_by(GROUP_NAME) |>
    dplyr::arrange(SPEC_SCORE, dplyr::desc(TARGET_EXPR), dplyr::desc(LOG2FC), FEATURE_ID, .by_group = TRUE) |>
    dplyr::mutate(RANKING = dplyr::row_number()) |>
    dplyr::ungroup()

  for (g in group_names) {
    subdf <- ranked_df |> dplyr::filter(GROUP_NAME == g)
    if (nrow(subdf) == 0) {
      markers.list[[paste0(g, "_markers")]] <- character(0)
    } else {
      markers.list[[paste0(g, "_markers")]] <- paste0(
        subdf$FEATURE_ID, " : ",
        format(subdf$SPEC_SCORE, scientific = FALSE, trim = TRUE)
      )
    }
  }
  message("Proteomics exclusivity applied.")
}

print(names(markers.list))

# -----------------------
# Sanity: ensure no gene appears in >1 list
# -----------------------
gene_lists <- lapply(markers.list, function(x) sub(" :.*", "", x))
gene_with_source <- unlist(gene_lists, use.names = FALSE)
sources <- rep(names(gene_lists), lengths(gene_lists))

dups_df <- data.frame(gene = gene_with_source, source = sources, stringsAsFactors = FALSE) |>
  dplyr::distinct(gene, source) |>
  dplyr::group_by(gene) |>
  dplyr::summarise(n_sources = dplyr::n(), .groups = "drop") |>
  dplyr::filter(n_sources > 1)

print(dups_df)
if (nrow(dups_df) > 0) stop("Error: Duplicated genes found across gene marker lists.", call. = FALSE)

# -----------------------
# Filter markers by specificity cutoff
# -----------------------
markers.list <- lapply(markers.list, function(v) {
  v <- v[grepl(":", v)]
  scr <- suppressWarnings(as.numeric(sub(".*: ", "", v)))
  v[!is.na(scr) & scr <= SPECIFICITY_SCORE_CUTOFF]
})
total_markers <- sum(vapply(markers.list, length, integer(1)))
if (total_markers == 0) {
  stop("No markers found after applying specificity score cutoff. Consider lowering the cutoff value.", call. = FALSE)
}
print(markers.list)

# -----------------------
# Build marker table (per group)
# -----------------------
marker_tables <- list()
for (nm in names(markers.list)) {
  v <- markers.list[[nm]]
  if (length(v) < 1) next
  parts <- strsplit(v, " : ")
  GENE_ID <- vapply(parts, `[`, "", 1)
  SPEC <- suppressWarnings(as.numeric(vapply(parts, `[`, "", 2)))
  group_label <- sub("_markers$", "", nm)
  grp <- summary_df |> dplyr::filter(label == group_label) |> dplyr::pull(group)
  n_assays <- summary_df |> dplyr::filter(label == group_label) |> dplyr::pull(n_assays)

  df <- data.frame(
    ACCESSION        = accession,
    GROUP            = grp,
    GROUP_NAME       = group_label,
    GENE_ID          = GENE_ID,
    SPECIFICITY_SCORE= SPEC,
    RANKING          = seq_along(GENE_ID),
    METRIC           = metric,
    NUMBER_SAMPLES   = n_assays,
    EXPRESSION       = NA_real_,
    stringsAsFactors = FALSE
  )
  marker_tables[[nm]] <- df
}

all_markers_df <- do.call(rbind, marker_tables)
print(head(all_markers_df))
print(dim(all_markers_df))

# -----------------------
# Complete to all gene x tissue combos
# -----------------------
tissue_map <- summary_df |> dplyr::select(GROUP = group, GROUP_NAME = label, NUMBER_SAMPLES = n_assays)
unique_genes <- unique(all_markers_df$GENE_ID)

complete_combos <- merge(
  expand.grid(GENE_ID = unique_genes, GROUP_NAME = tissue_map$GROUP_NAME, stringsAsFactors = FALSE),
  tissue_map,
  by = "GROUP_NAME",
  all.x = TRUE
)

existing_keys <- all_markers_df |> dplyr::select(GENE_ID, GROUP_NAME, GROUP, NUMBER_SAMPLES)

missing_combos <- dplyr::anti_join(complete_combos, existing_keys,
                                   by = c("GENE_ID", "GROUP_NAME", "GROUP", "NUMBER_SAMPLES"))

filled_rows <- missing_combos |>
  dplyr::mutate(
    SPECIFICITY_SCORE = -1,
    RANKING           = -1,
    METRIC            = metric,
    EXPRESSION        = NA_real_,
    ACCESSION         = accession
  )

final_df <- dplyr::bind_rows(all_markers_df, filled_rows) |>
  dplyr::arrange(GENE_ID, GROUP_NAME)

# -----------------------
# Add gene names and per-group expression from decorated matrix
# -----------------------
expr_decorated_path <- args[["expr_decorated"]]

# Use fread for robustness; handle both TSV and CSV
if (grepl("\\.tsv(\\.|$)", expr_decorated_path, ignore.case = TRUE)) {
  expression_decorated <- data.table::fread(expr_decorated_path)
} else {
  expression_decorated <- data.table::fread(expr_decorated_path) # still fread; sep auto-detected
}

print(head(expression_decorated))

if (!is_rnaseq) {
  # Map group codes to labels and clean decorated column names
  name_map <- setNames(summary_df$label, summary_df$group)
  new_names <- colnames(expression_decorated)

  for (code in names(name_map)) {
    pattern <- paste0("^", code, "(?=\\.|$)")
    new_names <- sub(pattern, name_map[[code]], new_names, perl = TRUE)
  }
  new_names <- sub("\\.WithInSampleAbundance", "", new_names)
  colnames(expression_decorated) <- new_names
}

find_col <- function(df, candidates, label, required = TRUE) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) > 0) return(hit[[1]])

  if (required) {
    stop(
      paste0(
        "Decorated expression file must contain ", label,
        ". Tried: ", paste(candidates, collapse = ", "),
        ". Available columns: ", paste(names(df), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  NA_character_
}

gene_id_col <- find_col(
  expression_decorated,
  c("GeneID", "Gene ID", "Gene.ID", "gene_id", "gene.id"),
  "gene ID column"
)

gene_name_col <- find_col(
  expression_decorated,
  c("Gene.Name", "Gene Name", "GeneName", "gene_name", "gene.name"),
  "gene name column",
  required = FALSE
)

if (is.na(gene_name_col)) {
  expression_decorated$GENE_NAME_FALLBACK <- ""
  gene_name_col <- "GENE_NAME_FALLBACK"
}

final_df <- final_df |>
  dplyr::left_join(
    data.frame(
      GENE_ID = expression_decorated[[gene_id_col]],
      GENE_NAME = expression_decorated[[gene_name_col]],
      stringsAsFactors = FALSE
    ),
    by = "GENE_ID"
  )

# Compute expression for each row (vectorized where possible)
final_df$EXPRESSION <- NA_real_

if (is_rnaseq) {
  # For RNA-seq, expression_decorated columns named by GROUP (assay code), values as "x,y,z" strings; take 3rd value.
  # Build a fast lookup for GeneID rows
  gid_index <- match(final_df$GENE_ID, expression_decorated[[gene_id_col]])

  # Safe extraction helper
  extract_third <- function(x) {
    vals <- suppressWarnings(as.numeric(strsplit(as.character(x), ",", fixed = TRUE)[[1]]))
    if (!length(vals)) return(NA_real_)
    if (length(vals) >= 3) vals[3] else tail(vals, 1)
  }

  # Iterate by row for varying GROUP columns (still O(n))
  for (i in seq_len(nrow(final_df))) {
    row_idx <- gid_index[i]
    if (is.na(row_idx)) next
    group_col <- final_df$GROUP[i]
    if (!is.na(group_col) && group_col %in% colnames(expression_decorated)) {
      final_df$EXPRESSION[i] <- extract_third(expression_decorated[row_idx, ..group_col][[1]])
    }
  }
} else {
  # Proteomics: columns named by GROUP_NAME, numeric single values
  gid_index <- match(final_df$GENE_ID, expression_decorated[[gene_id_col]])
  for (i in seq_len(nrow(final_df))) {
    row_idx <- gid_index[i]
    if (is.na(row_idx)) next
    gname <- final_df$GROUP_NAME[i]
    if (!is.na(gname) && gname %in% colnames(expression_decorated)) {
      final_df$EXPRESSION[i] <- suppressWarnings(as.numeric(expression_decorated[row_idx, ..gname][[1]]))
    }
  }
}

# -----------------------
# Apply expression cutoff & rerank
# -----------------------
genes_to_remove <- final_df |>
  dplyr::filter(RANKING != -1, EXPRESSION <= EXPRESSION_CUTOFF) |>
  dplyr::distinct(GENE_ID) |>
  dplyr::pull(GENE_ID)

filtered_df <- final_df |> dplyr::filter(!GENE_ID %in% genes_to_remove)

ranked_markers <- filtered_df |>
  dplyr::filter(RANKING != -1) |>
  dplyr::group_by(GROUP_NAME) |>
  dplyr::arrange(SPECIFICITY_SCORE, dplyr::desc(EXPRESSION), GENE_ID, .by_group = TRUE) |>
  dplyr::mutate(RANKING = dplyr::row_number()) |>
  dplyr::ungroup()

non_markers <- filtered_df |>
  dplyr::filter(RANKING == -1)

filtered_df <- dplyr::bind_rows(ranked_markers, non_markers)

ordered_genes <- filtered_df |>
  dplyr::filter(RANKING > 0) |>
  dplyr::group_by(GROUP_NAME) |>
  dplyr::arrange(RANKING, .by_group = TRUE) |>
  dplyr::reframe(GENE_ID = unique(GENE_ID)) |>
  dplyr::pull(GENE_ID)

filtered_df_sorted <- filtered_df |>
  dplyr::filter(GENE_ID %in% ordered_genes) |>
  dplyr::mutate(GENE_ID = factor(GENE_ID, levels = ordered_genes)) |>
  dplyr::arrange(GENE_ID, GROUP_NAME)

# Replace -1 with "NULL" strings for SQL loading
filtered_df_sorted <- filtered_df_sorted |>
  dplyr::mutate(
    RANKING = ifelse(RANKING == -1, "NULL", as.character(RANKING)),
    SPECIFICITY_SCORE = ifelse(SPECIFICITY_SCORE == -1, "NULL", as.character(SPECIFICITY_SCORE))
  )

# Final column order & names
colnames(filtered_df_sorted) <- c(
  "experiment_accession", "assay_id", "assay", "gene_id",
  "specificity_score", "marker_gene_rank", "expression_unit",
  "number_assays", "expression_level", "gene_name"
)

# Drop number_assays (per original intent)
filtered_df_sorted <- filtered_df_sorted[, !(names(filtered_df_sorted) %in% c("number_assays"))]

if (!is_rnaseq) {
  # Replace missing proteomics expression with a large sentinel as in original
  na_idx <- is.na(filtered_df_sorted$expression_level)
  if (any(na_idx)) filtered_df_sorted$expression_level[na_idx] <- 1000

  # Move assay_id to the end
  filtered_df_sorted <- filtered_df_sorted[c(setdiff(names(filtered_df_sorted), "assay_id"), "assay_id")]

  # If gene_name entirely missing, swap with gene_id
  if (all(is.na(filtered_df_sorted$gene_name))) {
    message("WARNING: Gene names missing - swapping gene_id/gene_name columns.")
    names(filtered_df_sorted)[names(filtered_df_sorted) == "gene_id"]   <- "temp_col_swap"
    names(filtered_df_sorted)[names(filtered_df_sorted) == "gene_name"] <- "gene_id"
    names(filtered_df_sorted)[names(filtered_df_sorted) == "temp_col_swap"] <- "gene_name"
  }
}

print(head(filtered_df_sorted))

# Ensure assay_id last (both branches)
filtered_df_sorted <- filtered_df_sorted[c(setdiff(names(filtered_df_sorted), "assay_id"), "assay_id")]
print(head(filtered_df_sorted))

# -----------------------
# Write output
# -----------------------
data.table::fwrite(
  filtered_df_sorted,
  file = args[["output_file"]],
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

message("Done. Wrote: ", args[["output_file"]])
