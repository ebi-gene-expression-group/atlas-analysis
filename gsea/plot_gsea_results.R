#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))

usage <- paste(
  "Usage:",
  "plot_gsea_results.R <gsea.tsv> <plot.png> <plot.svg> <title> <gene_set_type> [top_n] [gene_set_file] [significant_gene_file]",
  sep = "\n"
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop(usage, call. = FALSE)
}

gsea_file <- args[[1]]
png_file <- args[[2]]
svg_file <- args[[3]]
plot_title <- args[[4]]
gene_set_type <- args[[5]]
top_n <- if (length(args) >= 6) as.integer(args[[6]]) else 10
gene_set_file <- if (length(args) >= 7) args[[7]] else ""
significant_gene_file <- if (length(args) >= 8) args[[8]] else ""

if (is.na(top_n) || top_n < 1) {
  top_n <- 10
}

clean_title <- function(x) {
  lines <- trimws(strsplit(x, "\n", fixed = TRUE)[[1]])
  lines <- lines[nzchar(lines)]
  paste(strwrap(paste(lines, collapse = " "), width = 62), collapse = "\n")
}

gene_set_label <- function(x) {
  labels <- c(
    go = "GO terms",
    reactome = "Reactome pathways",
    interpro = "InterPro domains"
  )

  key <- tolower(x)
  if (key %in% names(labels)) {
    labels[[key]]
  } else {
    paste(x, "terms")
  }
}

wrap_label <- function(x, width = 48) {
  vapply(
    x,
    function(label) paste(strwrap(label, width = width), collapse = "\n"),
    character(1)
  )
}

parse_number <- function(x) {
  suppressWarnings(as.numeric(gsub(",", "", as.character(x), fixed = TRUE)))
}

find_column <- function(data, exact = character(), pattern = NULL) {
  column_names <- names(data)
  normalized <- tolower(gsub("[^a-z0-9]+", "", column_names))
  exact_normalized <- tolower(gsub("[^a-z0-9]+", "", exact))

  matched <- match(exact_normalized, normalized, nomatch = 0)
  if (any(matched > 0)) {
    return(column_names[matched[matched > 0][[1]]])
  }

  if (!is.null(pattern)) {
    pattern_match <- grep(pattern, normalized, perl = TRUE)
    if (length(pattern_match) > 0) {
      return(column_names[pattern_match[[1]]])
    }
  }

  NA_character_
}

looks_like_accession <- function(x) {
  grepl("^(GO:|R-[A-Z]+-|IPR[0-9]+)", as.character(x))
}

read_significant_genes <- function(path) {
  if (!nzchar(path) || !file.exists(path) || file.info(path)$size == 0) {
    return(character())
  }

  tokens <- unlist(strsplit(readLines(path, warn = FALSE), "\t", fixed = TRUE))
  tokens <- trimws(tokens)
  tokens <- tokens[nzchar(tokens)]
  tokens <- tokens[!tolower(tokens) %in% c("x", "v1", "gene", "genes")]
  unique(tokens)
}

read_gene_sets <- function(path, significant_genes, term_ids) {
  if (!nzchar(path) || !file.exists(path) || file.info(path)$size == 0) {
    return(list())
  }
  if (length(significant_genes) == 0 || length(term_ids) == 0) {
    return(list())
  }

  mapping <- tryCatch(
    read.table(path, sep = "\t", header = FALSE, quote = "\"", comment.char = "!", fill = TRUE, stringsAsFactors = FALSE),
    error = function(e) data.frame()
  )
  if (ncol(mapping) < 2 || nrow(mapping) == 0) {
    return(list())
  }

  mapping <- mapping[, 1:2, drop = FALSE]

  first_term_matches <- sum(as.character(mapping[[1]]) %in% term_ids, na.rm = TRUE)
  second_term_matches <- sum(as.character(mapping[[2]]) %in% term_ids, na.rm = TRUE)
  if (first_term_matches > second_term_matches) {
    mapping <- mapping[, c(2, 1), drop = FALSE]
  }

  colnames(mapping) <- c("gene", "term")
  mapping$gene <- as.character(mapping$gene)
  mapping$term <- as.character(mapping$term)
  mapping <- mapping[mapping$term %in% term_ids & mapping$gene %in% significant_genes, , drop = FALSE]

  split(mapping$gene, mapping$term)
}

find_overlap_groups <- function(plot_data, gene_sets, min_shared = 3, min_jaccard = 0.20, min_overlap = 0.50) {
  n_terms <- nrow(plot_data)
  if (n_terms < 2 || length(gene_sets) == 0) {
    return(rep(NA_integer_, n_terms))
  }

  adjacency <- matrix(FALSE, n_terms, n_terms)
  for (i in seq_len(n_terms - 1)) {
    for (j in (i + 1):n_terms) {
      genes_i <- unique(gene_sets[[plot_data$accession[[i]]]])
      genes_j <- unique(gene_sets[[plot_data$accession[[j]]]])
      if (length(genes_i) == 0 || length(genes_j) == 0) {
        next
      }

      shared <- length(intersect(genes_i, genes_j))
      union_size <- length(union(genes_i, genes_j))
      smaller_size <- min(length(genes_i), length(genes_j))
      jaccard <- if (union_size > 0) shared / union_size else 0
      overlap <- if (smaller_size > 0) shared / smaller_size else 0

      if (shared >= min_shared && (jaccard >= min_jaccard || overlap >= min_overlap)) {
        adjacency[i, j] <- TRUE
        adjacency[j, i] <- TRUE
      }
    }
  }

  groups <- rep(NA_integer_, n_terms)
  group_id <- 0L
  for (start in seq_len(n_terms)) {
    if (!is.na(groups[[start]]) || !any(adjacency[start, ])) {
      next
    }

    group_id <- group_id + 1L
    queue <- start
    groups[[start]] <- group_id
    while (length(queue) > 0) {
      current <- queue[[1]]
      queue <- queue[-1]
      neighbours <- which(adjacency[current, ] & is.na(groups))
      if (length(neighbours) > 0) {
        groups[neighbours] <- group_id
        queue <- c(queue, neighbours)
      }
    }
  }

  groups
}

make_placeholder_plot <- function(message) {
  ggplot() +
    annotate("text", x = 0, y = 0, label = message, size = 5) +
    xlim(-1, 1) +
    ylim(-1, 1) +
    labs(title = clean_title(plot_title)) +
    theme_void(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, margin = margin(b = 18)),
      plot.margin = margin(20, 20, 20, 20)
    )
}

save_plot <- function(plot, rows = top_n) {
  height <- max(4.5, 2.4 + 0.45 * rows)
  width <- 9

  ggsave(png_file, plot = plot, width = width, height = height, dpi = 300, bg = "white")

  svg_result <- tryCatch({
    if (requireNamespace("svglite", quietly = TRUE)) {
      ggsave(svg_file, plot = plot, width = width, height = height, device = svglite::svglite, bg = "white")
    } else {
      FALSE
    }
    file.exists(svg_file) && file.info(svg_file)$size > 0
  }, error = function(e) {
    FALSE
  })

  if (!isTRUE(svg_result)) {
    writeLines(
      c(
        '<svg xmlns="http://www.w3.org/2000/svg" width="900" height="450" viewBox="0 0 900 450">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<text x="450" y="210" text-anchor="middle" font-family="sans-serif" font-size="24">SVG export unavailable</text>',
        '<text x="450" y="250" text-anchor="middle" font-family="sans-serif" font-size="16">See the PNG dot plot generated alongside this file.</text>',
        '</svg>'
      ),
      con = svg_file
    )
  }
}

if (!file.exists(gsea_file) || file.info(gsea_file)$size == 0) {
  save_plot(make_placeholder_plot("No enrichment results available"), rows = 1)
  quit(save = "no", status = 0)
}

gsea <- tryCatch(
  read.delim(gsea_file, check.names = FALSE, stringsAsFactors = FALSE),
  error = function(e) data.frame()
)

if (nrow(gsea) == 0 || ncol(gsea) < 2) {
  save_plot(make_placeholder_plot("No enrichment results available"), rows = 1)
  quit(save = "no", status = 0)
}

if (ncol(gsea) == 1 && any(grepl("no results", gsea[[1]], ignore.case = TRUE))) {
  save_plot(make_placeholder_plot("No enrichment results passed the selected thresholds"), rows = 1)
  quit(save = "no", status = 0)
}

term_column <- find_column(gsea, exact = c("Term", "Name"))
accession_column <- find_column(gsea, exact = c("Accession", "ID"))

if (is.na(term_column)) {
  term_column <- names(gsea)[[1]]
}
if (is.na(accession_column) && ncol(gsea) >= 2) {
  accession_column <- names(gsea)[[2]]
}

term_values <- as.character(gsea[[term_column]])
accession_values <- if (!is.na(accession_column)) as.character(gsea[[accession_column]]) else term_values

first_is_accession <- mean(looks_like_accession(term_values), na.rm = TRUE) > 0.5
second_is_accession <- mean(looks_like_accession(accession_values), na.rm = TRUE) > 0.5

if (first_is_accession && !second_is_accession) {
  accessions <- term_values
  terms <- accession_values
} else {
  terms <- term_values
  accessions <- accession_values
}

terms[!nzchar(terms) | is.na(terms)] <- accessions[!nzchar(terms) | is.na(terms)]

effect_column <- find_column(
  gsea,
  exact = c("effect.size", "effect size", "Observed/Expected", "Observed Expected"),
  pattern = "(effect|observed.*expected)"
)
padj_column <- find_column(
  gsea,
  exact = c("p adj (non-dir.)", "Adjusted.p.value", "Adjusted p-value", "FDR", "padj"),
  pattern = "(padj|adjusted.*p|fdr)"
)
sig_column <- find_column(
  gsea,
  exact = c("Significant..in.gene.set.", "Significant (in gene set)", "Significant in gene set"),
  pattern = "significant.*in.*gene.*set"
)
genes_total_column <- find_column(
  gsea,
  exact = c("Genes (tot)", "Genes tot", "Genes"),
  pattern = "genes.*tot"
)

if (is.na(effect_column)) {
  save_plot(make_placeholder_plot("No observed/expected enrichment column found"), rows = 1)
  quit(save = "no", status = 0)
}

effect_size <- parse_number(gsea[[effect_column]])
padj <- if (!is.na(padj_column)) parse_number(gsea[[padj_column]]) else NA_real_
sig_genes <- if (!is.na(sig_column)) parse_number(gsea[[sig_column]]) else NA_real_
genes_total <- if (!is.na(genes_total_column)) parse_number(gsea[[genes_total_column]]) else NA_real_

plot_data <- data.frame(
  accession = accessions,
  term = terms,
  effect_size = effect_size,
  padj = padj,
  sig_genes = sig_genes,
  genes_total = genes_total,
  stringsAsFactors = FALSE
)

plot_data <- plot_data[is.finite(plot_data$effect_size), , drop = FALSE]
if (nrow(plot_data) == 0) {
  save_plot(make_placeholder_plot("No plottable enrichment results available"), rows = 1)
  quit(save = "no", status = 0)
}

if (all(is.na(plot_data$padj))) {
  plot_data$padj <- 1
}
if (all(is.na(plot_data$sig_genes))) {
  plot_data$sig_genes <- ifelse(is.na(plot_data$genes_total), 1, plot_data$genes_total)
}

plot_data <- plot_data[order(plot_data$padj, -plot_data$effect_size), , drop = FALSE]
plot_data <- head(plot_data, top_n)
plot_data$original_rank <- seq_len(nrow(plot_data))

significant_genes <- read_significant_genes(significant_gene_file)
gene_sets <- read_gene_sets(gene_set_file, significant_genes, plot_data$accession)
plot_data$overlap_group <- find_overlap_groups(plot_data, gene_sets)
has_overlap_groups <- any(!is.na(plot_data$overlap_group))

if (has_overlap_groups) {
  group_order <- aggregate(
    plot_data$padj,
    by = list(group = ifelse(is.na(plot_data$overlap_group), paste0("single_", plot_data$original_rank), paste0("group_", plot_data$overlap_group))),
    FUN = function(x) {
      x <- x[is.finite(x)]
      if (length(x) == 0) Inf else min(x)
    }
  )
  names(group_order)[[2]] <- "group_min_padj"
  plot_data$sort_group <- ifelse(is.na(plot_data$overlap_group), paste0("single_", plot_data$original_rank), paste0("group_", plot_data$overlap_group))
  plot_data <- merge(plot_data, group_order, by.x = "sort_group", by.y = "group", all.x = TRUE, sort = FALSE)
  plot_data <- plot_data[order(plot_data$group_min_padj, plot_data$padj, -plot_data$effect_size), , drop = FALSE]
}

plot_data$label <- wrap_label(plot_data$term)
duplicated_labels <- duplicated(plot_data$label) | duplicated(plot_data$label, fromLast = TRUE)
plot_data$label[duplicated_labels] <- paste0(
  plot_data$label[duplicated_labels],
  "\n",
  plot_data$accession[duplicated_labels]
)
plot_data$label <- factor(plot_data$label, levels = rev(plot_data$label))

strip_data <- plot_data[!is.na(plot_data$overlap_group), , drop = FALSE]
if (nrow(strip_data) > 0) {
  strip_data$overlap_group <- factor(strip_data$overlap_group, levels = sort(unique(strip_data$overlap_group)))
}

x_min <- min(c(0.9, plot_data$effect_size), na.rm = TRUE)
x_max <- max(c(1, plot_data$effect_size), na.rm = TRUE)
x_range <- max(x_max - x_min, 0.1)
x_strip <- x_min + (x_range * 0.025)
x_strip_width <- x_range * 0.012

overlap_palette <- c(
  "#0072B2", "#D55E00", "#009E73", "#CC79A7",
  "#E69F00", "#56B4E9", "#7E57C2", "#666666"
)
if (nrow(strip_data) > length(overlap_palette)) {
  overlap_palette <- grDevices::colorRampPalette(overlap_palette)(nrow(strip_data))
}
overlap_colors <- overlap_palette[seq_len(max(1, length(levels(strip_data$overlap_group))))]

subtitle <- paste0(
  "Top ",
  nrow(plot_data),
  " ",
  gene_set_label(gene_set_type),
  " by adjusted p-value"
)

dotplot <- ggplot(plot_data, aes(x = effect_size, y = label)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey65", linewidth = 0.4) +
  geom_segment(
    aes(x = 1, xend = effect_size, y = label, yend = label),
    color = "grey82",
    linewidth = 1.1
  ) +
  geom_point(aes(size = sig_genes, color = padj), alpha = 0.9) +
  scale_color_gradient(low = "#3b4cc0", high = "#f6b26b", name = "FDR") +
  scale_size_area(max_size = 9, name = "Significant genes") +
  labs(
    title = clean_title(plot_title),
    subtitle = subtitle,
    x = "Observed / expected",
    y = NULL,
    caption = if (has_overlap_groups) "Matching coloured strips mark terms with meaningful significant-gene overlap." else NULL
  ) +
  scale_x_continuous(limits = c(x_min, x_max + x_range * 0.06)) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 10, color = "grey15"),
    axis.text.x = element_text(color = "grey25"),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
    plot.subtitle = element_text(hjust = 0.5, color = "grey35", margin = margin(b = 12)),
    plot.caption = element_text(color = "grey40", size = 9),
    legend.position = "right",
    plot.margin = margin(18, 24, 18, 18)
  )

if (nrow(strip_data) > 0) {
  dotplot <- dotplot +
    geom_tile(
      data = strip_data,
      aes(x = x_strip, y = label, fill = overlap_group),
      width = x_strip_width,
      height = 0.72,
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    scale_fill_manual(values = overlap_colors)
}

save_plot(dotplot, rows = nrow(plot_data))
