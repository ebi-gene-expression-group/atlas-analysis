#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))

usage <- paste(
  "Usage:",
  "plot_gsea_results.R <gsea.tsv> <plot.png> <plot.svg> <title> <gene_set_type> [top_n]",
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

if (is.na(top_n) || top_n < 1) {
  top_n <- 10
}

clean_title <- function(x) {
  lines <- trimws(strsplit(x, "\n", fixed = TRUE)[[1]])
  lines <- lines[nzchar(lines)]
  paste(strwrap(paste(lines, collapse = " "), width = 75), collapse = "\n")
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
      ggsave(svg_file, plot = plot, width = width, height = height, device = grDevices::svg, bg = "white")
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

plot_data$label <- wrap_label(plot_data$term)
duplicated_labels <- duplicated(plot_data$label) | duplicated(plot_data$label, fromLast = TRUE)
plot_data$label[duplicated_labels] <- paste0(
  plot_data$label[duplicated_labels],
  "\n",
  plot_data$accession[duplicated_labels]
)
plot_data$label <- factor(plot_data$label, levels = rev(plot_data$label))

subtitle <- paste0(
  "Top ",
  nrow(plot_data),
  " ",
  gene_set_type,
  " terms by adjusted p-value"
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
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 10, color = "grey15"),
    axis.text.x = element_text(color = "grey25"),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, color = "grey35", margin = margin(b = 12)),
    legend.position = "right",
    plot.margin = margin(18, 24, 18, 18)
  )

save_plot(dotplot, rows = nrow(plot_data))
