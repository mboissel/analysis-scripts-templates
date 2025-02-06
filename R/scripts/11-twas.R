message(timestamp(quiet = TRUE))
### Project Setup ==================================================================================
library(here)
project_name <- sub("(.*)_*\\..*", "\\1", list.files(here(), pattern = ".Rproj$"))
output_directory <- here("outputs", "11-twas")
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE, mode = "0775")

phenotype <- here("docs", "phenotype.xlsx")

data_directory <- file.path("/disks/DATA/Projects", project_name, "QC")
ethnicity <- file.path(data_directory, "Omni2.5", paste0(project_name, "_ethnicity.csv"))
exclusion <- file.path(data_directory, "Omni2.5", paste0(project_name, "_QC_exclusion.xlsx"))

run_directory <- "/disks/RUN/Run_XXX/Output/RSEM"


### Load Packages ==================================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(tximport)
  library(DESeq2)
  library(biomaRt)
  library(httr)
  library(flashpcaR)
  library(ggplot2)
  library(ragg)
  library(patchwork)
  library(ggtext)
  library(scales)
  library(ggrepel)
})


### Tables and Figures Theme =======================================================================
options(
  ggplot2.discrete.colour = function(...) scale_colour_viridis_d(..., begin = 0.15, end = 0.85),
  ggplot2.discrete.fill = function(...) scale_fill_viridis_d(..., begin = 0.15, end = 0.85),
  ggplot2.continuous.colour = function(...) scale_colour_viridis_c(..., begin = 0.15, end = 0.85),
  ggplot2.continuous.fill = function(...) scale_fill_viridis_c(..., begin = 0.15, end = 0.85)
)
theme_set(theme_minimal(base_family = "Tex Gyre Termes"))
theme_update(
  plot.title.position = "plot",
  plot.caption.position = "plot",
  plot.title = element_markdown(),
  plot.subtitle = element_markdown(face = "italic", size = rel(0.80)),
  plot.caption = element_markdown(face = "italic", size = rel(0.65)),
  axis.title.x = element_markdown(),
  axis.text.x = element_markdown(),
  axis.title.y = element_markdown(),
  axis.text.y = element_markdown()
)


### Functions ======================================================================================
find_tissue <- function(
  files,
  type = "rsem", # or "matrix"
  tpm_threshold = 100,
  # url = "https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz" # ancien chemin
  url = "https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"
) {

  gtex_raw_data <- data.table::fread(url)[, -c("Description")]

  gtex_markers <- data.table::melt(
    data = gtex_raw_data[j = "Name" := gsub("\\..*$", "", Name)],
    id.vars = "Name",
    variable.name = "Tissue",
    value.name = "TPM"
  )[TPM >= tpm_threshold][j = Tissue := as.character(Tissue)]

  gtex_markers <- gtex_markers[
    i = gtex_markers[j = .SD[j = data.table::uniqueN(Tissue)], by = "Name"][V1 %in% 1:3],
    j = .SD,
    on = c("Name")
  ]

  if (type == "rsem") {
    samples_tissue <- Reduce(
      f = function(x, y) {
        merge(x, y, by = "Tissue", all = TRUE)[
          j =  `:=`(gene_id = sub("\\..*$", "", gene_id))
        ]
      },
      x = lapply(
        X = files,
        .gtex = gtex_markers,
        FUN = function(.file, .gtex) {
          out <- merge(
            x = .gtex,
            y = data.table::fread(
              file = .file,
              select = c("gene_id", "TPM"),
              col.names = c("Name", "sample_tpm")
            ),
            by = "Name",
            all = FALSE
          )[
            i = sample_tpm >= tpm_threshold & TPM >= sample_tpm,
            j = .N,
            by = "Tissue"
          ]
          data.table::setnames(out, old = "N", new = gsub(".genes.results$", "", basename(.file)))
        }
      )
    )

    class(samples_tissue) <- c("tissue", class(samples_tissue))

    samples_tissue <- data.table::setnafill(
      x = samples_tissue,
      fill = 0,
      cols = setdiff(names(samples_tissue), "Tissue")
    )[
      j = .(
        Tissue = Tissue[order(- rowMeans(.SD) / sum(rowMeans(.SD)))],
        .SD[order(- rowMeans(.SD) / sum(rowMeans(.SD)))]
      ),
      .SDcols = gsub(".genes.results$", "", basename(files))
    ]

    return(samples_tissue)

  }

  if (type == "matrix") {
    mat <- data.table::fread(files)
    names(mat)[1] <- "gene_id"
    mat$gene_id <- sub("\\..*$", "", mat$gene_id)
    samples_names <- names(mat)[-1]

    samples_tissue_tmp <- merge(gtex_markers, mat, by.x = "Name", by.y = "gene_id", all.y = TRUE)
    samples_tissue_tmp <- samples_tissue_tmp[!is.na(Tissue), ]
    samples_tissue_tmp[1:5, 1:5]

    samples_tissue <- Reduce(
      f = function(x, y) {
        merge(x, y, by = "Tissue", all = TRUE)
      },
      x = lapply(
        X = samples_names,
        .gtex = gtex_markers,
        FUN = function(.id, .gtex) {
          out <- samples_tissue_tmp[, .SD, .SDcols = c("Tissue", "TPM", .id)]
          data.table::setnames(out, old = .id, new = "sample_tpm")
          out <- out[
            i = sample_tpm >= tpm_threshold & TPM >= sample_tpm,
            j = .N,
            by = "Tissue"
          ]
          data.table::setnames(out, old = "N", new = .id)
          return(out)
        }
      )
    )
    class(samples_tissue) <- c("tissue", class(samples_tissue))

    samples_tissue <- data.table::setnafill(
      x = samples_tissue,
      fill = 0,
      cols = setdiff(names(samples_tissue), "Tissue")
    )[
      j = .(
        Tissue = Tissue[order(- rowMeans(.SD) / sum(rowMeans(.SD)))],
        .SD[order(- rowMeans(.SD) / sum(rowMeans(.SD)))]
      ),
      .SDcols = samples_names
    ]
    samples_tissue[1:5, 1:10]

    return(samples_tissue)
  }
}

is.unique <- function(object) {
  length(unique(summary.tissue(object)[["best_tissue"]])) == 1
}

summary.tissue <- function(object, ...) {
  object[
    j = data.table::data.table(
      t(sapply(.SD, function(x) c(best_tissue = Tissue[which.max(x)], n_gene = x[which.max(x)]))),
      keep.rownames = "sample_id"
    ),
    .SDcols = setdiff(names(object), "Tissue")
  ]
}

pval_trans <- function(alpha = NULL, md = FALSE, prefix = FALSE, colour = "#ee2c2c") {
  scales::trans_new(
    name = "pval",
    domain = c(0, 1),
    transform = function(x) {
      x[x < .Machine$double.xmin] <- .Machine$double.xmin
      -log(x, 10)
    },
    inverse = function(x) 10^-x,
    breaks = (function(n = 5) {
      function(x) {
        max <- floor(-log(min(c(x, alpha), na.rm = TRUE), base = 10))
        if (max == 0) 1 else sort(unique(c(10^-seq(0, max, by = floor(max / n) + 1), alpha)))
      }
    })(),
    format = (function(x) {
      if (md & nchar(system.file(package = "ggtext")) != 0) {
        prefix_text <- if (prefix) "&alpha; = " else ""
        x_fmt <- gsub(
          "^(.*)e[+]*([-]*)0*(.*)$",
          "\\1 &times; 10<sup>\\2\\3</sup>",
          format(x, scientific = TRUE)
        )
        x_fmt[x %in% c(0, 1)] <- x[x %in% c(0, 1)]
        x_fmt <- gsub("^1 &times; ", "", x_fmt)
        alpha_idx <- format(x, scientific = TRUE) == format(alpha, scientific = TRUE)
        x_fmt[alpha_idx] <- paste0("<b style='color:", colour, ";'>", prefix_text, x_fmt[alpha_idx], "</b>")
        x_fmt
      } else {
        prefix_text <- if (prefix) "alpha == " else ""
        x_fmt <- gsub(
          "^(.*)e[+]*([-]*)0*(.*)$",
          "\\1 %*% 10^\\2\\3",
          format(x, scientific = TRUE)
        )
        x_fmt[x %in% c(0, 1)] <- x[x %in% c(0, 1)]
        x_fmt <- gsub("^1 \\%\\*\\% ", "", x_fmt)
        alpha_idx <- format(x, scientific = TRUE) == format(alpha, scientific = TRUE)
        x_fmt[alpha_idx] <- paste0(prefix_text, x_fmt[alpha_idx])
        parse(text = x_fmt)
      }
    })
  )
}

draw_volcano <- function(
  data, x, y,
  label_x = "Fold-Change (log<sub>2</sub>)",
  label_y = "P-value",
  label_colour = label_x,
  alpha = 0.05,
  max_p = 1
) {
  data <- data.table::as.data.table(data)#[, .SD, .SDcols = c(x, y)]
  data.table::setnames(data, y, "pvalue", skip_absent = TRUE)
  if (!"contrast" %in% names(data)) data[j = contrast := "Volcano plot"]

  if ("OR" %in% x) {
    label_x <- gsub("Fold Change", "OR", label_x)
    data[j = c(x) := log2(.SD), .SDcols = x]
  }

  out <- lapply(unique(data[["contrast"]]), function(.contrast) {
    ggplot2::ggplot(data[contrast %in% .contrast][pvalue <= max_p]) +
      ggplot2::aes(x = .data[[x]], y = .data[["pvalue"]], colour = abs(.data[[x]])) +
      ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "black") +
      ggplot2::geom_point(size = 0.50) +
      ggrepel::geom_label_repel(
        data = ~ .x[!is.na(external_gene_name) & fdr < 0.05][order(pvalue)][1:10], # max 10
        mapping = aes(label = .data[["external_gene_name"]]),
        min.segment.length = unit(0, "lines"),
        size = 1.5,
        show.legend = FALSE,
        na.rm = TRUE,
        nudge_y = 0.5,
        fontface = "italic",
        segment.colour = "black",
        colour = "black"
      ) +
      ggplot2::scale_colour_viridis_c(
        trans = "sqrt",
        limits = c(0, data[pvalue <= 0.05, max(abs(.SD)), .SDcols = x])
      ) +
      ggplot2::labs(x = label_x, y = label_y, colour = label_colour) +
      ggplot2::theme(legend.position = "none") +
      ggplot2::annotate(
        geom = "rect",
        xmin = -Inf, xmax = Inf, ymin = 1, ymax = alpha,
        fill = "#ee2c2c", alpha = 0.2, colour = NA
      ) +
      ggplot2::geom_hline(yintercept = alpha, linetype = 2, colour = "#ee2c2c") +
      ggplot2::scale_y_continuous(
        trans = pval_trans(alpha = alpha, md = TRUE, colour = "#ee2c2c"),
        expand = ggplot2::expansion(mult = c(0, 0.2)),
        limits = c(max_p, NA)
      )
  })
  names(out) <- unique(data[["contrast"]])
  out
}

plot_planes <- function(pca_dfxy, ivar, fig_n_comp) {
  apply(
    X = combn(sprintf("PC%02d", seq_len(fig_n_comp)), 2),
    MARGIN = 2,
    FUN = function(.x) {
      ggplot(data = pca_dfxy[, .SD, .SDcols = c(ivar, .x)]) +
        aes(x = .data[[.x[1]]], y = .data[[.x[2]]], colour = .data[[ivar]]) +
        geom_hline(yintercept = 0, linetype = 2, na.rm = TRUE) +
        geom_vline(xintercept = 0, linetype = 2, na.rm = TRUE) +
        geom_point(size = 0.75, na.rm = TRUE) +
        {
          if (is.numeric(pca_dfxy[[ivar]])) {
            scale_colour_viridis_c(
              name = NULL,
              begin = 0,
              end = 0.75
            )
          } else {
            list(
              stat_ellipse(type = "norm", na.rm = TRUE, show.legend = FALSE),
              scale_colour_viridis_d(
                name = NULL,
                begin = if (pca_dfxy[, uniqueN(.SD), .SDcols = ivar] == 2) 0.25 else 0,
                end = 0.75,
                guide = guide_legend(override.aes = list(size = 4))
              ),
              if (length(unique(pca_dfxy[[ivar]])) > 10) {
                theme(legend.position = "none")
              } else {
                NULL
              }
            )
          }
        } +
        theme_minimal(base_size = 10, base_family = "Tex Gyre Termes") +
        theme(
          plot.title.position = "plot",
          plot.caption.position = "plot",
          plot.title = element_markdown(),
          plot.subtitle = element_markdown(face = "italic", size = rel(0.80)),
          plot.caption = element_markdown(face = "italic", size = rel(0.65)),
          axis.title.x = element_markdown(),
          axis.text.x = element_markdown(),
          axis.title.y = element_markdown(),
          axis.text.y = element_markdown(),
          panel.grid.minor = element_blank()
        )
    }
  )
}


### Setup biomaRt ==================================================================================
set_config(config(ssl_verifypeer = FALSE)) # Fix SSL check error

species <- "hsapiens_gene_ensembl"
build <- 38
version <- 103

get_mart <- quote(useEnsembl(
  biomart = "ensembl",
  dataset = species,
  version = version,
  GRCh = if (build == 37) build else NULL
))

mart <- try(eval(get_mart), silent = TRUE)
if (inherits(mart, "try-error")) mart <- eval(get_mart)
ensembl_version <- sprintf("GRCh%d-%d", build, version)


### Get QCed Phenotypes ============================================================================
sample_sheet_qc <- merge(
  x = setDT(read_excel(phenotype)),
  y = fread(ethnicity),
  by = "iid"
)[j = `:=`(
  "Sample_ID" = iid,
  "#IID" = iid,
  "vcf_id" = iid,
  "group" = as.integer(groupe) + 1L,
  "sex" = c("M" = 1L, "F" = 2L)[sex],
  "age" = age,
  "bmi" = bmi,
  "rnaseq_path" = file.path(run_directory, "...")
)]

relatedness <- setDT(read_excel(exclusion, 2))[gsub("_.*", "", IID1) != gsub("_.*", "", IID2)]
relatedness[j = paste0("vcf_id", 1:2) := list(paste0(FID1, "_", IID1), paste0(FID2, "_", IID2))]

sample_sheet_qc <- sample_sheet_qc[ # not related
  !vcf_id %in% unique(unlist(relatedness[j = c("vcf_id1", "vcf_id2")], use.names = FALSE))
][ # European descent
  super_pop_closest %in% "EUR"
]


### Traits to be analysed ==========================================================================
traits <- c(
  "Control/Case" = "group"
)

covariates <- list(
  "group" = c("sex", "age", "bmi", sprintf("PC%02d", 1:2))
)


### Check Tissue ===================================================================================
sample_tissue <- "Liver"
res_tissue <- summary(find_tissue(sample_sheet_qc[["rnaseq_path"]]))
if (!all(grepl(sample_tissue, unique(res_tissue[["best_tissue"]])))) {
  message(
    paste0('Tissue found using GTEx does not match "', sample_tissue, '" for all samples!'),
    appendLF = TRUE
  )
  fwrite(
    x = res_tissue[!best_tissue %in% sample_tissue, gsub("\\.", "-", sample_id)],
    file = file.path(output_directory, "tissue_issue.csv")
  )
  sample_sheet_qc <- sample_sheet_qc[
    i = Sample_ID %in% res_tissue[best_tissue %in% sample_tissue, gsub("\\.", "-", sample_id)]
  ]
}


### Analysis =======================================================================================
for (rna_level in c("genes", "isoforms")) {
  rna_level_name <- unname(c("genes" = "gene", "isoforms" = "transcript")[rna_level])

  ensembl_id <- sprintf("ensembl_%s_id", rna_level_name)

  rsem_files <- setNames(sample_sheet_qc[["rnaseq_path"]], sample_sheet_qc[["Sample_ID"]])
  txi_counts <- tximport(
    files = gsub("\\.genes\\.results$", paste0(".", rna_level, ".results"), rsem_files),
    type = "rsem",
    txIn = rna_level == "isoforms",
    txOut = rna_level == "isoforms",
    countsFromAbundance = "no"
  )
  txi_counts$length[txi_counts$length == 0] <- 1

  ## PCA -------------------------------------------------------------------------------------------
  n_comp <- 3
  fig_n_comp <- 3
  keep_technical <- c("Status" = "status", "Sex" = "sex", "Age" = "age", "BMI" = "bmi")

  pca_data <- txi_counts$counts
  pca_phenotypes <- sample_sheet_qc[Sample_ID %in% colnames(pca_data)]
  pca_res <- flashpca(X = t(pca_data), stand = "sd", ndim = n_comp)
  pca_dfxy <- as.data.table(pca_res[["vectors"]], keep.rownames = "Sample_ID")
  setnames(
    x = pca_dfxy,
    old = setdiff(names(pca_dfxy), "Sample_ID"),
    new = sprintf("PC%02d", as.numeric(gsub("V", "", setdiff(names(pca_dfxy), "Sample_ID"))))
  )
  pca_dfxy <- merge(x = pca_dfxy, y = pca_phenotypes, by = "Sample_ID")
  p_inertia <- ggplot(
    data = data.table(
      y = pca_res[["pve"]],
      x = sprintf("PC%02d", seq_along(pca_res[["pve"]]))
    )[x %in% sprintf("PC%02d", seq_len(fig_n_comp))]
  ) +
    aes(
      x = paste0(x, "<br><i style='font-size:6pt;'>(", percent_format(accuracy = 0.01, suffix = " %")(y), ")</i>"),
      y = y
    ) +
    geom_col(
      width = 1,
      colour = "white",
      fill = viridis_pal(begin = 0.5, end = 0.5)(1),
      na.rm = TRUE
    ) +
    scale_y_continuous(
      labels = percent_format(accuracy = 0.1, suffix = " %"),
      expand = expansion(mult = c(0, 0.05))
    ) +
    labs(
      x = "Principal Components",
      y = "Contribution"
    ) +
    theme_minimal(base_size = 10, base_family = "Tex Gyre Termes") +
    theme(
      plot.title.position = "plot",
      plot.caption.position = "plot",
      plot.title = element_markdown(),
      plot.subtitle = element_markdown(face = "italic", size = rel(0.80)),
      plot.caption = element_markdown(face = "italic", size = rel(0.65)),
      axis.title.x = element_markdown(),
      axis.text.x = element_markdown(),
      axis.title.y = element_markdown(),
      axis.text.y = element_markdown(),
      panel.grid.minor = element_blank()
    ) +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())

  asso_dt <- melt(
    data = pca_dfxy,
    measure.vars = grep("^PC[0-9]+$", names(pca_dfxy), value = TRUE),
    variable.name = "pc",
    value.name = "values"
  )[pc %in% sprintf("PC%02d", seq_len(n_comp))][
    j = {
      m <- model.matrix(
        object = as.formula(
          object = paste0("values ~ ", paste(sprintf("`%s`", keep_technical), collapse = " + "))
        ),
        data = .SD
      )

      if (qr(m)$rank == ncol(m)) {
        out <- as.data.table(
          anova(
            lm(
              formula = as.formula(
                object = paste0("values ~ ", paste(sprintf("`%s`", keep_technical), collapse = " + "))
              ),
              data = .SD
            )
          ),
          keep.rownames = "term"
        )[term != "Residuals"]
      } else {
        out <- rbindlist(
          lapply(X = sprintf("`%s`", keep_technical), .data = .SD, FUN = function(.x, .data) {
            as.data.table(
              anova(
                lm(
                  formula = tats::as.formula(paste0("values ~ ", .x)),
                  data = .SD
                )
              ),
              keep.rownames = "term"
            )[term != "Residuals"]
          })
        )
      }
      out[, full_rank := qr(m)$rank == ncol(m)][j = term := gsub("`", "", term)]
    },
    by = "pc"
  ]

  pca_asso <- ggplot(data = asso_dt) +
    aes(
      x = factor(.data[["pc"]]),
      y = factor(
        x = .data[["term"]],
        levels = setorderv(
          x = dcast(
            data = asso_dt[j = list(pc, term, `Pr(>F)` = fifelse(`Pr(>F)` <= 0.1, `Pr(>F)`, NA_real_))],
            formula = term ~ pc,
            value.var = "Pr(>F)"
          ),
          cols = levels(asso_dt[["pc"]])[seq_len(n_comp)],
          order = -1
        )[["term"]]
      ),
      fill = .data[["Pr(>F)"]]
    ) +
    geom_tile(colour = "white", na.rm = TRUE) +
    geom_richtext(
      mapping = aes(
        label = (function(x) {
          x_fmt <- x
          x_fmt[as.numeric(x) < 0.01] <- gsub(
            pattern = "(.*)e([-+]*)0*(.*)",
            replacement = "\\1<br>&times;<br>10<sup>\\2\\3</sup>",
            x = format(x_fmt[as.numeric(x) < 0.01], digits = 2, nsmall = 2, scientific = TRUE)
          )
          x_fmt[as.numeric(x) >= 0.01] <- format(
            x = as.numeric(x_fmt[as.numeric(x) >= 0.01]),
            digits = 2, nsmall = 2
          )
          x_fmt
        })(.data[["Pr(>F)"]])
      ),
      colour = "white",
      fill = NA,
      label.colour = NA,
      size = 2.5,
      na.rm = TRUE
    ) +
    scale_fill_viridis_c(end = 0.75, limits = c(0, 0.1)) +
    theme(panel.grid = element_blank()) +
    scale_x_discrete(
      expand = c(0, 0),
      labels = function(x) {
        paste0(
          x, "<br><i style='font-size:5pt;'>(",
          format(
            x = pca_res[["pve"]][as.numeric(gsub("PC", "", x))] * 100,
            digits = 2,
            nsmall = 3
          ),
          " %)</i>"
        )
      }
    ) +
    scale_y_discrete(
      expand = c(0, 0),
      labels = function(x) names(keep_technical[match(gsub("`", "", x), keep_technical)])
    ) +
    labs(
      x = "Principal Components",
      y = "Variables",
      title = "Association Tests Between Variables And Principal Components",
      caption = ifelse(
        test = all(asso_dt[["full_rank"]]),
        yes = "Variables are tested against principal components using ANOVA.",
        no = paste(
          "Variables are independently tested against principal components using ANOVA",
          "(*i.e.*, model matrix is not full rank)."
        )
      ),
      fill = "P-Value"
    ) +
    theme_minimal(base_size = 10, base_family = "Tex Gyre Termes") +
    theme(
      plot.title.position = "plot",
      plot.caption.position = "plot",
      plot.title = element_markdown(),
      plot.subtitle = element_markdown(face = "italic", size = rel(0.80)),
      plot.caption = element_markdown(face = "italic", size = rel(0.65)),
      axis.title.x = element_markdown(),
      axis.text.x = element_markdown(),
      axis.title.y = element_markdown(),
      axis.text.y = element_markdown(),
      panel.grid.minor = element_blank()
    ) +
    theme(plot.caption = element_markdown())

  agg_png(
    filename = file.path(output_directory, "pca", "gene_pca_asso.png"),
    width = 16, height = 12, units = "cm", res = 300, scaling = 1
  )
    print(pca_asso)
  invisible(dev.off())

  for (ivar in keep_technical) {
    p <- wrap_plots(
      c(
        plot_planes(pca_dfxy, ivar, fig_n_comp),
        list(p_inertia)
      ),
      guides = "collect"
    ) +
      plot_annotation(
        title = paste0(
          "Structure Detection For: '<i>",
          names(keep_technical[match(ivar, keep_technical)]),
          "</i>'"
        ),
        tag_levels = "A",
        theme = theme_minimal(base_size = 10, base_family = "Tex Gyre Termes") +
        theme(
          plot.title.position = "plot",
          plot.caption.position = "plot",
          plot.title = element_markdown(),
          plot.subtitle = element_markdown(face = "italic", size = rel(0.80)),
          plot.caption = element_markdown(face = "italic", size = rel(0.65)),
          axis.title.x = element_markdown(),
          axis.text.x = element_markdown(),
          axis.title.y = element_markdown(),
          axis.text.y = element_markdown(),
          panel.grid.minor = element_blank()
        )
      )

    agg_png(
      filename = sprintf("%s/pca/gene_pca_planes_%s.png",
        output_directory,
        tolower(gsub("[.]+$|^[.]+|^X+", "", gsub("[.]+", ".", make.names(ivar))))
      ),
      width = 16, height = 12, units = "cm", res = 300, scaling = 1
    )
      print(p)
    invisible(dev.off())
  }

  ## DE --------------------------------------------------------------------------------------------
  ensembl_dt <- setDT(getBM(
    attributes = c(
      ensembl_id,
      "entrezgene_id",
      "uniprotswissprot",
      "chromosome_name",
      "start_position",
      "end_position",
      "external_gene_name"
    ),
    filters = ensembl_id,
    values = list(rownames(txi_counts[["counts"]])),
    mart = mart
  ))[
    j = lapply(.SD, function(x) {
      out <- paste(setdiff(unique(x), ""), collapse = ";")
      fifelse(out == "", NA_character_, out)
    }),
    by = c(ensembl_id)
  ][j = ensembl_version := ensembl_version]


  entrez_dt <- setDT(getBM(
    attributes = c(
      ensembl_id,
      "entrezgene_id"
    ),
    filters = ensembl_id,
    values = list(rownames(txi_counts[["counts"]])),
    mart = mart
  ))[
    j = lapply(.SD, function(x) {
      out <- paste(setdiff(unique(x), ""), collapse = ";")
      fifelse(out == "", NA_character_, out)
    }),
    by = c(ensembl_id)
  ]

  uniprot_dt <- setDT(getBM(
    attributes = c(
      ensembl_id,
      "uniprotswissprot"
    ),
    filters = ensembl_id,
    values = list(rownames(txi_counts[["counts"]])),
    mart = mart
  ))[
    j = lapply(.SD, function(x) {
      out <- paste(setdiff(unique(x), ""), collapse = ";")
      fifelse(out == "", NA_character_, out)
    }),
    by = c(ensembl_id)
  ]

  annot_dt <- merge(
    x = ensembl_dt,
    y = merge(x = entrez_dt, y = uniprot_dt, by = ensembl_id, all = TRUE),
    by = ensembl_id,
    all.x = TRUE
  )

  if (length(setdiff(rownames(txi_counts[["counts"]]), annot_dt[[ensembl_id]])) != 0) {
    stop(sprintf('Not all "%s" have been found! Please check Ensembl version!', ensembl_id))
  }

  for (ntrait in names(traits)) {
    local({
      trait <- traits[[ntrait]]
      base_form <- as.formula(paste0("~ ", trait))

      dds_pheno <- na.exclude(sample_sheet_qc[j = .SD, .SDcols = c("Sample_ID", all.vars(base_form))])
      dds_counts <- lapply(
        X = txi_counts,
        .sample = as.character(dds_pheno[["Sample_ID"]]),
        FUN = function(.l, .samples) if (is.matrix(.l)) .l[, .samples] else .l
      )

      dds <- DESeqDataSetFromTximport(txi = dds_counts, colData = dds_pheno, design = base_form)
      dds <- dds[rowVars(counts(dds)) != 0 & rowMeans(counts(dds)) > 1 & rowMedians(counts(dds)) > 0, ]
      stats_dds <- replaceOutliers(
        object = nbinomWaldTest(estimateDispersions(estimateSizeFactors(dds)), maxit = 1000),
        minReplicates = ncol(dds)
      )

      is_issue <- data.table(
        x = rownames(dds),
        converge = mcols(stats_dds)$betaConv
      )
      setnames(x = is_issue, old = "x", new = ensembl_id)

      if (is.factor(sample_sheet_qc[[trait]])) {
        results_dt <- rbindlist(apply(
          X = combn(levels(sample_sheet_qc[[trait]]), 2),
          MARGIN = 2,
          .trait = trait,
          .dds = stats_dds,
          FUN = function(.contrast, .trait, .dds) {
            results_dds <- results(
              object = .dds,
              contrast = c(.trait, .contrast[2], .contrast[1]),
              pAdjustMethod = "BH",
              independentFiltering = FALSE,
              cooksCutoff = FALSE
            )
            results_dt <- as.data.table(results_dds, keep.rownames = ensembl_id)
            setnames(results_dt, old = "padj", new = "fdr")
            results_dt[j = contrast := paste0(.trait, ": ", .contrast[2], " Vs. ", .contrast[1], " (ref)")]
            results_dt
          }
        ))
      } else {
        results_dds <- results(
          object = stats_dds,
          name = trait,
          pAdjustMethod = "BH",
          independentFiltering = FALSE,
          cooksCutoff = FALSE
        )
        results_dt <- as.data.table(results_dds, keep.rownames = ensembl_id)
        setnames(results_dt, old = "padj", new = "fdr")
        results_dt[j = contrast := trait]
      }

      results_dt[j = Trait := ntrait]

      results_annot_dt <- merge(
        x = merge(x = results_dt, y = is_issue, by = ensembl_id),
        y = annot_dt,
        by = ensembl_id,
        all.x = TRUE
      )

      fwrite(
        x = results_annot_dt,
        file = file.path(
          output_directory,
          paste0(
            paste(
              project_name,
              "DE",
              trait,
              rna_level_name,
              sep = "_"
            ),
            ".csv.gz"
          )
        )
      )

      lvp <- draw_volcano(
        data = results_annot_dt, x = "log2FoldChange", y = "pvalue",
        label_x = "Fold-Change (log<sub>2</sub>)",
        label_y = "P-value",
        alpha = 0.05,
        max_p = 1
      )

      p_volcano <- wrap_plots(lvp) +
        labs(
          title = paste("Volcano Plot:", unique(results_annot_dt[["Trait"]])),
          subtitle = sprintf("Model: %s", paste(labels(terms(base_form)), collapse = " + ")),
          caption = paste(
            sapply(
              X = c("pvalue", "fdr"),
              lp = lvp,
              alpha = 0.05,
              FUN = function(lp, p, alpha) {
                out <- sprintf("<b>%s</b>) %s",
                  LETTERS[seq_along(lp)],
                  sapply(lp, function(.gg) {
                    format(sum(.gg$data[[p]] < alpha, na.rm = TRUE), big.mark = ",", digits = 1L, nsmall = 0L, trim = TRUE)
                  })
                )
                sprintf(
                  fmt = "%s < %s: %s and %s.",
                  c("pvalue" = "P-value", "fdr" = "FDR")[p],
                  0.05,
                  paste(out[-length(out)], collapse = ", "),
                  out[length(out)]
                )
              }
            ),
            collapse = "<br>"
          )
        ) +
        theme(plot.caption = element_markdown(face = "italic", size = rel(0.75)))

      agg_png(
        filename = file.path(
          output_directory,
          paste0(
            paste(
              project_name,
              "DE",
              trait,
              rna_level_name,
              sep = "_"
            ),
            ".png"
          )
        ),
        width = 16, height = 12, units = "cm", res = 300, scaling = 1
      )
        print(p_volcano)
      invisible(dev.off())
    })
  }
}


### Archive ========================================================================================
# if (!Sys.getenv("USER") %in% c("root", "") && file.exists("~/.fex/id")) {
#   local({
#     owd <- getwd()
#     setwd(normalizePath(output_directory))
#     archive_name <- file.path(
#       normalizePath(output_directory),
#       paste0(
#         format(Sys.Date(), format = "%Y%m%d"), "_",
#         project_name, "_",
#         gsub("[0-9]+\\-", "", basename(output_directory)), ".zip"
#       )
#     )
#     zip(archive_name, files = list.files())
#     fex_out <- system(paste("fexsend", archive_name, "."), intern = TRUE)
#     unlink(archive_name)
#     setwd(owd)
#   })
# }


### Set chmod ======================================================================================
Sys.chmod(
  list.files(output_directory, full.names = TRUE),
  mode = "0775", use_umask = FALSE
)
Sys.chmod(
  list.files(output_directory, full.names = TRUE, recursive = TRUE, all.files = TRUE),
  mode = "0775", use_umask = FALSE
)
invisible(system(paste("chgrp -R staff", output_directory), intern = TRUE))


### Complete =======================================================================================
message("Success!", appendLF = TRUE)
message(timestamp(quiet = TRUE))
