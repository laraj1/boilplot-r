#' Plot single-cell gene expression distributions
#'
#' Plot expression distributions for one or more genes across categories
#' stored in a Seurat object. Only non-zero expression values are shown
#' in the violin and jitter plots. The proportion of cells with non-zero
#' expression is shown as a bar below zero.
#'
#' Supports Seurat v4 and v5 objects that are compatible with the
#' installed SeuratObject version.
#'
#' @param object A Seurat object.
#' @param genes A character vector of gene names. A single gene may be
#' supplied as a character string.
#' @param category A column name in the Seurat object's metadata used to
#' define groups.
#' @param assay Assay from which expression should be extracted. Defaults to
#' the object's default assay.
#' @param layer Expression layer to use. Defaults to "data". In Seurat v4,
#' this is mapped to the corresponding assay slot.
#' @param palette Optional named character vector mapping category names to
#' colours. If NULL, colours are generated automatically.
#' @param category_order Optional character vector specifying the order of
#' categories on the x-axis.
#' @param min_nonzero Minimum number of non-zero observations required for a
#' violin to be displayed for a gene/category combination.
#' @param title Optional plot title.
#' @param point_size Size of jitter points.
#' @param point_alpha Transparency of jitter points.
#' @param jitter_width Width of horizontal jitter.
#'
#' @return A ggplot object.
#'
#' @export
boil_plot <- function(
    object,
    genes,
    category,
    assay = NULL,
    layer = "data",
    palette = NULL,
    category_order = NULL,
    min_nonzero = 10,
    title = NULL,
    point_size = 0.7,
    point_alpha = 0.6,
    jitter_width = 0.1
) {
  assay <- assay %||% SeuratObject::DefaultAssay(object)
  genes <- validate_inputs(
    object = object,
    genes = genes,
    category = category,
    assay = assay,
    layer = layer,
    category_order = category_order,
    palette = palette,
    min_nonzero = min_nonzero
  )
  
  categories <- get_categories(
    object = object,
    category = category,
    category_order = category_order
  )
  
  resolved_palette <- get_palette(
    palette = palette,
    categories = categories
  )
  
  expression <- get_expression(
    object = object,
    genes = genes,
    assay = assay,
    layer = layer
  )
  
  plot_data <- prepare_plot_data(
    object = object,
    genes = genes,
    category = category,
    categories = categories,
    expression = expression,
    min_nonzero = min_nonzero
  )
  
  build_plot(
    plot_data = plot_data,
    category = category,
    palette = resolved_palette,
    title = title,
    point_size = point_size,
    point_alpha = point_alpha,
    jitter_width = jitter_width
  )
}

# General utilities
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
nonnegative_breaks <- function(limits) {
  if (is.null(limits) || length(limits) != 2) {
    return(NULL)
  }
  
  breaks <- scales::extended_breaks(n = 5)(limits)
  breaks[breaks >= 0]
}

# Input validation
validate_assay_compatibility <- function(
    object,
    assay
) {
  assay_object <- object[[assay]]
  seurat_object_version <- utils::packageVersion(
    "SeuratObject"
  )
  
  if (
    seurat_object_version < "5.0.0" &&
    inherits(assay_object, "Assay5")
  ) {
    stop(
      paste0(
        "The assay '",
        assay,
        "' is a Seurat v5 Assay5 object, but SeuratObject ",
        seurat_object_version,
        " is installed. ",
        "Please install SeuratObject >= 5.0.0 to use this object."
      ),
      call. = FALSE
    )
  }
}

validate_inputs <- function(
    object,
    genes,
    category,
    assay,
    layer,
    category_order,
    palette,
    min_nonzero
) {
  if (!inherits(object, "Seurat")) {
    stop(
      "object must be a Seurat object.",
      call. = FALSE
    )
  }
  
  if (!is.character(genes)) {
    stop(
      "genes must be a character vector.",
      call. = FALSE
    )
  }
  
  if (length(genes) == 0) {
    stop(
      "genes must contain at least one gene.",
      call. = FALSE
    )
  }
  
  genes <- unique(genes)
  
  if (!is.character(category) || length(category) != 1) {
    stop(
      "category must be a single metadata column name.",
      call. = FALSE
    )
  }
  
  metadata <- object[[]]
  
  if (!category %in% colnames(metadata)) {
    stop(
      sprintf(
        "'%s' was not found in the Seurat object's metadata.",
        category
      ),
      call. = FALSE
    )
  }
  
  if (!is.character(assay) || length(assay) != 1) {
    stop(
      "assay must be a single assay name.",
      call. = FALSE
    )
  }
  
  available_assays <- names(object)
  
  if (!assay %in% available_assays) {
    stop(
      sprintf(
        "Assay '%s' was not found. Available assays: %s",
        assay,
        paste(available_assays, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  validate_assay_compatibility(
    object = object,
    assay = assay
  )
  
  if (!is.character(layer) || length(layer) != 1) {
    stop(
      "layer must be a single layer name.",
      call. = FALSE
    )
  }
  
  if (
    !is.numeric(min_nonzero) ||
    length(min_nonzero) != 1 ||
    is.na(min_nonzero) ||
    min_nonzero < 1 ||
    min_nonzero != floor(min_nonzero)
  ) {
    stop(
      "min_nonzero must be a positive integer.",
      call. = FALSE
    )
  }
  
  if (!is.null(category_order)) {
    if (!is.character(category_order)) {
      stop(
        "category_order must be a character vector.",
        call. = FALSE
      )
    }
    
    if (length(category_order) == 0) {
      stop(
        "`category_order` cannot be empty.",
        call. = FALSE
      )
    }
    
    observed_categories <- unique(
      as.character(
        stats::na.omit(metadata[[category]])
      )
    )
    
    missing_categories <- setdiff(
      category_order,
      observed_categories
    )
    
    if (length(missing_categories) > 0) {
      stop(
        sprintf(
          paste0(
            "The following categories in `category_order` ",
            "were not found in `%s`: %s"
          ),
          category,
          paste(missing_categories, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    
    if (
      length(unique(category_order)) !=
      length(category_order)
    ) {
      stop(
        "`category_order` must not contain duplicate categories.",
        call. = FALSE
      )
    }
    
  }
  if (!is.null(palette)) {
    if (!is.character(palette)) {
      stop(
        "palette must be a named character vector or NULL.",
        call. = FALSE
      )
    }
    
    if (
      is.null(names(palette)) ||
      any(names(palette) == "")
    ) {
      stop(
        "`palette` must be a named character vector.",
        call. = FALSE
      )
    }
    
  }
  genes
}

# Category handling
get_categories <- function(
    object,
    category,
    category_order = NULL
) {
  values <- object[[category, drop = TRUE]]
  values_character <- as.character(values)
  
  if (!is.null(category_order)) {
    return(category_order)
  }
  
  if (
    is.factor(values) &&
    is.ordered(values)
  ) {
    observed <- unique(
      values_character[
        !is.na(values_character)
      ]
    )
    
    return(
      levels(values)[
        levels(values) %in% observed
      ]
    )
    
  }
  unique(
    values_character[
      !is.na(values_character)
    ]
  )
}

# Palette handling
get_palette <- function(
    palette,
    categories
) {
  if (!is.null(palette)) {
    missing <- setdiff(
      categories,
      names(palette)
    )
    if (length(missing) > 0) {
      stop(
        sprintf(
          paste0(
            "The supplied palette is missing colours for ",
            "the following categories: %s"
          ),
          paste(missing, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    
    return(palette[categories])
    
  }
  colours <- scales::hue_pal()(
    length(categories)
  )
  
  names(colours) <- categories
  
  colours
}

# Expression extraction
get_expression <- function(
    object,
    genes,
    assay,
    layer
) {
  assay_object <- object[[assay]]
  seurat_object_version <- utils::packageVersion(
    "SeuratObject"
  )
  
  if (seurat_object_version >= "5.0.0") {
    available_layers <- SeuratObject::Layers(
      object = assay_object
    )
    
    if (!layer %in% available_layers) {
      stop(
        sprintf(
          paste0(
            "Layer '%s' was not found in assay '%s'. ",
            "Available layers: %s"
          ),
          layer,
          assay,
          paste(available_layers, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    
    expression <- SeuratObject::LayerData(
      object = object,
      assay = assay,
      layer = layer
    )
    
  } else {
    available_slots <- c(
      "counts",
      "data",
      "scale.data"
    )
    if (!layer %in% available_slots) {
      stop(
        sprintf(
          "For SeuratObject < 5, `layer` must be one of: %s",
          paste(available_slots, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    
    expression <- SeuratObject::GetAssayData(
      object = object,
      assay = assay,
      slot = layer
    )
    
  }
  available_genes <- rownames(expression)
  
  missing_genes <- setdiff(
    genes,
    available_genes
  )
  
  if (length(missing_genes) > 0) {
    stop(
      sprintf(
        paste0(
          "The following genes were not found in assay '%s', ",
          "layer '%s': %s"
        ),
        assay,
        layer,
        paste(missing_genes, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  # Return cells x genes for downstream data preparation
  # Return cells x genes to match the downstream data preparation.
  expression <- as.matrix(
    expression[
      genes,
      ,
      drop = FALSE
    ]
  )
  
  t(expression)
}

# Plot data preparation
prepare_plot_data <- function(
    object,
    genes,
    category,
    categories,
    expression,
    min_nonzero
) {
  metadata <- object[[]]
  expression_df <- as.data.frame(
    as.matrix(expression),
    check.names = FALSE
  )
  
  expression_df$cell <- rownames(expression_df)
  
  expression_df$category <- as.character(
    metadata[
      expression_df$cell,
      category
    ]
  )
  
  expression_df <- expression_df[
    !is.na(expression_df$category) &
      expression_df$category %in% categories,
    ,
    drop = FALSE
  ]
  
  expression_df$category <- factor(
    expression_df$category,
    levels = categories,
    ordered = TRUE
  )
  
  expression_long <- tidyr::pivot_longer(
    expression_df,
    cols = dplyr::all_of(genes),
    names_to = "gene",
    values_to = "expression"
  )
  
  # Preserve the order in which genes were supplied
  expression_long$gene <- factor(
    expression_long$gene,
    levels = genes,
    ordered = TRUE
  )
  category_positions <- stats::setNames(
    seq_along(categories) - 1,
    categories
  )
  
  expression_long$x <- unname(
    category_positions[
      as.character(expression_long$category)
    ]
  )
  
  nonzero_counts <- (
    expression_long
    |> dplyr::group_by(gene, category)
    |> dplyr::summarise(
      nonzero_count = sum(
        expression > 0,
        na.rm = TRUE
      ),
      .groups = "drop"
    )
  )
  
  expression_long <- dplyr::left_join(
    expression_long,
    nonzero_counts,
    by = c("gene", "category")
  )
  
  detection_data <- (
    expression_long
    |> dplyr::group_by(gene, category)
    |> dplyr::summarise(
      nonzero_prop = mean(
        expression > 0,
        na.rm = TRUE
      ),
      .groups = "drop"
    )
  )
  
  detection_data$x <- unname(
    category_positions[
      as.character(detection_data$category)
    ]
  )
  
  max_expression <- (
    expression_long
    |> dplyr::group_by(gene)
    |> dplyr::summarise(
      max_expression = max(
        expression,
        na.rm = TRUE
      ),
      .groups = "drop"
    )
  )
  
  detection_data <- dplyr::left_join(
    detection_data,
    max_expression,
    by = "gene"
  )
  
  bar_width <- 0.65
  
  detection_data$bar_height <-
    detection_data$max_expression * 0.15
  
  detection_data$bar_ymin <-
    -detection_data$bar_height
  
  detection_data$bar_ymax <-
    -detection_data$max_expression * 0.05
  
  detection_data$bar_xmin <-
    detection_data$x - bar_width / 2
  
  detection_data$bar_xmax <-
    detection_data$x + bar_width / 2
  
  detection_data$detection_xmax <-
    detection_data$bar_xmin +
    bar_width * detection_data$nonzero_prop
  
  detection_data$label <- paste0(
    round(
      detection_data$nonzero_prop * 100
    ),
    "%"
  )
  
  detection_data$label_y <-
    -detection_data$bar_height * 1.5
  
  strip_data <- (
    expression_long
    |> dplyr::filter(expression > 0)
  )
  
  violin_data <- (
    expression_long
    |> dplyr::filter(
      expression > 0,
      nonzero_count >= min_nonzero
    )
  )
  
  list(
    expression_long = expression_long,
    strip_data = strip_data,
    violin_data = violin_data,
    detection_data = detection_data,
    category_positions = category_positions
  )
}

# Plot construction
build_plot <- function(
    plot_data,
    category,
    palette,
    title,
    point_size,
    point_alpha,
    jitter_width
) {
  expression_long <- plot_data$expression_long
  strip_data <- plot_data$strip_data
  violin_data <- plot_data$violin_data
  detection_data <- plot_data$detection_data
  category_positions <- plot_data$category_positions
  categories <- names(category_positions)
  
  p <- ggplot2::ggplot(
    expression_long,
    ggplot2::aes(
      x = x,
      y = expression,
      fill = category
    )
  ) +
    ggplot2::geom_violin(
      data = violin_data,
      trim = TRUE,
      scale = "width",
      width = 0.9
    ) +
    ggplot2::geom_jitter(
      data = strip_data,
      alpha = point_alpha,
      shape = 21,
      size = point_size,
      width = jitter_width
    ) +
    ggplot2::geom_rect(
      data = detection_data,
      mapping = ggplot2::aes(
        xmin = bar_xmin,
        xmax = bar_xmax,
        ymin = bar_ymin,
        ymax = bar_ymax
      ),
      fill = "#bdbdbd",
      inherit.aes = FALSE
    ) +
    ggplot2::geom_rect(
      data = detection_data,
      mapping = ggplot2::aes(
        xmin = bar_xmin,
        xmax = detection_xmax,
        ymin = bar_ymin,
        ymax = bar_ymax,
        fill = category
      ),
      inherit.aes = FALSE
    ) +
    ggplot2::geom_text(
      data = detection_data,
      mapping = ggplot2::aes(
        x = x,
        y = label_y,
        label = label
      ),
      size = 3,
      colour = "#333333",
      inherit.aes = FALSE
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      colour = "#444444",
      linewidth = 0.3
    ) +
    ggplot2::facet_wrap(
      ~gene,
      scales = "free_y",
      ncol = 1
    ) +
    ggplot2::labs(
      title = title,
      fill = category,
      x = category,
      y = "Expression"
    ) +
    ggplot2::scale_fill_manual(
      values = palette,
      limits = categories,
      guide = ggplot2::guide_legend(
        override.aes = list(
          colour = NA
        )
      )
    ) +
    ggplot2::scale_x_continuous(
      breaks = unname(category_positions),
      labels = categories,
      expand = ggplot2::expansion(
        mult = c(0.1, 0.1)
      )
    ) +
    ggplot2::scale_y_continuous(
      breaks = nonnegative_breaks
    ) +
    ggplot2::theme_classic()
  
  p
}