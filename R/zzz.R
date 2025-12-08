# This loads the ggplot theme used for the figures.

.onLoad <- function(libname, pkgname) {

    .element_text <- ggplot2::element_text
    if (requireNamespace("ggtext", quietly = TRUE)) {
        .element_text <- ggtext::element_markdown
    }
    # ggplot theme:
    ggplot2::theme_set(
        ggplot2::theme_classic() +
        ggplot2::theme(strip.background = ggplot2::element_blank(),
                       strip.text = .element_text(size = 10),
                       axis.title = .element_text(color = "black", size = 11),
                       axis.title.x = .element_text(color = "black", size = 11),
                       axis.title.y = .element_text(color = "black", size = 11),
                       axis.text = .element_text(color = "black", size = 9),
                       axis.text.x = .element_text(color = "black", size = 9),
                       axis.text.y = .element_text(color = "black", size = 9),
                       legend.title = .element_text(color = "black", size = 11),
                       axis.ticks = ggplot2::element_line(color = "black"),
                       legend.background = ggplot2::element_blank(),
                       plot.title = .element_text(size = 14, hjust = 0.5)))
}
