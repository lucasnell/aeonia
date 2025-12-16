#' Save a plot to a PDF file using `cairo_pdf` or `svglite`
#'
#' @param filename Filename of plot. Extension should be `".pdf"`, `".svg"`, or `".png"`
#'     and will determine what type of plot is created.
#' @param plot Plot object or function to create plot
#' @param width Width in inches
#' @param height Height in inches
#' @param seed Integer to seed RNG for consistent jittering.
#' @param dpi Resolution to use if creating a png file.
#' @param fun_args List containing arguments to `plot` if it's a function and has
#'     required arguments.
#' @param bg Default background color for the plot. Defaults to `"transparent"`.
#' @param ... Other arguments to `cairo_pdf` or `svglite`.
#'
#' @importFrom svglite svglite
#'
#' @return `NULL`
#' @export
#'
save_plot <- function(filename,
                      plot,
                      width,
                      height,
                      seed = NULL,
                      dpi = 300,
                      fun_args = list(),
                      bg = "transparent",
                      ...) {
    stopifnot(is.list(fun_args))
    ext <- tail(strsplit(filename, "\\.")[[1]], 1)
    if (!ext %in% c("pdf", "svg", "png")) {
        stop("ERROR: file name extension must be \".pdf\", \".svg\", or \".png\"")
    }
    fn_dir <- dirname(filename)
    if (!dir.exists(fn_dir)) stop("ERROR: `", fn_dir, "` doesn't exist")
    if (!is.null(seed)) set.seed(seed)
    if (ext == "pdf") {
        cairo_pdf(filename = filename, width = width, height = height, bg = bg, ...)
        if (is.function(plot)) {
            do.call(plot, fun_args)
        } else {
            plot(plot)
        }
        dev.off()
    } else if (ext == "svg") {
        svglite(filename = filename, width = width, height = height, bg = bg, ...)
        if (is.function(plot)) {
            do.call(plot, fun_args)
        } else {
            plot(plot)
        }
        dev.off()
    } else {
        png(filename = filename, width = width, height = height, bg = bg,
            units = "in", res = dpi, ...)
        if (is.function(plot)) {
            do.call(plot, fun_args)
        } else {
            plot(plot)
        }
        dev.off()
    }
    invisible(NULL)
}
