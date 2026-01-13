# For avoiding warnings for comparing nonsensible inputs
comparable <- function(x) {
    return(!is.null(x) && (is.vector(x) || is.array(x)) && !any(is.na(x)))
}

# Check for a single, whole number, perhaps in range
single_integer <- function(x, par, .min = NULL, .max = NULL) {
    all_good <- TRUE
    if (!comparable(x)) {
        all_good <- FALSE
    } else {
        all_good <- is.numeric(x) && length(x) == 1 && x %% 1 == 0
        if (all_good && !is.null(.min)) all_good <- x >= .min
        if (all_good && !is.null(.max)) all_good <- x <= .max
    }
    if (!all_good) {
        err_suff <- "a single number"
        if (!is.null(.min)) err_suff <- paste0(err_suff, ", min = ", .min)
        if (!is.null(.max)) err_suff <- paste0(err_suff, ", max = ", .max)
        err_msg(par, err_suff)
    }
    invisible(NULL)
}
# Check for a single number, perhaps in range
single_number <- function(x, par, .min = NULL, .max = NULL) {
    all_good <- TRUE
    if (!comparable(x)) {
        all_good <- FALSE
    } else {
        all_good <- is.numeric(x) && length(x) == 1
        if (all_good && !is.null(.min)) all_good <- x >= .min
        if (all_good && !is.null(.max)) all_good <- x <= .max
    }
    if (!all_good) {
        err_suff <- "a single number"
        if (!is.null(.min)) err_suff <- paste0(err_suff, ", min = ", .min)
        if (!is.null(.max)) err_suff <- paste0(err_suff, ", max = ", .max)
        err_msg(par, err_suff)
    }
    invisible(NULL)
}
# Check for an array of integers from 0 to 3 (used for landscapes)
# .dims is the exact dimensions, and should be length 3 for most cases,
# but will be length 2 in `big_plantscape`
is_landscape_array <- function(x, par, .dims) {
    all_good <- TRUE
    if (!comparable(x)) {
        all_good <- FALSE
    } else {
        all_good <- is.numeric(x) && all(x >= 0) && all(x <= 3) && all(x %% 1 == 0)
        if (all_good) all_good <- is.array(x)
        if (!is.null(.dims) && all_good) all_good <- length(dim(x)) == length(.dims)
        if (!is.null(.dims) && all_good) all_good <- all(dim(x) == .dims)
    }
    if (!all_good) {
        err_msg(par, sprintf("a %s array containing integers from 0 to 3",
                             paste(.dims, collapse = " x ")))
    }
    invisible(NULL)
}
is_type <- function(x, par, type, L = NULL, .min = NULL, .max = NULL) {
    all_good <- TRUE
    if (!comparable(x)) {
        all_good <- FALSE
    } else {
        if (is.character(type)) all_good <- inherits(x, type)
        if (is.function(type)) all_good <- type(x)
        if (all_good && !is.null(L)) all_good <- length(x) %in% L
        if (all_good && !is.null(.min)) all_good <- all(x >= .min)
        if (all_good && !is.null(.max)) all_good <- all(x <= .max)
    }
    if (!all_good) {
        err_suff <- paste("an object of type", type)
        if (!is.null(L)) err_suff <- paste(err_suff, "with possible length(s) =",
                                           paste(L, collapse = ", "))
        if (!is.null(.min)) err_suff <- paste0(err_suff, "; min value(s) = ", .min)
        if (!is.null(.max)) err_suff <- paste0(err_suff, "; max value(s) = ", .max)
        err_msg(par, err_suff)
    }
    invisible(NULL)
}

# Standard way to show error messages (also to make input-checking less verbose).
#
#
err_msg <- function(par, ...) {
    stop(sprintf("\nArgument `%s` must be %s.",
                 par, paste(...)), call. = FALSE)
}







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
