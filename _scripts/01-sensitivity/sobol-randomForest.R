
#'
#' Small-scale sensitivity via Sobol indices and randomForest
#'
#'

library(randomForest)
library(iml)


source("_scripts/00-preamble.R")
source("_scripts/01-sensitivity/00-sobol-preamble.R")


# Output file names:
out_fn <- list(rf = c(ob = "-",
                      ob0 = "-0",
                      dob = "-diff_"),
               part = c(ob = "Parts-",
                         dob = "Parts-diff_"),
               inter = c(ob = "Int-",
                          dob = "Int-diff_"),
               imp = c(ob = "Imp-",
                        dob = "Imp-diff_")) |>
    map(\(x) {
        z <- sprintf("_scripts/interm-data/randomforest%soutbreak_size.rds", x)
        names(z) <- names(x)
        return(z)
    })



# =============================================================================*
# Generating datasets ----
# =============================================================================*

# Summarize outputs from each set of simulations:
sobol_summs <- read_rds("_scripts/interm-data/sobol-sims-summs.rds") |>
    mutate(spat_config = factor(spat_config))

# Same thing but looking at differences between with and without Pseudo:
diff_sobol_summs <- sobol_summs |>
    group_by(combo, alate_dens, across(all_of(names(vary_pars)))) |>
    summarize(across(all_of(c(yvars, "p_outbreak", "outbreak_size2",
                              "sd_outbreak_size")),
                     \(x) x[n_pseudo > 0] - x[n_pseudo == 0]),
              .groups = "drop") |>
    mutate(across(starts_with("outbreak_size"), \(x) round(x, 2)))

# Dataset for outbreak size (n_pseudo > 0):
set.seed(99240974) # used for shuffling dataset
ob_sim_df <- sobol_summs |>
    filter(alate_dens == 1, n_pseudo > 0) |>
    select(outbreak_size, all_of(names(vary_pars))) |>
    slice_sample(prop = 1)
# Dataset for outbreak size (n_pseudo == 0):
set.seed(558448128) # used for shuffling dataset
ob_sim_df0 <- sobol_summs |>
    filter(alate_dens == 1, n_pseudo == 0) |>
    select(outbreak_size, all_of(names(vary_pars))) |>
    slice_sample(prop = 1)

# Dataset for effect of pseudo on outbreak size:
set.seed(1976042860) # used for shuffling dataset
dob_sim_df <- diff_sobol_summs |>
    filter(alate_dens == 1) |>
    select(outbreak_size, all_of(names(vary_pars))) |>
    slice_sample(prop = 1)

n_train <- as.integer(round(nrow(ob_sim_df) * 0.75))
n_test <- nrow(ob_sim_df) - n_train

ob_train_df <- ob_sim_df[1:n_train,]
ob_test_df <- ob_sim_df[(n_train+1):nrow(ob_sim_df),]
ob_train_df0 <- ob_sim_df0[1:n_train,]
ob_test_df0 <- ob_sim_df0[(n_train+1):nrow(ob_sim_df),]
dob_train_df <- dob_sim_df[1:n_train,]
dob_test_df <- dob_sim_df[(n_train+1):nrow(dob_sim_df),]





# # Which value of `mtry` results in most explained variance?
# # Each step takes ~3 min
# set.seed(45181283)
# ob_m_df <- tibble(m = 2:length(vary_pars)) |>
#     mutate(v = map_dbl(m, \(.m) {
#         rf <- randomForest(ob_sim_df[1:n_train,names(vary_pars)],
#                            ob_sim_df[["outbreak_size"]][1:n_train],
#                            corr.bias = TRUE, mtry = .m)
#         return(tail(rf$rsq, 1))
#     }))
# set.seed(1547476225)
# dob_m_df <- tibble(m = 2:length(vary_pars)) |>
#     mutate(v = map_dbl(m, \(.m) {
#         rf <- randomForest(dob_sim_df[1:n_train,names(vary_pars)],
#                            dob_sim_df[["outbreak_size"]][1:n_train],
#                            corr.bias = TRUE, mtry = .m)
#         return(tail(rf$rsq, 1))
#     }))
# ob_m_df |>
#     ggplot(aes(m, v)) +
#     geom_line() +
#     geom_point() +
#     scale_x_continuous(breaks = 2:10) +
#     theme(panel.grid.major = element_line(color = "gray80"))
# dob_m_df |>
#     ggplot(aes(m, v)) +
#     geom_line() +
#     geom_point() +
#     scale_x_continuous(breaks = 2:10) +
#     theme(panel.grid.major = element_line(color = "gray80"))
# # The best for both was mtry = 5


# =============================================================================*
# Fit randomForests ----
# =============================================================================*

if (!all(file.exists(out_fn$rf))) {
    set.seed(359013909)
    ob_rf <- randomForest(ob_train_df[,names(vary_pars)],
                          ob_train_df[["outbreak_size"]],
                          importance = TRUE, corr.bias = TRUE, mtry = 5)
    set.seed(794345405)
    ob_rf0 <- randomForest(ob_train_df0[,names(vary_pars)],
                          ob_train_df0[["outbreak_size"]],
                          importance = TRUE, corr.bias = TRUE, mtry = 5)
    set.seed(84490683)
    dob_rf <- randomForest(dob_train_df[,names(vary_pars)],
                           dob_train_df[["outbreak_size"]],
                           importance = TRUE, corr.bias = TRUE, mtry = 5)
    write_rds(ob_rf, out_fn$rf[["ob"]], compress = "gz")
    write_rds(ob_rf0, out_fn$rf[["ob0"]], compress = "gz")
    write_rds(dob_rf, out_fn$rf[["dob"]], compress = "gz")
} else {
    ob_rf <- read_rds(out_fn$rf[["ob"]])
    ob_rf0 <- read_rds(out_fn$rf[["ob0"]])
    dob_rf <- read_rds(out_fn$rf[["dob"]])
}



# variances explained (partial r^2):
(ob_rf_rsq <- tail(ob_rf$rsq, 1))
# [1] 0.9630947
(ob_rf_rsq0 <- tail(ob_rf0$rsq, 1))
# [1] 0.9841875
(dob_rf_rsq <- tail(dob_rf$rsq, 1))
# [1] 0.8367532


#' #' Does predicting from separate randomForests for with and without Pseudomonas,
#' #' then taking difference improve compared to a single randomForest on
#' #' the difference?
#' #' It does, but it's much more annoying to analyze!
#'
#' pred <- predict(ob_rf, dob_test_df[,names(vary_pars)]) -
#'     predict(ob_rf0, dob_test_df[,names(vary_pars)])
#' pred0 <- predict(dob_rf, dob_test_df[,names(vary_pars)])
#'
#' cor(pred, dob_test_df[["outbreak_size"]])^2
#' # [1] 0.8809224
#' cor(pred0, dob_test_df[["outbreak_size"]])^2
#' # [1] 0.812233
#'
#' # Takes ~
#' t0 <- Sys.time()
#' set.seed(2106841574)
#' ob_mse <- lapply(1:100, \(i) {
#'     xm0 <- dob_train_df[,names(vary_pars)]
#'     y <- dob_train_df[["outbreak_size"]]
#'     mse0 <- Metrics::mse(y, predict(ob_rf, xm0) - predict(ob_rf0, xm0))
#'     out <- rep(list(0), length(names(vary_pars))) |>
#'         set_names(names(vary_pars)) |>
#'         as_tibble()
#'     for (n in names(vary_pars)) {
#'         xm <- xm0
#'         xm[[n]] <- sample(xm[[n]])
#'         predicted2 <- predict(ob_rf, xm) - predict(ob_rf0, xm)
#'         mse2 <- Metrics::mse(y, predicted2)
#'         out[[n]] <- mse2 / mse0
#'     }
#'     return(out)
#' }) |>
#'     list_rbind()
#' t1 <- Sys.time()
#' t1 - t0
#'
#' apply(as.matrix(ob_mse), 2, median)
#' dob_rf_imp$results



# =============================================================================*
# Calc. part. dependence ----
# =============================================================================*

# Calculate data for partial dependence plots (takes a while to run, so this
# should be done separately from plots)
partial_calc <- function(x, pred_data) {

    # x = ob_rf; pred_data = ob_train_df
    # rm(x, pred_data, n)

    suppressPackageStartupMessages(library(future.apply))
    with(plan(multisession, workers = options()[["mc.cores"]], gc = TRUE), local = TRUE)

    n <- nrow(pred_data)

    out <- x$forest$xlevels |>
        names() |>
        future_lapply(\(x_name) {
            # x_name = "spat_config"
            # rm(x_name, xv, n_pt, x_pt, y_pt, x_data, i)
            xv <- pred_data[[x_name]]
            n_pt <- min(length(unique(xv)), 51)
            if (is.numeric(xv)) {
                x_pt <- seq(min(xv), max(xv), length = n_pt)
            } else if (is.factor(xv)) {
                x_pt <- sort(unique(xv))
            } else stop("Only programmed for numeric and factors")
            y_pt <- numeric(n_pt)
            x_data <- pred_data
            for (i in 1:n_pt) {
                x_data[, x_name] <- rep(x_pt[i], n)
                y_pt[i] <- mean(predict(x, x_data), na.rm = TRUE)
            }
            out <- tibble(param = x_name, x = as.numeric(x_pt), y = y_pt)
            # Because spat_config goes from 0 to 4, not 1 to 5:
            if (x_name == "spat_config") out$x <- out$x - 1
            return(out)
        }, future.packages = "randomForest") |>
        list_rbind()

    return(out)
}


if (!all(file.exists(out_fn$part))) {
    # Each step takes ~ 1.3 min
    ob_rf_parts <- partial_calc(ob_rf, ob_train_df)
    dob_rf_parts <- partial_calc(dob_rf, dob_train_df)
    write_rds(ob_rf_parts, out_fn$part[["ob"]], compress = "gz")
    write_rds(dob_rf_parts, out_fn$part[["dob"]], compress = "gz")
} else {
    ob_rf_parts <- read_rds(out_fn$part[["ob"]])
    dob_rf_parts <- read_rds(out_fn$part[["dob"]])
}




# =============================================================================*
# Fit interactions and predictors ----
# =============================================================================*

#' For interactions:
#' The interaction is measured by Friedman's H-statistic (square root of the
#' H-squared test statistic) and takes on values between 0 (no interaction)
#' to 1 (100% of standard deviation of f(x) du to interaction).


#' Create `iml` Predictor objects.
ob_rf_pred <- Predictor$new(ob_rf, data = ob_train_df[,names(vary_pars)],
                            y = ob_train_df[["outbreak_size"]])
dob_rf_pred <- Predictor$new(dob_rf, data = dob_train_df[,names(vary_pars)],
                             y = dob_train_df[["outbreak_size"]])
if (!all(file.exists(out_fn$inter))) {

    #' Create `iml` Interaction objects.
    #' This function creates an Interaction option and temporarily allows
    #' the future package to access large objects:
    make_inter <- \(rf_obj) {
        suppressPackageStartupMessages(library(future.apply))
        with(plan(multisession, workers = options()[["mc.cores"]], gc = TRUE),
             local = TRUE)
        oopts <- options(future.globals.maxSize = 1.5e9)  ## 1.5 GB
        on.exit(options(oopts))
        Interaction$new(rf_obj)
    }
    # Takes ~1.5 min each
    ob_rf_int <- make_inter(ob_rf_pred)
    dob_rf_int <- make_inter(dob_rf_pred)
    write_rds(ob_rf_int, out_fn$inter[["ob"]], compress = "gz")
    write_rds(dob_rf_int, out_fn$inter[["dob"]], compress = "gz")
} else {
    ob_rf_int <- read_rds(out_fn$inter[["ob"]])
    dob_rf_int <- read_rds(out_fn$inter[["dob"]])
}



# =============================================================================*
# Fit importance ----
# =============================================================================*

if (!all(file.exists(out_fn$imp))) {

    #' Create `iml` FeatureImp objects.
    #' This function creates an FeatureImp option and temporarily allows
    #' the future package to access large objects:
    make_imp <- \(rf_obj) {
        suppressPackageStartupMessages(library(future.apply))
        with(plan(multisession, workers = options()[["mc.cores"]], gc = TRUE),
             local = TRUE)
        oopts <- options(future.globals.maxSize = 1.5e9)  ## 1.5 GB
        on.exit(options(oopts))
        FeatureImp$new(rf_obj, loss = "mse", n.repetitions = 100)
    }
    # Take ~4 min each
    ob_rf_imp <- make_imp(ob_rf_pred)
    dob_rf_imp <- make_imp(dob_rf_pred)
    write_rds(ob_rf_imp, out_fn$imp[["ob"]], compress = "gz")
    write_rds(dob_rf_imp, out_fn$imp[["dob"]], compress = "gz")
} else {
    ob_rf_imp <- read_rds(out_fn$imp[["ob"]])
    dob_rf_imp <- read_rds(out_fn$imp[["dob"]])
}






# =============================================================================*
# Plot functions ----
# =============================================================================*


partial_plot <- function(x, .ylab) {
    # x = ob_rf_parts; .ylab = "Outbreak size"
    # rm(x, .ylab)
    if ("K" %in% x$param) x$x <- ifelse(x$param == "K", x$x / 1e3, x$x)
    x |>
        mutate(param = pretty_params(param) |>
                   factor(levels = pretty_params(names(vary_pars)))) |>
        ggplot(aes(x, y)) +
        geom_line() +
        geom_point() +
        facet_wrap(~ param, scales = "free_x", strip.position = "bottom") +
        labs(y = .ylab) +
        theme(plot.title = element_markdown(),
              axis.title.y = element_markdown(),
              strip.text = element_markdown(),
              strip.placement = "outside",
              axis.title.x = element_blank())
}




# Plot predicted vs observed:
pred_plot <- function(x, test_x, test_y, .title = "") {
    # x = ob_rf
    # test_x = ob_test_df[,names(vary_pars)]
    # test_y = ob_test_df[["outbreak_size"]]
    # .title = ""
    # rm(x, test_x, test_y, .title, pred_rf, rf_r2)
    stopifnot(inherits(x, "randomForest"))
    stopifnot(inherits(test_x, "data.frame"))
    stopifnot(inherits(test_y, "numeric") || inherits(test_y, "factor"))
    stopifnot(length(test_y) == nrow(test_y))
    stopifnot(inherits(.title, "character") && length(.title) == 1L)
    pred_rf <- predict(x, test_x)
    rf_r2 <- sprintf("r^2 == %.3f", cor(pred_rf, test_y)^2)
    tibble(Predicted = pred_rf, Observed = test_y) |>
        ggplot(aes(Observed, Predicted)) +
        geom_point(shape = 1) +
        ggtitle(.title) +
        geom_text(data = tibble(lab = rf_r2, x = min(pred_rf), y = max(pred_rf)),
                  aes(x, y, label = lab), hjust = 0, vjust = 1, parse = TRUE) +
        geom_abline(slope = 1, intercept = 0, linetype = 2, color = "red") +
        coord_equal() +
        theme(plot.title = element_markdown())
}

# Plot variable importance:
imp_plot <- function(x, .title = "") {
    # x = ob_rf; .title = ""; .sort = TRUE
    # rm(x, .title, .sort)
    stopifnot(inherits(x, c("randomForest", "FeatureImp")))
    stopifnot(inherits(.title, "character") && length(.title) == 1L)
    if (inherits(x, "randomForest")) {
        stopifnot("%IncMSE" %in% colnames(x$importance))
        importance(x) |>
            as.data.frame() |>
            rownames_to_column("par") |>
            select(par, `%IncMSE`) |>
            rename(inc_mse = `%IncMSE`) |>
            mutate(par = pretty_params(par),
                   par = fct_reorder(par, inc_mse, .na_rm = TRUE)) |>
            ggplot(aes(inc_mse, par)) +
            geom_vline(xintercept = 0, linewidth = 1, color = "gray70") +
            geom_point(size = 3) +
            geom_segment(aes(xend = 0, yend = par), linewidth = 1) +
            ggtitle(.title) +
            xlab("Mean increase in MSE") +
            theme(axis.title.y = element_blank(),
                  axis.text.y = element_markdown(size = 11, color = "black"),
                  plot.title = element_markdown())
    } else {
        x$results |>
            as_tibble() |>
            mutate(feature = pretty_params(feature),
                   feature = fct_reorder(feature, importance, .na_rm = TRUE)) |>
            ggplot(aes(importance, feature)) +
            geom_vline(xintercept = 1, linewidth = 1, color = "gray70") +
            geom_segment(aes(xend = 1, yend = feature), linewidth = 0.5, color = "gray60") +
            geom_pointrange(aes(xmin = `importance.05`, xmax = `importance.95`),
                            linewidth = 1.5) +
            ggtitle(.title) +
            xlab("Mean relative MSE") +
            theme(axis.title.y = element_blank(),
                  axis.text.y = element_markdown(size = 11, color = "black"),
                  plot.title = element_markdown())
    }
}


# Plot variable interaction strength:
inter_plot <- function(x, .title = "") {
    # x = ob_rf_int; .title = ""
    # rm(x, .title)
    stopifnot(inherits(x, "Interaction"))
    stopifnot(inherits(.title, "character") && length(.title) == 1L)

    p <- x$results |>
        as_tibble() |>
        filter(!is.na(.interaction)) |>
        mutate(.feature = pretty_params(.feature),
               .feature = fct_reorder(.feature, .interaction)) |>
        ggplot(aes(.interaction, .feature)) +
        geom_vline(xintercept = 0, linewidth = 1, color = "gray70") +
        geom_segment(aes(xend = 0, yend = .feature), linewidth = 0.5, color = "gray60") +
        geom_point(size = 3) +
        ggtitle(.title) +
        labs(x = "Interaction strength (Friedman's <i>H</i>)") +
        theme(axis.title.y = element_blank(),
              axis.title.x = element_markdown(),
              axis.text.y = element_markdown(size = 11, color = "black"),
              plot.title = element_markdown())
    if (".class" %in% colnames(x$results)) {
        p <- p + facet_wrap(~ .class)
    }
    return(p)
}


imp_inter_plot <- function(rf, rf_int, .sort = NA, .title = "") {

    # rf = ob_rf; rf_int = ob_rf_int; .title = ""; .sort = "inter"
    # rm(rf, rf_int, .sort, .title, .df, mse_max, int_max, trans, itrans, .mse_col, .int_col)
    stopifnot(inherits(rf, "randomForest"))
    stopifnot("%IncMSE" %in% colnames(rf$importance))
    stopifnot(inherits(rf_int, "Interaction"))
    stopifnot(inherits(.title, "character") && length(.title) == 1L)
    stopifnot(is.na(.sort) || (inherits(.sort, "character") && length(.sort) == 1L))
    stopifnot(is.na(.sort) || (.sort %in% c("inter", "imp")))

    .df <- importance(rf) |>
        as.data.frame() |>
        rownames_to_column("par") |>
        select(par, `%IncMSE`) |>
        rename(inc_mse = `%IncMSE`) |>
        left_join(rf_int[["results"]] |> as_tibble() |>
                      rename(par = .feature, inter = .interaction),
                  by = "par") |>
        mutate(par = pretty_params(par) |>
                   factor(levels = rev(pretty_params(names(vary_pars)))))
    if (!is.na(.sort)) {
        if (.sort == "imp") .sort <- "inc_mse"
        .df$par <- fct_reorder(.df$par, .df[[.sort]])
    }
    .df <- .df |>
        pivot_longer(-par, names_to = "measure")

    mse_max <- max(.df$value[.df$measure == "inc_mse"])
    int_max <- max(.df$value[.df$measure == "inter"])
    trans <- \(x) x * mse_max / int_max
    itrans <- \(x) x * int_max / mse_max

    .mse_col <- "dodgerblue"
    .int_col <- "goldenrod"

    .df |>
        mutate(value = ifelse(measure == "inter", trans(value), value)) |>
        ggplot(aes(value, par, fill = measure)) +
        geom_vline(xintercept = 0, linewidth = 1, color = "gray70") +
        geom_col(position = position_dodge(0.5), color = NA, width = 0.4) +
        ggtitle(.title) +
        scale_x_continuous("Mean increase in MSE",
                           sec.axis = sec_axis(itrans, "Overall interaction strength")) +
        scale_fill_manual(values = c(inc_mse = .mse_col, inter = .int_col), guide = "none") +
        theme(axis.title.y = element_blank(),
              axis.text.y = element_markdown(size = 11, color = "black"),
              axis.title.x.top = element_markdown(color = .int_col, face = "bold"),
              axis.text.x.top = element_markdown(color = .int_col),
              axis.ticks.x.top = element_line(color = .int_col),
              axis.title.x.bottom = element_markdown(color = .mse_col, face = "bold"),
              axis.text.x.bottom = element_markdown(color = .mse_col),
              axis.ticks.x.bottom = element_line(color = .mse_col),
              plot.title = element_markdown())

}






# =============================================================================*
# Plots ----
# =============================================================================*



# (partial_plot(ob_rf_parts, "Outbreak size") |
#         partial_plot(dob_rf_parts, "Effect of *Pseudomonas* on outbreak size")) &
#     theme(strip.text = element_markdown(size = 8.5))






# rf_pred_p <- (pred_plot(ob_rf, ob_test_df[,names(vary_pars)],
#            ob_test_df[["outbreak_size"]],
#            "Outbreak size") |
#         pred_plot(dob_rf, dob_test_df[,names(vary_pars)],
#                   dob_test_df[["outbreak_size"]],
#                   "Effect of *Pseudomonas* on outbreak size")) +
#     plot_layout(nrow = 1) +
#     plot_annotation(title = "Observed vs predicted",
#                     theme = theme(plot.title = element_markdown(size = 18)))
# rf_pred_p

# save_plot("_plots/rf-pred.pdf", rf_pred_p, width = 8, height = 5)



# (imp_plot(ob_rf, "Outbreak size") | imp_plot(dob_rf, "Effect of *Pseudomonas* on outbreak size")) +
#     plot_layout(nrow = 1) +
#     plot_annotation(title = "Variable importance",
#                     theme = theme(plot.title = element_markdown(size = 18)))
#
# (inter_plot(ob_rf_int, "Outbreak size") | inter_plot(dob_rf_int, "Effect of *Pseudomonas* on outbreak size")) +
#     plot_layout(nrow = 1) +
#     plot_annotation(title = "Variable interaction strength",
#                     theme = theme(plot.title = element_markdown(size = 18)))


# rf_imp_inter_p <- (imp_inter_plot(ob_rf, ob_rf_int, "imp", .title = "Outbreak size") |
#         imp_inter_plot(dob_rf, dob_rf_int, "imp", .title = "Effect of *Pseudomonas* on<br>outbreak size")) +
#     plot_layout(nrow = 1) +
#     plot_annotation(title = "Variable importance and interaction strength",
#                     theme = theme(plot.title = element_markdown(size = 18)))
# rf_imp_inter_p
# save_plot("_plots/rf-imp_inter.pdf", rf_imp_inter_p, width = 10, height = 5)



# (imp_plot(dob_rf) | imp_plot(dob_rf_imp)) +
#     plot_annotation(title = "Effect of *Pseudomonas* on outbreak size",
#                     theme = theme(plot.title = element_markdown()))


# (imp_plot(dob_rf_imp) | inter_plot(dob_rf_int)) +
#     plot_annotation(title = "Effect of *Pseudomonas* on outbreak size",
#                     theme = theme(plot.title = element_markdown()))





# =============================================================================*
# Predicting ----
# =============================================================================*





#' Create a dataframe containing some of the predictor variables that vary
#' (`par_x`, `par_y`), then add predictions to it (where all variables not
#' in `c(par_x, par_y)` are at their means):
make_preds <- function(par_x, par_y, n_lvls = 101) {
    # par_x = "Y0"; par_y = "mean_N"; n_lvls = 101
    # rm(par_x, par_y, n_lvls, d, nd, n, z)
    stopifnot(length(par_x) == 1L && par_x %in% names(vary_pars))
    stopifnot(length(par_y) == 1L && par_y %in% names(vary_pars))

    d <- c(par_x, par_y) |>
        set_names() |>
        map(\(x) {
            if (x == "spat_config") factor(0:4, levels = 0:4)
            else seq(vary_pars[[x]][1], vary_pars[[x]][2], length.out = n_lvls)
        }) |>
        do.call(what = crossing)
    nd <- d
    for (n in names(vary_pars)[!names(vary_pars) %in% c(par_x, par_y)]) {
        if (n == "spat_config") nd[[n]] <- factor(1L, levels = 0:4)
        else nd[[n]] <- mean(vary_pars[[n]])
    }
    z <- predict(dob_rf, nd)
    d$pred <- z
    return(d)
}


vary_pars |> as_tibble()

#' Variables of interest:
#'
#' - `Y0`
#' - `mean_N`
#' - `K`
#' - `spat_config`
#' - `zeta`
#'

make_preds("Y0", "mean_N") |>
    ggplot(aes(Y0, mean_N)) +
    geom_raster(aes(fill = pred)) +
    geom_contour(aes(z = pred), color = "gray70") +
    scale_fill_viridis_c() +
    labs(x = pretty_params("Y0"), y = pretty_params("mean_N")) +
    theme(axis.title.x = element_markdown(),
          axis.title.y = element_markdown())

make_preds("K", "zeta") |>
    ggplot(aes(K, zeta)) +
    geom_raster(aes(fill = pred)) +
    geom_contour(aes(z = pred), color = "gray70") +
    scale_fill_viridis_c() +
    labs(x = pretty_params("K"), y = pretty_params("zeta")) +
    theme(axis.title.x = element_markdown(),
          axis.title.y = element_markdown())


predict(dob_rf, tibble(Y0 = 9, mean_N = 25, sd_N = 25, K = 12500,
                       virus_attract = 1.5, pseudo_repel = 1.5,
                       pseudo_surv = 0.9, zeta = 0.1,
                       spat_config = factor(1L, levels = 0:4)))



# sim_args <- list(alate_dens = 1, Y0 = 9, mean_N = 25,
#                  sd_N = 25, K = 12500,
#                  virus_attract = 1.5, pseudo_repel = 1.5,
#                  pseudo_surv = 0.9, zeta = 0.1,
#                  spat_config = 1L, n_sims = 1000)
#
# sim_df <- do.call(one_combo, c(sim_args, list(n_pseudo = 3L))) |>
#     select(rep:outbreak_size)
# sim0_df <- do.call(one_combo, c(sim_args, list(n_pseudo = 0L))) |>
#     select(rep:outbreak_size)
#
# sim_df$outbreak_size |> mean()
# sim0_df$outbreak_size |> mean()

# Set parameter names in order for objective function:
par_names <- c("K", "virus_attract", "pseudo_repel", "pseudo_surv", "zeta",
               "Y0", "mean_N", "sd_N")
# Set bounds:
lower_bounds <- map_dbl(head(vary_pars, -1), \(x) x[1])[par_names]
upper_bounds <- map_dbl(head(vary_pars, -1), \(x) x[2])[par_names]
# Force parameters inside bounds:
trans_pars <- function(p) {
    if (length(p) == 8L) {
        x <- exp(-exp(-p)) * (upper_bounds - lower_bounds) + lower_bounds
    } else if (!is.null(names(p))) {
        nms <- names(p)
        stopifnot(all(nms %in% names(lower_bounds)))
        x <- exp(-exp(-p)) * (upper_bounds[nms] - lower_bounds[nms]) +
            lower_bounds[nms]
    } else stop("length(p) == 8L || !is.null(names(p))")
    return(x)
}
# Convert to version used in objective function:
inv_trans_pars <- function(x) {
    if (length(x) == 8L) {
        p <- -log(-log((x - lower_bounds) / (upper_bounds - lower_bounds)))
    } else if (!is.null(names(x))) {
        nms <- names(x)
        stopifnot(all(nms %in% names(lower_bounds)))
        p <- -log(-log((x - lower_bounds[nms]) / (upper_bounds[nms] -
                                                      lower_bounds[nms])))
    } else stop("length(x) == 8L || !is.null(names(x))")
    return(p)
}




dob_pred_f <- function(pars, .Y0 = 5, .mean_N = 55, .sd_N = 25) {

    stopifnot(length(pars) %in% c(5L, 8L))

    if (length(pars) == 5L) {
        pars <- c(pars, inv_trans_pars(c(Y0 = .Y0, mean_N = .mean_N, sd_N = .sd_N)))
    }

    pars <- trans_pars(pars)

    -1 * predict(dob_rf, tibble(K = pars[1],
                           virus_attract = pars[2],
                           pseudo_repel = pars[3],
                           pseudo_surv = pars[4],
                           zeta = pars[5],
                           Y0 = pars[6],
                           mean_N = pars[7],
                           sd_N = pars[8],
                           spat_config = factor(1L, levels = 0:4)))
}



x0 <- c(K = 24999, virus_attract = 3, pseudo_repel = 3,
        pseudo_surv = 0.99, zeta = 0.9, Y0 = 8.9, mean_N = 25, sd_N = 25)

op <- nloptr::neldermead(inv_trans_pars(x0), dob_pred_f,
                         lower = rep(-10, length(x0)),
                         upper = rep(100, length(x0)))

op$value
tibble(par = names(x0), value = trans_pars(op$par))

op_sim_args <- op$par |>
    set_names(names(x0)) |>
    trans_pars() |>
    as.list() |>
    c(list(spat_config = 1, alate_dens = 1, n_sims = 1000))

op_sim_df <- do.call(one_combo, c(op_sim_args, list(n_pseudo = 3L))) |>
    select(rep:outbreak_size)
op_sim0_df <- do.call(one_combo, c(op_sim_args, list(n_pseudo = 0L))) |>
    select(rep:outbreak_size)

op_sim_df$outbreak_size |> mean()
op_sim0_df$outbreak_size |> mean()


library(estimatePMR)
library(mirai)

if (!file.exists("_scripts/interm-data/randomforest-optim.rds")) {
    # Takes ~ 10.5 min
    wop <- winnowing_optim(dob_pred_f,
                           rep(-10, 8) |> set_names(par_names),
                           rep(100, 8) |> set_names(par_names),
                           n_bevals = 10L, n_boxes = 100L,
                           fn_args = list(), n_outputs = c(75L, 50L, 25L),
                           optimizers = rep(list(nloptr::bobyqa), 3),
                           multithread = FALSE)
    write_rds(wop, "_scripts/interm-data/randomforest-optim.rds")
} else {
    wop <- read_rds("_scripts/interm-data/randomforest-optim.rds")
}


# # Same thing, but fix starting conditions (Y0, mean_N, and sd_N):
# lower_bounds2 <- head(lower_bounds, -3)
# upper_bounds2 <- head(upper_bounds, -3)

if (!file.exists("_scripts/interm-data/randomforest-optim2.rds")) {
    # Takes ~6.5 min
    wop2 <- winnowing_optim(dob_pred_f, rep(-10, 5), rep(100, 5),
                            n_bevals = 10L, n_boxes = 100L,
                           fn_args = list(), n_outputs = c(75L, 50L, 25L),
                           optimizers = rep(list(nloptr::bobyqa), 3))
    write_rds(wop2, "_scripts/interm-data/randomforest-optim2.rds")
} else {
    wop2 <- read_rds("_scripts/interm-data/randomforest-optim2.rds")
}
if (!file.exists("_scripts/interm-data/randomforest-optim3.rds")) {
    # Takes ~5.5 min
    wop3 <- winnowing_optim(dob_pred_f, rep(-10, 5), rep(100, 5),
                            n_bevals = 10L, n_boxes = 100L,
                           fn_args = list(.Y0 = 9, .mean_N = 55, .sd_N = 25),
                           n_outputs = c(75L, 50L, 25L),
                           optimizers = rep(list(nloptr::bobyqa), 3))
    write_rds(wop3, "_scripts/interm-data/randomforest-optim3.rds")
} else {
    wop3 <- read_rds("_scripts/interm-data/randomforest-optim3.rds")
}
if (!file.exists("_scripts/interm-data/randomforest-optim4.rds")) {
    # Takes ~5.5 min
    wop4 <- winnowing_optim(dob_pred_f, rep(-10, 5), rep(100, 5),
                            n_bevals = 10L, n_boxes = 100L,
                           fn_args = list(.Y0 = 5, .mean_N = 100, .sd_N = 25),
                           n_outputs = c(75L, 50L, 25L),
                           optimizers = rep(list(nloptr::bobyqa), 3))
    write_rds(wop4, "_scripts/interm-data/randomforest-optim4.rds")
} else {
    wop4 <- read_rds("_scripts/interm-data/randomforest-optim4.rds")
}
if (!file.exists("_scripts/interm-data/randomforest-optim5.rds")) {
    # Takes ~5.5 min
    wop5 <- winnowing_optim(dob_pred_f, rep(-10, 5), rep(100, 5),
                            n_bevals = 10L, n_boxes = 100L,
                            fn_args = list(.Y0 = 9, .mean_N = 100, .sd_N = 25),
                            n_outputs = c(75L, 50L, 25L),
                            optimizers = rep(list(nloptr::bobyqa), 3))
    write_rds(wop5, "_scripts/interm-data/randomforest-optim5.rds")
} else {
    wop5 <- read_rds("_scripts/interm-data/randomforest-optim5.rds")
}


wop_plotter <- function(.wop, .title = NULL) {
    .wop |>
        imap(\(x, i) {
            xp <- x$par
            names(xp) <- par_names[1:length(xp)]
            tibble(rep = factor(i, levels = 1:length(.wop)),
                   par = par_names[1:length(xp)],
                   value = trans_pars(xp),
                   outbreak_size = -1 * x$value)
        }) |>
        list_rbind() |>
        mutate(par = factor(pretty_params(par),
                            levels = pretty_params(names(vary_pars))) |>
                   droplevels(),
               pseudo_effect = ifelse(outbreak_size > 0, "hurts", "helps") |>
                   factor(levels = c("helps", "hurts"))) |>
        ggplot(aes(par, value)) +
        geom_jitter(aes(size = outbreak_size, color = pseudo_effect),
                    shape = 1, stroke = 1) +
        facet_wrap(~ par, scales = "free") +
        scale_color_manual("Effect of<br>*Pseudomonas*<br>on plants",
                           values = c(hurts = "firebrick", helps = "gray60")) +
        scale_size("Effect of<br>*Pseudomonas*<br>on outbreak size") +
        labs(title = .title) +
        theme(strip.text = element_blank(),
              strip.background = element_blank(),
              legend.title = element_markdown(),
              plot.title = element_markdown(),
              axis.ticks.x = element_blank(),
              axis.text.x = element_markdown(color = "black", size = 10),
              axis.title.x = element_blank())
}



wop |> wop_plotter()

wop2 |> wop_plotter(sprintf("Y<sub>0</sub> = %s, &mu;<sub>N</sub> = %s", 5, 55))


wop3 |> wop_plotter(sprintf("Y<sub>0</sub> = %s, &mu;<sub>N</sub> = %s", 9, 55))

wop[[1]][["value"]]
wop2[[1]][["value"]]
wop3[[1]][["value"]]
wop4[[1]][["value"]]
wop5[[1]][["value"]]



wop_sim_args <- wop[[1]][["par"]] |>
    set_names(par_names) |>
    trans_pars() |>
    c(spat_config = 1, alate_dens = 1, n_sims = 1000)

wop_sim_df <- do.call(one_combo, c(wop_sim_args, list(n_pseudo = 3L))) |>
    select(rep:outbreak_size)
wop_sim0_df <- do.call(one_combo, c(wop_sim_args, list(n_pseudo = 0L))) |>
    select(rep:outbreak_size)

op_sim_df$outbreak_size |> mean()
op_sim0_df$outbreak_size |> mean()





wop_sim_args2 <- wop2[[1]][["par"]] |>
    set_names(par_names[1:5]) |>
    trans_pars() |>
    c(Y0 = 5, mean_N = 55, sd_N = 25,
      spat_config = 1, alate_dens = 1, n_sims = 1000)

wop_sim_df2 <- do.call(one_combo, c(wop_sim_args2, list(n_pseudo = 3L))) |>
    select(rep:outbreak_size)
wop_sim0_df2 <- do.call(one_combo, c(wop_sim_args2, list(n_pseudo = 0L))) |>
    select(rep:outbreak_size)

wop_sim_df2$outbreak_size |> mean()
wop_sim0_df2$outbreak_size |> mean()



wop_sim_args3 <- wop3[[1]][["par"]] |>
    set_names(par_names[1:5]) |>
    trans_pars() |>
    c(Y0 = 9, mean_N = 55, sd_N = 25,
      spat_config = 1, alate_dens = 1, n_sims = 1000)

wop_sim_df3 <- do.call(one_combo, c(wop_sim_args3, list(n_pseudo = 3L))) |>
    select(rep:outbreak_size)
wop_sim0_df3 <- do.call(one_combo, c(wop_sim_args3, list(n_pseudo = 0L))) |>
    select(rep:outbreak_size)

wop_sim_df3$outbreak_size |> mean()
wop_sim0_df3$outbreak_size |> mean()
