#' Create a list in which to store insect population info.
#'
#' This object is used as an argument to [sim_plantscape()].
#'
#' To smooth aphid stage structure out over time,
#' not all 4th instar aphids move to adulthood immediately.
#' If adulthood starts on day `t`, then half of aphids at age `t-2` move to
#' adulthood, and half at age `t-1` do, too.
#' I adjusted age `t-2`, too, to avoid this affecting the growth rate too much.
#'
#' @param pseudo_surv Single numeric indicating aphid survival from *Pseudomonas*.
#' @param zeta Constant between 0 and 1 that affects the extent to which
#'     parasitoids respond to aphid density, where 0 results in an even
#'     distribution of parasitoids, and 1 results in a linear relationship
#'     between aphid density and parasitoids.
#' @param extinct_N Single numeric indicating the extinction threshold.
#'     Defaults to `0`.
#' @param demog_error Single logical for whether to include demographic
#'     stochasticity. Defaults to `FALSE`.
#' @param environ_error Single logical for whether to include environmental
#'     stochasticity. Defaults to `FALSE`.
#' @param K A single number indicating aphid density dependence.
#'     Defaults to `12.5e3` because this caused simulations to
#'     approximately match experiments.
#' @param K_p_mult A single number that is multiplied by `K` to get density
#'     dependence for parasitized aphids.
#'     Defaults to `1 / 1.57`, which is from Meisner et al. (2014).
#' @param pred_surv Daily survival from predation for aphids and mummies.
#'     Defaults to `0.9`.
#' @param alate_b0 The proportion of offspring from apterous aphids is
#'     `inv_logit(alate_b0 + alate_b1 * N)` where `N` is the total number of
#'     aphids in that field.
#'     Defaults to `-5`.
#' @param alate_b1 The proportion of offspring from apterous aphids is
#'     `inv_logit(alate_b0 + alate_b1 * N)` where `N` is the total number of
#'     aphids in that field.
#'     Defaults to `0.0022`, which makes alate production only mildly
#'     density dependent.
#' @param fly_p Single numeric indicating the proportion of alates that fly
#'     off plants each day. Defaults to `0.05`.
#' @param distr_0 A 2-column matrix indicating the starting stage distribution.
#'     (The total abundance of aphids is adjusted in the simulation function.)
#'     If `NULL`, there will be no starting alates and apterous densities
#'     will follow the stable age distribution.
#'     If a matrix, it must have 5 rows (rows indicate instar)
#'     or a row per aphid age (rows indicate age in days).
#'     Matrix column indicates apterous vs alate.
#'     The matrix will be coerced to sum to 1.
#'     Defaults to `NULL`.
#' @param resistant Logical or vector of survivals of
#'     singly attacked and multiply attacked aphids.
#'     If a logical, `FALSE` is equivalent to `c(0,0)` and results in no
#'     resistance.
#'     `TRUE` results in the resistance values for a resistance line
#'     from unpublished work by Anthony Ives.
#'     Defaults to `FALSE`.
#' @param surv_juv_apterous A single number for the juvenile survival rate for
#'     apterous aphids. Defaults to `NULL`, which results in estimates
#'     from a medium-reproduction line.
#' @param surv_adult_apterous A vector of adult survival probabilities for
#'     apterous aphids. Defaults to `NULL`, which results in estimates
#'     from a medium-reproduction line.
#' @param repro_apterous A vector of fecundities for for apterous aphids.
#'     Defaults to `NULL`, which results in estimates from a
#'     medium-reproduction line.
#' @param surv_juv_alates A single number for the juvenile survival rate for
#'     alates aphids. Defaults to `"low"`, which results in estimates
#'     from a low-reproduction line.
#' @param surv_adult_alates A vector of adult survival probabilities for
#'     alates aphids. Defaults to `"low"`, which results in estimates
#'     from a low-reproduction line.
#' @param repro_alates A vector of fecundities for for alates aphids.
#'     Defaults to `"low"`, which results in estimates from a low-reproduction
#'     line.
#' @param surv_paras A single number for the juvenile survival rate for
#'     parasitized aphids.
#'     Because parasitized aphids don't make it to adulthood, this is the only
#'     survival rate necessary.
#'     Defaults to `"low"`, which results in estimates from a
#'     low-reproduction line.
#' @param temp Single string specifying `"low"` (20º C) or `"high"` (27º C)
#'     temperature. Defaults to `"low"`.
#' @param p_instar_smooth A value greater than zero here makes not all 4th
#'     instar aphids go to the adult stage at the same time.
#'     Of the 4th instars that would've moved to adulthood, `p_instar_smooth`
#'     remain as 4th instars.
#'     To keep the growth rate approximately equal, for the aphids that
#'     would've transitioned to 4th instar, `p_instar_smooth` move to
#'     adulthood instead.
#'     Defaults to `0.5`.
#' @param sigma_x Standard deviation of environmental stochasticity for aphids.
#'     This argument has no effect (is changed to `0`) when
#'     `environ_error = FALSE`.
#'     Defaults to the internal object `environ$sigma_x`,
#'     which is from Meisner et al. (2014).
#' @param rho Environmental correlation among instars.
#'     This argument has no effect if used in simulations with no environmental
#'     error (i.e., when `environ_error = FALSE`).
#'     Defaults to `environ$rho`.
#' @param s_y Daily survival rate of adult wasps.
#'     Defaults to `populations$s_y`, which is from Meisner et al. (2014).
#' @param sex_ratio Sex ratio of adult wasps. Defaults to `0.5`.
#' @param a Parasitoid attack rate. Defaults to the internal object
#'     `wasp_attack$a`, which is from Meisner et al. (2014).
#' @param k Aggregation parameter of the negative binomial distribution.
#'     Defaults to the internal object `wasp_attack$k`,
#'     which is from Meisner et al. (2014).
#' @param h Parasitoid handling time. Defaults to the internal object
#'     `wasp_attack$h`, which is from Meisner et al. (2014).
#' @param rel_attack Relative parasitoid attack rate among instars.
#'     Defaults to `wasp_attack$rel_attack`, which is from
#'     Meisner et al. (2014).
#' @param sigma_y Standard deviation of environmental stochasticity for wasps.
#'     This argument has no effect (is changed to `0`) if
#'     `environ_error = FALSE`.
#'     Defaults to the internal object `environ$sigma_y`,
#'     which is from Meisner et al. (2014).
#' @param mum_smooth Proportion of mummies that will NOT take exactly 3 days
#'     to develop. As this value approaches 2/3, it will provide greater
#'     smoothing of wasp numbers through time.
#'     Defaults to `0.4`.
#'
#'
#' @return An `externalptr` object that points to a C++ object that can
#' be passed to [sim_plantscape()].
#'
#' @export
#'
make_insect_ptr <- function(pseudo_surv,
                            zeta,
                            extinct_N = 0,
                            demog_error = FALSE,
                            environ_error = FALSE,
                            K = populations$K,
                            K_p_mult = populations$K_p_mult,
                            pred_surv = 0.9,
                            alate_b0 = -5,
                            alate_b1 = 0.0022,
                            fly_p = 0.05,
                            distr_0 = NULL,
                            resistant = FALSE,
                            surv_juv_apterous = NULL,
                            surv_adult_apterous = NULL,
                            repro_apterous = NULL,
                            surv_juv_alates = "low",
                            surv_adult_alates = "low",
                            repro_alates = "low",
                            surv_paras = "low",
                            temp = "low",
                            p_instar_smooth = 0.5,
                            sigma_x = environ$sigma_x,
                            rho = environ$rho,
                            s_y = populations$s_y,
                            sex_ratio = populations$sex_ratio,
                            a = wasp_attack$a,
                            k = wasp_attack$k,
                            h = wasp_attack$h,
                            rel_attack = wasp_attack$rel_attack,
                            sigma_y = environ$sigma_y,
                            mummy_smooth = 0.4) {

    single_logical(environ_error, "environ_error")
    if (!environ_error) {
        sigma_x <- 0
        sigma_y <- 0
    }

    aphid_ptr <- aphid_pop(pseudo_surv = pseudo_surv,
                           name = "aphid",
                           K = K,
                           K_p_mult = K_p_mult,
                           pred_surv = pred_surv,
                           alate_b0 = alate_b0,
                           alate_b1 = alate_b1,
                           fly_p = fly_p,
                           distr_0 = distr_0,
                           resistant = resistant,
                           surv_juv_apterous = surv_juv_apterous,
                           surv_adult_apterous = surv_adult_apterous,
                           repro_apterous = repro_apterous,
                           surv_juv_alates = surv_juv_alates,
                           surv_adult_alates = surv_adult_alates,
                           repro_alates = repro_alates,
                           surv_paras = surv_paras,
                           temp = temp,
                           p_instar_smooth = p_instar_smooth,
                           sigma_x = sigma_x,
                           rho = rho,
                           extinct_N = extinct_N,
                           demog_error = demog_error)

    wasp_ptr <- wasp_pop(zeta = zeta,
                         s_y = s_y,
                         sex_ratio = sex_ratio,
                         a = a,
                         k = k,
                         h = h,
                         rel_attack = rel_attack,
                         sigma_y = sigma_y,
                         demog_error = demog_error,
                         temp = temp,
                         n_aphid_age_stages = get_aphid_n_age_stages(aphid_ptr))

    mummy_ptr <- mummy_pop(pred_surv = pred_surv,
                           mummy_smooth = mummy_smooth,
                           extinct_N = extinct_N)

    ptr <- make_insect_ptr_cpp(aphid_ptr, wasp_ptr, mummy_ptr)

    return(ptr)

}




aphid_pop <- function(pseudo_surv,
                      name,
                      K,
                      K_p_mult,
                      pred_surv,
                      alate_b0,
                      alate_b1,
                      fly_p,
                      distr_0,
                      resistant,
                      surv_juv_apterous,
                      surv_adult_apterous,
                      repro_apterous,
                      surv_juv_alates,
                      surv_adult_alates,
                      repro_alates,
                      surv_paras,
                      temp,
                      p_instar_smooth,
                      sigma_x,
                      rho,
                      extinct_N,
                      demog_error) {

    single_number(pseudo_surv, "pseudo_surv", .min = 0, .max = 1)
    single_string(name, "name")
    single_number(K, "K", .min = 0)
    single_number(K_p_mult, "K_p_mult", .min = 0)
    single_number(pred_surv, "pred_surv", .min = 0, .max = 1)
    single_number(alate_b0, "alate_b0")
    single_number(alate_b1, "alate_b1")
    single_number(fly_p, "fly_p", .min = 0, .max = 1)
    single_number(surv_paras, "surv_paras", .min = 0, .max = 1)
    single_string(temp, "temp")
    single_number(p_instar_smooth, "p_instar_smooth", .min = 0, .max = 1)
    single_number(sigma_x, "sigma_x", .min = 0)
    single_number(rho, "rho", .min = 0)
    single_number(extinct_N, "extinct_N", .min = 0)
    single_logical(demog_error, "demog_error")
    # The remaining arguments (distr_0, resistant, surv_juv_apterous,
    # surv_adult_apterous, repro_apterous, surv_juv_alates, surv_adult_alates,
    # repro_alates) are checked below

    temp <- match.arg(temp, c("low", "high"))
    temp <- paste0(temp, "T")


    # --------------*
    # Construct Leslie matrices
    # --------------*

    leslie <- list(apterous = NA,
                   alates = NA,
                   paras = NA)
    # In `leslie_mat` below, items in vector are aphid lines, slices are
    # apterous/alate/parasitized.

    # Set with default values for Leslie matrix calculations
    def_L_args <- list(instar_days = dev_times$instar_days[[temp]],
                       surv_juv = mean(do.call(c, populations$surv_juv)),
                       surv_adult = colMeans(do.call(rbind,
                                                     populations$surv_adult)),
                       repro = colMeans(do.call(rbind, populations$repro)))

    # Note that we don't have values for parasitized aphid adult survival
    # or fecundity because they don't reach adulthood anyway.
    # Below, the defaults will be filled in for these, then the downstream
    # code will ignore that info.
    inputs <- list("surv_juv_apterous" = surv_juv_apterous,
                   "surv_juv_alates" = surv_juv_alates,
                   "surv_juv_paras" = surv_paras,
                   "surv_adult_apterous" = surv_adult_apterous,
                   "surv_adult_alates" = surv_adult_alates,
                   "repro_apterous" = repro_apterous,
                   "repro_alates" = repro_alates)

    for (x in c("apterous", "alates", "paras")) {
        leslie_args <- def_L_args
        for (y in names(inputs)[grepl(paste0(x, "$"), names(inputs))]) {
            arg_name <- gsub(paste0("_", x), "", y)
            if (!is.null(inputs[[y]])) {
                if (is.numeric(inputs[[y]])) {
                    leslie_args[[arg_name]] <- inputs[[y]]
                } else if (inputs[[y]] %in% c("low", "high")) {
                    leslie_args[[arg_name]] <- populations[[arg_name]][[
                        inputs[[y]]]]
                } else {
                    msg <- paste0("\ninput argument `", y,
                                  "` to the `clonal_line` function should be ",
                                  "NULL, a numeric vector, \"low\", or ",
                                  "\"high\"")
                    stop(msg)
                }
            }
        }
        leslie[[x]] <- do.call(leslie_matrix, leslie_args)
    }

    if (!identical(dim(leslie[[1]]), dim(leslie[[2]])) ||
        !identical(dim(leslie[[1]]), dim(leslie[[3]]))) {
        stop("\nLeslie matrices for apterous, alates, and parasitized",
             " aphids must all be of the same dimensions\n")
    }

    leslie_array <- array(do.call(c, leslie), dim = c(dim(leslie[[1]]), 3))
    ns <- nrow(leslie_array)  # number of stages; used later

    # To make 4th instars not always move to adulthood at the same time:
    if (p_instar_smooth > 0) {
        .adult <- sum(head(dev_times$instar_days[[temp]], -1)) + 1
        # I included `- 1` bc we don't to adjust the Leslie matrix for the
        # parasitized aphids
        for (j in 1:(dim(leslie_array)[3] - 1)) {
            # Of the aphids that would've moved to adulthood, make
            # `p_instar_smooth` remain as 4th instars of age `.adult-1` instead:
            .t <- .adult - 1
            surv_t <- leslie_array[.adult, .t, j]
            leslie_array[.t, .t, j] <- surv_t * p_instar_smooth
            leslie_array[.adult, .t, j] <- surv_t * (1 - p_instar_smooth)
            # Of the aphids that would've moved to age `.adult-1`, make
            # `p_instar_smooth` move to adulthood instead:
            .t <- .adult - 2
            surv_t <- leslie_array[.t+1, .t, j]
            leslie_array[.adult, .t, j] <- surv_t * p_instar_smooth
            leslie_array[.t+1, .t, j] <- surv_t * (1 - p_instar_smooth)
        }
    }



    # --------------*
    # Fill starting age distributions
    # --------------*

    if (is.null(distr_0)) {
        # If provided with just one number, assume stable age distribution:
        sad <- as.numeric(sad_leslie(leslie_array[,,1]))
        distr_0 <- cbind(sad, rep(0, length(sad)))
    } else if (is.numeric(distr_0) && is.matrix(distr_0) &&
               identical(dim(distr_0), c(5L, 2L))) {
        d0 <- distr_0
        distr_0 <- matrix(0, ns, 2)
        # Going from instar to days old, using stable age distribution to
        # calculate the proportion for each age within each instar.
        sads <- cbind(as.numeric(sad_leslie(leslie_array[,,1])),
                      as.numeric(sad_leslie(leslie_array[,,2])))
        idays <- dev_times$instar_days[[temp]]
        d0_inds <- cbind(c(1, head(cumsum(idays), -1) + 1),
                         c(head(cumsum(idays), -1), nrow(sads)))
        for (i in 1:nrow(d0)) {
            irange <- d0_inds[i,1]:d0_inds[i,2]
            for (j in 1:2) {
                sad_ij = sads[irange, j] / sum(sads[irange,j])
                distr_0[irange,j] <- sad_ij * d0[i,j]
            }
        }
    } else if (! (is.numeric(distr_0) && is.matrix(distr_0) &&
                  identical(dim(distr_0), c(ns, 2L)))) {
        stop("\nThe `distr_0` arg to the ",
             "`clonal_line` function must be a single number or a 5x2 or ",
             ns, "x2 numeric matrix")
    }
    if (any(distr_0 < 0)) {
        stop("\nThe `distr_0` arg to the `clonal_line` function ",
             "cannot contain values < 0.")
    }
    if (sum(distr_0) != 1) distr_0 <- distr_0 / sum(distr_0)


    # --------------*
    # Fill other info
    # --------------*

    attack_surv <- c(0, 0)
    if (!(is.logical(resistant) && length(resistant) == 1) &&
        !(is.numeric(resistant) && all(resistant >= 0 & resistant <= 1))) {
        stop("`resistant` must be a single logical or a numeric vector with ",
             "all elements >= 0 and <= 1.")
    }
    if (is.logical(resistant) && resistant) attack_surv <- wasp_attack$attack_surv
    if (is.numeric(resistant)) attack_surv <- resistant

    # assume only adult alates can disperse across fields:
    field_disp_start <- sum(head(dev_times$instar_days[[temp]], -1))

    living_days <- dev_times$mum_days[[1]]


    K_p <- K * K_p_mult

    # --------------*
    # Make final output
    # --------------*

    aphid_ptr <- make_aphid_ptr(K = K,
                                K_p = K_p,
                                pseudo_surv = pseudo_surv,
                                pred_surv = pred_surv,
                                extinct_N = extinct_N,
                                demog_error = demog_error,
                                sigma_x = sigma_x,
                                rho = rho,
                                fly_p = fly_p,
                                attack_surv = attack_surv,
                                aphid_name = name,
                                leslie_mat = leslie_array,
                                aphid_density_0 = distr_0,
                                alate_b0 = alate_b0,
                                alate_b1 = alate_b1,
                                fly_start = field_disp_start,
                                living_days = living_days)

    return(aphid_ptr)

    # output <- structure(list(name = name,
    #                          distr_0 = distr_0,
    #                          attack_surv = attack_surv,
    #                          leslie = leslie_array,
    #                          temp = temp,
    #                          ptr = aphid_ptr),
    #                     class =  "AphidPop")
    #
    # return(output)

}


# #'
# #' @export
# #' @noRd
# #'
# print.AphidPop <- function(x, ...) {
#
#     cat("< Aphid clonal line >\n")
#     cat("Name: ", x$name, "\n", sep = "")
#     cat("Fields:\n")
#     cat("  * name <string>\n")
#     cat("  * distr_0 <matrix>\n")
#     cat("  * attack_surv <vector>\n")
#     cat("  * leslie <3D array>\n")
#     cat("  * temp <string>\n")
#     cat("  * ptr <externalptr>\n")
#
#     invisible(x)
#
# }






wasp_pop <- function(zeta,
                     s_y,
                     sex_ratio,
                     a,
                     k,
                     h,
                     rel_attack,
                     sigma_y,
                     demog_error,
                     temp,
                     n_aphid_age_stages) {

    single_number(zeta, "zeta", .min = 0, .max = 1)
    single_number(s_y, "s_y", .min = 0, .max = 1)
    single_number(sex_ratio, "sex_ratio", .min = 0, .max = 1)
    single_number(a, "a", .min = 0)
    single_number(k, "k", .min = 0)
    single_number(h, "h", .min = 0)
    is_type(rel_attack, "rel_attack", is.numeric, .min = 0, L = 5L)
    single_number(sigma_y, "sigma_y", .min = 0)
    single_logical(demog_error, "demog_error")

    rel_attack <- rel_attack / sum(rel_attack)
    dt <- dev_times$instar_days[[temp]]
    n_adult_days <- n_aphid_age_stages - sum(head(dt, -1))
    stopifnot(n_adult_days >= 0)
    if (tail(dt, 1) != n_adult_days) dt[length(dt)] <- n_adult_days
    # Commented version isn't used bc it wasn't done this way when fitting
    # the model.
    # rel_attack__ <- mapply(function(.x, .y) rep(.x, .y) / .y,
    #                        rel_attack,  dt)
    rel_attack__ <- mapply(rep, rel_attack, dt)
    rel_attack <- do.call(c, rel_attack__)

    ptr <- make_wasps_ptr(rel_attack = rel_attack,
                          a = a,
                          k = k,
                          h = h,
                          sex_ratio = sex_ratio,
                          s_y = s_y,
                          zeta = zeta,
                          sigma_y = sigma_y,
                          demog_error = demog_error)

    return(ptr)

}


mummy_pop <- function(pred_surv,
                      mummy_smooth,
                      extinct_N) {

    single_number(pred_surv, "pred_surv", .min = 0, .max = 1)
    single_number(mummy_smooth, "mummy_smooth", .min = 0, .max = 2/3)
    single_number(extinct_N, "extinct_N", .min = 0)

    mummy_dev_time <- dev_times$mum_days[2]
    single_integer(mummy_dev_time, "mummy_dev_time", .min = 1)

    ptr <- make_mummies_ptr(pred_surv = pred_surv,
                            mummy_smooth = mummy_smooth,
                            mummy_dev_time = mummy_dev_time,
                            extinct_N = extinct_N)

    return(ptr)
}






