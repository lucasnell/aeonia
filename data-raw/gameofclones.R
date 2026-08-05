
# This creates similar datasets to what's available in `gameofclones`, but with
# a couple of changes specific to `aeonia`.
#
# Note that I've changed the names of these objects for this package (added
# prefix "an_") to avoid conflicts when I load gameofclones for testing.


# ** Uncomment these to run this script **
# if (!require("gameofclones", quietly = TRUE)) remotes::install_github("lucasnell/gameofclones@v1.0.2")
# library(gameofclones)  # v1.0.2

# Development times for aphids and wasps.
an_dev_times <- dev_times
# Save to package:
usethis::use_data(an_dev_times, overwrite = TRUE)



# Population rates and starting values for aphids and wasps.
an_populations <- populations
# Rename and remove some items:
an_populations[["K"]] <- 12.5e3
an_populations[["K_p_mult"]] <- 1 / 1.57
an_populations[["K_y"]] <- NULL
an_populations[["aphids_0"]] <- NULL
an_populations[["wasps_0"]] <- NULL
an_populations[["prop_resist"]] <- NULL
# reorder:
an_populations <- an_populations[c("surv_juv", "surv_adult", "repro", "K",
                                   "K_p_mult", "s_y", "sex_ratio")]
# Save to package:
usethis::use_data(an_populations, overwrite = TRUE)



# Wasp attack rate parameters.
an_wasp_attack <- wasp_attack
# Save to package:
usethis::use_data(an_wasp_attack, overwrite = TRUE)



# Parameters associated with environmental effects and stochasticity.
an_environ <- environ
# Remove some items:
an_environ[["harvest_surv"]] <- NULL
an_environ[["disp_aphid"]] <- NULL
an_environ[["disp_wasp"]] <- NULL
an_environ[["pred_rate"]] <- NULL
an_environ[["cycle_length"]] <- NULL
#
# From "The Role of Aphid Behaviour in the Epidemiology of Potato Virus Y:
# a Simulation Study" by Thomas Nemecek (1993; p. 72;
# doi: 10.3929/ethz-a-000909462), dispersal distances
# follow a Weibull distribution with shape = 0.6569 and scale = 9.613.
#
# For larger landscapes, we use a `radius` argument that is the median of this
# distribution.
# I'm dividing by 0.75 to convert from meters to plant locations that are
# 0.75 meters apart (typical spacing for pea):
#
an_environ[["radius"]] <- qweibull(0.5, 0.6569, 9.613) / 0.75
# Save to package:
usethis::use_data(an_environ, overwrite = TRUE)





