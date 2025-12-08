
get_maxes <- function(dd, .alate_max, .aphid_max, .wasp_max) {

    dens_maxes <- c("aphids", "alates", "wasps") |>
        set_names() |>
        map(\(spp) {
            dd$density[dd$species == spp] |>
                max(na.rm = TRUE)
        })
    if (!is.null(.alate_max)) {
        if (.alate_max < dens_maxes$alates) {
            stop("\nIncrease .alate_max to at least ", dens_maxes$alates)
        }
        dens_maxes$alates <- .alate_max
    }
    if (!is.null(.aphid_max)) {
        if (.aphid_max < dens_maxes$aphids) {
            stop("\nIncrease .aphid_max to at least ", dens_maxes$aphids)
        }
        dens_maxes$aphids <- .aphid_max
    }
    if (!is.null(.wasp_max)) {
        if (.wasp_max < dens_maxes$wasps) {
            stop("\nIncrease .wasp_max to at least ", dens_maxes$wasps)
        }
        dens_maxes$wasps <- .wasp_max
    }

    return(dens_maxes)

}






sim_plotter <- function(sims,
                        zeta = NULL,
                        .title  = waiver(),
                        .tag = waiver(),
                        .alate_max = NULL,
                        .aphid_max = NULL,
                        .wasp_max = NULL) {
    # zeta = 0; .title = waiver(); .tag = waiver(); .alate_max = NULL; .aphid_max = NULL; .wasp_max = NULL
    # rm(zeta, .title, .tag, .alate_max, .aphid_max, .wasp_max, dd, dens_maxes)
    # rm(trans, itrans, aphid_breaks, aphid_labels, aphid_ylab, p)
    # aphids <- sims |>
    #     filter(is.na(wasps)) |>
    #     mutate(plant = interaction(x, y), rep = factor(rep)) |>
    #     mutate(aphids = aphids + parasitized) |>
    #     select(rep, plant, time, aphids, alates, mummies)
    # virus <- sims |>
    #     filter(is.na(wasps)) |>
    #     mutate(plant = interaction(x, y), rep = factor(rep)) |>
    #     group_by(plant) |>
    #     filter(virus == 1) |>
    #     filter(time == min(time)) |>
    #     ungroup() |>
    #     select(rep, plant, time, virus)
    if (!is.null(zeta)) {
        dd <- sims |>
            group_by(rep, time) |>
            mutate(z_hat = (aphids + alates + parasitized) / sum(aphids[is.na(x)] + alates[is.na(x)] + parasitized[is.na(x)]),
                   wasps = ifelse(is.na(wasps), wasps[!is.na(wasps)] * ((1 - zeta) / 9 + zeta * z_hat),
                                  wasps)) |>
            ungroup() |>
            filter(!is.na(x)) |>
            mutate(plant = interaction(x, y), rep = factor(rep),
                   aphids = aphids + parasitized) |>
            select(rep, plant, time, aphids, alates, wasps) |>
            pivot_longer(aphids:wasps, names_to = "species", values_to = "density")
    } else {
        dd <- sims |>
            filter(is.na(x)) |>
            mutate(rep = factor(rep),
                   aphids = aphids + parasitized) |>
            select(rep, time, aphids, alates, wasps) |>
            pivot_longer(aphids:wasps, names_to = "species", values_to = "density")
    }

    dens_maxes <- get_maxes(dd, .alate_max, .aphid_max, .wasp_max)


    # convert from wasps or alates --> aphids:
    trans <- c("wasps", "alates") |>
        set_names() |>
        map(\(n) {
            denom <- paste0(dens_maxes[[n]])
            eval(parse(text = paste0("function(x) { return(" ,
                                     "x * ", dens_maxes$aphids, " / ", denom, ") }")))
        })
    # convert from aphids --> wasps or alates:
    itrans <- c("wasps", "alates") |>
        set_names() |>
        map(\(n) {
            numer <- paste0(dens_maxes[[n]])
            eval(parse(text = paste0("function(x) { return(" ,
                                     "x * ", numer, "/ ", dens_maxes$aphids, ") }")))
        })

    aphid_breaks <- scales::breaks_extended(n = 4)(c(0, dens_maxes$aphids))
    aphid_labels <- sprintf(paste0("%s (<span style=\"color: ", spp_pal[["alates"]],
                                   ";\">%.1f</span>)"), aphid_breaks,
                            itrans$alates(aphid_breaks))
    aphid_ylab <- paste0("Aphid density (<span style=\"color: ", spp_pal[["alates"]],
                         ";\">alate density</span>)")

    p <- dd |>
        mutate(density = case_when(species == "aphids" ~ density,
                                   species == "alates" ~ trans$alates(density),
                                   species == "wasps" ~ trans$wasps(density),
                                   .default = NA)) |>
        ggplot(aes(time, density)) +
        geom_line(aes(color = species), linewidth = 1) +
        facet_wrap( ~ plant, nrow = 3) +
        scale_y_continuous(aphid_ylab,
                           breaks = aphid_breaks,
                           labels = aphid_labels,
                           sec.axis = sec_axis(itrans$wasps, "Wasp density")) +
        scale_color_manual(values = spp_pal, guide = "none") +
        labs(x = "Time (days)",
             tag = .tag,
             title = .title) +
        coord_cartesian(ylim = c(0, dens_maxes$aphids)) +
        theme(plot.title = element_markdown(hjust = 0.5),
              plot.tag = element_markdown(size = 16),
              plot.tag.location = "plot",
              panel.grid.major.y = element_line(color = "gray80"),
              axis.title.y.left = element_markdown(color = spp_pal[["aphids"]],
                                                   face = "bold"),
              axis.title.y.right = element_markdown(color = spp_pal[["wasps"]],
                                                    face = "bold"),
              axis.text.y.left = element_markdown(color = spp_pal[["aphids"]]),
              axis.text.y.right = element_markdown(color = spp_pal[["wasps"]]))


    # p <- aphids |>
    #     ggplot(aes(time)) +
    #     geom_line(aes(y = aphids, color = plant)) +
    #     geom_line(data = wasps, aes(y = trans(wasps)), linewidth = 1) +
    #     scale_y_continuous("Aphid density",
    #                        sec.axis = sec_axis(itrans, "Wasp density")) +
    #     labs(title = .title) +
    #     theme(plot.title = element_markdown(hjust = 0, family = ""))
    #
    # if (!is.null(zeta)) {
    #     p <- p +
    #         # geom_vline(data = virus, aes(xintercept = time), color = "#EC008C",
    #         #            linetype = "22") +
    #         facet_wrap( ~ plant, nrow = 3) +
    #         scale_color_viridis_d(begin = 0.1, end = 0.9, guide = "none")
    # } else p <- p + scale_color_viridis_d(begin = 0.1, end = 0.9)

    return(p)

}
