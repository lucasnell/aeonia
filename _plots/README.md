
# `_plots`



## Folder contents (after scripts run):

```
├── main
│   ├── densities-large.pdf
│   ├── extremes
│   │   ├── disease-n_infected.pdf
│   │   └── timeseries.pdf
│   ├── extremes-manip-lines-pseudo_surv.pdf
│   ├── large-baselines.pdf
│   ├── large-outcomes.pdf
│   ├── lines-zeta.pdf
│   └── z-p_alates.pdf
├── README.md
└── supp
    ├── dens-dep-dispersal.pdf
    ├── extremes-disease-outcomes.pdf
    ├── extremes-manip-lines-K-all-outcomes.pdf
    ├── extremes-manip-lines-pseudo_surv-all-outcomes.pdf
    ├── extremes-manip-lines-zeta-all-outcomes.pdf
    ├── extremes-manip-thresholds.pdf
    ├── extremes-manip-Y0-N0-heatmaps.pdf
    ├── large-baseline-wt_pp.pdf
    ├── large-empirical-zeta.pdf
    ├── large-manips-interior.pdf
    ├── large-outcomes-small-outbreaks.pdf
    ├── large-small-manips.pdf
    ├── large-zeta-n_infected.pdf
    ├── sd_N-densities.pdf
    └── virus-attract.pdf
```


## Folder descriptions

The `main` folder has figures found in the main text, while the `supp`
folder contains figures for the supplement.
The main text figures are bare bones (e.g., no axis titles) because they were
imported into Illustrator for the final versions.
No changes were made to the panels themselves in Illustrator.



## File descriptions:

- `main` folder:
    - `densities-large.pdf`: landscape-wide densities through time 
      for total aphids, alates, and parasitoids for large landscapes
      and across occupancies of *Pseudomonas*
      (created in `_scripts/03-large/01b-large-results.R`)
    - `disease-n_infected.pdf`: frequencies of the numbers of peak infected
      plants by *Pseudomonas* presence in small landscapes
      (created in `_scripts/01-small/01-extremes.R`)
    - `timeseries.pdf`: densities of total aphids, alates, parasitoids, and
      infected plants through time for small landscapes and with or without
      *Pseudomonas* (created in `_scripts/01-small/01-extremes.R`)
    - `extremes-manip-lines-pseudo_surv.pdf`: number of peak infected plants 
      versus aphid *Pseudomonas* infection in small landscapes
      (created in `_scripts/01-small/02b-extremes-manip.R`)
    - `large-outcomes.pdf`: number of peak infected plants in large landscapes 
      across occupancies of *Pseudomonas* and across whether 
      the virus starts on a *Pseudomonas* plant
      (created in `_scripts/03-large/01b-large-results.R`)
    - `lines-zeta.pdf`: number of peak infected plants in small landscapes 
      versus parasitoid responsiveness to aphid densities
      (created in `_scripts/01-small/02b-extremes-manip.R`)
    - `z-p_alates.pdf`: proportion of offspring that are alates versus the total
      density of aphids on a plant
      (created in `_scripts/01-small/01-extremes.R`)
- `README.md`: this file
- `supp` folder:
    - `dens-dep-dispersal.pdf`: virus outcomes in small landscapes versus
      parasitoid responsiveness to aphid densities with density dependent and
      fixed alate production
      (created in `_scripts/01-small/01-extremes.R`)
    - `extremes-disease-outcomes.pdf`: probabilities of emergence and
      frequencies of the outbreak sizes versus *Pseudomonas* presence in 
      small landscapes
      (created in `_scripts/01-small/01-extremes.R`)
    - `extremes-manip-lines-K-all-outcomes.pdf`: virus outcomes in small 
      landscapes versus aphid carrying capacities
      (created in `_scripts/01-small/02b-extremes-manip.R`)
    - `extremes-manip-lines-pseudo_surv-all-outcomes.pdf`: probabilities of 
      emergence and mean outbreak sizes in small landscapes versus
      aphid survival from *Pseudomonas*
      (created in `_scripts/01-small/02b-extremes-manip.R`)
    - `extremes-manip-lines-zeta-all-outcomes.pdf`: probabilities of 
      emergence and mean outbreak sizes in small landscapes versus
      parasitoid responsiveness to aphid densities
      with and without *Pseudomonas*
      (created in `_scripts/01-small/02b-extremes-manip.R`)
    - `extremes-manip-thresholds.pdf`: time series of densities when the ratio
      of starting parasitoids to aphids is too high or too low
      (created in `_scripts/01-small/02b-extremes-manip.R`)
    - `extremes-manip-Y0-N0-heatmaps.pdf`: heatmaps showing how starting 
      parasitoid and aphid densities in small landscapes change the effects
      of *Pseudomonas* presence on emergence probabilities, outbreak sizes,
      and numbers of peak infected plants
      (created in `_scripts/01-small/02b-extremes-manip.R`)
    - `large-baseline-wt_pp.pdf`: numbers of peak infected plants versus
      *Pseudomonas* occupancy when *Pseudomonas* locations are clustered
      (created in `_scripts/03-large/01b-large-results.R`)
    - `large-baselines.pdf`: number of peak infected plants in large landscapes 
      across occupancies of *Pseudomonas* and across whether 
      the virus starts on a *Pseudomonas* plant,
      *Pseudomonas* repels alates, and virus infection attracts alates
      (created in `_scripts/03-large/01b-large-results.R`)
    - `large-empirical-zeta.pdf`: observed / predicted parasitoid abundances
      across parasitoid responsiveness to aphid densities, where predictions
      are based on empirical data
      (created in `_scripts/03-large/03-large-emp-zeta-plots.R`)
    - `large-manips-interior.pdf`: number of peak infected plants in large
      landscapes where the virus starts on an interior plant,
      across occupancies of *Pseudomonas* and across whether 
      the virus starts on a *Pseudomonas* plant,
      *Pseudomonas* repels alates, and virus infection attracts alates
      (created in `_scripts/03-large/01b-large-results.R`)
    - `large-outcomes-small-outbreaks.pdf`: virus outcomes in large landscapes 
      when loading probabilities are low
      across occupancies of *Pseudomonas*
      (created in `_scripts/03-large/01b-large-results.R`)
    - `large-small-manips.pdf`: virus outcomes in large landscapes 
      when loading probabilities are low
      across occupancies of *Pseudomonas* and across whether 
      the virus starts on a *Pseudomonas* plant,
      *Pseudomonas* repels alates, and virus infection attracts alates
      (created in `_scripts/03-large/01b-large-results.R`)
    - `large-zeta-n_infected.pdf`: number of peak infected plants versus
      parasitoid responsiveness to aphid densities, for varying levels
      of *Pseudomonas* occupancy; also present on this figure are points 
      indicating the levels of parasitoid responsiveness corresponding
      to empirical estimates
      (created in `_scripts/03-large/03-large-emp-zeta-plots.R`)
    - `sd_N-densities.pdf`: total aphid density through time, with and 
      without variation among plants in initial aphid abundances
      (created in `_scripts/03-large/01b-large-results.R`)
    - `virus-attract.pdf`: for alate virus attraction being absent or present,
      numbers of incoming alates for infected plants across the landscape
      and time series of newly infected plants and numbers of adult alates
      (created in `_scripts/03-large/04-virus-attract.R`)
