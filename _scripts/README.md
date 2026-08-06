
# `_scripts`

## Folder contents:

```
├── 00-preamble.R
├── 01-small
│   ├── 01-extremes.R
│   ├── 02a-heatmap-funs.R
│   └── 02b-extremes-manip.R
├── 02-large-hpc
│   ├── 00-large-hpc-preamble.R
│   ├── 01-large-main.sh
│   ├── 02-large-interior.sh
│   ├── 03-large-densities.R
│   └── 04-large-emp-zeta.R
├── 03-large
│   ├── 01a-large-plot-funs.R
│   ├── 01b-large-results.R
│   ├── 02-large-densities-plots.R
│   ├── 03-large-emp-zeta-plots.R
│   └── 04-virus-attract.R
├── interm-data
│   └── ...
└── README.md
```



## Folder descriptions:

- `01-small`: scripts for small landscape simulations
- `02-large-hpc`: scripts for large landscape simulations that had to be run
  on the high-performance cluster
- `03-large`: scripts for large landscape simulations
- `interm-data`: folder containing intermediate data from simulations;
  this folder has it's own `README.md` file describing its contents.


## File descriptions:

- `01-small` files:
    - `01-extremes.R`: small landscape simulations of extremes (weak and strong
      parasitoid responsiveness to aphid densities) with respect to
      *Pseudomonas* effects on virus transmission
    - `02a-heatmap-funs.R`: functions used to generate heatmaps of virus 
      outcomes (or effects of *Pseudomonas* on those outcomes) versus two
      different parameters
    - `02b-extremes-manip.R`: extremes in small landscapes where we manipulate
      one or two parameters at a time
- `02-large-hpc` files:
    - `00-large-hpc-preamble.R`: preamble used for all large-hpc scripts;
      the biggest job of this script is to create the landscapes of
      *Pseudomonas* locations
    - `01-large-main.sh`: primary large landscape simulations where we
      manipulate parasitoid reponsiveness to aphid densities,
      loading probabilities, variability in initial aphid densities,
      alate attraction to virus-infected plants, and
      alate repellence from *Pseudomonas*-inhabited plants
    - `02-large-interior.sh`: same as `01-large-main.sh`, but initially
      virus-infected plants are located on the interior of the landscape
    - `03-large-densities.R`: simulations for densities of insects and viruses
      across varying *Pseudomonas* occupancies and parasitoid responsiveness
    - `04-large-emp-zeta.R`: empirical zeta estimates in larger landscapes
- `03-large` files:
    - `01a-large-plot-funs.R`: miscellaneous large landscape helper (mostly 
      plotting) functions used in `01b-large-results.R`
    - `01b-large-results.R`: primary large landscape plots
    - `02-large-densities-plots.R`: plots of densities through time in large 
      landscapes
    - `03-large-emp-zeta-plots.R`: plots for empirical estimates of zeta in 
      large landscapes
    - `04-virus-attract.R`: simulations and plots for how alate attraction to 
      virus-infected plants affects virus disease dynamics through time and 
      space
- `README.md`: this file
