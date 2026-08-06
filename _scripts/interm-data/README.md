
# `interm-data`

## Folder contents (after scripts run):

```
├── dd_disp-sims.csv
├── density-sims.csv
├── extremes-manip-sims.rds
├── extremes-manip2-sims.rds
├── large-interior
│   ├── large-interior-01.rds
│   ├── ...
│   └── large-interior-32.rds
├── large-landscapes.rds
├── large-main
│   ├── large-main-01.rds
│   ├── ...
│   └── large-main-16.rds
├── large-zeta-ninf-summs.csv
├── large-zeta-nopseudo-summs.csv
├── large-zeta-summs.csv
└── README.md
```



## File / folder descriptions:

- `dd_disp-sims.csv`: simulations of density-dependent and fixed alate 
  production in small landscapes (created in `_scripts/01-small/01-extremes.R`)
- `density-sims.csv`: large-landscape simulations for densities of insects and 
  viruses across varying *Pseudomonas* occupancies and parasitoid responsiveness
  (created in `_scripts/02-large-hpc/03-large-densities.R`)
- `extremes-manip-sims.rds`: small-landscape simulations where we manipulate
  one parameter at a time for weak and strong parasitoid responsiveness
  (created in `_scripts/01-small/02b-extremes-manip.R`)
- `extremes-manip2-sims.rds`: small-landscape simulations where we manipulate
  two parameters at a time (only starting aphid and parasitoid densities)
  for weak and strong parasitoid responsiveness
  (created in `_scripts/01-small/02b-extremes-manip.R`)
- `large-interior`: same as `large-main`, but initially
  virus-infected plants are located on the interior of the landscape
  this folder contains 32 files, `large-interior-01.rds` -- `large-interior-32`
  (created in `_scripts/02-large-hpc/02-large-interior.sh`)
- `large-landscapes.rds`: dataframe containing simulated large landscapes,
  where we varied the percent of plants with *Pseudomonas*, 
  whether the virus started on a *Pseudomonas*-inhabited plant, and
  clustering of *Pseudomonas* locations
  (created in `_scripts/02-large-hpc/00-large-hpc-preamble.R`)
- `large-main`: primary large landscape simulations where we
  manipulate parasitoid reponsiveness to aphid densities,
  loading probabilities, variability in initial aphid densities,
  alate attraction to virus-infected plants, and
  alate repellence from *Pseudomonas*-inhabited plants; 
  this folder contains 16 files, `large-main-01.rds` -- `large-main-16.rds`
  (created in `_scripts/02-large-hpc/01-large-main.sh`)
- `large-zeta-ninf-summs.csv`: In larger landscape simulations, number of
  peak infected plants summarized by parasitoid responsiveness and
  number of *Pseudomonas* plants on landscape
  (created in `_scripts/02-large-hpc/04-large-emp-zeta.R`)
- `large-zeta-nopseudo-summs.csv`: In larger landscape simulations, number of
  peak infected plants when there are no *Pseudomonas* plants on landscape
  (created in `_scripts/02-large-hpc/04-large-emp-zeta.R`)
- `large-zeta-summs.csv`: In larger landscape simulations, observed divided
  by predicted (empirically) parasitoid density summarized by parasitoid 
  responsiveness, number of *Pseudomonas* plants on landscape, and 
  whether plant has *Pseudomonas*
  (created in `_scripts/02-large-hpc/04-large-emp-zeta.R`)



- `README.md`: this file

