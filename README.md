

# Alate Effects On Nonpersistent Infection Advancement (aeonia)


<!-- [![DOI](https://zenodo.org/badge/XXXXXXXXX.svg)](https://zenodo.org/badge/latestdoi/XXXXXXXXX) -->


This repo contains the R package `aeonia` that simulates aphid--parasitoid
dynamics and aphid-borne virus transmission.
It also contains scripts (many of which use the R package) that are used to
create figures in the manuscript "Natural enemies can mediate the impact of 
plant microbiota on insect-borne virus transmission".


# Organization

This repo contains the following top-level files/folders:

```
├── _plots
├── _scripts
├── aeonia.Rproj
├── data
├── data-raw
├── DESCRIPTION
├── inst
├── LICENSE.md
├── man
├── NAMESPACE
├── R
├── README.md
└── src
```

Here are descriptions of each:

* `_plots`: Folder used to store figures created using R scripts.
    None are present in this repo, but they can be created using
    the provided scripts.
* `_scripts`: Folder containing R scripts that create the plots used
    in the manuscript. They are ordered by use in the paper, and
    each file contains a description of what it contains.
* `aeonia.Rproj`: RStudio project file for this project.
* `data`: Folder containing data used in the R package.
* `data-raw`: Folder containing a script (`gameofclones.R`) showing how I
  created the objects used in `data`.
* `DESCRIPTION`: Description file for the `aeonia` package.
* `inst`: folder containing necessary files to use the PCG family of random
  number generators in the C++ code for the `aeonia` package.
* `LICENSE.md`: License file for the `aeonia` package.
* `man`: Folder containing the documentation files for the `aeonia` package.
* `NAMESPACE`: File defining imports and export for the `aeonia` package.
* `R`: Folder containing the R files for the `aeonia` package.
* `README.md`: This file, the top-level README.
* `src`: Folder containing the C++ files for the `aeonia` package.






# Replicating R environment

I used R version 4.6.1 (platform: aarch64-apple-darwin23) for all my scripts.

This project uses the `renv` package, so if you want to use this, you must
first install it:

```r
install.packages("renv")
```

Then to install all the packages I used for these analyses, you can
simply run the following while having this project’s main directory as
your working directory:

``` r
renv::restore()
```

If you’d rather avoid `renv`, then you can install all the packages (in
the versions I used) this way:

``` r
if (!requireNamespace("jsonlite", quietly = TRUE)) install.package("jsonlite")
pkg_list <- jsonlite::read_json("renv.lock")[["Packages"]]
pkg_vers <- lapply(pkg_list, \(x) paste0(x[["Package"]], "@", x[["Version"]]))
pkg_vers[["aeonia"]] <- NULL  # bc we're getting this one from github

# Most packages from CRAN:
install.packages(do.call(c, pkg_vers))

# Now aeonia from GitHub:
remotes::install_github(paste0("lucasnell/aeonia@",
                               pkg_list[["aeonia"]][["RemoteSha"]]))
```
