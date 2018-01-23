[![GPLv3](https://img.shields.io/badge/license-GPLv3-blue.svg)](https://github.com/lifewatch/marinespeed/blob/master/LICENSE.md)
[![Build Status](https://travis-ci.org/lifewatch/marinespeed.svg?branch=master)](https://travis-ci.org/lifewatch/marinespeed)
[![Coverage Status](http://codecov.io/github/lifewatch/marinespeed/coverage.svg?branch=master)](http://codecov.io/github/lifewatch/marinespeed?branch=master)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/marinespeed)](https://CRAN.R-project.org/package=marinespeed)

# MarineSPEED
R code for downloading and working with the Marine SPEcies and Environmental Data ([MarineSPEED.org](http://MarineSPEED.org)) benchmark dataset.

Installation:

```R
    install.packages("marinespeed")
    # or 
    devtools::install_github("lifewatch/marinespeed")
```

Example usage:

```R
    library(marinespeed)
    
    ## get species list
    species <- list_species()
    View(species)
    
    ## count number of occurrences for all species 
    get_occ_count <- function(speciesname, occ) {
      nrow(occ)
    }
    record_counts <- lapply_species(get_occ_count)
    print(sum(unlist(record_counts)))
    
    ## plot first 2 folds for the first 10 species
    plot_occurrences <- function(speciesname, data, k) {
       title <- paste0(speciesname, " (fold = ", k, ")")
       plot(data$occurrence_train[,c("longitude", "latitude")], pch=".", 
            col="blue", main = title)
       points(data$occurrence_test[,c("longitude", "latitude")], pch=".", 
              col="red")
    }
    
    # plot training (blue) and test (red) occurrences of the first 2 disc folds 
    # for the first 10 species
    species <- list_species()
    lapply_kfold_species(plot_occurrences, species=species[1:10,],
                         fold_type = "disc", k = 1:2)
```
