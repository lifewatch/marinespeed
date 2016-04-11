# "species.csv.gz"
# filename <- "occurrences_raw.zip"
# "occurrences.zip"
# "pseudodisc_background_5cv_folds.csv.gz"
# "pseudodisc_species_5cv_folds.csv.gz"
# "random_background.csv.gz"
# "random_background_5cv_folds.csv.gz"
# "random_species_5cv_folds.csv.gz"
# "targetgroup_background.csv.gz"
# "targetgroup_background_5cv_folds.csv.gz"


get_version <- function() {
  "v1"
}


get_datadir <- function() {
  datadir <- getOption("marinespeed_datadir")
  if(is.null(datadir)) {
    datadir <- normalizePath(file.path(path.expand("~"), "R", "marinespeed", get_version()))
    if(!dir.exists(datadir)) {
      dir.create(datadir, recursive = TRUE)
    }
  }
  datadir
}

get_file <- function(filename) {
  datadir <- get_datadir()
  outfile <- file.path(datadir, filename)
  outfile_nozip <- file.path(datadir,sub("[.]zip$", "", filename))
  if(!file.exists(outfile) && !dir.exists(outfile_nozip)) {
    root <- paste0("http://www.phycology.ugent.be/research/marinespeed/", get_version(), "/")
    tryCatch({
      download.file(paste0(root, filename), outfile)
    }, error = function(e) { remove.file(outfile)})
    if(grepl("[.]zip$", filename)) {
      unzip(outfile, exdir = outfile_nozip)
      file.remove(outfile)
    }
  }
  outfile_nozip
}

#' List species
#'
#' \code{list_species} returns a dataframe with all the species names and their
#' WoRMS identifier (\code{aphia_id}).
#'
#' @usage list_species()
#'
#' @details If the file with the list of species is not present on the hard disk
#'   then it is downloaded and stored in the data directory. The data directory
#'   can be set with \code{options(marinespeed_datadir = ".")}.
#'
#' @examples
#' species <- list_species()
#' species$species
#' species$aphia_id
#' @export
list_species <- function() {
  read.csv(get_file("species.csv.gz"))
}

#' Get filenames of the occurrences files for some or all species
#'
#' \code{get_occurrence_files} returns a character vector with the full paths to
#' the occurences files for all or a selection of species
#'
#' @usage get_occurrence_file(species = NULL, raw = FALSE)
#'
#' @param species dataframe or character vector. Dataframe like returned by \code{\link{list_species}} or the names of the species. If \code{NULL} then all species occurrence files are returned.
#' @param raw logical. If \code{TRUE} the occurrence records without 25km² grid filtering and outlier removal. Default is \code{FALSE}.
#'
#' @examples
#' \dontrun{
#' all_paths <- get_occurrence_files()
#'
#' species <- list_species()
#' first10 <- get_occurrence_files(species[1:10,])
#'
#' abalistes_stellatus <- get_occurrence_files("Abalistes stellatus")
#' }
#' @export
get_occurrence_files <- function(species = NULL, raw = FALSE) {
  if(raw) {
    dir <- get_file("occurrences_raw.zip")
  } else {
    dir <- get_file("occurrences.zip")
  }
  paths <- list.files(dir, pattern = "[.]csv[.]gz", full.names = TRUE)
  if(is.null(species)) {
    species <- list_species()
  }
  if(NCOL(species) == 2) {
    #paste()
    ## TODO match() or something else ???
  } else {

  }
  stop("TODO implement this")
  paths
}


#' Load folds
#'
#' @details The different supported \code{fold_type} are:
#'
#' \code{"pseudodisc"}: 5-fold disc partitioning of occurrences with pairwise
#' distance sampled and buffer filtered random background points, equivalent to
#' calling \code{\link{kfold_occurrence_background}} with
#' \code{occurrence_fold_type = "disc", k = 5, pwd_sample = TRUE,
#' background_buffer = 200*1000}
#'
#' \code{"random"}:
#'
#' \code{"targetgroup"}:
#' @export
load_folds <- function(fold_type = "pseudodisc") {
  fold_name <- as.character(fold_type)
  if(fold_type == "pseudodisc") {
    bg <- "pseudodisc_background_5cv_folds.csv.gz"
    species <- "pseudodisc_background_5cv_folds.csv.gz"
  } else if(fold_type == "random") {
    bg <- "random_background_5cv_folds.csv.gz"
    species <- "random_species_5cv_folds.csv.gz"
  } else if(fold_type == "targetgroup") {
    bg <- "targetgroup_background_5cv_folds.csv.gz"
    species <- "random_species_5cv_folds.csv.gz"
  } else {
    stop("fold_type not supported")
  }
  list(bg_folds = get_file(bg), species_folds = get_file(species))
}

#' @export
lapply_species <- function(fun, ..., raw = FALSE) {
  fun <- match.fun(fun)
  result <- list()
  species <- list_species()
  for(i in 1:NROW(species)) {
    row <- species[i,]
    file <- paste0(get_occurrences_dir(raw = raw), row$species, " ", row$aphia_id, ".csv.gz")
    if(file.exists(file)) {
      data <- read.csv(file, stringsAsFactors = FALSE)

      result[[i]] <- fun(row, data, ...)
    } else {
      print(paste("OCCURRENCES FILE NOT FOUND", file))
    }
  }
  result
}

#' @export
lapply_kfold_species <- function(fun, ..., kfolds = 1:5, fold_type = "pseudodisc",  raw = FALSE) {
  stop("TODO implement this")
#   fun <- match.fun(fun)
#   result <- list()
#   species <- list_species()
#   for(i in 1:NROW(species)) {
#     row <- species[i,]
#     file <- paste0(get_occurrences_dir(raw = raw), row$species, " ", row$aphia_id, ".csv.gz")
#     if(file.exists(file)) {
#       data <- read.csv(file, stringsAsFactors = FALSE)
#
#       result[[i]] <- fun(row, data, ...)
#     } else {
#       print(paste("OCCURRENCES FILE NOT FOUND", file))
#     }
#   }
#   result
}

# "occurrences_raw.zip"
# "occurrences.zip"
# "random_background.csv.gz"
# "targetgroup_background.csv.gz"

