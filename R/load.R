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
#' @seealso \code{\link{lapply_kfold_species}} \code{\link{lapply_species}}
#'   \code{\link{get_fold_data}} \code{\link{get_occurrences}}
#'   \code{\link{species_traits}}
#' @examples
#' species <- list_species()
#' species$species
#' species$aphia_id
#' @export
list_species <- function() {
  csv2rds(get_file("species.csv.gz"))
}

#' Species traits
#'
#' \code{species_traits} returns a dataframe with species traits information.
#'
#' @usage species_traits()
#'
#' @details Traits information includes information about the taxonomy from the
#'   World Register of Marine Species (WoRMS) and habitat traits from the
#'   Encyclopedia of Life (EOL).
#'
#' @seealso \code{\link{list_species}}
#'
#' @examples
#' traits <- species_traits()
#' traits$kingdom
#' @export
species_traits <- function() {
  csv2rds(get_file("traits.csv.gz"))
}

#' Get occurrence records
#'
#' \code{get_occurrences} returns a data.frame with all occurrence records for
#' one or more species, first columns are species, longitude and latitude
#' followed by environmental data columns.
#'
#' @usage get_occurrences(species, raw = FALSE)
#'
#' @param species dataframe or character vector. Dataframe like returned by
#'   \code{\link{list_species}} or the names of the species.
#' @param raw logical. If \code{FALSE} then 25 square kilometer grid and manual
#'   outlier filtered occurrence records are returned.
#'
#' @return Dataframe with as columns: species, longitude, latitude and the
#'   environmental variable columns.
#'
#' @seealso \code{\link{list_species}} \code{\link{get_fold_data}}
#'   \code{\link{get_background}}
#'
#' @examples \dontrun{
#' abalistes_stellatus <- get_occurrences("Abalistes stellatus")
#'
#' species <- list_species()
#' first10 <- get_occurrences(species[1:10,])
#' }
#' @export
get_occurrences <- function(species = NULL, raw = FALSE) {
  if(raw) {
    dir <- get_file("occurrences_raw.zip")
  } else {
    dir <- get_file("occurrences.zip")
  }
  paths <- list.files(dir, pattern = "[.]csv[.]gz", full.names = TRUE)
  if(!is.null(species)) {
    species <- get_species_names(species)

    filenames <- basename(paths)
    filespecies <- sub("\\s?[0-9]*[.]csv[.]gz", "", filenames)
    paths <- paths[filespecies %in% species]

    if(length(paths) == 0) {
      warning("No occurrence files found for the provided species")
    }
  }
  do.call("rbind", lapply(paths, function(p) { csv2rds(p) }))
}

#' Get fold data
#'
#' \code{get_fold_data} returns a list of training and test occurrence and
#' background data fold(s) for one or more species.
#'
#' @usage get_fold_data(species, fold_type, k)
#'
#' @param species dataframe or character vector. Row from the dataframe returned
#'   by \code{\link{list_species}} or the name of the species.
#' @param fold_type character. Type of partitioning you want to use, default is
#'   \code{"disc"}.
#' @param k integer vector. Numbers of the folds you want to get data for, if
#'   you want all 5-folds pass use \code{1:5}, which is the default.
#'
#' @details The different \code{fold_type} are:
#'
#'   \code{"disc"}: 5-fold disc partitioning of occurrences with pairwise
#'   distance sampled and buffer filtered random background points, equivalent
#'   to calling \code{\link{kfold_occurrence_background}} with
#'   \code{occurrence_fold_type = "disc", k = 5, pwd_sample = TRUE,
#'   background_buffer = 200*1000}
#'
#'   \code{"grid_4"} and \code{"grid_9"}: 4-fold and 9-fold grid partitioning of
#'   occurrences with pairwise distance sampled and buffer filtered random
#'   background points, equivalent to calling
#'   \code{\link{kfold_occurrence_background}} with \code{occurrence_fold_type =
#'   "grid", k = 4, pwd_sample = TRUE, background_buffer = 200*1000}
#'
#'   \code{"random"}: 5-fold random partitioning of occurrences and random
#'   background points, equivalent to calling
#'   \code{\link{kfold_occurrence_background}} with \code{occurrence_fold_type =
#'   "random", k = 5, pwd_sample = FALSE, background_buffer = 0}
#'
#'   \code{"targetgroup"}: same way of partitioning as the \code{"random"} folds
#'   but instead of random background points, a random subset of all occurrences
#'   points was used creating a targetgroup background points set which has the
#'   same sampling bias as the entire dataset.
#'
#' @return A 5 element list with fold data filled in for all \code{k}. Fold data
#'   consists of a list with 4 elements: \code{occurrence_training},
#'   \code{occurrence_test}, \code{background_training} and
#'   \code{background_test}.
#'
#' @seealso \code{\link{list_species}} \code{\link{lapply_kfold_species}}
#'   \code{\link{lapply_species}} \code{\link{kfold_data}}
#'
#' @examples \dontrun{
#' aba_folds <- get_fold_data("Abalistes stellatus", "random", k = 1:5)
#' k1 <- aba_folds[[1]]
#' k1$occurrence_training
#' k1$occurrence_test
#' k1$background_training
#' k1$background_test
#' }
#' @export
get_fold_data <- function(species, fold_type, k) {
  if(NROW(species) > 1) {
    ## otherwise might get into memory problems
    stop("get_fold_data expects only 1 species")
  }
  max_k <- 5
  max_k <- ifelse(fold_type == "grid_4", 4, max_k)
  max_k <- ifelse(fold_type == "grid_9", 9, max_k)
  if(min(k) < 1 || max(k) > max_k) {
    stop(paste0("k should be between 1 and ", max_k))
  }
  species <- get_species_names(species)
  occurrences <- get_occurrences(species)
  folds <- get_folds(fold_type)
  if(fold_type == "targetgroup") {
    bg <- get_background("targetgroup")
  } else {
    bg <- get_background("random")
  }
#   specieslist <- list()
#   for (si in 1:length(species)) {
  klist <- list(NULL, NULL, NULL, NULL, NULL)
  for (fold in k) {
    occ_train <- kfold_data(species, occurrences, folds$species, fold, training = TRUE)
    occ_test <- kfold_data(species, occurrences, folds$species, fold, training = FALSE)
    bg_train <- kfold_data(species, bg, folds$background, fold, training = TRUE)
    bg_test <- kfold_data(species, bg, folds$background, fold, training = FALSE)
    bg_train$species <- rep("background", nrow(bg_train))
    bg_test$species <- rep("background", nrow(bg_test))
    klist[[fold]] <- list(occurrence_training=occ_train, occurrence_test=occ_test,
                        background_training=bg_train, background_test=bg_test)
  }
  # }
  klist
}

#' Get MarineSPEED version
#'
#' \code{get_version} returns the currently used MarineSPEED version, this can
#' be changed by setting \code{options(marinespeed_version="<version
#' information>")}.
#'
#' @usage get_version()
#'
#' @return Character with the current version ("V1") or another version if the
#'   \code{marinespeed_version} has been set.
#'
#' @examples
#' print(get_version())
#'
#' @export
get_version <- function() {
  v <- getOption("marinespeed_version")
  if(is.null(v)) {
    v <- "v1"
  }
  v
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
    root <- paste0("http://marinespeed.samuelbosch.com/", get_version(), "/")
    tryCatch({
      download.file(paste0(root, filename), outfile)
    }, error = function(e) { file.remove(outfile)})
    if(grepl("[.]zip$", filename)) {
      unzip(outfile, exdir = outfile_nozip)
      file.remove(outfile)
    }
  }
  normalizePath(outfile_nozip)
}

get_species_names <- function(species) {
  if(NCOL(species) == 2 && colnames(species) == c("species", "aphia_id")) {
    species <- species[,"species"]
  }
  if(!is.character(species) && !is.factor(species)){
    stop("invalid species input")
  }
  as.character(species)
}

csv2rds <- function(file, extension = ".rds") {
  rds_file <- sub("[.]csv[.]gz$", extension, file)
  data <- NULL
  if(file.exists(rds_file)) { ## cache as rds (faster)
    data <- readRDS(rds_file)
  }
  if(extension == ".rds" && (is.null(data) || getOption("stringsAsFactors") != is.factor(data[1,1]))) {
    data <- read.csv(file)
    saveRDS(data, rds_file)
  } else if (is.null(data) && grepl("cv_folds_bit[.]rds$", rds_file)) {
    folds <- read.csv(file)
    data <- structure(list(), class = "marinespeed_folds")

    species <- folds[,1]
    for(sp in unique(species)) {
      w <- as.which(species == sp)
      data[["species"]][[sp]] <- c(min(w), max(w))
    }
    for(ki in 2:ncol(folds)) {
      if(any(is.na(folds[,ki]))) {
        data[[paste0(colnames(folds)[ki], "_NOTNA")]] <- !as.bit(is.na(folds[,ki]))
      }
      data[[colnames(folds)[ki]]] <- as.bit(folds[,ki])
    }
    saveRDS(data, rds_file)
  } else if(is.null(data)) {
    stop("unsupported")
  }
  data
}

#' Get folds
#'
#' \code{get_folds} returns the different pre-generated folds information. To
#' get the fold data for a species see also \code{\link{get_fold_data}}.
#'
#' @usage get_folds(type = "disc")
#'
#' @param type character. The type of partitioning you want to load.
#'
#' @details The different supported \code{type} are:
#'
#'   \code{"disc"}: 5-fold disc partitioning of occurrences with pairwise
#'   distance sampled and buffer filtered random background points, equivalent
#'   to calling \code{\link{kfold_occurrence_background}} with
#'   \code{occurrence_fold_type = "disc", k = 5, pwd_sample = TRUE,
#'   background_buffer = 200*1000}
#'
#'   \code{"grid_4"} and \code{"grid_9"}: 4-fold and 9-fold grid partitioning of
#'   occurrences with pairwise distance sampled and buffer filtered random
#'   background points, equivalent to calling
#'   \code{\link{kfold_occurrence_background}} with \code{occurrence_fold_type =
#'   "grid", k = 4, pwd_sample = TRUE, background_buffer = 200*1000}
#'
#'   \code{"random"}: 5-fold random partitioning of occurrences and random
#'   background points, equivalent to calling
#'   \code{\link{kfold_occurrence_background}} with \code{occurrence_fold_type =
#'   "random", k = 5, pwd_sample = FALSE, background_buffer = 0}
#'
#'   \code{"targetgroup"}: same way of partitioning as the \code{"random"} folds
#'   but instead of random background points, a random subset of all occurrences
#'   points was used creating a targetgroup background points set which has the
#'   same sampling bias as the entire dataset.
#'
#' @return A list with two entries \code{"background"} and \code{"species"},
#'   each entry is a dataframe with species name column and 5 fold columns as
#'   created by \code{\link{kfold_occurrence_background}}
#'
#' @seealso \code{\link{lapply_kfold_species}} \code{\link{get_fold_data}}
#'   \code{\link{get_occurrences}} \code{\link{get_background}}
#'   \code{\link{kfold_data}}
#'
#' @examples \dontrun{
#' folds <- get_folds("random")
#'
#' abalistes <- "Abalistes stellatus"
#' occ <- get_occurrences(abalistes)
#' bg <- get_background("random")
#'
#' occ_train <- kfold_data(abalistes, occ, folds$species, k=1, training=TRUE)
#' occ_test <- kfold_data(abalistes, occ, folds$species, k=1, training=FALSE)
#' bg_train <- kfold_data(abalistes, bg, folds$background, k=1, training=TRUE)
#' bg_test <- kfold_data(abalistes, bg, folds$background, k=1, training=FALSE)
#' }
#' @export
get_folds <- function(type = "disc") {
  type <- as.character(type)
  if(type == "disc") {
    bg <- "pseudodisc_background_5cv_folds"
    species <- "pseudodisc_species_5cv_folds"
  } else if(type == "grid_4") {
    bg <- "grid_background_4cv_folds"
    species <- "grid_species_4cv_folds"
  } else if(type == "grid_9") {
    bg <- "grid_background_9cv_folds"
    species <- "grid_species_9cv_folds"
  } else if(type == "random") {
    bg <- "random_background_5cv_folds"
    species <- "random_species_5cv_folds"
  } else if(type == "targetgroup") {
    bg <- "targetgroup_background_5cv_folds"
    species <- "random_species_5cv_folds"
  } else {
    stop("fold_type not supported")
  }
  extension <- getOption("marinespeed_folds_extension")
  extension <- ifelse(is.null(extension), "_bit.rds", extension)
  bg <- get_file(paste0(bg, extension))
  species <- get_file(paste0(species, extension))
  if(extension == ".csv.gz") {
    bg_folds <- csv2rds(bg)
    species_folds <- csv2rds(bg)
  } else {
    bg_folds <- readRDS(bg)
    species_folds <- readRDS(species)
  }

  list(background = bg_folds, species = species_folds)
}

#' Get background data
#'
#' \code{get_background} returns pre-generated background data.
#'
#' @usage get_background(type)
#'
#' @param type character. Either \code{"random"} or \code{"targetgroup"}.
#'
#' @details The targetgroup background was created by subsampling an average of
#'   37 occurrence records (20000 in total) from each species in the dataset
#'   providing in essence the same sampling bias as the entire dataset.
#'
#' @return A dataframe with all background points.
#'
#' @seealso \code{\link{get_fold_data}} \code{\link{get_occurrences}}
#'   \code{\link{get_folds}}
#'
#' @examples \dontrun{
#' random_bg <- get_background("random")
#' plot(random_bg[,2:3], pch=".", main="random background")
#'
#' targetgroup_bg <- get_background("targetgroup")
#' plot(targetgroup_bg[,2:3], pch=".", main="targetgroup background")
#' }
#' @export
get_background <- function(type) {
  if(!(as.character(type) %in% c("random", "targetgroup"))) {
    stop("background type not recognized")
  }
  read.csv(get_file(paste0(as.character(type), "_background.csv.gz")))
}
