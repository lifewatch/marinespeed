#' Apply a function over a set of species
#'
#' \code{lapply_species} returns a list where each element is the result of
#' applying \code{fun} to all species or the provided subset of species.
#'
#' @usage lapply_species(fun, ..., species = NULL, raw = FALSE)
#'
#' @param fun function. The function to be applied to the occurrence records of
#'   each species. Parameters are the species name and a dataframe with the
#'   occurrence records.
#' @param ... optional arguments to \code{fun}.
#' @param species dataframe or character vector. Dataframe like returned by
#'   \code{\link{list_species}} or the names of the species. If \code{NULL}
#'   (default) then \code{fun} is applied for all species.
#' @param raw logical. If \code{FALSE} then 25 square kilometer grid and manual
#'   outlier filtered occurrence records are returned.
#'
#' @details The parameters passed to \code{fun} are \code{speciesname} and
#'   \code{data}, which is a dataframe with the occurrence records and their
#'   environmental data.
#'
#' @return A list with one named entry for every species provided or for all
#'   species.
#'
#' @seealso \code{\link{list_species}} \code{\link{lapply_kfold_species}}
#'   \code{\link{get_occurrences}}
#'
#' @examples \dontrun{
#' get_occ_count <- function(speciesname, occ) {
#'  nrow(occ)
#' }
#'
#' record_counts <- lapply_species(get_occ_count)
#' sum(unlist(record_counts))
#'
#' # count first 10
#' species <- list_species()
#' lapply_species(get_occ_count, species=species[1:10,])
#' }
#' @export
#' @seealso \code{\link{lapply_kfold_species}}, \code{\link{get_occurrences}},
#'   \code{\link{list_species}}
lapply_species <- function(fun, ..., species = NULL, raw = FALSE) {
  fun <- match.fun(fun)
  result <- list()
  if(is.null(species)) {
    species <- list_species()
  }
  species <- get_species_names(species)
  for(i in 1:NROW(species)) {
    speciesname <- species[i]
    data <- get_occurrences(speciesname, raw = raw)
    result[[speciesname]] <- fun(speciesname, data, ...)
  }
  result
}

#' Apply a function over the folds of a set of species
#'
#' \code{lapply_kfold_species} returns a list of lists where each element is the
#' result of applying \code{fun} to all species or the provided subset of
#' species for the specified folds.
#'
#' @usage lapply_kfold_species(fun, ..., species = NULL, fold_type = "disc", k =
#'   1:5)
#'
#' @param fun function. The function to be applied to the occurrence records of
#'   each species. Parameters are the species name, a list with the occurrence
#'   and background training and test records and a fold number.
#' @param ... optional arguments to \code{fun}.
#' @param species dataframe or character vector. Dataframe like returned by
#'   \code{\link{list_species}} or the names of the species. If \code{NULL}
#'   (default) then \code{fun} is applied for all species.
#' @param fold_type character. Type of partitioning you want to use, default is
#'   \code{"disc"}.
#' @param k integer vector. Numbers of the folds you want to get data for, if
#'   you want all 5-folds pass use \code{1:5}, which is the default.
#'
#' @details The parameters passed to \code{fun} are \code{speciesname},
#'   \code{data} where \code{data} is a list with 4 elements
#'   (\code{occurrence_training}, \code{occurrence_test},
#'   \code{background_training} and \code{background_test}) and a parameter
#'   \code{fold} which contains the fold number.
#'
#'   The different \code{fold_type} are:
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
#' @return  A list with one named entry for every species provided or for all
#'   species. Every list entry is a list with \code{k} as names and the result
#'   of \code{fun} as value.
#'
#' @seealso \code{\link{list_species}} \code{\link{lapply_species}}
#'   \code{\link{get_fold_data}}
#'
#' @examples \dontrun{
#' plot_occurrences <- function(speciesname, data, fold) {
#'    title <- paste0(speciesname, " (fold = ", fold, ")")
#'    plot(data$occurrence_train[,c("longitude", "latitude")], pch=".",
#'         col="blue", main = title)
#'    points(data$occurrence_test[,c("longitude", "latitude")], pch=".",
#'         col="red")
#' }
#'
#' # plot training (blue) and test (red) occurrences
#' # of the first 2 folds for the first 10 species
#' species <- list_species()
#' lapply_kfold_species(plot_occurrences, species=species[1:5,],
#'                      fold_type = "disc", k = 1:2)
#' }
#' @export
#' @seealso \code{\link{lapply_species}}, \code{\link{get_fold_data}},
#'   \code{\link{list_species}}
lapply_kfold_species <- function(fun, ..., species = NULL, fold_type = "disc", k = 1:5) {
  fun <- match.fun(fun)
  result <- list()
  if(is.null(species)) {
    species <- list_species()
  }
  species <- get_species_names(species)

  folds <- get_folds(fold_type)
  if(fold_type == "targetgroup") {
    bg <- get_background("targetgroup")
  } else {
    bg <- get_background("random")
  }

  for(i in 1:NROW(species)) {
    speciesname <- species[i]
    occurrences <- get_occurrences(speciesname)

    klist <- list(NULL, NULL, NULL, NULL, NULL)
    for (fold in k) {
      occ_train <- kfold_data(speciesname, occurrences, folds$species, fold, training = TRUE)
      occ_test <- kfold_data(speciesname, occurrences, folds$species, fold, training = FALSE)
      bg_train <- kfold_data(speciesname, bg, folds$background, fold, training = TRUE)
      bg_test <- kfold_data(speciesname, bg, folds$background, fold, training = FALSE)
      bg_train$species <- rep("background", nrow(bg_train))
      bg_test$species <- rep("background", nrow(bg_test))
      data <- list(occurrence_training=occ_train, occurrence_test=occ_test,
                   background_training=bg_train, background_test=bg_test)
      klist[[fold]] <- fun(speciesname, data, fold, ...)
    }
    result[[speciesname]] <- klist
  }
  result
}
