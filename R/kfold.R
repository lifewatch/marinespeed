#' Get the data associated with a fold
#'
#' \code{kfold_data} returns rows from data for a specific species for the
#' kth fold
#'
#' @usage kfold_data(species_name, data, folds, k, training)
#'
#' @param species_name Character vector. The name of the species you want get
#'   the fold data for.
#' @param data Dataframe. The occurrence or background data you want to get the
#'   fold data for.
#' @param folds Dataframe. The occurrence or background folds as created by
#'   \code{kfold_species_background}. Essentially a data.frame with a species
#'   column and the k-fold cross-validation logical vectors.
#' @param k Integer. Indicate which kth fold you want to return data for.
#' @param training Logical. If \code{TRUE} then training data is returned else
#'   if \code{FALSE} then test data is returned.
#'
#' @return Filtered version of data.
#'
#' @examples
#' set.seed(42)
#' occurrence_data <- data.frame(species = rep("Abalistes stellatus", 50), longitude = runif(50, -180, 180), latitude = runif(50, -90, 90))
#' background_data <- data.frame(species = rep("background", 10000), longitude = runif(1000, -180, 180), latitude = runif(1000, -90, 90))
#' folds <- kfold_occurrence_background(occurrence_data, background_data)
#'
#' ## TODO replace above with download of the MarineSPEED occurrence data and fold data
#'
#' ## get training and test data for the first fold
#' occ_training <- kfold_data("Abalistes stellatus", occurrence_data, folds, k = 1, training = TRUE)
#' occ_test <- kfold_data("Abalistes stellatus", occurrence_data, folds, k = 1, training = FALSE)
#'
#' background_training <- kfold_data("Abalistes stellatus", background_data, folds, k = 1, training = TRUE)
#' background_test <- kfold_data("Abalistes stellatus", background_data, folds, k = 1, training = FALSE)
kfold_data <- function(species_name, data, folds, k, training) {
  ## TODO refactor/change to make it easier/less typing for the user (higher level API)
  # TODO export ? @export
  kcol <- paste0("k",k)
  if(!(kcol %in% colnames(folds))) {
    stop("k not found in folds")
  }
  filter <- c(as.character(species_name), "background")
  f <- folds[as.character(folds$species) %in% filter, kcol] == training
  d <- data[f & !is.na(f),] ## handle NA's (from e.g. pseudo-disc background)
  d
}

#' Create k folds of occurrence and background data for cross-validation
#'
#' \code{kfold_occurrence_background} creates a k-fold partitioning of
#' occurrence and background data for cross-validation using random and
#' stratified folds. Returns a list with the occurrence folds and the background
#' folds, folds are represented as TRUE/FALSE/NA columns of a dataframe, 1
#' column for each fold.
#'
#' @usage kfold_occurrence_background(occurrence_data, background_data,
#' fold_type = "disc", k = 5, pwd_sample = TRUE, background_buffer = 200*1000,
#' distfun = geosphere::distGeo)
#'
#' @param occurrence_data Dataframe. Occurrence points of the species, the first
#'   column should be the scientific name of the species followed by two columns
#'   representing the longitude and latitude (or x,y coordinates if \code{lonlat
#'   = FALSE}).
#' @param background_data Dataframe. Background data points, the first column is
#'   a dummy column followed by two columns representing the longitude and
#'   latitude (or x,y coordinates if \code{lonlat = FALSE}).
#' @param occurrence_fold_type Character vector. How occurrence folds should be
#'   generated, currently \code{"disc"} (see \code{\link{kfold_disc}}) and
#'   \code{"random"} are supported.
#' @param k Integer. The number of folds (partitions) that have to be created.
#'   By default 5 folds are created.
#' @param pwd_sample Logical. Whether backgound points should be picked by doing
#'   pair-wise distance sampling (see \code{\link[dismo]{pwdSample}}). It is
#'   recommended to install the FNN package if you want to do pair-wise distance
#'   sampling.
#' @param lonlat Logical. If \code{TRUE} (default) then Great Circle distances
#'   are calculated else if \code{FALSE} Euclidean (planar) distances are
#'   calculated.
#' @param background_buffer Numeric. Distance in meters around species test
#'   points where background data should be excluded from.
#'
#' @return A list with 2 dataframes, \code{occurrence} and \code{background},
#'   with as first column the scientifc name or \code{"background"} and k
#'   columns containing \code{TRUE}, \code{FALSE} or \code{NA}.
#'
#' @details Note that which and how many background points get selected in each
#' fold depends on the \code{fold_type}, \code{pwd_sample} and the
#' \code{background_buffer} and whether \code{pwd_sample} is \code{TRUE} or
#' \code{FALSE}, even leading in some cases to the selection of no background
#' data. Background points that are neither selected for the training fold nor
#' for the test fold are set to \code{NA} in the background folds. Random
#' assignment of background points to the folds can be achieved by setting
#' \code{pwd_sample} to \code{FALSE} and \code{background_buffer} to 0. Note
#' also that when \code{pwd_sample} is \code{TRUE}, the same background point
#' might be assigned to different folds.
#'
#' @references Hijmans, R. J. (2012). Cross-validation of species distribution
#' models: removing spatial sorting bias and calibration with a null model.
#' Ecology, 93(3), 679–688. doi:10.1890/11-0826.1 Radosavljevic, A., & Anderson,
#' R. P. (2013). Making better Maxent models of species distributions:
#' complexity, overfitting and evaluation. Journal of Biogeography.
#' doi:10.1111/jbi.12227
#'
#' @seealso \code{\link{kfold_disc}} \code{\link{geographic_filter}}
#'   \code{\link[dismo]{pwdSample}}
#' @examples
#' set.seed(42)
#' occurrence_data <- data.frame(species = rep("Abalistes stellatus", 50), longitude = runif(50, -180, 180), latitude = runif(50, -90, 90))
#' background_data <- data.frame(species = rep("background", 10000), longitude = runif(1000, -180, 180), latitude = runif(1000, -90, 90))
#' disc_folds <- kfold_occurrence_background(occurrence_data, background_data, "disc")
#' random_folds <- kfold_occurrence_background(occurrence_data, background_data, "random", pwd_sample = FALSE, background_buffer = 0)
#'
#' @export
kfold_occurrence_background <- function(occurrence_data, background_data, occurrence_fold_type = "disc", k = 5, pwd_sample = TRUE, lonlat = TRUE, background_buffer = 200*1000) {
  # 1) partition species presences with pseudo-discs (1st fold = real disc, other folds = )
  # 2) select testing background points with pairwise distance sampling (pwdSample) to reduce spatial sorting bias
  # 3) select training background points by filtering background points that are within background_buffer distance of the points in the training fold
  # 4) make sure training and testing background are different
  # 5) set NA background points that are not in the training and not in the testing set
  if(!pwd_sample) {
    background_partitions <- dismo::kfold(background_data, k)
  }
  if(occurrence_fold_type == "disc") {
    occurrence_partitions <- kfold_disc(occurrence_data[,2:3], k, lonlat)
  } else if (occurrence_fold_type == "random") {
    occurrence_partitions <- dismo::kfold(occurrence_data, k)
  } else {
    stop("Unknown fold type")
  }
  occurrence_folds <- data.frame(species = occurrence_data[,1], stringsAsFactors = is.factor(occurrence_data[,1])) ## first column with scientific name
  background_folds <- data.frame(species = rep("background", NROW(background_data)), stringsAsFactors = is.factor(occurrence_data[,1]))
  for(ki in 1:k) {
    ## training folds
    occurrence_folds[,paste0("k",ki)] <- occurrence_partitions != ki

    occurrence_test <- occurrence_data[occurrence_partitions == ki, 2:3]
    occurrence_training <- occurrence_data[occurrence_partitions != ki, 2:3]

    if(pwd_sample) {
      test_sample <- pwdSample(occurrence_test, background_data[,2:3], occurrence_training, n=5) ## try to get 5 background points for each testing presence point, you'll get less than that
    } else {
      test_sample <- which(background_partitions==ki) ## use randomly generated partitions
    }
    background_test_i <- unique(na.omit(as.vector(test_sample)))
    background_training_i <- base::setdiff(1:NROW(background_data), background_test_i)

    if(background_buffer > 0) {
      filtered <- geographic_filter(background_data[background_training_i,2:3], occurrence_test, lonlat, background_buffer)
      background_training_i <- (1:NROW(background_data))[background_training_i][filtered]
    }

    na_i <- setdiff(1:NROW(background_data), c(background_training_i,background_test_i))

    ## "training" folds
    background_folds[,paste0("k",ki)] <- (1:NROW(background_data) %in% background_training_i)
    background_folds[na_i,paste0("k",ki)] <- NA
  }
  list(occurrence=occurrence_folds, background=background_folds)
}

#' Create k disc based folds for cross-validation
#'
#' \code{kfold_disc} creates a k-fold partitioning of geographical data for
#' cross-validation based on the distance between points. The n points nearest
#' to a selected point are put into a group. Returns a vector with fold numbers
#' ranging from 1 to k.
#'
#' @usage kfold_disc(data, k = 5, distfun = geosphere::distGeo)
#'
#' @param data Matrix or dataframe. The first two columns should represent the
#'   longitude and latitude (or x,y coordinates if \code{lonlat = FALSE}).
#' @param k Integer. The number of folds (partitions) that have to be created.
#'   By default 5 folds are created.
#' @param lonlat Logical. If \code{TRUE} (default) then Great Circle distances
#'   are calculated else if \code{FALSE} Euclidean (planar) distances are
#'   calculated.
#'
#' @return A vector with fold numbers ranging from 1 to k.
#'
#' @examples
#' set.seed(42)
#' lonlat_data <- cbind(runif(11, -180, 180), runif(11, -90, 90))
#' folds <- kfold_disc(lonlat_data, k = 5)
#' plot_folds(lonlat_data, folds)
#'
#' # use the euclidean distance
#' xy_data <- cbind(runif(11, 0, 100), runif(11, 0, 100))
#' folds <- kfold_disc(xy_data, k = 5, lonlat = FALSE)
#' plot_folds(xy_data, folds)
#'
#' @export
kfold_disc <- function(data, k = 5, lonlat = TRUE) {
  distfun <- get_distfun(lonlat)
  k <- as.integer(k)
  if(is.na(k) || k < 1) {
    stop("k should at least be 1")
  } else if(k > (NROW(data)/2)) {
    stop("k should be less then or equal to half the number of rows in data")
  }
  d <- as.data.frame(data[,1:2])
  d[,3] <- 1:NROW(d)
  colnames(d) <- c("x", "y", "index")

  npoints <- floor(NROW(d) / k)
  remainder <- NROW(d) - (npoints * k)

  partitions <- integer(NROW(d))
  itest <- .Machine$integer.max ## intialize to max value to prevent rows being removed from d in the 1st fold

  for(ki in 1:k) {
    if(ki == 1) {
      i <- sample(1:NROW(d), 1)
    } else {
      i <- sample(order(dists, decreasing = TRUE)[1:npoints], 1)
    }
    center <- d[i,1:2] ## next center: 1 of the n points furthest away from the previous cluster
    d <- d[-itest,]
    dists <- distfun(center, d[, 1:2])
    stopifnot(all(!is.na(dists)))
    if(remainder >= ki) { ## balance the fold lengths
      itest <- order(dists, decreasing = FALSE)[1:(npoints+1)]
    } else {
      itest <- order(dists, decreasing = FALSE)[1:npoints]
    }
    partitions[d[itest,"index"]] <- ki
  }
  partitions
}

#' plot folds
#'
#' \code{plot_folds} makes a rudimentary plot of the data and the folds created
#' with e.g. \code{\link{kfold_disc}} or \code{\link[dismo]{kfold}}.
#'
#' @usage plot_folds(data, folds)
#'
#' @param data Matrix or dataframe. Data for which the folds where created. The
#'   first two columns should represent the longitude and latitude (or x,y
#'   coordinates).
#' @param folds NUmeric vector with group assignments from e.g.
#'   \code{\link{kfold_disc}} or \code{\link[dismo]{kfold}}.
#'
#' @examples
#' set.seed(42)
#' lonlat_data <- cbind(runif(11, -180, 180), runif(11, -90, 90))
#' folds <- kfold_disc(lonlat_data, k = 5)
#' plot_folds(lonlat_data, folds)
#'
#' @export
plot_folds <- function(data, folds) {
  plot(data, pch=".")
  k <- max(folds)
  cols <- rainbow(k)
  for(i in 1:k) {
    text(data[folds==i,], labels = i, pch=20, col=cols[i])
  }
}
