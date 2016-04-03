#' Create k folds of species and background data for cross-validation
#'
#' \code{kfold_species_background} creates a k-fold partitioning of species and background data for cross-validation using random and stratified folds.
#' Returns a list with the species folds and the background folds, folds are represented as TRUE/FALSE/NA columns of a dataframe, 1 column for each fold.
#'
#' @usage
#' kfold_species_background(species_data, background_data, fold_type = "disc", k = 5, pwd_sample = TRUE, background_buffer = 200*1000, distfun = geosphere::distGeo, seed = 42)
#'
#' @param species_data matrix or dataframe. The first two columns should represent the longitude and latitude (or x,y coordinates if distfun is provided).
#' @param background_data matrix or dataframe. The first two columns should represent the longitude and latitude (or x,y coordinates if distfun is provided).
#' @param species_fold_type character. How species folds should be generated, currently \code{"disc"} (see \code{\link{kfold_disc}}) and \code{"random"} are supported.
#' @param k integer. The number of folds (partitions) that have to be created. By default 5 folds are created.
#' @param pwd_sample logical. Whether backgound points should be picked by doing pair-wise distance sampling (see \code{\link[dismo]{pwdSample}}). It is recommended to install the FNN package if you want to do pair-wise distance sampling.
#' @param background_buffer numeric. Distance in meters around species test points where background data should be excluded from.
#' @param distfun function. Distance function to calculate distances between a point and an other matrix or dataframe of points. By default the \code{geosphere::distGeo} function is used.
#' @param seed a single value, interpreted as an integer or NULL (see \code{\link{base::set.seed}}).
#'
#' @return A list with 2 dataframes, species and background, each with k columns containing \code{TRUE}, \code{FALSE} or \code{NA}.
#'
#' @details
#' Note that which and how many background points get selected in each fold depends on the \code{fold_type}, \code{pwd_sample} and the \code{background_buffer} and
#' whether \code{pwd_sample} is \code{TRUE} or \code{FALSE}, even leading in some cases to the selection of no background data.
#' Background points that are neither selected for the training fold nor for the test fold are set to \code{NA} in the background folds.
#' Random assignment of background points to the folds can be achieved by setting \code{pwd_sample} to \code{FALSE} and \code{background_buffer} to 0.
#' Note also that when \code{pwd_sample} is \code{TRUE}, the same background point might be assigned to different folds.
#'
#' @references
#' Hijmans, R. J. (2012). Cross-validation of species distribution models: removing spatial sorting bias and calibration with a null model. Ecology, 93(3), 679–688. doi:10.1890/11-0826.1
#' Radosavljevic, A., & Anderson, R. P. (2013). Making better Maxent models of species distributions: complexity, overfitting and evaluation. Journal of Biogeography. doi:10.1111/jbi.12227
#'
#' @seealso \code{\link{kfold_disc}} \code{\link{geographic_filter}} \code{\link[dismo]{pwdSample}}
#' @examples
#'
#' ## TODO
#'
#' @export
kfold_species_background <- function(species_data, background_data, species_fold_type = "disc", k = 5, pwd_sample = TRUE, lonlat = TRUE, background_buffer = 200*1000, seed = 42) {
  set.seed(seed)
  #' 1) partition species presences with pseudo-discs (1st fold = real disc, other folds = )
  #' 2) select testing background points with pairwise distance sampling (pwdSample) to reduce spatial sorting bias
  #' 3) select training background points by filtering background points that are within background_buffer distance of the points in the training fold
  #' 4) make sure training and testing background are different
  #' 5) set NA background points that are not in the training and not in the testing set
  if(fold_type == "disc") {
    species_partitions <- kfold_disc(species_data, k, lonlat, seed = seed)
  } else if (fold_type == "random") {
    set.seed(seed)
    species_partitions <- dismo::kfold(species_data, k)
  } else {
    stop("Unknown fold type")
  }
  if(!pwd_sample) {
    background_partitions <- dismo::kfold(background_data, k)
  }
  species_folds <- data.frame(row.names = 1:NROW(species_data))
  background_folds <- data.frame(row.names = 1:NROW(background_data))
  for(ki in 1:k) {
    ## training folds
    species_folds[,paste0("k",ki)] <- species_partitions != ki

    species_test <- species_data[species_partitions == ki, 1:2]
    species_training <- species_data[species_partitions != ki, 1:2]

    if(pwd_sample) {
      test_sample <- pwdSample(species_test, background_data[,1:2], species_training, n=5) ## try to get 5 background points for each testing presence point, you'll get less than that
    } else {
      test_sample <- which(background_partitions==k) ## use randomly generated partitions
    }
    background_test_i <- unique(na.omit(as.vector(test_sample)))
    background_training_i <- base::setdiff(background_training_i, background_test_i)

    if(background_buffer > 0) {
      background_training_i <- geographic_filter(background_data[background_training_i,], species_test, lonlat, background_buffer)
    }

    na_i <- setdiff(1:NROW(background_data), c(background_training_i,background_test_i))

    ## "training" folds
    background_folds[,paste0("k",ki)] <- (1:NROW(background_data) %in% background_training_i)
    background_folds[na_i,paste0("k",ki)] <- NA
  }
  list(species=species_folds, background=background_folds)
}


#' Create k disc based folds for cross-validation
#'
#' \code{kfold_disc} creates a k-fold partitioning of geographical data for cross-validation based on the distance between points. The n points nearest to a selected point are put into a group.
#' Returns a vector with fold numbers ranging from 1 to k.
#'
#' @usage
#' kfold_disc(data, k = 5, distfun = geosphere::distGeo, seed = 42)
#'
#' @param data dataframe or matrix. First two columns should contain the longitude and latitude values (or x,y if a different distance function is supplied).
#' @param k integer. The number of folds (partitions) that have to be created. By default 5 folds are created.
#' @param lonlat logical. TODO: document this
#' @param seed a single value, interpreted as an integer or NULL (see \code{\link{base::set.seed}}).
#'
#' @return A vector with fold numbers ranging from 1 to k.
#'
#' @examples
#' set.seed(42)
#' lonlat_data <- cbind(runif(11, -180, 180), runif(11, -90, 90))
#' folds <- kfold_disc(lonlat_data, k = 5, seed = 1)
#' plot_folds(lonlat_data, folds)
#'
#' # use the euclidean distance
#' xy_data <- cbind(runif(11, 0, 100), runif(11, 0, 100))
#' folds <- kfold_disc(xy_data, k = 5, lonlat = FALSE)
#' plot_folds(xy_data, folds)
#'
#' @export
kfold_disc <- function(data, k = 5, lonlat = TRUE, seed = 42) {
  set.seed(seed)
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

#' @export
plot_folds <- function(data, folds) {
  plot(data, pch=".")
  k <- max(folds)
  cols <- rainbow(k)
  for(i in 1:k) {
    text(data[folds==i,], labels = i, pch=20, col=cols[i])
  }
}
