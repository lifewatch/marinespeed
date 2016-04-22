## TODO grid_filter

#' Filter points based on distance to other points
#'
#' \code{geographic_filter} returns the indexes of all points in \code{data}
#' that are not within \code{bufferdistance} of the \code{filter_data}.
#'
#' @usage geographic_filter(data, filter_data, buffer_distance = 200*1000,
#'   lonlat = TRUE)
#'
#' @param data Matrix or dataframe. The first two columns should represent the
#'   longitude and latitude (or x,y coordinates if \code{lonlat = FALSE}).
#' @param filter_data Matrix or dataframe. The first two columns should
#'   represent the longitude and latitude (or x,y coordinates if \code{lonlat =
#'   FALSE}).
#' @param buffer_distance Positive number. The minimal distance a point in
#'   \code{data} should be from a point in \code{filter_data}.
#' @param lonlat Logical. If \code{TRUE} (default) then Great Circle distances
#'   are calculated else Euclidean (planar) distances are calculated.
#'
#' @return Vector of integer with the indexes of the rows in data that are not
#'   within bufferdistance of the filter_data.
#'
#' @examples
#' set.seed(42)
#' data <- cbind(runif(10, -180, 180), runif(10, -90, 90))
#' filter_data <- cbind(runif(10, -180, 180), runif(10, -90, 90))
#' ## remove points from data data are within a 1000km buffer around the points in filter_data
#' filtered <- geographic_filter(data, filter_data, buffer_distance = 1000*1000, lonlat = TRUE)
#'
#' data_filtered <- data[filtered,]
#' data_removed <- data[-filtered,]
#'
#' @export
geographic_filter <- function(data, filter_data, buffer_distance=200*1000, lonlat = TRUE) {
  if(buffer_distance <= 0) {
    stop("buffer_distance should be positive")
  }
  distfun <- get_distfun(lonlat)
  data <- as.matrix(data[,1:2])
  filter_data <- as.matrix(filter_data[,1:2])
  removed <- c()
  for(i in 1:NROW(filter_data)) {
    if(NROW(data) > 0) {
      dists <- distfun(filter_data[i,,drop = FALSE], data)
      new_removed <- which(dists < buffer_distance)
      if(length(new_removed) > 0) {
        data <- data[-new_removed,]
        removed <- c(removed, new_removed)
      }
    }
  }
  setdiff(1:NROW(data), removed)
}
