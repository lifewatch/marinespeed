## TODO grid_filter

#' Filter points based on distance to other points
#'
#'
#' @param data
#' @param filter_data
#' @param bufferdistance
#' @param distfun
#'
#' @return Vector of integer with the indexes of the rows in data that are not within bufferdistance of the filter_data.
#'
#' @details
#' @examples
#' # TODO
#'
#' @export
geographic_filter <- function(data, filter_data, bufferdistance=200*1000, distfun = geosphere::distGeo) {
  data <- as.matrix(data[,1:2])
  filter_data <- as.matrix(filter_data[,1:2])
  removed <- c()
  for(i in 1:NROW(filter_data)) {
    dists <- distfun(filter_data[i,], data)
    new_removed <- which(dists < bufferdistance)
    if(length(new_removed) > 0) {
      data <- data[-new_removed,]
      removed <- c(removed, new_removed)
    }
  }
  setdiff(1:NROW(data), removed)
}
