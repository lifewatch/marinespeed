#' @importFrom geosphere distGeo
get_distfun <- function(lonlat, dismo = FALSE) {
  if(lonlat) {
    if(dismo) {
      dismo_distGeo
    } else {
      geosphere::distGeo
    }
  } else {
    if(dismo) {
      dismo_distPlane
    } else {
      dist_plane
    }
  }
}

lonlat_xyz <- function(data, lonlat_cols) {
  rad_coordinates <- data[,lonlat_cols] * (pi/180)
  lon <- rad_coordinates[,lonlat_cols[1]]
  lat <- rad_coordinates[,lonlat_cols[2]]
  r <- cbind(cos(lat)*cos(lon), cos(lat)*sin(lon), sin(lat))
  colnames(r) <- c("x", "y", "z")
  return(r)
}

dist_plane <- function(a, b) {
  ## Euclidean distance from point to data.frame/matrix of points
  sqrt((a[, 1] - b[, 1])^2 + (a[, 2] - b[, 2])^2)
}

dismo_distHaversine <- function(p1, p2) {
  r <- 6378137
  toRad <- pi/180
  p1 <- p1 * toRad
  p2 <- p2 * toRad
  p <- cbind(p1[, 1], p1[, 2], p2[, 1], p2[, 2])
  dLat <- (p[, 4] - p[, 2])
  dLon <- (p[, 3] - p[, 1])
  a <- sin(dLat/2) * sin(dLat/2) + cos(p[, 2]) * cos(p[,
                                                       4]) * sin(dLon/2) * sin(dLon/2)
  dist <- 2 * atan2(sqrt(a), sqrt(abs(1 - a))) * r
  as.vector(dist)
}
dismo_distGeo <- function(x, y) {
  n <- nrow(x)
  m <- nrow(y)
  dm <- matrix(ncol = m, nrow = n)
  for (i in 1:n) {
    dm[i, ] <- distHaversine(x[i, , drop = FALSE], y)
  }
  return(dm)
}
dismo_distPlane <- function(x, y) {
  n = nrow(x)
  m = nrow(y)
  dm = matrix(ncol = m, nrow = n)
  for (i in 1:n) {
    dm[i, ] = dist_plane(x[i, , drop = FALSE], y)
  }
  return(dm)
}
