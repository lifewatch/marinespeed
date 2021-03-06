lonlat_xyz <- function(data, lonlat_cols) {
  rad_coordinates <- data[,lonlat_cols] * (pi/180)
  lon <- rad_coordinates[,lonlat_cols[1]]
  lat <- rad_coordinates[,lonlat_cols[2]]
  r <- cbind(cos(lat)*cos(lon), cos(lat)*sin(lon), sin(lat))
  colnames(r) <- c("x", "y", "z")
  return(r)
}

dist_plane_point_matrix <- function(a, b) {
  ## Euclidean distance from point to data.frame/matrix of points
  sqrt((a[, 1] - b[, 1])^2 + (a[, 2] - b[, 2])^2)
}

dist_haversine <- function(p1, p2) {
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
dist_geo <- function(x, y) {
  n <- nrow(x)
  m <- nrow(y)
  dm <- matrix(ncol = m, nrow = n)
  for (i in 1:n) {
    dm[i, ] <- dist_haversine(x[i, , drop = FALSE], y)
  }
  return(dm)
}
dist_plane <- function(x, y) {
  n = nrow(x)
  m = nrow(y)
  dm = matrix(ncol = m, nrow = n)
  for (i in 1:n) {
    dm[i, ] = dist_plane_point_matrix(x[i, , drop = FALSE], y)
  }
  return(dm)
}

#' @importFrom geosphere distGeo
get_distfun <- function(lonlat, dismo_like = FALSE) {
  if(lonlat) {
    if(dismo_like) {
      dist_geo
    } else {
      geosphere::distGeo
    }
  } else {
    if(dismo_like) {
      dist_plane
    } else {
      dist_plane_point_matrix
    }
  }
}

mindist <- function(distfun, a, b, lonlat) {
  if(requireNamespace("FNN")) {
    if(lonlat) {
      az <- lonlat_xyz(a,1:2)
      bz <- lonlat_xyz(b,1:2)
      nn <- FNN::get.knnx(bz, az, k=1)
      sapply(1:(nrow(a)), function(i) distfun(a[i,,drop=F], b[nn$nn.index[i],,drop=F]))
    } else {
      nn <- FNN::get.knnx(b, a, k=1)
      nn$nn.dist
    }
  } else if(requireNamespace("dismo")) {
    warning("package FNN not found, using a slower approach")
    partition_count <- (1 %/% (1000 / NROW(b))) + 1
    parts <- dismo::kfold(x=b, k=partition_count)
    r <- c()
    for(i in 1:partition_count) {
      mind <- apply(distfun(a, b[parts==i,]), 1, min)
      r <- cbind(r, mind)
    }
    apply(r, 1, min)
  } else {
    stop("Either FNN (preferably) or dimo should be installed for this to work")
  }
}
