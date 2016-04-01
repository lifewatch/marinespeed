#' @export
ssb <- function(p, a, reference, lonlat = TRUE, avg = TRUE) {
  if (lonlat) {
    distfun <- distGeo
  }
  else {
    distfun <- distPlane
  }
  if (inherits(p, "SpatialPoints"))
    p <- sp::coordinates(p)
  if (inherits(a, "SpatialPoints"))
    a <- sp::coordinates(a)
  if (inherits(reference, "SpatialPoints"))
    reference <- sp::coordinates(reference)
  p <- as.matrix(p)[, 1:2]
  a <- as.matrix(a)[, 1:2]
  reference <- as.matrix(reference)[, 1:2]

  mindist <- function(distfun, a, b) {
    if(require(FNN)) {
      if(lonlat) {
        az <- lonlat_xyz(a,1:2)
        bz <- lonlat_xyz(b,1:2)
        nn <- FNN::get.knnx(bz, az, k=1)
        sapply(1:(nrow(a)), function(i) distfun(a[i,,drop=F], b[nn$nn.index[i],,drop=F]))
      } else {
        nn <- FNN::get.knnx(b, a, k=1)
        nn$nn.dist
      }
    } else {
      warning("package FNN not found, using a slower approach")
      partition_count <- (1 %/% (1000 / NROW(b))) + 1
      parts <- dismo::kfold(x=b, k=partition_count)
      r <- c()
      for(i in 1:partition_count) {
        mind <- apply(distfun(a, b[parts==i,]), 1, min)
        r <- cbind(r, mind)
      }
      apply(r, 1, min)
    }
  }
  pd <- mindist(distfun, p, reference) #   pd <- apply(distfun(p, reference), 1, min)
  ad <- mindist(distfun, a, reference) #   ad <- apply(distfun(a, reference), 1, min)
  if (avg) {
    res <- cbind(mean(pd), mean(ad))
    colnames(res) <- c("p", "a")
    return(res)
  }
  else {
    return(list(pd, ad))
  }
}

#' @export
pwdSample <- function(fixed, sample, reference, tr = 0.33, nearest= TRUE, n=1, lonlat = TRUE, warn = TRUE) {
  if (lonlat) {
    distfun <- distGeo
  }
  else {
    distfun <- distPlane
  }
  stopifnot(tr > 0)
  n <- round(n)
  stopifnot(n >= 1)
  if (inherits(fixed, "SpatialPoints"))
    fixed <- sp::coordinates(fixed)
  if (inherits(sample, "SpatialPoints"))
    sample <- sp::coordinates(sample)
  if (inherits(reference, "SpatialPoints"))
    reference <- sp::coordinates(reference)
  fixed <- as.matrix(fixed)[, 1:2]
  sample <- as.matrix(sample)[, 1:2]
  reference <- as.matrix(reference)[, 1:2]
  if (warn) {
    if (nrow(sample) < nrow(fixed)) {
      warning("nrow(sample) < nrow(fixed)")
    }
  }
  mindist <- function(distfun, a, b) {
    if(require(FNN)) {
      if(lonlat) {
        az <- lonlat_xyz(a,1:2)
        bz <- lonlat_xyz(b,1:2)
        nn <- FNN::get.knnx(bz, az, k=1)
        sapply(1:(nrow(a)), function(i) distfun(a[i,,drop=F], b[nn$nn.index[i],,drop=F]))
      } else {
        nn <- FNN::get.knnx(b, a, k=1)
        nn$nn.dist
      }
    } else {
      warning("package FNN not found, using a slower approach")
      partition_count <- (1 %/% (1000 / NROW(b))) + 1
      parts <- dismo::kfold(x=b, k=partition_count)
      r <- c()
      for(i in 1:partition_count) {
        mind <- apply(distfun(a, b[parts==i,]), 1, min)
        r <- cbind(r, mind)
      }
      apply(r, 1, min)
    }
  }

  fromd <- mindist(distfun, fixed, reference) ##apply(distfun(fixed, reference), 1, min)
  tod <- mindist(distfun, sample, reference) ##apply(distfun(sample, reference), 1, min)
  ngb <- matrix(NA, nrow = length(fromd), ncol = n)
  iter <- sample(1:nrow(fixed))
  for (j in 1:n) {
    for (i in iter) {
      d <- abs(tod - fromd[i])
      if (min(d) < (tr * fromd[i])) {
        if (nearest) {
          x <- which.min(d)
        }
        else {
          x <- sample(which(d < (tr * fromd[i])), size = 1)
        }
        ngb[i, j] <- x
        tod[x] <- Inf
      }
    }
  }
  return(ngb)
}

lonlat_xyz <- function(data, lonlat_cols) {
  rad_coordinates <- data[,lonlat_cols] * (pi/180)
  lon <- rad_coordinates[,lonlat_cols[1]]
  lat <- rad_coordinates[,lonlat_cols[2]]
  r <- cbind(cos(lat)*cos(lon), cos(lat)*sin(lon), sin(lat))
  colnames(r) <- c("x", "y", "z")
  return(r)
}
distHaversine <- function(p1, p2) {
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
distGeo <- function(x, y) {
  n <- nrow(x)
  m <- nrow(y)
  dm <- matrix(ncol = m, nrow = n)
  for (i in 1:n) {
    dm[i, ] <- distHaversine(x[i, , drop = FALSE], y)
  }
  return(dm)
}
distPlane <- function(x, y) {
  dfun <- function(x, y) {
    sqrt((x[, 1] - y[, 1])^2 + (x[, 2] - y[, 2])^2)
  }
  n = nrow(x)
  m = nrow(y)
  dm = matrix(ncol = m, nrow = n)
  for (i in 1:n) {
    dm[i, ] = dfun(x[i, , drop = FALSE], y)
  }
  return(dm)
}
