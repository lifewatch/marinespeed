#' Spatial sorting bias
#'
#' For more information see \code{\link[dismo]{ssb}}
ssb <- function(p, a, reference, lonlat = TRUE, avg = TRUE) {
  distfun <- get_distfun(lonlat, dismo = TRUE)
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
    if(requireNamespace(FNN)) {
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

#' Pair-wise distance sampling
#'
#' For more information see \code{\link[dismo]{pwdSample}}
pwdSample <- function(fixed, sample, reference, tr = 0.33, nearest= TRUE, n=1, lonlat = TRUE, warn = TRUE) {
  distfun <- get_distfun(lonlat, dismo = TRUE)
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
