#' @importFrom sp coordinates
get_coordinates <- function(x) {
  if (inherits(x, "SpatialPoints")) {
    x <- sp::coordinates(x)
  }
  as.matrix(x)[, 1:2]
}

#' @importFrom sp coordinates
ssb <- function(p, a, reference, lonlat = TRUE, avg = TRUE) {
  distfun <- get_distfun(lonlat, dismo_like = TRUE)
  p <- get_coordinates(p)
  a <- get_coordinates(a)
  reference <- get_coordinates(a)

  pd <- mindist(distfun, p, reference, lonlat)
  ad <- mindist(distfun, a, reference, lonlat)
  if (avg) {
    res <- cbind(mean(pd), mean(ad))
    colnames(res) <- c("p", "a")
    return(res)
  }
  else {
    return(list(pd, ad))
  }
}

#' @importFrom sp coordinates
pwd_sample <- function(fixed, sample, reference, tr = 0.33, nearest= TRUE, n=1, lonlat = TRUE, warn = TRUE) {
  distfun <- get_distfun(lonlat, dismo_like = TRUE)
  stopifnot(tr > 0)
  n <- round(n)
  stopifnot(n >= 1)
  fixed <- get_coordinates(fixed)
  sample <- get_coordinates(sample)
  reference <- get_coordinates(reference)
  if (warn & nrow(sample) < nrow(fixed)) {
    warning("nrow(sample) < nrow(fixed)")
  }
  fromd <- mindist(distfun, fixed, reference, lonlat)
  tod <- mindist(distfun, sample, reference, lonlat)
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
