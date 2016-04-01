library(marinespeed)

context("kfold")

lonlat_data <- function(nrow=50, seed=42) {
  set.seed(seed)
  lonlat <- cbind(lon=runif(nrow, -180, 180), lat=runif(nrow, -90, 90))
  lonlat
}


test_that("kfold disc returns a vector of valid folds", {
  fold_lengths <- function(folds, k=5) {
    sapply(1:5, function(k) sum(folds==k))
  }
  d <- lonlat_data(nrow=50)
  folds <- kfold_disc(d, k=5)
  lengths <- fold_lengths(folds, k=5)
  expect_identical(sum(lengths), 50)
  expect_identical(max(lengths) - min(lengths), 0)

  d <- lonlat_data(nrow=49)
  folds <- kfold_disc(d, k=5)
  lengths <- fold_lengths(folds, k=5)
  expect_identical(sum(lengths), 49)
  expect_identical(max(lengths) - min(lengths), 1)
})

test_that("kfold disc works with alternative distance function", {
  euclidean <- function(a, b) sqrt(sum((rep(a,NROW(b)) - b) ^ 2))

  d <- lonlat_data(nrow=50)
  folds <- kfold_disc(d, k=5, distfun = euclidean)
  lengths <- fold_lengths(folds, k=5)
  expect_identical(sum(lengths), 50)
  expect_identical(max(lengths) - min(lengths), 0)
})

test_that("kfold_disc different seed, different result, no seed, different result", {
  d <- lonlat_data(nrow=20)
  fold42a <- kfold_disc(d, k=5)
  fold42b <- kfold_disc(d, k=5)
  fold1a <- kfold_disc(d, k=5, seed=1)
  fold1b <- kfold_disc(d, k=5, seed=1)
  foldNULL <- kfold_disc(d, k=5, seed=NULL)
  foldNA1 <- kfold_disc(d, k=5, seed=NA)
  foldNA2 <- kfold_disc(d, k=5, seed=NA)
  expect_identical(fold42a, fold42b)
  expect_false(identical(fold42a, fold1a))
  expect_identical(fold1a, fold1b)
  expect_false(identical(foldNULL, foldNA1))
  expect_false(identical(foldNA1, foldNA2))
})

test_that("kfold_species_background with default settings works", {
  species <- lonlat_data(nrow=20)
  background <- lonlat_data(nrow=200)

  folds <- kfold_species_background(species, background)
  ## TODO add expect...
})


plot_folds <- function(d, folds) {
  plot(d, pch=".")
  cols <- rainbow(5)
  for(i in 1:5) {
    text(d[folds==i,], labels = i, pch=20, col=cols[i])
  }
}
