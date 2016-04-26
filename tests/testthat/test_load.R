library(marinespeed)

context("load")

setup_load <- function() {
  skip_on_cran()
  # skip_on_travis()
}


test_that("get datadir", {
  tmp <- tempdir()
  options(marinespeed_datadir=tmp)
  expect_equal(marinespeed:::get_datadir(), tmp)
  options(marinespeed_datadir=NULL)
  expect_true(dir.exists(marinespeed:::get_datadir()))
})

test_that("get version", {
  expect_equal(get_version(), "v1")
  options(marinespeed_version="V.test")
  expect_equal(get_version(), "V.test")
  options(marinespeed_version=NULL)
  expect_equal(get_version(), "v1")
})

test_that("list species", {
  setup_load()

  species <- list_species()
  expect_more_than(NROW(species), 500)
  expect_equal(NCOL(species), 2)
  expect_equal(colnames(species), c("species", "aphia_id"))
})

test_that("get occurrences works", {
  setup_load()

  abalistes_stellatus <- get_occurrences("Abalistes stellatus")
  expect_more_than(NROW(abalistes_stellatus), 10)
  expect_more_than(NCOL(abalistes_stellatus), 50)

  species <- list_species()
  occ <- get_occurrences(species[1:3,])
  expect_more_than(NROW(occ), 10)
  expect_more_than(NCOL(occ), 50)
  expect_equal(length(unique(occ$species)), 3)
})

test_that("get_fold_data random", {
  setup_load()

  check_fold <- function(fold) {
    expect_more_than(NROW(fold$occurrence_training), 10)
    expect_more_than(NROW(fold$occurrence_test), 10)
    expect_more_than(NROW(fold$background_training), 1000)
    expect_more_than(NROW(fold$background_test), 100)
    expect_more_than(NCOL(fold$occurrence_training), 50)
    expect_equal(NCOL(fold$occurrence_test), NCOL(fold$occurrence_training))
    expect_equal(NCOL(fold$background_training), NCOL(fold$occurrence_training))
    expect_equal(NCOL(fold$background_test), NCOL(fold$occurrence_test))
  }

  folds <- get_fold_data("Abalistes stellatus", fold_type = "random", k=c(2,4))

  expect_null(folds[[1]])
  check_fold(folds[[2]])
  expect_null(folds[[3]])
  check_fold(folds[[4]])
  expect_null(folds[[5]])
  ## fold records should be different
  expect_more_than(length(setdiff(folds[[2]]$occurrence_training[,"longitude"], folds[[4]]$occurrence_training[,"longitude"])), nrow(folds[[2]]$occurrence_training) / 6)
})

test_that("get_file works", {
  setup_load()
  tmp <- tempdir()
  options(marinespeed_datadir=tmp)
  file <- get_file("species.csv.gz")
  expect_true(file.exists(file))
  species <- read.csv(file)
  expect_more_than(NROW(species), 500)
  file.remove(file)
  expect_false(file.exists(file))
  options(marinespeed_datadir=NULL)
})

test_that("get_folds works", {
  setup_load()
  check_folds <- function(folds) {
    expect_equal(NCOL(folds$background),6)
    expect_equal(NCOL(folds$species),6)
    expect_more_than(NROW(folds$background),1000)
    expect_more_than(NROW(folds$species),10000)
  }
  check_folds(get_folds("disc"))
  check_folds(get_folds("random"))
  check_folds(get_folds("targetgroup"))
})

test_that("get_background works", {
  setup_load()
  check_background <- function(bg) {
    expect_more_than(NCOL(bg), 50)
    expect_more_than(NROW(bg), 1000)
  }
  check_background(get_background("random"))
  check_background(get_background("random"))
})
