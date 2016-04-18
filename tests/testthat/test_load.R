library(marinespeed)

context("load")

set_datadir_tmp <- function() {
  tmp <- tempdir()
  options(marinespeed_datadir=tmp)
  tmp
}
test_that("get datadir", {
  tmp <- set_datadir_tmp()
  expect_equal(marinespeed:::get_datadir(), tmp)
  options(marinespeed_datadir=NULL)
  expect_true(dir.exists(marinespeed:::get_datadir()))
})

test_that("get version", {
  expect_true(get_version(), "V1")
  options(marinespeed_version="V.test")
  expect_true(get_version(), "V.test")
  options(marinespeed_version=NULL)
  expect_true(get_version(), "V1")
})

test_that("list species", {
  tmp <- set_datadir_tmp()
  species <- list_species()
  expect_more_than(NROW(species), 500)
  expect_equal(NCOL(species), 2)
  expect_equal(colnames(species), c("species", "aphia_id"))
})

test_that("get occurrences files works", {
  all_paths <- get_occurrence_files()

  species <- list_species()
  first10 <- get_occurrence_files(species[1:10,])

  abalistes_stellatus <- get_occurrence_files("Abalistes stellatus")
})

test_that("get_fold random", {
  species <- list_species()
  fold <- get_folds(species[sample(1:nrow(species), 1),], fold_type = "random", k=1)
  fold$occurrence_training
  fold$occurrence_test
  fold$background_training
  fold$background_test

})

t <- function() {
  d <- marinespeed::get_occurrences_dir(raw = FALSE)
  s <- marinespeed::list_species()
  marinespeed::load_folds()
  marinespeed::lapply_species(function(speciesinfo, data) { print(paste(speciesinfo, NROW(data)))}, raw = FALSE)


  marinespeed:::get_file("test.csv.gz")
  marinespeed:::get_file("test.zip")
  marinespeed:::get_datadir()
}
