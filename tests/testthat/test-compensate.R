test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

library(tibble)

tb_real <- tibble(
    "marker1" = c(3, 5, 17, 3, 17, 2), 
    "barcode" = rep("none", 6),
    "type" = rep("real cells", 6)
    )

tb_bead <- tibble(
    "marker1" = c(2, 3, 2),
    "barcode" = rep("marker2", 3),
    "type" = rep("beads", 3)
    )

target_marker <- "marker1"
spillover_markers <- "marker2"

compensate(tb_real, tb_bead, target_marker, spillover_markers, runmed_k = 1)
