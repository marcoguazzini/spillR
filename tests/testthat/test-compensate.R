test_that("spillover probability works", {
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

    suppressWarnings({
        res <- spillR:::compensate(
            tb_real, tb_bead, target_marker, spillover_markers, n_iter = 2
            )
    })
    
    y_support <- sort(unique(tb_real$marker1))
    p_pkg <- res$tb_spill_prob |>
        dplyr::filter(marker1 %in% y_support) |>
        dplyr::pull(spill_prob) |>
        round(digits = 3)

    # from step-by-step calculation in the paper (appendix A)
    p_paper <- c(0.631, 0.449, 0, 0)

    expect_equal(p_pkg[1], p_paper[1])
    expect_equal(p_pkg[2], p_paper[2])
    expect_equal(p_pkg[3], p_paper[3])
    expect_equal(p_pkg[4], p_paper[4])
})
