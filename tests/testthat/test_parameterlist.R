context("STGARMAGARCH testing.")


test_that("parameter vector and list", {
          expect_equal(class(parvector2list()), "list")
          expect_equal(length(parvector2list()), 6)
}
)
