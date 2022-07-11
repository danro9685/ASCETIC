Sys.setenv("R_TESTS" = "")

library("testthat")
library("ASCETIC")

test_check("ASCETIC")
