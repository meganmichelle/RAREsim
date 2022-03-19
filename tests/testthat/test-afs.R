test_that("afs function gives correct oupts and warnings", {
    mac_test <- data.frame(Lower = c(1), Upper = c(1))
    expect_equal(afs(b = 100, beta = 10, alpha = 0, mac = mac_test)$Prop, 100)
    expect_error(afs(b = 1, beta = 2, alpha = 3))
    expect_error(afs(b = 1, beta = 2, mac = mac_test))
    
    expect_error(fit_afs(Observed_bin_props = data.frame(Lower = c(2,1),
                                              Upper = c(1,2), Prop = c(.5,.2))))
    expect_error(fit_afs(Observed_bin_props = data.frame(Lower = c(1,2),
                                              Upper = c(1,2), Prop = c(NA,.3))))
})