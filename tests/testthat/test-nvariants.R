test_that("nvariant gives correct outputs and warnings", {
    expect_equal(nvariant(N = 2000, phi = 1, omega = 1),2000)
    expect_warning(nvariant(N = 1000, phi = 1, omega = 1))
    expect_error(nvariant(N = 1000))
    
    expect_error(fit_nvariant(Observed_variants_per_kb = data.frame(n = c(2,1),
                                                         per_kb = c(1,2))))
    expect_error(fit_nvariant(Observed_variants_per_kb = data.frame(n = c(1,2),
                                                         per_kb = c(NA,2))))
})