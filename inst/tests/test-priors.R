context("binomial prior\n")
test_that("prior.binomial works", {
  expect_that(prior.binomial(-1,100,3) , throws_error())
  expect_that(prior.binomial(500,100,3) , throws_error())
  expect_that(prior.binomial(0:10,100,3), is_equivalent_to(dbinom(0:10, size=100, prob=3/100)))
})

context("beta binomial prior\n")
test_that("prior.betabinomial works", {
  expect_that(prior.betabinomial(-1,100,3) , throws_error())
  expect_that(prior.betabinomial(500,100,3) , throws_error())
  expect_that(prior.betabinomial(0:10,100,3,0.9) , throws_error())
  expect_that(prior.betabinomial(0:10,100,3,1), is_equivalent_to(dbinom(0:10, size=100, prob=3/100)))
  expect_that(prior.betabinomial(0:10,100,3,1.2), is_equivalent_to(dbetabinom(0:10, size=100, prob=3/100, rho=0.2/99))) 
})
