context("binomial prior\n")
test_that("prior.binomial works", {
  expect_that(prior.binomial(-1,100,3) , throws_error())
  expect_that(prior.binomial(500,100,3) , throws_error())
  expect_that(prior.binomial(x=0:10,n=100,expected=3,value="prob"), is_equivalent_to(dbinom(0:10, size=100, prob=3/100)))
})

context("beta binomial prior\n")
test_that("prior.betabinomial works", {
  expect_that(prior.betabinomial(-1,100,3) , throws_error())
  expect_that(prior.betabinomial(500,100,3) , throws_error())
  expect_that(prior.betabinomial(0:10,100,3,0.9) , throws_error())
  expect_that(prior.betabinomial(0:10,100,3,1,value="prob"), is_equivalent_to(dbinom(0:10, size=100, prob=3/100)))
  expect_that(prior.betabinomial(0:10,100,3,1.2,value="prob"), is_equivalent_to(dbetabinom(0:10, size=100, prob=3/100, rho=0.2/99))) 
})

context("value argument\n")
test_that("value argument works", {
  expect_that(prior.binomial(1,1,1,value="jjj"), throws_error())
  expect_that(prior.binomial(x=1:10,n=100,expected=3,value="odds",pi0=0.1),
              is_equivalent_to(prior.binomial(x=1:10,n=100,expected=3, value="prob")/0.1))
  expect_that(prior.betabinomial(1:10,100,3,2,value="odds",pi0=0.1),
              is_equivalent_to(prior.betabinomial(x=1:10,n=100,expected=3,overdispersion=2, value="prob")/0.1))
})

test_that("prior odds for 0 is 1", {
  expect_that(prior.binomial(x=0,n=100,expected=3,value="odds"),
              is_equivalent_to(1))
  expect_that(prior.betabinomial(x=0,n=100,expected=3,overdispersion=2,value="odds"),
              is_equivalent_to(1))
})
