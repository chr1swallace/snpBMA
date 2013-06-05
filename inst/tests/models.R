## example data

snps<-letters[1:12]
groups <- list(letters[1:3], letters[4:5], letters[6:8],letters[9], letters[10], letters[11], letters[12])
n.use <- 2

single.groups <- make.models.single(snps[9:12],n.use,quiet=TRUE)
multi.groups <- make.models.multi(snps[1:8],n.use,groups[1:3],quiet=TRUE)
both.groups <- make.models(snps,n.use,groups,quiet=TRUE)

context("make.models run\n")

test_that("make.groups produces dgCMatrix-es", {
  expect_that(single.groups, is_a("dgCMatrix")) 
  expect_that(multi.groups, is_a("dgCMatrix")) 
  expect_that(both.groups, is_a("dgCMatrix"))
})

test_that("make.groups nrows are correct", {
  expect_that(nrow(single.groups), equals(6))
  expect_that(nrow(multi.groups), equals(21))
  expect_that(nrow(both.groups), equals(59))
})

test_that("make.groups have correct number of 1s per row", {
  expect_that(all(apply(single.groups==1,1,sum)==n.use), is_true())
  expect_that(all(apply(multi.groups==1,1,sum)==n.use), is_true())
  expect_that(all(apply(both.groups==1,1,sum)==n.use), is_true())
})

test_that("mdrop keeps correct rows", {
  dropped <- mdrop(cbind2(single.groups,multi.groups[1:nrow(single.groups),]),
                   single.groups[1:3,], quiet=TRUE)
  expect_identical(dropped[,1:ncol(single.groups)], single.groups[-c(1:3),])
})
          
context("make.models deals with exception: too many SNPs to model\n")

test_that("make.groups with more modelled SNPs than exist returns nrow(0) Matrices", {
  n.use <- length(groups)
  single.groups <- make.models.single(snps=snps[9:12],n.use=n.use,quiet=TRUE)
  multi.groups <- make.models.multi(snps=snps[1:8],n.use=n.use,groups=groups[1:3],quiet=TRUE)
  both.groups <- make.models(snps,n.use=n.use,groups,quiet=TRUE)  

  expect_that(nrow(single.groups), equals(0))
  expect_that(nrow(multi.groups), equals(0))
  expect_that(nrow(both.groups), equals(prod(sapply(groups,length))))
})

