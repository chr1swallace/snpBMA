## example data

snps<-letters[1:12]
groups <- list(letters[1:3], letters[4:5], letters[6:8],letters[9], letters[10], letters[11], letters[12])
n.use <- 2
snps.single <- snps[9:12]
snps.multi <- snps[1:8]

context("make.models and make.groups\n")

single.groups <- make.models.single(snps.single,n.use,quiet=TRUE)
multi.groups <- make.models.multi(snps.multi,n.use,groups[1:3],quiet=TRUE)
both.groups <- make.models(snps,n.use,groups,quiet=TRUE)

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

################################################################################

context("max.models\n")
test_that("max.models works\n", {
  expect_equal( nrow(single.groups), max.models.single(n=length(snps.single), n.use=n.use))
  expect_equal( nrow(multi.groups), max.models.multi(groups=groups[1:3], n.use=n.use ))
  expect_equal( nrow(both.groups), max.models(snps=snps, n.use=n.use, groups=groups))
})    
    
################################################################################
context("mdrop\n")
test_that("mdrop keeps correct rows", {
  dropped <- mdrop(cbind2(single.groups,multi.groups[1:nrow(single.groups),]),
                   single.groups[1:3,], quiet=TRUE)
  expect_identical(dropped[,1:ncol(single.groups)], single.groups[-c(1:3),])
})

################################################################################

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


################################################################################

context("growing models\n")

model.string <- function(m) {
  sort(apply(m,1,paste,collapse=""))
}

test_that("mgrow works on singles", {
  single.3 <- make.models.single( snps.single, n.use+1, quiet=TRUE )
  bma <- new("snpBMA",
             snps=snps.single,
             models=single.groups,
             groups=list(),
             bf=matrix(0,nrow(single.groups),1),
             nsnps=2,
             nmodels=0)
  single.grow <- mgrow(bma, quiet=TRUE)
  expect_identical(  model.string(single.3), model.string(single.grow) )
})

test_that("mgrow works on groups", {
  multi.3 <- make.models.multi( snps.multi, n.use+1, groups=groups[1:3], quiet=TRUE )
  bma <- new("snpBMA",
             snps=snps.multi,
             models=multi.groups,
             groups=groups,
             bf=matrix(0,nrow(multi.groups),1),
             nsnps=2,
             nmodels=0)
  multi.grow <- mgrow(bma)

  expect_identical( model.string(multi.3), model.string(multi.grow) )
})


test_that("mgrow works on mixed groups and singles", {
  both.3 <- make.models( snps, n.use+1, groups=groups, quiet=TRUE )
  bma <- new("snpBMA",
             snps=snps,
             models=both.groups,
             groups=groups,
             bf=matrix(0,nrow(both.groups),1),
             nsnps=2,
             nmodels=0)
  both.grow <- mgrow(bma)

  expect_identical( model.string(both.3), model.string(both.grow) )
})

## models.group.collapse
