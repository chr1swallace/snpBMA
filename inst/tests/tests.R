## load sample data
##context("prep\n")
data(for.exercise, package="snpStats")
test.data <- snps.10[,11:20]

## tagging
context("tagging\n")
test_that("tag finds the right number of tags", {
  expect_equal(length(unique(tag(X=test.data,tag.threshold=0.5))), 3)
})

## make.group
## test_that("make.group checks parameters", {
##   expect_that( make.group(-1), throws_error())
##   expect_that( make.group(1.1), throws_error())
## })

## do.merge
## groups <- list(a=c("a","b","c"),
##                d=c("d","e","f"),
##                g=c("g","h"))
## test_that("do.merge works", {
##   expect_equal( length(do.merge(groups, c("a","d"))), 2)
##   expect_that( do.merge(groups, c("a","e")), throws_error())
## })
          
################################################################################

## l.rbind/ l.c
context("l.rbind/l.c\n")

test_that("l.c works", {
  test.list <- list(list(A=list(a=1:2,b=3:5)),
                    list(A=list(a=11:15,b=21:28)))
  expect_identical(l.c( test.list, "A", "b" ), 
                   c(3:5,21:28))
})

test_that("l.rbind works", {
  test.list2 <- list(list(A=list(a=matrix(1:4,2,2))),
                     list(A=list(a=matrix(5:10,3,2))))
  expect_identical(l.rbind( test.list2, "A", "a" ), 
                   rbind(matrix(1:4,2,2), matrix(5:10,3,2)))
})

################################################################################

## models

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

