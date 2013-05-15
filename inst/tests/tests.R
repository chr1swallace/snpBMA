## load sample data
data(for.exercise, package="snpStats")
test.data <- snps.10[,11:20]

## tagging
test_that("tag finds the right number of tags", {
  expect_equal(length(unique(tag(snps=colnames(test.data),tag.threshold=0.5,X=test.data))), 3)
})

## make.group
test_that("make.group checks parameters", {
  expect_that( make.group(-1), throws_error())
  expect_that( make.group(1.1), throws_error())
})

## do.merge
groups <- list(a=c("a","b","c"),
               d=c("d","e","f"),
               g=c("g","h"))
test_that("do.merge works", {
  expect_equal( length(do.merge(groups, c("a","d"))), 2)
  expect_that( do.merge(groups, c("a","e")), throws_error())
})
          

## l.rbind/ l.c
test_that("l.rbind works", {
  test.list <- list(list(A=list(a=1:2,b=3:5)),
                    list(A=list(a=11:15,b=21:28)))
  expect_identical(l.c( test.list, "A", "b" ), 
                   c(3:5,21:28))
})

test_that("l.c works", {
  test.list2 <- list(list(A=list(a=matrix(1:4,2,2))),
                     list(A=list(a=matrix(5:10,3,2))))
  expect_identical(l.rbind( test.list2, "A", "a" ), 
                   rbind(matrix(1:4,2,2), matrix(5:10,3,2)))
})


