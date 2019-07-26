library(RUnit)
source("mtcarsR.R")
#-------------------------------------------------------------------------------
run_tests <- function()
{
  test_one_predictor()
  test_two_predictors()
}#testing all
#----------------------------------------------------------------------------------------------------
test_one_predictor <- function()
{
  message(sprintf("--- test_one_predictor"))
  model.cyl1 <- lm(formula = mtcars$mpg ~ mtcars$cyl)
  checkEquals(one_predictor_rsquared("mpg", "cyl"),summary(model.cyl1)$r.squared)
  model.qsec1 <- lm(formula = mtcars$mpg ~ mtcars$qsec)
  checkEquals(one_predictor_rsquared("mpg", "qsec"),summary(model.qsec1)$r.squared)
}#test_one_predictor
#---------------------------------------------------------------------------
test_two_predictors <- function()
{
  message(sprintf("--- test_two_predictors"))
  model.cylwt <- lm(formula = mtcars$mpg ~ mtcars$cyl + mtcars$wt)
  checkEquals(two_predictors_rsquared("mpg", "cyl", "wt"),summary(model.cylwt)$r.squared)
  model.hpam <- lm(formula = mtcars$mpg ~ mtcars$hp + mtcars$am)
  checkEquals(two_predictors_rsquared("mpg","hp", "am"),summary(model.hpam)$r.squared)
}#test_two_predictors
#---------------------------------------------------------------------------
if(!interactive())
  run_tests()
