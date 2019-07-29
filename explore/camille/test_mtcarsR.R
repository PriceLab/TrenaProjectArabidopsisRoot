library(RUnit)
source("mtcarsR.R")
#-------------------------------------------------------------------------------
runTests <- function()
{
  test_one_predictor_rsquared()
  test_two_predictors_rsquared()
  test_multiple_predictors_model()
  #test_two_predictors_residuals
}#testing all
#----------------------------------------------------------------------------------------------------
test_one_predictor_rsquared <- function()
{
  message(sprintf("--- test_one_predictor"))
  model.cyl1 <- lm(formula = mtcars$mpg ~ mtcars$cyl)
  checkEquals(one_predictor_rsquared("mpg", "cyl"),summary(model.cyl1)$r.squared)
  model.qsec1 <- lm(formula = mtcars$mpg ~ mtcars$qsec)
  checkEquals(one_predictor_rsquared("mpg", "qsec"),summary(model.qsec1)$r.squared)
  
}#test_one_predictor
#---------------------------------------------------------------------------
test_two_predictors_rsquared <- function()
{
  message(sprintf("--- test_two_predictors"))
  model.cylwt <- lm(formula = mtcars$mpg ~ mtcars$cyl + mtcars$wt)
  checkEquals(two_predictors_rsquared("mpg", "cyl", "wt"),summary(model.cylwt)$r.squared)
  model.hpam <- lm(formula = mtcars$mpg ~ mtcars$hp + mtcars$am)
  checkEquals(two_predictors_rsquared("mpg","hp", "am"),summary(model.hpam)$r.squared)
}#test_two_predictors
#---------------------------------------------------------------------------
test_multiple_predictors_model <- function()
{
  message(sprintf("--- test_multiple_predictors_model"))
  model.amwt <- lm(formula = mtcars$mpg ~ mtcars$am + mtcars$wt)
  checkEquals(multiple_predictors_model(mtcars, "mpg", c("am","wt")),summary(model.amwt)$r.squared)
  model.cylhpvs <- lm(formula = mtcars$mpg ~ mtcars$cyl + mtcars$hp + mtcars$vs)
  checkEquals(multiple_predictors_model(mtcars, "mpg", c("cyl","hp", "vs")),summary(model.cylhpvs)$r.squared)
}#test_multiple_predictors_model
#---------------------------------------------------------------------------


#test_two_predictors_residuals
#---------------------------------------------------------------------------
if(!interactive())
  runTests()
