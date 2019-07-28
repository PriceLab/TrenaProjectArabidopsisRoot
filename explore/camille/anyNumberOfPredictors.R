library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
fitModel <- function(tbl, target, predictors)
{
   variables.proposed <- c(target, predictors)
   stopifnot(all(variables.proposed %in% colnames(tbl)))

   part.1 <- sprintf("%s ~ ", target)
   part.2 <- paste(predictors, collapse=" + ")
   formula.as.string <- paste(part.1, part.2, sep="")

   model <- lm(formula.as.string, data=tbl)

   return(model)

} # fitModel
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_fitModelAnyNumberOfPredictors()
}
#------------------------------------------------------------------------------------------------------------------------
# for now, this test ensures only that any number of predictors can be used by the function defined abovd, "fitModel"
test_fitModelAnyNumberOfPredictors <- function()
{
   m0 <- fitModel(mtcars, target="mpg", "cyl")
   expected <- c("(Intercept)","cyl")
   checkEquals(names(coefficients(m0)), expected)

   m1 <- fitModel(mtcars, target="mpg", c("cyl", "hp"))
   expected <- c("(Intercept)","cyl", "hp")
   checkEquals(names(coefficients(m1)), expected)

   all.predictors <- setdiff(colnames(mtcars), "mpg")  # everything but the target
   m2 <- fitModel(mtcars, target="mpg", all.predictors)
   expected <- c("(Intercept)", all.predictors)
   checkEquals(expected, expected)


} # tets_fitModelAnyNumberOfPredictors
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
