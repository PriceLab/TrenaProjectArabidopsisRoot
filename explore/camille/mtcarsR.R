#calculates r squared for each predictor of mpg
model.cyl <- lm(formula = mpg ~ cyl, data=mtcars)
summary(model.cyl)$r.squared # [1] 0.72618

model.disp <- lm(formula = mpg ~ disp, data=mtcars)
summary(model.disp)$r.squared  # [1] 0.7183433

model.hp <- lm(formula = mpg ~ hp, data=mtcars)
summary(model.hp)$r.squared # [1] 0.6024373

model.drat <- lm(formula = mpg ~ drat, data=mtcars)
summary(model.drat)$r.squared # [1] 0.4639952

model.wt <- lm(formula = mpg ~ wt, data=mtcars)# highest r-squared
summary(model.wt)$r.squared # [1] 0.7528328

model.qsec <- lm(formula = mpg ~ qsec, data=mtcars)# lowest r-squared
summary(model.qsec)$r.squared # [1] 0.1752963

model.vs <- lm(formula = mpg ~ vs, data=mtcars)
summary(model.vs)$r.squared # [1] 0.4409477

model.am <- lm(formula = mpg ~ am, data=mtcars)
summary(model.am)$r.squared # [1] 0.3597989

model.gear <- lm(formula = mpg ~ gear, data=mtcars)
summary(model.gear)$r.squared # [1] 0.2306734

model.carb <- lm(formula = mpg ~ carb, data=mtcars)
summary(model.carb)$r.squared # [1] 0.3035184

library(RUnit)

#----------------------------------------------------------------------------------------------------
one_predictor_rsquared <- function(target, predictor) #input with quotes
{
   model.x <- lm(formula = mtcars[[target]] ~ mtcars[[predictor]])
   summary(model.x)$r.squared
}#function for r squared values (mpg vs one predictor)
#--------------------------------------------------------------------------
two_predictors_rsquared <- function(target, predictor1, predictor2) #input with quotes
{
  model.b <-lm(formula = mtcars[[target]] ~ 
                 mtcars[[predictor1]] + mtcars[[predictor2]])
  summary(model.b)$r.squared
}#function for r squared values (mpg vs two predictors)
#---------------------------------------------------------------------------
multiple_predictors_model <- function(tbl, target, predictors) #input with quotes
{
    variables.proposed <- c(target, predictors)
    stopifnot(all(variables.proposed %in% colnames(tbl)))
    
    part.1 <- sprintf("%s ~ ", target)
    part.2 <- paste(predictors, collapse=" + ")
    formula.as.string <- paste(part.1, part.2, sep="")
    
    model <- lm(formula.as.string, data=tbl)
    
    return(summary(model)$r.squared)
}#multiple_predictors_model
#------------------------------------------------------------------------
r.squared <- c(0.72618,0.7183433,0.6024373,0.4639952,0.7528328,0.1752963,0.4409477,0.3597989,0.2306734,0.3035184)
names(r.squared) <- c("cyl","disp","hp","drat","wt","qsec","vs","am","gear","carb")
r_squared_data_frame <- data.frame(rSquared = r.squared)
#makes data frame of r squared values
#--------------------------------------------------------------------------
histogram_rsquared <- function()
{
par(las=2,mar=c(4,4,2,1))
with(r_squared_data_frame, barplot(rSquared, names.arg = names(r.squared), main = "R Squared Values", xlab = "Predictors", ylab = "R-Squared"))
}#uses data frame to make bar graph
#---------------------------------------------------------------------------
for(i in 2:11)
{
  model.x <- lm(formula = mtcars$mpg ~ mtcars[[i]])
  print(summary(model.x)$r.squared)
}#r squared calculation repeated for each predictor
#----------------------------------------------------------------------------
two_predictors_residuals <- function(target, predictor1, predictor2) #input with quotes
{
  model.n <-lm(formula = mtcars[[target]] ~ 
                 mtcars[[predictor1]] + mtcars[[predictor2]])
  summary(model.n)$residuals
}#model for residuals with multiple predictors
#----------------------------------------------------------------------------
rsquared <- function()
{
for(i in 2:11)
{
  for(c in 2:11)
  {
    model.b <-lm(formula = mtcars$mpg ~ 
                   mtcars[[i]] + mtcars[[c]])
    vector.b <- c(vector.b, i)
    vector.c <-c(vector.c, summary(model.b)$r.squared)
  }
}
  print(vector.c)
  data_frame_rsquared <- data.frame(predictor1 = vector.b, 
                                    predictor2 = c(rep(colnames(mtcars[,-1]), 10)), 
                                    rSquared = vector.c)#data frame doesnt work
}
#r squared calculated for each combo of 2 predictors
#----------------------------------------------------------------------------
rsquared2 <- function()
{
  all.predictors <- c(setdiff(colnames(mtcars),"mpg"))
  list1 <- as.list(c(rep(all.predictors,10)))
  list2 <- as.list(c(rep("cyl",10),rep("disp",10),rep("hp",10),rep("drat",10),rep("wt",10),
                     rep("qsec",10),rep("vs",10),rep("am",10),rep("gear",10),rep("carb",10)))
  two.combo.list <- cbind(list1,list2)
  apply(two.combo.list,1, multiple_predictors_model, target = "mpg", tbl = mtcars)
}#r squared calculated for each combo of 2 predictors
#----------------------------------------------------------------------------
