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

#loading RUnit package
library(RUnit)

#function for r squared values (mpg vs one predictor)
#----------------------------------------------------------------------------------------------------
model <- function(target, predictors)
{
   model.x <- lm(formula = mtcars[[target]] ~ mtcars[[predictor]])
   model.x <-lm(formula = mpg ~ x, data=mtcars)
   summary(model.x)$r.squared
}
#----------------------------------------------------------------------------------------------------


#test model function
test_model <- function(){
  model.cyl1 <- lm(formula = mpg ~ cyl, data=mtcars)
  checkEquals(model(cyl),summary(model.cyl1)$r.squared)
  model.qsec1 <- lm(formula = mpg ~ qsec, data=mtcars)
  checkEquals(model(qsec),summary(model.qsec1)$r.squared)
}

#function for r squared values (mpg vs two predictors)
model1 <- function(x,y)
{model.b <-lm(formula = mpg ~ x + y, data=mtcars)
summary(model.b)$r.squared}

#test function
test_model1 <- function(){
  model.cylwt <- lm(formula = mpg ~ cyl + wt, data=mtcars)
  checkEquals(model1(cyl, wt),summary(model.cylwt)$r.squared)
  model.hpam <- lm(formula = mpg ~ hp + am, data=mtcars)
  checkEquals(model1(hp, am),summary(model.hpam)$r.squared)
}

#makes data frame of r squared values
r.squared <- c(0.72618,0.7183433,0.6024373,0.4639952,0.7528328,0.1752963,0.4409477,0.3597989,0.2306734,0.3035184)
names(r.squared) <- c("cyl","disp","hp","drat","wt","qsec","vs","am","gear","carb")
r_squared_data_frame <- data.frame(rSquared = r.squared)

#uses data frame to make bar graph
par(las=2,mar=c(4,4,2,1))
with(r_squared_data_frame, barplot(rSquared, names.arg = names(r.squared), main = "R Squared Values", xlab = "Predictors", ylab = "R-Squared"))

#model repeated for each predictor
for(i in 2:11)
{
  b <- as.numeric(mtcars[,i])
  model(b)
}


#par(las=1,mar=c(3,2,1,1))      slightly cut off left side
#par(las=2,mar=c(3,3,2,1))      vertical labels, no axis labels

#barplot(r.squared, main = "R Squared Values", xlab = "Predictors", ylab = "R-Squared")
#works using vector instead of data frame to make bar graph
