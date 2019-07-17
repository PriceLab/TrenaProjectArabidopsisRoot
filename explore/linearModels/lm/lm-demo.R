model.disp <- lm(formula = mpg ~ disp, data=mtcars)
plot(mtcars$disp, mtcars$mpg)
abline(model.disp, col="red")

# this will give you four successive plots.
# you will br prompted to "Hit <return> to see the next plot"
# I admit to not fully understanding these plots, but perhaps you can
# study up, and then explain them to us?

plot(model.disp)
