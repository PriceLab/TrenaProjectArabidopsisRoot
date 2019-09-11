#get matrix (set path to "/Users/bioadmin/github/TrenaProjectArabidopsisRoot/inst/extdata/expression")
load("aluru.18617x1938.orfRowNames.RData")

# Get high expression data
wbc19.vector <- mtx[c("AT3G55130", "AT3G59060"),]
plotWBC19 <- plot((mtx[c("AT3G55130"),]), (mtx[c("AT3G59060"),]), main= "WBC19 Expression in Roots")
new.wbc.vector <- t(wbc19.vector)
number1 <- new.wbc.vector[(new.wbc.vector[,1] > 2),]
number2 <- number1[(number1[,2] > 2),]
plot(number2, main="High WBC19 Expression in Roots")
highExpressionWBC19 <- t(number2)
rownames(highExpressionWBC19)
mtx2 <- highExpressionWBC19
mtx2.colnames <- colnames(mtx2)
mtx2 <-  mtx[,mtx2.colnames]

# created matrix of high expression data
# mtx3 has 169 samples
ordered.number2<- number2[order(number2[,1], number2[,2], decreasing=TRUE),]
# allows you to view sample names with highest expression
ordered.number2.50 <- ordered.number2[1:50,]
evenHigherExpressionWBC19 <- t(ordered.number2.50)
rownames(evenHigherExpressionWBC19)
mtx3.colnames <- colnames(evenHigherExpressionWBC19)
mtx3 <- mtx[,mtx3.colnames]plotWBC19 <- plot((mtx[c("AT3G55130"),]), (mtx[c("AT3G59060"),]), main= "WBC19 Expression in Roots")

#now build low expression
new.wbc.vector1 <- t(wbc19.vector)
number3 <- new.wbc.vector[(new.wbc.vector[,1] < -2),]
number4 <- number3[(number3[,2] < -2),]
plot(number4, main="low WBC19 Expression in Roots")
lowExpressionWBC19 <- t(number4)
rownames(lowExpressionWBC19)
[1] "AT3G55130" "AT3G59060"
mtx4 <- lowExpressionWBC19
mtx4.colnames <- colnames(mtx4)
mtx4 <-  mtx[,mtx4.colnames]
length(mtx4)
[1] 120
ordered.number4<- number4[order(number4[,1], number4[,2], decreasing=FALSE),]
# allows you to view sample names with lowest expression first
ordered.number4.50 <- ordered.number4[1:50,]
evenLowerExpressionWBC19 <- t(ordered.number4.50)
rownames(evenLowerExpressionWBC19)
[1] "AT3G55130" "AT3G59060"
mtx5.colnames <- colnames(evenLowerExpressionWBC19)
mtx5 <- mtx[,mtx5.colnames]
highAndLowMtx <- cbind(mtx3, mtx5)



