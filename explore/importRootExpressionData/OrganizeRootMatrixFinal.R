library(RUnit)
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_splitNames()
   test_eliminateDuplicateRows()

} # runTests
#----------------------------------------------------------------------------------------------------
if(!exists("test.tbl"))
   load("test.tbl.RDAta")

#----------------------------------------------------------------------------------------------------
eliminateDuplicateRows <- function(tbl)
{
   last.column <- ncol(tbl)
   mtx <- as.matrix(tbl[, 3:last.column])
   variances.of.just.the.number.columns.after.making.a.matrix <- apply(mtx, 1, var)
   tbl$variance <- variances.of.just.the.number.columns.after.making.a.matrix
   # adding the column "variance" to the data frame "tbl"
   # used the apply function to get the variance for each row in "tbl"
   tbl.sorted <- tbl[order(tbl$Alias, tbl$variance, decreasing= TRUE), ]
   # sorted the table in decreasing order of the "Alias" & "variance" columns.
   # sorting the data frame this way allows us see the duplicate Alias names and the variance associated
   # By having it in decreasing order, we can then delete the duplicates with the lowest variance
   # Get rid of duplicates
   tbl.sorted.nodupes <- tbl.sorted[-which(duplicated(tbl.sorted$Alias)), ]
   # deleting the duplicates (with the lower variance) from the Alias column
   # setting these changes in another version of the data frame, named "tbl.sorted.nodupes"
   rownames(tbl.sorted.nodupes) <- tbl.sorted.nodupes$Alias
   # Setting the Alias column names to now be the rownames of the data frame (it was previously in numerical order)
   tbl.sorted.nodupes <- tbl.sorted.nodupes[, -c(1:2)]
   #deleting columns 1 & 2 from the data frame: the "id" column & "Alias" column.
   #In the previous step, the Alias column was used for the rownames, so it is not unnecessary to have the column
   variance.column <- grep("variance", colnames(tbl.sorted.nodupes))
   #search for the column name with "variance". set it as variance.column
   tbl.sorted.nodupes <- tbl.sorted.nodupes[, -variance.column]
   #delete variance column
   return(tbl.sorted.nodupes)
   #result is an organized table with no duplicate rows

} # eliminateDuplicateRows
#----------------------------------------------------------------------------------------------------
test_eliminateDuplicateRows <- function()
   #create test for eliminateDuplicateRows
{
   message(sprintf("--- test_eliminateDuplicateRows"))
   #print message stated
   there.are.duplicates <- any(duplicated(test.tbl$Alias))
   #set a variable stating that there are duplicates in the Alias column on the test.tbl
   checkTrue(there.are.duplicates)
   #check if there are duplicates
   tbl.fixed <- eliminateDuplicateRows(test.tbl)
   #set a variable showing the result table
   checkTrue(nrow(tbl.fixed) < nrow(test.tbl))
   #check if there are less rows in the result table than than the original table
   checkTrue(ncol(tbl.fixed) < ncol(test.tbl))
   #check if there are less columns in the result table than the original table
   checkEquals((length(grep(";", test.tbl$Alias))), (length(grep(";", rownames(tbl.fixed)))))
   #check if the # of times ";" is found in the orginal table is the same as the result table
   #in our first original table (the smaller table before the entire root.tbl is used), there should be 15

} # test_loadData
#----------------------------------------------------------------------------------------------------
# create function splitNames
splitNames <- function(string)
{
   if(!grepl(";", string))
   {return(string)}
   #if there isnt a ";" in the parameter, just simply return the the parameter
   else if(grepl(";", string)){
      #if there is a ";" in our input, do the following:

      singleOrfNames <- unlist(strsplit(string, ";"))
      #split the rownames by ";", this alone will print a list, so "unlist" is used to print out characters/strings
      return(singleOrfNames)}
}
#---------------------------------------------------------------------------------------------------
splitNames.2 <- function(names)
{
   strsplit(names, ";")

} # splitNames.2
#---------------------------------------------------------------------------------------------------
test_splitNames.2 <- function()
{
   message(sprintf("--- test_splitNames.2"))
   checkEquals(splitNames.2("abc"), list("abc"))

   checkEquals(splitNames.2("abcd;efgh"), list(c("abcd", "efgh")))

   checkEquals(splitNames.2(c("abc", "abcd;efgh")), list("abc", c("abcd", "efgh")))

} # test_splitNames.2
#---------------------------------------------------------------------------------------------------
test_splitNames <- function()
{
   message(sprintf("--- test_splitNames"))
   chesplitNames.2("abc")

   checkEquals(splitNames("abc"), "abc")
   checkEquals(splitNames(""), "")
   #checking for no splits

   result <- splitNames.2("abcd;efgh")
   checkEquals(result, c("abcd", "efgh"))
   #checking for 1 splits

   result.2 <- splitNames("abcd;efgh;ijkl")
   checkEquals(result.2, c("abcd", "efgh", "ijkl"))
   #checking for 2 splits

   result.3 <- splitNames("abcd;efgh;ijkl;mnop")
   checkEquals(result.3, c("abcd", "efgh", "ijkl", "mnop"))
   #checking for 3 splits

   result.4 <- splitNames("abcd;efgh;ijkl;mnop;qrst")
   checkEquals(result.4, c("abcd", "efgh", "ijkl", "mnop", "qrst"))
   #checking for 4 splits

} # test_splitNames
#----------------------------------------------------------------------------------------------------
#create function splitNamesDoubleRows by using previously made splitNames function

tbl <- eliminateDuplicateRows(test.tbl)
splitNamesRepeatRows <- function(tbl)
{
   mtx <- as.matrix(numRepeatRows)
   numRepeatRows <- list(apply(mtx, 1, length))
   x <- lapply(1:nrow(tbl), function(i) splitNames(tbl[i,]))
   # (1:nrow(tbl)) creates a sequence/list to loop over by making a vector of integers calculated from the # of rows in "tbl"
   # in this case, tbl 400 rows so integers 1-25
   # "i" allows any variable to be used in the function (mainly loops); is simply used as an index
   # lapply takes the vector of integers calculated from the # of rows in tbl & the function as inputs so the function is applied to each element in the vector
   # sets it to a variable "x" to make life easier
   tbl.combined <- do.call(rbind, x)
}
#----------------------------------------------------------------------------------------------------
test_splitNamesRepeatRows <- function()
{
   message(sprintf("--- test_splitNamesRepeatRows"))

   stopifnot(exists("test.tbl"))
   tbl.fixed <- eliminateDuplicateRows(test.tbl)
   tbl.tiny <- tbl.fixed[1:6, 1:5]


} # test_splitNamesRepeatRows
#----------------------------------------------------------------------------------------------------
replicate.rows.when.needed <- function(row)
{
  browser()
  new.row.names <- unlist(splitNames.2(rownames(row)))
  total.row.count <- length(new.row.names)

  if(total.row.count == 1)
     return(row)

    # we get here only if total.row.count > 1

  tbl.expanded <- matrix(rep(t(row), total.row.count), ncol = ncol(row) , byrow = TRUE )
  colnames(tbl.expanded) <- colnames(row)
  rownames(tbl.expanded) <- new.row.names

  return(tbl.expanded)

} # replicate.rows.when.needed
#----------------------------------------------------------------------------------------------------
test_replicate.rows.when.needed <- function()
{
   stopifnot(exists("test.tbl"))
   tbl.fixed <- eliminateDuplicateRows(test.tbl)
   tbl.tiny <- as.matrix(tbl.fixed[1:6, 1:5])

   tbl.expanded2 <- replicate.rows.when.needed(tbl.tiny[2,])

   tbl.expanded <- replicate.rows.when.needed(tbl.tiny[1,])
     # do a bunch of checks here

   list.of.tables <- apply(tbl.tiny, 1, replicate.rows.when.needed)
   tbl.fixed <- do.call(rbind, list.of.tables)
   checkEquals(dim(tbl.fixed), c(8, 5))

} # test_replicate.rows.when.needed
#----------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
#---------------------------------------------------------------------------------------------------
