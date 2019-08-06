library(RUnit)
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_splitNames()
   test_eliminateDuplicateRows()
   
} # runTests
#----------------------------------------------------------------------------------------------------
if(!exists("test.tbl"))
   get(load("test.tbl.RDAta"))

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
   
} # test_loadData
#----------------------------------------------------------------------------------------------------
# create function splitNames
splitNames <- function(string)
{
   if(!grepl(";", string))
   {return(string)}
   else if(grepl(";", string)){
      #if there isnt a ";" in the parameter, just simply the parameter
      
      singleOrfNames <- unlist(strsplit(string, ";"))
      #split the rownames by ";", this alone will print a list, so "unlist" is used
      return(singleOrfNames)}
}
#---------------------------------------------------------------------------------------------------
test_splitNames <- function()
{
   message(sprintf("--- test_splitNames"))
   checkEquals(splitNames("abc"), "abc")   
   checkEquals(splitNames(""), "")
   #checking for no splits

   result <- splitNames("abcd;efgh")
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
   #checking for 3 splits

} # test_splitNames
#----------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
#---------------------------------------------------------------------------------------------------
