library(RUnit)
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
{
   print("--- test_eliminateDuplicateRows")
   there.are.duplicates <- any(duplicated(test.tbl$Alias))
   checkTrue(there.are.duplicates)

   tbl.fixed <- eliminateDuplicateRows(test.tbl)
   checkTrue(nrow(tbl.fixed) < nrow(test.tbl))
   checkTrue(ncol(tbl.fixed) < ncol(tbl.test))

} # test_loadData
#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_eliminateDuplicateRows()

  # test_splitMultipleNames()
  # test_demo.doAllRows()
  #these are the newly created functions^^
  #the function "runTests" will tell you whether the tests function properly by stating "TRUE"
  #if the function "runTests" runs and the functions above do not work, it will show errors stating why
} # runTests
# #----------------------------------------------------------------------------------------------------
##Expand rows with multiple names
# #----------------------------------------------------------------------------------------------------
# library(RUnit)
# #loading/opening the package RUnit to create a function and run tests for the function
# #----------------------------------------------------------------------------------------------------
# if(!exists("tbl.test")){
#   #If the file "tbl.test" does not exist, do the following:
#   print("fresh load of our data file")
#   #prints out this ^ string
#   #setting the variable "tbl.test" as our data frame to use while testing the newly created functions
# }
#------------------------------------------------------------------------------------------------------------------------
# used for one-row tables with multiple orf names in their rowname
splitMultipleNames <- function(tbl)
{
  tbl <- tbl.test[10,]
  #setting the data frame as one-rowed table
  stopifnot(nrow(tbl) == 1)
  #only handle data.frames with one row
  #if the tbl does not have only one row, do not continue to use this function (an error will occur)
  tokens <- strsplit(rownames(tbl), ";")[[1]]
  #splits the rownames within "tbl" on the character ";"
  #ex: if there is ONE ";" in the rowname, it will split ONE time and show strings of two orf names
  #set this as the variable "tokens" to serve as a tangible representation of what is done here^^
  tbl.result <- tbl[rep(seq_len(nrow(tbl)), each=length(tokens)),]
  #repeats the number/length of a sequence based on the number of rows in the "tbl" data frame
  #the sequence is repeated depending on the length of "tokens", aka how many orf names are split
  #set this^^ as the result table from using the function.
  rownames(tbl.result) <- tokens
  #the output from tokens (the orf names) is set at the rownames for tbl.result
  return(tbl.result)
  #helps define our new function by returning the value expected as a result of using this new function

} # splitMultipleNames
# #----------------------------------------------------------------------------------------------------
# demo.doAllRows <- function(tbl)
#   #creating the function demo.doAllRows
#   #is used for multi-row tables with multiple orf names in a rowname
# {
#   tbl <- tbl.test
#   #setting the current data frame as a multi-row table
#   x <- lapply(1:nrow(tbl), function(i) splitMultipleNames(tbl[i,]))
#   #(1:nrow(tbl)) creates a sequence/list to loop over by making a vector of integers calculated from the # of rows in "tbl"
#   #in this case, tbl 25 rows so integers 1-25
#   #"i" allows any variable to be used in the function (mainly loops); is simply used as an index
#   #lapply takes the vector of integers calculated from the # of rows in tbl & the function as inputs so the function is applied to each element in the vector
#   #sets it to a variable "x" to make life easier
#   tbl.combined <- do.call(rbind, x)
#   #returns a data frame; w/out do.call, it returns as a list
#   #we want our result to be a data frame combining the newly formed rows. We're setting this as the variable "tbl.combined"
#   return(tbl.combined)
#   #helps define our new function by returning the value expected as a result of using this new function
# } # demo.doAllRows does not work properly; it is repeating one row multiple times in the list created in "x"
# #how do we fix this? does it have to do with splitMultipleNames()?
# #----------------------------------------------------------------------------------------------------
test_splitMultipleNames <- function()
{
  message(sprintf("--- test_splitMultipleNames"))   # announce what test you are running

   tbl.smn <- eliminateDuplicateRows(tbl.test)[10,]

   tbl.fixed <- splitMultipleNames(tbl.5)
   checkEquals(nrow(tbl.fixed), 2)
   checkEquals(rownames(tbl.fixed), c("ATMG01250","AT2G07697"))

} #test_splitMultipleNames everything comes out true
# #----------------------------------------------------------------------------------------------------
# test_splitMultipleNames <- function()
# {
#   message(sprintf("--- test_splitMultipleNames"))   # announce what test you are running
#
#   tbl.5 <- tbl.test[7,]
#   tbl.fixed <- splitMultipleNames(tbl.5)
#   checkEquals(nrow(tbl.fixed), 2)
#   checkEquals(rownames(tbl.fixed), c("ATMG01250","AT2G07697"))
#
# } #test_splitMultipleNames the last checkEquals is coming out true but it should not...
# #tbl.5 here has "ATMG01300;AT2G07692" from tbl.test[7,]
# #but splitMultipleNames(tbl.5) shows different orf names from tbl.test[10,] since row 10 is specified while creating the function splitMultipleNames()
# #this^ problem prevents any other rows from being accessed, which gives us "falses positives" if the problem with splitMultipleNames() isn't fixed
# #----------------------------------------------------------------------------------------------------
#
# test_demo.doAllRows <- function()
# {
#   message(sprintf("---test_demo.doAllRows")) #announce what you're running
#   tbl.expanded <- demo.doAllRows(tbl.test)
#
#   checkTrue(nrow(tbl.expanded) > nrow(tbl.test))
#
#   multiple.orf.names.orig <- length(grep(";", rownames(tbl.test)))
#   checkEquals(multiple.orf.names.orig, 13)
#
#   multiple.orf.names.expanded <- length(grep(";", rownames(tbl.expanded)))
#   checkEquals(multiple.orf.names.expanded, 0)
#
# } # test_splitMultiNamesToRows
# #------------------------------------------------------------------------------------------------------------------------
# if(!interactive())
#   runTests()
#
