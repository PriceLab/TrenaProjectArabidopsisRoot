library(RUnit)
#----------------------------------------------------------------------------------------------------
if(!exists("tbl.test")){
    tbl.test <- get(load("10by10rootTable.RData"))
    }
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
  test_splitMultipleNamesIntoTheirOwnRows()

} # runTests
#----------------------------------------------------------------------------------------------------
splitMultipleNamesIntoTheirOwnRows <- function(tbl)
{
  stopifnot(nrow(tbl) == 1)    # we only handle data.frames with one row

  if(!grepl(";", rownames(tbl))){
     return(tbl)  # just give it back to the caller
     }

  tokens <- strsplit(rownames(tbl), ";")[[1]]

  if(length(tokens) == 2)
     tbl.result <- rbind(tbl, tbl)

  if(length(tokens) == 3)
     tbl.result <- rbind(tbl, tbl, tbl)

  if(length(tokens) == 4)
     tbl.result <- rbind(tbl, tbl, tbl, tbl)

  rownames(tbl.result) <- tokens


  stop("no code yet to handle multiple orf names embedded in a row name")

  # your code to split the rowname, and duplicate the row values, building
  # up new data.frame, of two or more rows,  will happen here

} # splitMultipleNamesIntoTheirOwnRows
#----------------------------------------------------------------------------------------------------
test_splitMultipleNamesIntoTheirOwnRows <- function()
{
  message(sprintf("--- test_splitMultipleNamesIntoTheirOwnRows"))   # announce what test you are running

     #------------------------------------------------------------------------
     #  do the simplest test first: a row with a simple name (no semicolons)
     #------------------------------------------------------------------------

  tbl.1 <- tbl.test[1,]
  tbl.fixed <- splitMultipleNamesIntoTheirOwnRows(tbl.1)
  checkEquals(nrow(tbl.fixed), 1)

  tbl.2 <- tbl.test[10,]
  tbl.fixed <- splitMultipleNamesIntoTheirOwnRows(tbl.2)
  checkEquals(nrow(tbl.fixed), 2)


} # test_splitMultiple
#----------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
