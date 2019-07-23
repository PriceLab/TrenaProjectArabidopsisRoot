library(RUnit)
#----------------------------------------------------------------------------------------------------
if(!exists("tbl.test")){
  print("fresh load of our data file")
    tbl.test <- get(load(file= "10by10rootTable.RData"))
    }
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
  test_obsoleteSplitMultipleNames()
  test_splitMultipleNames()
  test_demo.doAllRows()

} # runTests
#----------------------------------------------------------------------------------------------------
obsolete.splitMultipleNamesIntoTheirOwnRows <- function(tbl)
{
  stopifnot(nrow(tbl) == 1)    # we only handle data.frames with one row

  if(!grepl(";", rownames(tbl))){
     return(tbl)  # just give it back to the caller
     }

  tokens <- strsplit(rownames(tbl), ";")[[1]]
  tbl[rep(seq_len(nrow(tbl)), each=length(tokens)),]


  if(length(tokens) == 2)
     tbl.result <- rbind(tbl, tbl)

  if(length(tokens) == 3)
     tbl.result <- rbind(tbl, tbl, tbl)

  if(length(tokens) == 4)
     tbl.result <- rbind(tbl, tbl, tbl, tbl)

  rownames(tbl.result) <- tokens

  return(tbl.result)

} # obsolete.splitMultipleNamesIntoTheirOwnRows
#----------------------------------------------------------------------------------------------------
splitMultipleNames <- function(tbl)
{
  stopifnot(nrow(tbl) == 1)    # we only handle data.frames with one row
  tokens <- strsplit(rownames(tbl), ";")[[1]]
  tbl.result <- tbl[rep(seq_len(nrow(tbl)), each=length(tokens)),]
  rownames(tbl.result) <- tokens

  return(tbl.result)

} # splitMultipleNames
#----------------------------------------------------------------------------------------------------
test_obsoleteSplitMultipleNames <- function()
{
  message(sprintf("--- test_obsoleteSplitMultipleNames"))   # announce what test you are running

     #------------------------------------------------------------------------
     #  do the simplest test first: a row with a simple name (no semicolons)
     #------------------------------------------------------------------------

  tbl.1 <- tbl.test[1,]
  tbl.fixed <- splitMultipleNames(tbl.1)
  checkEquals(nrow(tbl.fixed), 1)
  checkEquals(rownames(tbl.fixed), rownames(tbl.1))

  tbl.2 <- tbl.test[10,]
  tbl.fixed <- splitMultipleNames(tbl.2)
  checkEquals(nrow(tbl.fixed), 2)
  checkEquals(rownames(tbl.fixed),c("ATMG01250", "AT2G07697"))

  tbl.3 <- tbl.1
  rownames(tbl.3) <- "aaa;bbb;ccc"
  tbl.fixed <- splitMultipleNames(tbl.3)
  checkEquals(nrow(tbl.fixed), 3)
  checkEquals(rownames(tbl.fixed),c("aaa", "bbb", "ccc"))

  tbl.4 <- tbl.1
  rownames(tbl.4) <- "abc;def;ghi;jkl"
  tbl.fixed <- splitMultipleNames(tbl.4)
  checkEquals(nrow(tbl.fixed), 4)
  checkEquals(rownames(tbl.fixed),c("abc", "def", "ghi", "jkl"))

} # test_obsoleteSplitMultipleNames
#----------------------------------------------------------------------------------------------------
test_splitMultipleNames <- function()
{
  message(sprintf("--- test_splitMultipleNames"))   # announce what test you are running

  tbl.5 <- tbl.test[10,]
  tbl.fixed <- splitMultipleNames(tbl.5)
  checkEquals(nrow(tbl.fixed), 2)
  checkEquals(rownames(tbl.fixed), c("ATMG01250","AT2G07697"))

} #test_splitMultipleNames
#----------------------------------------------------------------------------------------------------
demo.doAllRows <- function(tbl)
{
   tbl <- get(load("10by10rootTable.RData"))
   x <- lapply(1:nrow(tbl), function(i) splitMultipleNames(tbl[i,]))
   tbl.combined <- do.call(rbind, x)
   return(tbl.combined)

} # demo.doAllRows
#------------------------------------------------------------------------------------------------------------------------
test_demo.doAllRows <- function()
{
   message(sprintf("---test_demo.doAllRows")) #announce what you're running
   tbl.expanded <- demo.doAllRows(tbl.test)

   checkTrue(nrow(tbl.expanded) > nrow(tbl.test))

   multiple.orf.names.orig <- length(grep(";", rownames(tbl.test)))
   checkEquals(multiple.orf.names.orig, 3)

   multiple.orf.names.expanded <- length(grep(";", rownames(tbl.expanded)))
   checkEquals(multiple.orf.names.expanded, 0)

} # test_demo.doAllRows
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()

