library(RUnit)
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_splitNames()
   test_eliminateDuplicateRows()
   test_expand.single.row.when.needed()
   test_do.it.all()
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
old.splitNames <- function(string)
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
splitNames <- function(names)
{
   unlist(strsplit(names, ";"))

} # splitNames
#---------------------------------------------------------------------------------------------------
test_splitNames <- function()
{
   message(sprintf("--- test_splitNames"))
   checkEquals(splitNames("abc"), "abc")

   checkEquals(splitNames("abcd;efgh"), c("abcd", "efgh"))

   checkEquals(splitNames(c("abc", "abcd;efgh")), c("abc","abcd", "efgh"))

} # test_splitNames
#---------------------------------------------------------------------------------------------------
test_old.splitNames <- function()
{
   message(sprintf("--- test_old.splitNames"))

   checkEquals(old.splitNames("abc"), "abc")
   checkEquals(old.splitNames(""), "")
   #checking for no splits

   result <- old.splitNames("abcd;efgh")
   checkEquals(result, c("abcd", "efgh"))
   #checking for 1 splits

   result.2 <- old.splitNames("abcd;efgh;ijkl")
   checkEquals(result.2, c("abcd", "efgh", "ijkl"))
   #checking for 2 splits

   result.3 <- old.splitNames("abcd;efgh;ijkl;mnop")
   checkEquals(result.3, c("abcd", "efgh", "ijkl", "mnop"))
   #checking for 3 splits

   result.4 <- old.splitNames("abcd;efgh;ijkl;mnop;qrst")
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
expand.single.row.when.needed <- function(x)
{
   browser()
   stopifnot(nrow(x)==1)
   new.row.names <- splitNames(rownames(x[,,drop=FALSE]))
   total.row.count <- length(new.row.names)

   if(total.row.count == 1)
     return(x)

    # we get here only if total.row.count > 1

   tbl.expanded <- matrix(rep(t(x), total.row.count), ncol = ncol(x) , byrow = TRUE )
   colnames(tbl.expanded) <- colnames(x)
   rownames(tbl.expanded) <- new.row.names

   return(tbl.expanded)

} # expand.single.row.when.needed
#----------------------------------------------------------------------------------------------------
test_expand.single.row.when.needed <- function()
{
   message(sprintf("--- test_expand.single.row.when.needed"))
   stopifnot(exists("test.tbl"))
   tbl.fixed <- eliminateDuplicateRows(test.tbl)
   tbl.tiny <- as.matrix(tbl.fixed[1:6, 1:5])

       # row 2 of tbl.tiny has a simple single orf name, ATMG00840
   tbl.expanded2 <- expand.single.row.when.needed(tbl.tiny[2,,drop=FALSE])
   checkEquals(dim(tbl.expanded2), c(1, 5))
      # row 1 has a double orf name: ATMG00850;AT2G07682
   tbl.expanded.1 <- expand.single.row.when.needed(tbl.tiny[1,,drop=FALSE])
   checkEquals(dim(tbl.expanded.1), c(2, 5))
   checkEquals(rownames(tbl.expanded.1), c("ATMG00850", "AT2G07682"))

     # do a bunch of checks here
   checkEquals(rownames(tbl.expanded.1), splitNames(rownames(tbl.tiny[1,,drop=FALSE])))
   checkEquals(rownames(tbl.expanded2), splitNames(rownames(tbl.tiny[2,,drop=FALSE])))

     # now for a more complete and realistic test, iterate over the whole of the matrix,
     # row by row, ucalling expand.single.row.when.needed  on each row.  since apply (an
     # otherwise sensible solution) drops row names, let's  try an old-fashioned for loop

   result <- list()
   for(i in 1:nrow(tbl.tiny)){
      result[[i]] <- expand.single.row.when.needed(tbl.tiny[i,,drop=FALSE])
      }

   mtx.expanded <- do.call(rbind, result)
   checkEquals(dim(mtx.expanded), c(8, 5))
   expected.rownames <- unlist(strsplit(rownames(tbl.tiny), ";"))
   checkEquals(rownames(mtx.expanded), expected.rownames)

} # test_expand.single.row.when.needed
#----------------------------------------------------------------------------------------------------
do.it.all <- function(d)
{
   browser()
   organizedTbl <- eliminateDuplicateRows(test.tbl)
   tbl.Expanded <- apply(organizedTbl, 1, expand.single.row.when.needed)
   return(tbl.Expanded)
}
#----------------------------------------------------------------------------------------------------
test_do.it.all <- function()
{
   browser()
   message(sprintf("--- test_expand.single.row.when.needed"))
   
   small.tbl <- tbl[1:10, 1:10]
   rownamesOfTbl <- splitNames(rownames(small.tbl))
   
   checkTrue(nrow(small.tbl) < nrow(do.it.all(small.tbl)))
   checkEquals(dim(do.it.all(small.tbl)), c(12,10))
   checkTrue(!grep(";", rownames(do.it.all(small.tbl))))
   checkTrue(rownames(do.it.all(small.tbl)), rownamesOfTbl)
   
}
#----------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
#---------------------------------------------------------------------------------------------------
# expand.single.row.when.needed(tbl)
#error occurs. it states:
#Error in dimnames(x) <- dn : 
   #length of 'dimnames' [1] not equal to array extent

#I know the apply function can only take an array or a matrix if the dimension of the array is 2.
#So I do not understand the problem if the dimensions of tbl is 2 dimensions (shown below):
# > dimnames(tbl)
# [[1]]
# [1] "ATMG00850;AT2G07682" "ATMG00840"           "ATMG00820"           "ATMG00750;AT2G07686" "ATMG00740"           "ATMG00720"          
# [7] "ATMG00710"           "ATMG00690"           "ATMG00680"           "ATMG00670"           "ATMG00660"           "ATMG00650"          
# [13] "ATMG00640"           "ATMG00060"           "AT5G66230"           "AT5G65360"           "AT5G65110"           "AT5G63420"          
# [19] "AT5G63320"           "AT5G62440"           "AT5G62000"           "AT5G60310"           "AT5G57670"           "AT5G57260"          
# [25] "AT5G57150"           "AT5G54130"           "AT5G50740"           "AT5G49680"           "AT5G48250"           "AT5G47260"          
# [31] "AT5G46470"           "AT5G45260"           "AT5G44820"           "AT5G44650"           "AT5G43640"           "AT5G43230"          
# [37] "AT5G41240"           "AT5G41080"           "AT5G39785"           "AT5G39130"           "AT5G37590"           "AT5G32540"          
# [43] "AT5G30420"           "AT5G27410"           "AT5G24690"           "AT5G20540"           "AT5G19800"           "AT5G16980"          
# [49] "AT5G14180"           "AT5G12170"           "AT5G08230"           "AT5G07740"           "AT5G05860"           "AT5G04560"          
# [55] "AT5G04310"           "AT4G39480"           "AT4G39420"           "AT4G37950"           "AT4G34830"           "AT4G33520"          
# [61] "AT4G33330"           "AT4G33180"           "AT4G32700"           "AT4G31230"           "AT4G31030"           "AT4G28130"          
# [67] "AT4G26630"           "AT4G25450"           "AT4G24860"           "AT4G24350"           "AT4G20840"           "AT4G20760"          
# [73] "AT4G20130"           "AT4G16020"           "AT4G15810"           "AT4G15590"           "AT4G15230"           "AT4G15200"          
# [79] "AT4G15180"           "AT4G13750"           "AT4G13150"           "AT4G11830"           "AT4G10955;AT4G10960" "AT4G09900"          
# [85] "AT4G08410"           "AT4G07920"           "AT4G07600"           "AT4G00800"           "AT4G00630"           "AT4G00400"          
# [91] "AT4G00280"           "AT4G00040"           "AT4G00020"           "AT3G61010"           "AT3G60420"           "AT3G59570"          
# [97] "AT3G55890;AT3G55910" "AT3G54320"           "AT3G53130"           "AT3G51290"           "AT3G49080"           "AT3G47980"          
# [103] "AT3G47910"           "AT3G44510"           "AT3G44500;AT3G47260" "AT3G43020"           "AT3G42360"           "AT3G42220"          
# [109] "AT3G30260"           "AT3G29510"           "AT3G27810"           "AT3G27785"           "AT3G27120"           "AT3G26980"          
# [115] "AT3G26300"           "AT3G25820;AT3G25830" "AT3G24982"           "AT3G23510;AT3G23530" "AT3G22430"           "AT3G20960"          
# [121] "AT3G20300"           "AT3G20210"           "AT3G19720"           "AT3G18930"           "AT3G17998;AT3G18000" "AT3G17970"          
# [127] "AT3G15550"           "AT3G15520"           "AT3G15240"           "AT3G13230"           "AT3G13160"           "AT3G12800"          
# [133] "AT3G10180"           "AT3G09600"           "AT3G09550"           "AT3G09260"           "AT3G07250"           "AT3G06810"          
# [139] "AT3G06520"           "AT3G06433"           "AT3G04750"           "AT3G04710"           "AT2G48060"           "AT2G47730"          
# [145] "AT2G47300"           "AT2G46340"           "AT2G45960"           "AT2G40030"           "AT2G39140"           "AT2G38940"          
# [151] "AT2G38840"           "AT2G37840"           "AT2G35710"           "AT2G35570"           "AT2G35450"           "AT2G35390"          
# [157] "AT2G32560"           "AT2G32190;AT2G32210" "AT2G31340"           "AT2G31110"           "AT2G31010"           "AT2G30600"          
# [163] "AT2G30520"           "AT2G27900"           "AT2G27395"           "AT2G25870"           "AT2G22540"           "AT2G22010"          
# [169] "AT2G20500"           "AT2G20210"           "AT2G20100"           "AT2G20050"           "AT2G19360"           "AT2G18670;AT2G18680"
# [175] "AT2G13890"           "AT2G13790"           "AT2G13680"           "AT2G13230"           "AT2G10630"           "AT2G10400"          
# [181] "AT2G10140"           "AT2G10070"           "AT2G07783;ATMG00830" "AT2G07698"           "AT2G07674"           "AT2G06310"          
# [187] "AT2G06040"           "AT2G04970;AT2G06440" "AT2G03390"           "AT2G02410"           "AT2G01390"           "AT1G80320"          
# [193] "AT1G79170"           "AT1G79150"           "AT1G77960"           "AT1G77890"           "AT1G77080"           "AT1G76820"          
# [199] "AT1G74770"           "AT1G73805"           "AT1G73090"           "AT1G73060;AT1G73066" "AT1G72410"           "AT1G70920"          
# [205] "AT1G70250"           "AT1G68010"           "AT1G65680"           "AT1G65540"           "AT1G65450"           "AT1G65390"          
# [211] "AT1G64940"           "AT1G64100"           "AT1G61210"           "AT1G60200"           "AT1G56500"           "AT1G54470"          
# [217] "AT1G54180"           "AT1G52240"           "AT1G51270"           "AT1G51110"           "AT1G50240"           "AT1G47480"          
# [223] "AT1G42970"           "AT1G42650"           "AT1G42605;AT4G10580" "AT1G33680"           "AT1G31710"           "AT1G31690"          
# [229] "AT1G31050"           "AT1G29600"           "AT1G26300"           "AT1G23640"           "AT1G23440"           "AT1G23380"          
# [235] "AT1G22330"           "AT1G21830"           "AT1G21120"           "AT1G19520"           "AT1G19025"           "AT1G18390"          
# [241] "AT1G16040"           "AT1G14110"           "AT1G13370"           "AT1G11130"           "AT1G10090"           "AT1G07870"          
# [247] "AT1G07730;AT1G07740" "AT1G05490"           "AT1G05400"           "AT1G04380"           "AT1G04050"           "AT1G03630"          
# [253] "AT1G02205"           "AT1G01320"           "257901_at"           "246357_x_at"        
# 
# [[2]]
# [1] "E.GEOD.11532_1"  "E.GEOD.11532_2"  "E.GEOD.11532_3"  "E.GEOD.11532_4"  "E.GEOD.11532_5"  "E.GEOD.11532_6"  "E.GEOD.13783_1"  "E.GEOD.13783_2" 
# [9] "E.GEOD.13783_3"  "E.GEOD.13783_4"  "E.GEOD.13783_5"  "E.GEOD.13783_6"  "E.GEOD.13783_7"  "E.GEOD.13783_8"  "E.GEOD.13783_9"  "E.GEOD.13783_10"
# [17] "E.GEOD.13783_11" "E.GEOD.13783_12"

