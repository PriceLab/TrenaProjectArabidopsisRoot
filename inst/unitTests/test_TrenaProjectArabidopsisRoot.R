library(TrenaProjectArabidopsisRoot)
library(RUnit)
library(MotifDb)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("tp")) {
   message(sprintf("--- creating instance of TrenaProjectArabidopsisRoot"))
   tp <- TrenaProjectArabidopsisRoot();
   }

if(!exists("tbl.names"))
   tbl.names <- get(load("~/github/TrenaProjectArabidopsisRoot/inst/extdata/misc/geneIdMap.RData"))

if(!exists("tbl.geneInfo"))
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_constructor()
   test_getGeneNames()
   test_supportedGenes()
   test_getTranscriptionFactors()

   test_footprintDatabases()
   test_expressionMatrices()
   test_setTargetGene()
   test_canonicalizeName()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_constructor <- function()
{
   message(sprintf("--- test_constructor"))

   checkTrue(all(c("TrenaProjectArabidopsisRoot", "TrenaProject") %in% is(tp)))

} # test_constructor
#------------------------------------------------------------------------------------------------------------------------
test_getGeneNames <- function()
{
   message(sprintf("--- test_getGeneNames"))

   symbol <- "WBC19"
   orf <- "AT3G55130"

   x <- getGeneNames(tp, symbol)
   checkEquals(x$symbol, symbol)
   checkEquals(x$orf, orf)

   x <- getGeneNames(tp, orf)
   checkEquals(x$symbol, symbol)
   checkEquals(x$orf, orf)

   x <- getGeneNames(tp, "bogus")
   checkTrue(all(is.null(x)))

} # test_getGeneNames
#------------------------------------------------------------------------------------------------------------------------
test_supportedGenes <- function()
{
   message(sprintf("--- test_supportedGenes"))

   subset.expected <- c("WBC19")
   checkTrue(all(subset.expected %in% getSupportedGenes(tp)))

} # test_supportedGenes
#------------------------------------------------------------------------------------------------------------------------
test_footprintDatabases <- function()
{
   message(sprintf("--- test_footprintDatabases"))

   expected <- c()
   checkTrue(is.na(getFootprintDatabaseNames(tp)))
   checkTrue(is.na(getFootprintDatabaseHost(tp)))

} # test_footprintDatabases
#------------------------------------------------------------------------------------------------------------------------
test_expressionMatrices <- function()
{
   message(sprintf("--- test_expressionMatrices"))

   expected <- c("munoz.zinc.21201x42")
   checkTrue(all(expected %in% getExpressionMatrixNames(tp)))

   mtx <- getExpressionMatrix(tp, expected[1])
   checkEquals(dim(mtx), c(21201, 42))

} # test_expressionMatrices
#------------------------------------------------------------------------------------------------------------------------
# setting the target gene implies a few other assignements, all tested here:
#   geneInfo (temporarily also masquerading at tbl.transcripts
#   geneRegion
#   geneEnhancersRegion (when avaialable, defaults to geneRegion)
#
test_setTargetGene <- function()
{
   message(sprintf("--- test_setTargetGene"))

      # popularly called WBC19, with TAIR geneID AT3G55130, the tair10 gtf file calls it ABCG19
   setTargetGene(tp, "ABCG19")
   checkEquals(getTargetGene(tp), "ABCG19")

   message(sprintf("    transcripts"))
   tbl.transcripts <- getTranscriptsTable(tp)
   checkTrue(nrow(tbl.transcripts) >= 1)
   checkEquals(tbl.transcripts$chr, "3")
   checkEquals(tbl.transcripts$start, 20433839)
   checkEquals(tbl.transcripts$end , 20436453)
   checkEquals(tbl.transcripts$tss, 20436453)
   checkEquals(tbl.transcripts$strand, -1)

   message(sprintf("    geneRegion"))
   region <- getGeneRegion(tp, flankingPercent=0)
   checkTrue(all(c("chromLocString", "chrom", "start", "end") %in% names(region)))
   checkEquals(region$chromLocString, "3:20433839-20436453")

} # test_setTargetGene
#------------------------------------------------------------------------------------------------------------------------
# setting the target gene implies a few other assignements, all tested here:
#   geneInfo (temporarily also masquerading at tbl.transcripts
#   geneRegion
#   geneEnhancersRegion (when avaialable, defaults to geneRegion)
#
test_getTranscriptionFactors <- function()
{
   message(sprintf("--- test_getTranscriptionFactors"))

   tfs.go <- getAllTranscriptionFactors(tp, "Gene Ontology")
   tfs.mdb <- getAllTranscriptionFactors(tp, "MotifDb")

   checkTrue(length(tfs.go) >  1650)

   checkEquals(length(tfs.go), length(grep("^AT", tfs.go)))  # all should start with AT
   checkEquals(length(tfs.mdb), length(grep("^AT", tfs.mdb)))  # all should start with AT

   checkTrue(length(tfs.mdb) > 380)  # just from jaspar2018

   in.both <-  length(intersect(tfs.go, tfs.mdb))
   checkTrue(in.both > 350)
   checkTrue(in.both < 380)

} # test_getTranscriptionFactors
#------------------------------------------------------------------------------------------------------------------------
test_buildSingleGeneModel <- function()
{
   printf("--- test_buildSingleGeneModel")
   setTargetGene(tp, "AT3G55130")
   tbl.transcripts <- getTranscriptsTable(tp)
   getExpressionMatrixNames(tp)
   mtx <- getExpressionMatrix(tp, "munoz.zinc.19121x42.orfRowNames")
   tfs <- getAllTranscriptionFactors(tp)
   length(tfs)
   length(intersect(tfs, rownames(mtx))) # 382

} # test_buildSingleGeneModel
#------------------------------------------------------------------------------------------------------------------------
test_canonicalizeName <- function()
{
   message(sprintf("--- test_canonicalizeName"))
   mtx.munoz <- get(load(system.file(package="TrenaProjectArabidopsisRoot", "extdata", "expression", "munoz.zinc.21201x42.RData")))
   head(rownames(mtx.munoz))
      # some gene sybols from mtx.munoz, some (which start with "Arth") which fail to canonicalize
   checkEquals(canonicalizeName(tp, "ORF25"), "ATMG00640")
   checkEquals(canonicalizeName(tp, "NAD4L"), "ATMG00650")
   checkEquals(canonicalizeName(tp, "ATWBC19"), "AT3G55130")
   checkEquals(canonicalizeName(tp, "ArthMp060"),  "ArthMp060")

      # from MotifDb: head(mcols(query(MotifDb, c("thaliana", "jaspar2018")))$geneSymbol)

   checkEquals(canonicalizeName(tp, "AGL3"), "AT2G03710")
   checkEquals(canonicalizeName(tp, "AG"), "AT4G18960")
   checkEquals(canonicalizeName(tp, "HAT5"), "AT3G01470")
   checkEquals(canonicalizeName(tp, "ATHB-5"), "AT5G65310")
   checkEquals(canonicalizeName(tp, "ARR10"), "AT4G31920")
   checkEquals(canonicalizeName(tp, "AGL15"), "AT5G13790")

} # test_canonicalizeName
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
