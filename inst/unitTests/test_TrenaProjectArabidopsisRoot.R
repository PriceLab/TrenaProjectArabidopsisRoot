library(TrenaProjectArabidopsisRoot)
library(RUnit)
library(MotifDb)
library(trenaSGM)
library(org.At.tair.db)
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
   test_findCandidateTranscriptionFactorsByMotifInSequence()

   test_footprintDatabases()
   test_expressionMatrices()
   test_setTargetGene()
   test_canonicalizeName()
   test_buildSingleGeneModel_WBC19()

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
test_findCandidateTranscriptionFactorsByMotifInSequence <- function()
{
   message(sprintf("--- test_findCandidateTranscriptionFactorsByMotifInSequence"))

   start.loc <- 20438600
   end.loc   <- start.loc + 30
   tbl.region <- data.frame(chrom="Chr3", start=start.loc, end=end.loc, stringsAsFactors=FALSE)
   tbl.tfs <- findCandidateTranscriptionFactorsByMotifInSequence(tp, tbl.region, 95L)
   checkEquals(dim(tbl.tfs), c(2, 14))
   checkEquals(tbl.tfs$orf, c("AT5G62165", "AT1G80840"))

   tbl.tfs <- findCandidateTranscriptionFactorsByMotifInSequence(tp, tbl.region, 80L)
   checkEquals(dim(tbl.tfs), c(37, 14))
   checkEquals(length(unique(tbl.tfs$orf)), 30)
   checkTrue(all(c("AT5G62165", "AT1G80840") %in% tbl.tfs$orf))

} # test_findCandidateTranscriptionFactorsByMotifInSequence
#------------------------------------------------------------------------------------------------------------------------
test_buildSingleGeneModel_WBC19 <- function()
{
   message(sprintf("--- test_buildSingleGeneModel_WBC19"))

   genome <- "tair10"
   mtx <- getExpressionMatrix(tp, "munoz.zinc.18993x42.canonicalRowNames")

   targetGene <- canonicalizeName(tp, "WBC19")

     #--------------------------------------------------------------------------
     # first model: only 18 TFs with >= 95% motif match in short root-specific
     # ATAC-seq about 2k upstream of WBC19
     #--------------------------------------------------------------------------

   tbl.atac <- data.frame(chrom=rep("Chr3", 2),
                          start=c(20438126,20438543),
                          end=c(20438246,20438984),
                          stringsAsFactors=FALSE)
   tbl.tfs <- findCandidateTranscriptionFactorsByMotifInSequence(tp, tbl.atac, 95L)
   candidate.tfs <- sort(unique(tbl.tfs$orf))
   length(candidate.tfs)

   recipe <- list(title="WBC19",
                  type="noDNA.tfsSupplied",
                  matrix=mtx,
                  candidateTFs=candidate.tfs,
                  tfPool=getAllTranscriptionFactors(tp, "MotifDb"),
                  tfPrefilterCorrelation=0.2,
                  annotationDbFile=dbfile(org.At.tair.db),
                  orderModelByColumn="rfScore",
                  solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"),
                  quiet=TRUE)

   builder <- NoDnaModelBuilder(genome, targetGene, recipe, quiet=TRUE)
   x <- build(builder)
   tbl.model <- x$model[1:10,]
       # add geneSymbols, easier for humans to read
   tbl.model$geneSymbol <- unlist(lapply(tbl.model$gene, function(gene) getGeneNames(tp, gene)$symbol))

   checkEquals(x$regulatoryRegions, data.frame())

     #------------------------------------------------------------------------------
     # second model: 48 TFs with >= 95% motif match within traditional 2kb promoter.
     # no ATAC-seq information used.
     #------------------------------------------------------------------------------

   tbl.promoter <- data.frame(chrom="Chr3", start=20436000, end=20438000, stringsAsFactors=FALSE)
   tbl.tfs <- findCandidateTranscriptionFactorsByMotifInSequence(tp, tbl.promoter, 95L)
   candidate.tfs <- sort(unique(tbl.tfs$orf))
   length(candidate.tfs)

     #--------------------------------------------------------------------------
     # just one change needed: substitute in the new largner number of TFs
     #--------------------------------------------------------------------------

   recipe$candidateTFs=candidate.tfs
   builder <- NoDnaModelBuilder(genome, targetGene,  recipe, quiet=TRUE)
   x2 <- build(builder)
   tbl.model.2 <- x2$model
   tbl.model.2$geneSymbol <- unlist(lapply(tbl.model.2$gene, function(gene) getGeneNames(tp, gene)$symbol))

   dim(tbl.model.2)

} # test_buildSingleGeneModel_WBC19
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
