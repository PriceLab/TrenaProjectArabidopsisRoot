library(GSEABase)
library(GOstats)
library(GO.db)
library(Category)
library(trena)
library(trenaSGM)
library(TrenaProjectArabidopsisRoot)
library(org.At.tair.db)
#------------------------------------------------------------------------------------------------------------------------
geneOntology <- function(orf)
{
   tbl.go <- select(org.At.tair.db, keytype="TAIR", keys="AT3G59060", columns="GO")
   return(select(GO.db, keys=tbl.go$GO, columns=c("TERM"))$TERM)

} # geneOntology
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
  test_geneOntology()
  test_goEnrichment()
  
} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_geneOntology <- function()
{
   message(sprintf("--- test_geneOntology"))
   geneSymbol <- "PIF5"
   orfName <- getGeneNames(tp, geneSymbol)$orf
   tbl.go <- geneOntology(orfName)

} # test_geneOntology
#------------------------------------------------------------------------------------------------------------------------
goEnrichment <- function(orfNames, maxCategories)
{
  gene.universe = character(0)

  go.params <- new("GOHyperGParams", geneIds=orfNames,
                   universeGeneIds=gene.universe, annotation = "org.At.tair.db",
                   ontology = 'BP', pvalueCutoff = 0.01, conditional = FALSE,
                   testDirection = "over")
  message(sprintf("about to calculate GO enrichment"))

  go.bp.hgr <- hyperGTest(go.params)
  tbl.go <- summary(go.bp.hgr)
  if(nrow(tbl.go) > maxCategories)
     tbl.go <- head(tbl.go, n=maxCategories)

  message(sprintf("about to look up genes for each of %d GO categories", nrow(tbl.go)))

  suppressMessages(
    geneSymbols <- lapply(tbl.go$GOBPID,
                       function(goTerm){
                          if(!goTerm %in% keys(org.At.tair.db, keytype="GOALL")) {
                             printf("failed goTerm: %s", goTerm)
                             return ("")
                             }
                          keepers <- intersect(orfNames,select(org.At.tair.db, keys=goTerm, keytype="GOALL", columns="TAIR")$TAIR)
                          keeper.geneSymbols <- unlist(lapply(keepers, function(orf) getGeneNames(tp,orf)$symbol))
                          paste(keeper.geneSymbols, collapse=";")
                          }))

  tbl.go$genes <- unlist(geneSymbols)

  tbl.go

} # goEnrichment
#------------------------------------------------------------------------------------------------------------------------
test_goEnrichment <- function()
{
   message(sprintf("--- test_goEnrichment"))
   orfNames <- c("AT3G55130", "AT3G59060", "AT4G00050", "AT4G36540", "AT5G62165", "AT5G18830","AT5G60200",
                 "AT1G53160", "AT1G09530", "AT3G15270","AT1G80840")

   tbl.enrichment <- goEnrichment(orfNames, maxCategories=5)  # 10 is more typical

} # test_goEnrichment
#------------------------------------------------------------------------------------------------------------------------
if(!exists("tp")) {
   message(sprintf("--- creating instance of TrenaProjectArabidopsisRoot"))
   tp <- TrenaProjectArabidopsisRoot();
   }

genome <- "tair10"
getExpressionMatrixNames(tp)
mtx <- getExpressionMatrix(tp, "aluru.18617x1938.orfRowNames")
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
               tfPrefilterCorrelation=0.1,
               annotationDbFile=dbfile(org.At.tair.db),
               orderModelByColumn="rfScore",
               solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"),
               quiet=TRUE)

builder <- NoDnaModelBuilder(genome, targetGene, recipe, quiet=TRUE)
x <- build(builder)
lapply(x, dim)
tbl.model <- x$model[1:10,]
tf.symbols <- unlist(lapply(tbl.model$gene, function(tf) getGeneNames(tp, tf)$symbol))
tbl.model$symbol <- tf.symbols

#tbl.model.2:
# gene   betaLasso   lassoPValue pearsonCoeff   rfScore   betaRidge spearmanCoeff bindingSites symbol
# 5  AT3G59060  0.25196238 7.586023e-257    0.6738150 1341.7531  0.17726536     0.5077358           NA   PIF5
# 6  AT4G00050  0.25841300 1.790158e-129    0.5765029  752.1774  0.23942232     0.4916015           NA  UNE10
# 7  AT4G36540  0.20233725  1.016042e-88    0.5333334  695.9929  0.19719870     0.4806764           NA   BEE2
# 10 AT5G62165 -0.09225101  6.961638e-41   -0.4991596  477.2526 -0.11708188    -0.3849297           NA  AGL42
# 8  AT5G18830 -0.25105221  6.186848e-25   -0.3776184  336.8471 -0.28733313    -0.3622195           NA   SPL7
# 9  AT5G60200 -0.10839858  4.796574e-21   -0.3877858  332.8228 -0.14395389    -0.3551172           NA DOF5.3
# 2  AT1G53160  0.06427679  4.696022e-21    0.3917817  289.1844  0.09060847     0.3554472           NA   SPL4
# 1  AT1G09530  0.01668285  1.433186e-05    0.4333553  256.7569  0.06262309     0.3551013           NA   PIF3
# 4  AT3G15270  0.07912486  2.032701e-26    0.4072103  241.1361  0.10120961     0.3111859           NA   SPL5
# 3  AT1G80840  0.00000000  9.856132e-01    0.1509693  176.1332  0.01338915     0.1306107           NA WRKY40


          #------------------------------------------------------------------------------
          # second model: 16 TFs with >= 95% motif match within 500kb promoter.
          # ATAC-seq about 2k upstream of WBC19
          #------------------------------------------------------------------------------

tbl.atac <- data.frame(chrom=rep("Chr3", 2),
                       start=c(20437000,20434000),
                       end=c(20437500,20434500),
                       stringsAsFactors=FALSE)
tbl.tfs <- findCandidateTranscriptionFactorsByMotifInSequence(tp, tbl.atac, 95L)
candidate.tfs <- sort(unique(tbl.tfs$orf))
length(candidate.tfs)

recipe <- list(title="WBC19",
               type="noDNA.tfsSupplied",
               matrix=mtx,
               candidateTFs=candidate.tfs,
               tfPool=getAllTranscriptionFactors(tp, "MotifDb"),
               tfPrefilterCorrelation=0.1,
               annotationDbFile=dbfile(org.At.tair.db),
               orderModelByColumn="rfScore",
               solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"),
               quiet=TRUE)

builder <- NoDnaModelBuilder(genome, targetGene, recipe, quiet=TRUE)
x2 <- build(builder)
lapply(x2, dim)
tbl.model.2 <- x2$model[1:10,]
tf.symbols <- unlist(lapply(tbl.model.2$gene, function(tf) getGeneNames(tp, tf)$symbol))
tbl.model.2$symbol <- tf.symbols
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
  runTests()

