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
mtx <- load("aluru.18617x1938.orfRowNames.RData")
# set working directory path to "/Users/bioadmin/github/TrenaProjectArabidopsisRoot/inst/extdata/expression/"
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

#------------------------------------------------------------------------------------------------------------------------

if(!interactive())
  runTests()
