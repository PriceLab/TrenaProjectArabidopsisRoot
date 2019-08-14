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

# tbl.model:
#         gene   betaLasso   lassoPValue pearsonCoeff   rfScore   betaRidge spearmanCoeff bindingSites symbol
# 5  AT3G59060  0.25182969 7.586023e-257    0.6738150 1301.6266  0.17793969     0.5077358           NA   PIF5
# 6  AT4G00050  0.25848974 1.790158e-129    0.5765029  738.1581  0.23987441     0.4916015           NA  UNE10
# 7  AT4G36540  0.20257033  1.016042e-88    0.5333334  701.2759  0.19764529     0.4806764           NA   BEE2
# 10 AT5G62165 -0.09238406  6.961638e-41   -0.4991596  513.0297 -0.11703297    -0.3849297           NA  AGL42
# 9  AT5G60200 -0.10874331  4.796574e-21   -0.3877858  347.2974 -0.14396653    -0.3551172           NA DOF5.3
# 8  AT5G18830 -0.25172200  6.186848e-25   -0.3776184  331.5325 -0.28773203    -0.3622195           NA   SPL7
# 2  AT1G53160  0.06447171  4.696022e-21    0.3917817  298.7136  0.09052672     0.3554472           NA   SPL4
# 1  AT1G09530  0.01684827  1.433186e-05    0.4333553  256.6698  0.06233615     0.3551013           NA   PIF3
# 4  AT3G15270  0.07932286  2.032701e-26    0.4072103  243.6807  0.10121717     0.3111859           NA   SPL5
# 3  AT1G80840  0.00000000  9.856132e-01    0.1509693  175.8583  0.01323395     0.1306107           NA WRKY40


          #------------------------------------------------------------------------------
          # second model: 16 TFs with >= 95% motif match within 500kb (each) promoter.
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

# tbl.model.2
#         gene   betaLasso   lassoPValue pearsonCoeff   rfScore   betaRidge spearmanCoeff bindingSites    symbol
# 5  AT4G14465 -0.46401927 1.919310e-162   -0.5628267 1331.0327 -0.35798325    -0.3313744           NA     AHL20
# 10 AT5G62165 -0.26067798 8.382164e-104   -0.4991596  871.5585 -0.22696959    -0.3849297           NA     AGL42
# 9  AT5G60200 -0.15703209  1.020370e-22   -0.3877858  573.3732 -0.18761945    -0.3551172           NA    DOF5.3
# 6  AT5G13080 -0.13021371  1.794286e-21   -0.3233200  352.6772 -0.13125359    -0.2561005           NA    WRKY75
# 7  AT5G46350  0.26386017  1.382743e-16    0.2225107  333.4381  0.25837985     0.1982006           NA     WRKY8
# 4  AT2G46130 -0.06135211  2.192553e-12   -0.3191933  326.8847 -0.08656127    -0.2546452           NA    WRKY43
# 8  AT5G49520 -0.05820226  6.901810e-05   -0.2456283  322.7049 -0.11555394    -0.2372729           NA    WRKY48
# 1  AT1G64620  0.05390654  1.121712e-02   -0.1025043  269.3130  0.10217146    -0.1231658           NA AT1G64620
# 3  AT2G37590  0.00000000  1.613487e-01   -0.2068212  250.8909 -0.09618521    -0.1980168           NA    DOF2.4
# 2  AT1G80840  0.02794550  1.308692e-02    0.1509693  242.7149  0.08498833     0.1306107           NA    WRKY40


          #------------------------------------------------------------------------------
          # third model: 36 TFs with >= 95% motif match within 1000kb (each) promoter.
          #------------------------------------------------------------------------------
tbl.atac <- data.frame(chrom=rep("Chr3", 2),
                       start=c(20434500,20436000),
                       end=c(20435500,20437000),
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
x3 <- build(builder)
lapply(x3, dim)
tbl.model.3 <- x2$model[1:10,]
tf.symbols <- unlist(lapply(tbl.model.3$gene, function(tf) getGeneNames(tp, tf)$symbol))
tbl.model.3$symbol <- tf.symbols

# tbl.model.3
#         gene   betaLasso   lassoPValue pearsonCoeff   rfScore   betaRidge spearmanCoeff bindingSites    symbol
# 5  AT4G14465 -0.46401927 1.919310e-162   -0.5628267 1331.0327 -0.35798325    -0.3313744           NA     AHL20
# 10 AT5G62165 -0.26067798 8.382164e-104   -0.4991596  871.5585 -0.22696959    -0.3849297           NA     AGL42
# 9  AT5G60200 -0.15703209  1.020370e-22   -0.3877858  573.3732 -0.18761945    -0.3551172           NA    DOF5.3
# 6  AT5G13080 -0.13021371  1.794286e-21   -0.3233200  352.6772 -0.13125359    -0.2561005           NA    WRKY75
# 7  AT5G46350  0.26386017  1.382743e-16    0.2225107  333.4381  0.25837985     0.1982006           NA     WRKY8
# 4  AT2G46130 -0.06135211  2.192553e-12   -0.3191933  326.8847 -0.08656127    -0.2546452           NA    WRKY43
# 8  AT5G49520 -0.05820226  6.901810e-05   -0.2456283  322.7049 -0.11555394    -0.2372729           NA    WRKY48
# 1  AT1G64620  0.05390654  1.121712e-02   -0.1025043  269.3130  0.10217146    -0.1231658           NA AT1G64620
# 3  AT2G37590  0.00000000  1.613487e-01   -0.2068212  250.8909 -0.09618521    -0.1980168           NA    DOF2.4
# 2  AT1G80840  0.02794550  1.308692e-02    0.1509693  242.7149  0.08498833     0.1306107           NA    WRKY40

# side note: AGL42 is within the top 5 among all 3 models (repressor)
# another note: AGL40 is the first one among models 2 & 3


#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
  runTests()

