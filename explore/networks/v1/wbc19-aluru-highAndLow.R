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
#get matrix (set path to "/Users/bioadmin/github/TrenaProjectArabidopsisRoot/inst/extdata/expression")
mtx <- as.matrix(get(load("highAndLowMtx.RData")))
targetGene <- canonicalizeName(tp, "WBC19")

          #----------------------------------------------------------------------------------------
          # 1st model: 28 TFs with >= 95% motif match within 2000 upstream, 2000 downstream of TSS
          #----------------------------------------------------------------------------------------

tbl.atac <- data.frame(chrom=rep("Chr3", 2),
                       start=c(20432112,20434113),
                       end=c(20434112,20436113),
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
# Warning messages (occurred 7x):
#1: In Solver(mtx.assay = mtx.assay, targetGene = targetGene, candidateRegulators = candidateRegulators,  :
#               Target gene mean expression is in the bottom 10% of all genes in the assay matrix
lapply(x, dim)
tbl.model <- x$model[1:6,]
tf.symbols <- unlist(lapply(tbl.model$gene, function(tf) getGeneNames(tp, tf)$symbol))
tbl.model$symbol <- tf.symbols

# tbl.model
#         gene   betaLasso  lassoPValue pearsonCoeff  rfScore  betaRidge spearmanCoeff bindingSites  symbol
# 16 AT5G60200 -0.53112793 3.461610e-38   -0.9051483 504.5427 -0.3700104    -0.7535674           NA  DOF5.3
# 17 AT5G62165 -0.40599282 4.459571e-33   -0.8771186 290.4697 -0.2177964    -0.7197240           NA   AGL42
# 14 AT5G18830 -0.87901002 8.278885e-31   -0.8600045 261.4355 -0.6302244    -0.7223522           NA    SPL7
# 1  AT1G52150 -0.36776837 2.899774e-25   -0.8698721 245.1601 -0.4900214    -0.7418542           NA ATHB-15
# 11 AT4G14465 -0.50349276 3.488459e-37   -0.9004465 123.2128 -0.1872414    -0.7942754           NA   AHL20
# 2  AT1G53160  0.07126078 1.525717e-14    0.8273910 110.3970  0.2054830     0.7663726           NA    SPL4

          #----------------------------------------------------------------------------------------
          # 2nd model: 56 TFs with >= 95% motif match within 2000 upstream, 2000 downstream of TSS
          # AGL42 regulator
          #----------------------------------------------------------------------------------------

targetGene <- canonicalizeName(tp, "AGL42")
tbl.atac <- data.frame(chrom=rep("Chr5", 2),
                       start=c(24963074,24965075),
                       end=c(24965074,24967075),
                       stringsAsFactors=FALSE)
tbl.tfs <- findCandidateTranscriptionFactorsByMotifInSequence(tp, tbl.atac, 95L)
candidate.tfs <- sort(unique(tbl.tfs$orf))
length(candidate.tfs)

recipe <- list(title="AGL42",
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
tbl.model.2 <- x2$model[1:6,]
tf.symbols <- unlist(lapply(tbl.model.2$gene, function(tf) getGeneNames(tp, tf)$symbol))
tbl.model.2$symbol <- tf.symbols

# tbl.model.2
#         gene  betaLasso  lassoPValue pearsonCoeff   rfScore   betaRidge spearmanCoeff bindingSites  symbol
# 37 AT5G60200  0.3162493 6.952895e-17    0.8205313 102.63296  0.11654657     0.7632043           NA  DOF5.3
# 28 AT4G36540 -0.3257270 3.740934e-31   -0.8652954  70.56216 -0.08816573    -0.7547675           NA    BEE2
# 32 AT5G18830  0.0000000 9.037812e-01    0.7300341  66.18541  0.13912380     0.6837444           NA    SPL7
# 2  AT1G52150  0.0000000 6.389678e-01    0.7183179  43.26442  0.08459405     0.6831683           NA ATHB-15
# 3  AT1G53160  0.0000000 2.222368e-01   -0.7194983  32.94846 -0.05359661    -0.7004260           NA    SPL4
# 26 AT4G31920  0.4368234 2.349344e-17    0.8010088  21.53626  0.19902339     0.7508911           NA   ARR10


          #----------------------------------------------------------------------------------------
          # 3rd model: 44 TFs with >= 95% motif match within 2000 upstream, 2000 downstream of TSS
          # DOF5.3 regulator
          #----------------------------------------------------------------------------------------
targetGene <- canonicalizeName(tp, "DOF5.3")
tbl.atac <- data.frame(chrom=rep("Chr5", 2),
                       start=c(24239077,24241078),
                       end=c(24241077,24243078),
                       stringsAsFactors=FALSE)
tbl.tfs <- findCandidateTranscriptionFactorsByMotifInSequence(tp, tbl.atac, 95L)
candidate.tfs <- sort(unique(tbl.tfs$orf))
length(candidate.tfs)

recipe <- list(title="DOF5.3",
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
tbl.model.3 <- x3$model[1:6,]
tf.symbols <- unlist(lapply(tbl.model.3$gene, function(tf) getGeneNames(tp, tf)$symbol))
tbl.model.3$symbol <- tf.symbols
# tbl.model.3
#         gene   betaLasso  lassoPValue pearsonCoeff   rfScore   betaRidge spearmanCoeff bindingSites  symbol
# 23 AT4G36540 -0.12771897 2.758304e-32   -0.8746272 34.244405 -0.05377290    -0.7648965           NA    BEE2
# 3  AT1G52150  0.48330397 2.243489e-33    0.8795316 23.485641  0.14304595     0.8378278           NA ATHB-15
# 31 AT5G62165  0.06774890 7.402972e-19    0.8205313 22.529689  0.04361396     0.7632043           NA   AGL42
# 20 AT4G14465  0.05441007 9.279738e-18    0.8174130 16.886400  0.04017035     0.7240564           NA   AHL20
# 4  AT1G53160  0.00000000 2.202929e-01   -0.7424131  9.788317 -0.04145836    -0.7293969           NA    SPL4
# 27 AT5G12870  0.00000000 7.383427e-01    0.6637216  4.300935  0.03802582     0.6611581           NA   MYB46


          #----------------------------------------------------------------------------------------
          # 4th model: 55 TFs with >= 95% motif match within 2000 upstream, 2000 downstream of TSS
          # SPL7 regulator
          #----------------------------------------------------------------------------------------

targetGene <- canonicalizeName(tp, "SPL7")
tbl.atac <- data.frame(chrom=rep("Chr5", 2),
                       start=c(6274114,6276115),
                       end=c(6276114,6278115),
                       stringsAsFactors=FALSE)
tbl.tfs <- findCandidateTranscriptionFactorsByMotifInSequence(tp, tbl.atac, 95L)
candidate.tfs <- sort(unique(tbl.tfs$orf))
length(candidate.tfs)

recipe <- list(title="SPL7",
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
x4 <- build(builder)
lapply(x4, dim)
tbl.model.4 <- x4$model[1:6,]
tf.symbols <- unlist(lapply(tbl.model.4$gene, function(tf) getGeneNames(tp, tf)$symbol))
tbl.model.4$symbol <- tf.symbols

# tbl.model.4
#         gene   betaLasso  lassoPValue pearsonCoeff  rfScore   betaRidge spearmanCoeff bindingSites    symbol
# 33 AT5G65210  0.05362110 1.189573e-27    0.8390267 9.929997  0.01775844     0.7674287           NA      TGA1
# 24 AT4G36540 -0.02607314 5.176391e-16   -0.8141420 4.361003 -0.01816434    -0.7240924           NA      BEE2
# 30 AT5G60200  0.00000000 7.215740e-01    0.7579193 4.320408  0.01970143     0.7203720           NA    DOF5.3
# 18 AT3G57600  0.00000000 6.241473e-03   -0.7424516 4.078474 -0.01825572    -0.7310171           NA AT3G57600
# 31 AT5G62165  0.00000000 7.632455e-01    0.7300341 2.742569  0.01123549     0.6837444           NA     AGL42
# 3  AT1G52150  0.22609869 2.425951e-25    0.8272799 2.652548  0.05628861     0.7957876           NA   ATHB-15

          #----------------------------------------------------------------------------------------
          # 5th model: 43 TFs with >= 95% motif match within 2000 upstream, 500 downstream of TSS
          # ATHB-15 regulator
          #----------------------------------------------------------------------------------------

targetGene <- canonicalizeName(tp, "ATHB-15")
tbl.atac <- data.frame(chrom=rep("Chr1", 2),
                       start=c(19407912,19409913),
                       end=c(19409912,19410413),
                       stringsAsFactors=FALSE)
tbl.tfs <- findCandidateTranscriptionFactorsByMotifInSequence(tp, tbl.atac, 95L)
candidate.tfs <- sort(unique(tbl.tfs$orf))
length(candidate.tfs)

recipe <- list(title="ATHB-15",
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
x5 <- build(builder)
lapply(x5, dim)
tbl.model.5 <- x5$model[1:6,]
tf.symbols <- unlist(lapply(tbl.model.5$gene, function(tf) getGeneNames(tp, tf)$symbol))
tbl.model.5$symbol <- tf.symbols

# tbl.model.5
#         gene betaLasso  lassoPValue pearsonCoeff   rfScore   betaRidge spearmanCoeff bindingSites symbol
# 24 AT5G60200 0.3302012 2.243489e-33    0.8795316 17.506781  0.07546730     0.8378278           NA DOF5.3
# 25 AT5G62165 0.0000000 5.187714e-01    0.7183179  8.959486  0.02015833     0.6831683           NA  AGL42
# 23 AT5G59780 0.1094763 7.736923e-18    0.7672236  7.036857  0.05384100     0.7766937           NA  MYB59
# 2  AT1G53160 0.0000000 5.093381e-01   -0.7128332  6.654776 -0.03047347    -0.7441584           NA   SPL4
# 18 AT5G12870 0.0533393 6.776370e-11    0.7030407  4.649051  0.05135668     0.7241164           NA  MYB46
# 8  AT2G42200 0.0000000 5.205207e-01   -0.6571580  2.745416 -0.03253889    -0.6369757           NA   SPL9


          #----------------------------------------------------------------------------------------
          # 6th model: 58 TFs with >= 95% motif match within 2000 upstream, 500 downstream of TSS
          # AHL20 regulator
          #----------------------------------------------------------------------------------------

targetGene <- canonicalizeName(tp, "AHL20")
tbl.atac <- data.frame(chrom=rep("Chr4", 2),
                       start=c(8318970,8320971),
                       end=c(8320970,8322971),
                       stringsAsFactors=FALSE)
tbl.tfs <- findCandidateTranscriptionFactorsByMotifInSequence(tp, tbl.atac, 95L)
candidate.tfs <- sort(unique(tbl.tfs$orf))
length(candidate.tfs)

recipe <- list(title="AHL20",
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
x6 <- build(builder)
lapply(x6, dim)
tbl.model.6 <- x6$model[1:6,]
tf.symbols <- unlist(lapply(tbl.model.6$gene, function(tf) getGeneNames(tp, tf)$symbol))
tbl.model.6$symbol <- tf.symbols

# tbl.model.6
#         gene  betaLasso  lassoPValue pearsonCoeff   rfScore   betaRidge spearmanCoeff bindingSites  symbol
# 37 AT5G65210  0.2674667 8.574817e-38    0.9032854 110.18597  0.07082518     0.8475608           NA    TGA1
# 4  AT1G52150  0.0000000 5.531711e-01    0.7774755  68.49835  0.14513811     0.7320252           NA ATHB-15
# 36 AT5G60200  0.1381018 2.504228e-16    0.8174130  65.67806  0.11972797     0.7240564           NA  DOF5.3
# 19 AT3G01220  0.1404743 6.125165e-25    0.8662509  39.74159  0.05572179     0.8393039           NA  ATHB20
# 15 AT2G42200 -0.2081760 4.551701e-21   -0.8125636  27.10923 -0.09467280    -0.7916472           NA    SPL9
# 27 AT4G40060  0.2442242 1.110684e-13    0.7884207  18.99265  0.20938867     0.7775218           NA ATHB-16

          #----------------------------------------------------------------------------------------
          # 7th model: 58 TFs with >= 95% motif match within 2000 upstream, 500 downstream of TSS
          # use both matrices mtx and mtx3
          #----------------------------------------------------------------------------------------

targetGene <- canonicalizeName(tp, "SPL4")
tbl.atac <- data.frame(chrom=rep("Chr1", 2),
                       start=c(19804477,19806478),
                       end=c(19806477,19808478),
                       stringsAsFactors=FALSE)
tbl.tfs <- findCandidateTranscriptionFactorsByMotifInSequence(tp, tbl.atac, 95L)
candidate.tfs <- sort(unique(tbl.tfs$orf))
length(candidate.tfs)

recipe <- list(title="SPL4",
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
x7 <- build(builder)
# Warning messages:
#   1: In Solver(mtx.assay = mtx.assay, targetGene = targetGene, candidateRegulators = candidateRegulators,  :
#                  Target gene mean expression is in the bottom 10% of all genes in the assay matrix
lapply(x7, dim)
tbl.model.7 <- x7$model[1:6,]
tf.symbols <- unlist(lapply(tbl.model.7$gene, function(tf) getGeneNames(tp, tf)$symbol))
tbl.model.7$symbol <- tf.symbols

# tbl.model.7
#         gene   betaLasso  lassoPValue pearsonCoeff  rfScore   betaRidge spearmanCoeff bindingSites symbol
# 26 AT5G65210 -0.15594394 5.482777e-24   -0.8143540 62.05995 -0.05759788    -0.7727453           NA   TGA1
# 25 AT5G62165 -0.01743150 1.389467e-04   -0.7194983 37.21068 -0.04361129    -0.7004260           NA  AGL42
# 14 AT4G36540  0.08632200 2.273436e-10    0.7638330 28.78918  0.05588922     0.7307291           NA   BEE2
# 13 AT4G14465 -0.23710844 2.962709e-25   -0.8178321 24.92285 -0.05937281    -0.8397840           NA  AHL20
# 24 AT5G60200 -0.00865482 2.823662e-02   -0.7424131 21.20526 -0.07950806    -0.7293969           NA DOF5.3
# 9  AT3G01220 -0.05916779 1.825529e-09   -0.7785281 16.70120 -0.04286229    -0.7792379           NA ATHB20


