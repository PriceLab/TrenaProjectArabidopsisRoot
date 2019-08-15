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
# make sure to review the others as well

          #------------------------------------------------------------------------------
          # fourth model: 19 TFs with >= 95% motif match within 1200kb (each) promoter.
          #------------------------------------------------------------------------------


tbl.atac <- data.frame(chrom=rep("Chr3", 2),
                       start=c(20432000,20433210),
                       end=c(20433200,20434410),
                       stringsAsFactors=FALSE)
tbl.tfs <- findCandidateTranscriptionFactorsByMotifInSequence(tp, tbl.atac, 95L)
candidate.tfs <- sort(unique(tbl.tfs$orf))
length(candidate.tfs)
#side note: next model, go back to around the third model's region.

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
x4 <- build(builder)
lapply(x4, dim)
tbl.model.4 <- x4$model[1:10,]
tf.symbols <- unlist(lapply(tbl.model.4$gene, function(tf) getGeneNames(tp, tf)$symbol))
tbl.model.4$symbol <- tf.symbols

# tbl.model.4
#         gene   betaLasso   lassoPValue pearsonCoeff   rfScore   betaRidge spearmanCoeff bindingSites    symbol
# 7  AT4G14465 -0.63306716 1.919310e-162   -0.5628267 1372.8482 -0.44939720    -0.3313744           NA     AHL20
# 9  AT5G60200 -0.23040752  9.592281e-32   -0.3877858  673.3668 -0.24152133    -0.3551172           NA    DOF5.3
# 1  AT1G52150 -0.01352807  3.660479e-02   -0.2822424  493.6272 -0.10277135    -0.3354909           NA   ATHB-15
# 4  AT2G22430  0.53927909  3.009286e-43    0.3242114  481.4268  0.48990652     0.2671454           NA    ATHB-6
# 10 AT5G65590  0.19838070  5.449207e-39    0.3535122  454.8653  0.19894077     0.2939522           NA AT5G65590
# 6  AT3G01470  0.00000000  9.838639e-01   -0.2254802  327.7194 -0.10693530    -0.1848806           NA    ATHB-1
# 2  AT1G64620  0.04547039  1.688527e-02   -0.1025043  292.8133  0.10958571    -0.1231658           NA AT1G64620
# 5  AT2G37590  0.00000000  5.330965e-01   -0.2068212  261.5152 -0.08311396    -0.1980168           NA    DOF2.4
# 8  AT5G39660 -0.03144685  5.430603e-03   -0.2213106  254.0163 -0.12676544    -0.1426431           NA      CDF2
# 3  AT1G80840  0.00000000  3.945345e-01    0.1509693  252.7218  0.03454721     0.1306107           NA    WRKY40

          #------------------------------------------------------------------------------
          # fifth model: 30 TFs with >= 95% motif match within 1500 (each) promoter.
          #------------------------------------------------------------------------------

tbl.atac <- data.frame(chrom=rep("Chr3", 2),
                       start=c(22532400,22538210),
                       end=c(22533900,22539910),
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
x5 <- build(builder)
lapply(x5, dim)
tbl.model.5 <- x5$model[1:10,]
tf.symbols <- unlist(lapply(tbl.model.5$gene, function(tf) getGeneNames(tp, tf)$symbol))
tbl.model.5$symbol <- tf.symbols

# tbl.model.5
#         gene   betaLasso   lassoPValue pearsonCoeff   rfScore   betaRidge spearmanCoeff bindingSites    symbol
# 8  AT4G14465 -0.45126281 1.919310e-162   -0.5628267 1217.9599 -0.33094480    -0.3313744           NA     AHL20
# 4  AT2G01760  0.35079119  1.811735e-68    0.4259098  591.9681  0.31450068     0.3974638           NA     ARR14
# 9  AT4G31920 -0.41681921  2.657399e-51   -0.3963221  578.9558 -0.38395182    -0.3279602           NA     ARR10
# 12 AT5G60200 -0.21672612  2.056503e-33   -0.3877858  526.6647 -0.22047585    -0.3551172           NA    DOF5.3
# 2  AT1G67710 -0.04803401  3.653619e-14   -0.3991527  319.9506 -0.09684702    -0.3220879           NA     ARR11
# 13 AT5G65590  0.14140722  6.260262e-34    0.3535122  297.9758  0.15024399     0.2939522           NA AT5G65590
# 10 AT4G40060  0.00000000  2.658003e-01   -0.2832086  235.4189 -0.01086425    -0.2614751           NA   ATHB-16
# 1  AT1G64620  0.03707026  1.295114e-02   -0.1025043  207.6641  0.10161573    -0.1231658           NA AT1G64620
# 5  AT2G37590 -0.04281820  3.412511e-03   -0.2068212  200.8512 -0.14316299    -0.1980168           NA    DOF2.4
# 11 AT5G39660 -0.08074957  4.412580e-06   -0.2213106  191.9243 -0.14452532    -0.1426431           NA      CDF2

          #------------------------------------------------------------------------------
          # sixth model:  TFs with >= 95% motif match within 1000 (each) promoter.
          #------------------------------------------------------------------------------
tbl.atac <- data.frame(chrom=rep("Chr3", 2),
                       start=c(13532400,1353710),
                       end=c(13533400,13534710),
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
x6 <- build(builder)
lapply(x6, dim)
tbl.model.6 <- x6$model[1:10,]
tf.symbols <- unlist(lapply(tbl.model.6$gene, function(tf) getGeneNames(tp, tf)$symbol))
tbl.model.6$symbol <- tf.symbols



#Get high expression data
plotWBC19 <- plot((mtx[c("AT3G55130"),]), (mtx[c("AT3G59060"),]))
new.wbc.vector <- t(wbc19.vector)
number1 <- new.wbc.vector[(new.wbc.vector[,1] > 2),]
number2 <- number1[(number1[,2] > 2),]
plot(number2)
highExpressionWBC19 <- t(number2)
rownames(highExpressionWBC19)
mtx2 <- highExpressionWBC19
#created matrix of high expression data

tbl.atac <- data.frame(chrom=rep("Chr3", 2),
                       start=c(22632400,22638210),
                       end=c(22633900,22639910),
                       stringsAsFactors=FALSE)
tbl.tfs <- findCandidateTranscriptionFactorsByMotifInSequence(tp, tbl.atac, 95L)
candidate.tfs <- sort(unique(tbl.tfs$orf))
length(candidate.tfs)

recipe <- list(title="WBC19",
               type="noDNA.tfsSupplied",
               matrix=mtx2,
               candidateTFs=candidate.tfs,
               tfPool=getAllTranscriptionFactors(tp, "MotifDb"),
               tfPrefilterCorrelation=0.1,
               annotationDbFile=dbfile(org.At.tair.db),
               orderModelByColumn="rfScore",
               solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"),
               quiet=TRUE)

builder <- NoDnaModelBuilder(genome, targetGene, recipe, quiet=TRUE)
x7 <- build(builder)
lapply(x7, dim)
tbl.model.7 <- x7$model[1:10,]
tf.symbols <- unlist(lapply(tbl.model.7$gene, function(tf) getGeneNames(tp, tf)$symbol))
tbl.model.7$symbol <- tf.symbols

#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
  runTests()

