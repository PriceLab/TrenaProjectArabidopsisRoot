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
          # sixth model: 29 TFs with >= 95% motif match within 1000 (each) promoter.
          #------------------------------------------------------------------------------
tbl.atac <- data.frame(chrom=rep("Chr3", 2),
                       start=c(13032400,13037100),
                       end=c(13033400,13038100),
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
# tbl.model.6
#         gene   betaLasso   lassoPValue pearsonCoeff   rfScore   betaRidge spearmanCoeff bindingSites    symbol
# 8  AT4G14465 -0.38623478 1.919310e-162   -0.5628267 1072.8117 -0.29493570    -0.3313744           NA     AHL20
# 15 AT5G62165 -0.17274948 8.382164e-104   -0.4991596  660.4390 -0.16416343    -0.3849297           NA     AGL42
# 9  AT4G31920 -0.30060591  3.462460e-41   -0.3963221  465.5176 -0.30569923    -0.3279602           NA     ARR10
# 14 AT5G60200 -0.12411684  3.034644e-23   -0.3877858  402.2469 -0.16687658    -0.3551172           NA    DOF5.3
# 1  AT1G17950 -0.05301677  3.260175e-34   -0.4291068  352.9953 -0.08665296    -0.3200268           NA     MYB52
# 7  AT2G42200  0.17331619  1.024211e-34    0.3494224  256.2832  0.17743521     0.3171417           NA      SPL9
# 5  AT2G22430  0.34975419  6.445556e-38    0.3242114  246.3073  0.31453242     0.2671454           NA    ATHB-6
# 11 AT5G02320 -0.17835622  1.354184e-32   -0.3254427  218.6656 -0.17944006    -0.3079975           NA   MYB3R-5
# 16 AT5G65590  0.10574985  1.332595e-33    0.3535122  216.4876  0.12160481     0.2939522           NA AT5G65590
# 12 AT5G11510 -0.06323905  1.388600e-09   -0.2168312  177.7344 -0.09146526    -0.2150125           NA   MYB3R-4

#------------------------------------------------------------------------------------------------------------------------------
# Get high expression data
wbc19.vector <- mtx[c("AT3G55130", "AT3G59060"),]
plotWBC19 <- plot((mtx[c("AT3G55130"),]), (mtx[c("AT3G59060"),]))
new.wbc.vector <- t(wbc19.vector)
number1 <- new.wbc.vector[(new.wbc.vector[,1] > 2),]
number2 <- number1[(number1[,2] > 2),]
plot(number2, labels= rownames(number2))
highExpressionWBC19 <- t(number2)
rownames(highExpressionWBC19)
mtx2 <- highExpressionWBC19
mtx2.colnames <- colnames(mtx2)
mtx2 <-  mtx[,mtx2.colnames]

# created matrix of high expression data
# mtx3 has 169 samples
ordered.number2<- number2[order(number2[,1], number2[,2], decreasing=TRUE),]
# allows you to view sample names with highest expression
ordered.number2.50 <- ordered.number2[1:50,]
evenHigherExpressionWBC19 <- t(ordered.number2.50)
rownames(evenHigherExpressionWBC19)
mtx3.colnames <- colnames(evenHigherExpressionWBC19)
mtx3 <- mtx[,mtx3.colnames]
#build smaller matrix with even higher expression of WBC19

# top 1-3,7 samples are from E-GEOD-24348
# Experiment: Transcription profiling by array of Arabidopsis mutant for nas4x after iron deprivation

# top 4-6,9-10 samples are from experiment 489 from NASC arrays: 
# Experiment: Investigating the molecular response of Arabidopsis to combined 
# nematode and dehydration stress using the Affymetrix ATH1 microarray.

# top 8 sample is from Experiment #490 from NASC arrays:
# Experiment: Cell-specific nitrogen responses in the Arabidopsis root
#------------------------------------------------------------------------------------------------------------------------------
#Modeling with high expression data:

          #------------------------------------------------------------------------------
          # seventh model: 32 TFs with >= 95% motif match within 1500 (each) promoter.
          #------------------------------------------------------------------------------
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

# tbl.model.7
#         gene     betaLasso  lassoPValue pearsonCoeff  rfScore   betaRidge spearmanCoeff bindingSites symbol
# 10 AT3G61150  1.844280e-01 6.937774e-06    0.3380855 6.106435  0.13379306     0.3369922           NA   HDG1
# 11 AT5G16560  1.117630e-02 2.977540e-02    0.2091608 5.729552  0.03780677     0.2138357           NA    KAN
# 13 AT5G60200 -7.479631e-02 1.935798e-03   -0.2496830 5.132500 -0.08682439    -0.3265924           NA DOF5.3
# 7  AT2G37590 -1.357610e-01 5.089497e-04   -0.2755380 4.436605 -0.14352526    -0.3060042           NA DOF2.4
# 6  AT2G22430  4.522021e-02 2.389368e-02    0.2079384 4.399553  0.10152916     0.2029337           NA ATHB-6
# 4  AT1G73360  8.112228e-03 4.952906e-02    0.1015424 4.388088  0.05521869     0.2124932           NA   EDT1
# 12 AT5G25830  3.193348e-02 1.265326e-02    0.1590170 3.715073  0.04473851     0.1023321           NA GATA12
# 9  AT3G54810 -4.519279e-02 1.264134e-02   -0.1717809 3.307638 -0.07737611    -0.1361245           NA  GATA8
# 14 AT5G62940  0.000000e+00 8.643032e-01   -0.1254659 2.961352 -0.02415076    -0.1940754           NA DOF5.6
# 8  AT2G46130 -7.883023e-05 6.007031e-02   -0.1284293 2.786565 -0.02130100    -0.1615932           NA WRKY43
#------------------------------------------------------------------------------------------------------------------------

#use smaller mtx (mtx3) for even higher expression data

          #------------------------------------------------------------------------------
          # eighth model: 79 TFs with >= 95% motif match within 1500 (each) promoter.
          #------------------------------------------------------------------------------
tbl.atac <- data.frame(chrom=rep("Chr3", 2),
                       start=c(18632400,18638210),
                       end=c(18633900,18639710),
                       stringsAsFactors=FALSE)
tbl.tfs <- findCandidateTranscriptionFactorsByMotifInSequence(tp, tbl.atac, 95L)
candidate.tfs <- sort(unique(tbl.tfs$orf))
length(candidate.tfs)

recipe <- list(title="WBC19",
               type="noDNA.tfsSupplied",
               matrix=mtx3,
               candidateTFs=candidate.tfs,
               tfPool=getAllTranscriptionFactors(tp, "MotifDb"),
               tfPrefilterCorrelation=0.1,
               annotationDbFile=dbfile(org.At.tair.db),
               orderModelByColumn="rfScore",
               solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"),
               quiet=TRUE)

builder <- NoDnaModelBuilder(genome, targetGene, recipe, quiet=TRUE)
x8 <- build(builder)
lapply(x8, dim)
tbl.model.8 <- x8$model[1:10,]
tf.symbols <- unlist(lapply(tbl.model.8$gene, function(tf) getGeneNames(tp, tf)$symbol))
tbl.model.8$symbol <- tf.symbols

# tbl.model.8
#         gene    betaLasso lassoPValue pearsonCoeff   rfScore    betaRidge spearmanCoeff bindingSites    symbol
# 2  AT1G25550 -0.015443951  0.02172242   -0.3272428 0.9130473 -0.027804138   -0.29190876           NA AT1G25550
# 8  AT1G74930  0.000000000  0.83419203   -0.2764086 0.4672606 -0.007247037   -0.27875150           NA    ERF018
# 33 AT5G65210  0.000000000  0.13996296   -0.2391059 0.4449350 -0.025328820   -0.16869148           NA      TGA1
# 31 AT5G62165  0.000000000  0.24719669    0.2466379 0.3578841  0.017652083    0.25051621           NA     AGL42
# 4  AT1G61110  0.000000000  0.04399367    0.2745188 0.3098864  0.027149326    0.36998800           NA   anac025
# 3  AT1G53160  0.000000000  0.51715079    0.1885299 0.2887157  0.022863901    0.18703481           NA      SPL4
# 20 AT3G54810 -0.005432485  0.03930759   -0.2940763 0.2704288 -0.043393328   -0.20384154           NA     GATA8
# 29 AT5G39610  0.000000000  0.09783550   -0.2399027 0.2386434 -0.012143183   -0.27452581           NA    ATNAC2
# 18 AT3G06740  0.029844783  0.01895477    0.3307985 0.2139190  0.045933432    0.27352397           NA    GATA15
# 32 AT5G62940  0.000000000  0.20621717    0.1423515 0.2113919  0.018586458    0.08821128           NA    DOF5.6

          #------------------------------------------------------------------------------
          # ninth model: 24 TFs with >= 95% motif match within 1000 (each) promoter.
          #------------------------------------------------------------------------------

tbl.atac <- data.frame(chrom=rep("Chr3", 2),
                       start=c(13532400,13538210),
                       end=c(13533400,13539210),
                       stringsAsFactors=FALSE)
tbl.tfs <- findCandidateTranscriptionFactorsByMotifInSequence(tp, tbl.atac, 95L)
candidate.tfs <- sort(unique(tbl.tfs$orf))
length(candidate.tfs)

recipe <- list(title="WBC19",
               type="noDNA.tfsSupplied",
               matrix=mtx3,
               candidateTFs=candidate.tfs,
               tfPool=getAllTranscriptionFactors(tp, "MotifDb"),
               tfPrefilterCorrelation=0.1,
               annotationDbFile=dbfile(org.At.tair.db),
               orderModelByColumn="rfScore",
               solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"),
               quiet=TRUE)

builder <- NoDnaModelBuilder(genome, targetGene, recipe, quiet=TRUE)
x9 <- build(builder)
lapply(x9, dim)
tbl.model.9 <- x9$model[1:10,]
tf.symbols <- unlist(lapply(tbl.model.9$gene, function(tf) getGeneNames(tp, tf)$symbol))
tbl.model.9$symbol <- tf.symbols

# tbl.model.9
#         gene    betaLasso lassoPValue pearsonCoeff   rfScore    betaRidge spearmanCoeff bindingSites    symbol
# 6  AT3G54810 -0.060234503  0.03603976   -0.2940763 0.6984862 -0.063738119   -0.20384154           NA     GATA8
# 12 AT5G62940  0.000000000  0.24576433    0.1423515 0.6378016  0.028938147    0.08821128           NA    DOF5.6
# 7  AT4G32890  0.003008057  0.10549938    0.2345584 0.6144417  0.027627205    0.19164466           NA     GATA9
# 5  AT3G06740  0.101587678  0.01895477    0.3307985 0.6094221  0.073410073    0.27352397           NA    GATA15
# 2  AT1G53160  0.000000000  0.34003156    0.1885299 0.5776807  0.026310324    0.18703481           NA      SPL4
# 10 AT5G25830  0.000000000  0.25792875    0.1818998 0.5764563  0.015273616    0.07130852           NA    GATA12
# 9  AT5G16560  0.000000000  0.56268364    0.1300772 0.5381874  0.008683821    0.15668667           NA       KAN
# 8  AT5G02840  0.000000000  0.12336829   -0.1900216 0.3819700 -0.043429621   -0.11750300           NA      LCL1
# 13 AT5G65590  0.000000000  0.45605691   -0.1758766 0.3539472 -0.009495115   -0.30208884           NA AT5G65590
# 3  AT1G64620  0.000000000  0.15354621    0.2198560 0.3528297  0.049762578    0.21949580           NA AT1G64620

          #------------------------------------------------------------------------------
          # 10th model: 27 TFs with >= 95% motif match within 1000 (each) promoter.
          #------------------------------------------------------------------------------

tbl.atac <- data.frame(chrom=rep("Chr3", 2),
                       start=c(20434000,20435500),
                       end=c(20435000,20436500),
                       stringsAsFactors=FALSE)
tbl.tfs <- findCandidateTranscriptionFactorsByMotifInSequence(tp, tbl.atac, 95L)
candidate.tfs <- sort(unique(tbl.tfs$orf))
length(candidate.tfs)

recipe <- list(title="WBC19",
               type="noDNA.tfsSupplied",
               matrix=mtx3,
               candidateTFs=candidate.tfs,
               tfPool=getAllTranscriptionFactors(tp, "MotifDb"),
               tfPrefilterCorrelation=0.1,
               annotationDbFile=dbfile(org.At.tair.db),
               orderModelByColumn="rfScore",
               solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"),
               quiet=TRUE)

builder <- NoDnaModelBuilder(genome, targetGene, recipe, quiet=TRUE)
x10 <- build(builder)
lapply(x10, dim)
tbl.model.10 <- x10$model[1:10,]
tf.symbols <- unlist(lapply(tbl.model.10$gene, function(tf) getGeneNames(tp, tf)$symbol))
tbl.model.10$symbol <- tf.symbols

# tbl.model.10
#         gene    betaLasso lassoPValue pearsonCoeff   rfScore   betaRidge spearmanCoeff bindingSites    symbol
# 11 AT5G46760  0.000000000  0.72496950    0.1201494 0.5582605  0.02368710    0.18722689           NA AT5G46760
# 5  AT3G06740  0.099267109  0.01895477    0.3307985 0.5543676  0.06976204    0.27352397           NA    GATA15
# 12 AT5G62165  0.003440748  0.08922270    0.2466379 0.5489340  0.02126300    0.25051621           NA     AGL42
# 4  AT2G36270  0.000000000  0.21446488   -0.2409433 0.5174508 -0.01353082   -0.11510204           NA      ABI5
# 13 AT5G62940  0.000000000  0.23139773    0.1423515 0.5087755  0.02582457    0.08821128           NA    DOF5.6
# 7  AT3G54810 -0.054040991  0.03603976   -0.2940763 0.4856329 -0.04923703   -0.20384154           NA     GATA8
# 1  AT1G53160  0.000000000  0.26685569    0.1885299 0.4703625  0.02634593    0.18703481           NA      SPL4
# 6  AT3G47500  0.000000000  0.27552762    0.2092642 0.4636700  0.02766993    0.29901561           NA      CDF3
# 10 AT5G16560  0.000000000  0.60334479    0.1300772 0.4396909  0.01021718    0.15668667           NA       KAN
# 2  AT1G64620  0.000000000  0.13975644    0.2198560 0.4252454  0.06111624    0.21949580           NA AT1G64620


          #------------------------------------------------------------------------------
          # 11th model: 19 TFs with >= 95% motif match within 1000 (each) promoter.
          #------------------------------------------------------------------------------

tbl.atac <- data.frame(chrom=rep("Chr3", 2),
                       start=c(20433010,20439400),
                       end=c(20434010,20440400),
                       stringsAsFactors=FALSE)
tbl.tfs <- findCandidateTranscriptionFactorsByMotifInSequence(tp, tbl.atac, 95L)
candidate.tfs <- sort(unique(tbl.tfs$orf))
length(candidate.tfs)

recipe <- list(title="WBC19",
               type="noDNA.tfsSupplied",
               matrix=mtx3,
               candidateTFs=candidate.tfs,
               tfPool=getAllTranscriptionFactors(tp, "MotifDb"),
               tfPrefilterCorrelation=0.1,
               annotationDbFile=dbfile(org.At.tair.db),
               orderModelByColumn="rfScore",
               solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"),
               quiet=TRUE)

builder <- NoDnaModelBuilder(genome, targetGene, recipe, quiet=TRUE)
x11 <- build(builder)
lapply(x11, dim)
tbl.model.11 <- x11$model[1:10,]
tf.symbols <- unlist(lapply(tbl.model.11$gene, function(tf) getGeneNames(tp, tf)$symbol))
tbl.model.11$symbol <- tf.symbols

# tbl.model.11
#         gene    betaLasso lassoPValue pearsonCoeff   rfScore    betaRidge spearmanCoeff bindingSites    symbol
# 1  AT1G25550 -0.022746363  0.02172242   -0.3272428 1.2660687 -0.039665255   -0.29190876           NA AT1G25550
# 5  AT3G54810 -0.025411789  0.04902443   -0.2940763 0.6545343 -0.051632292   -0.20384154           NA     GATA8
# 10 AT5G62165  0.001002055  0.10201106    0.2466379 0.5926299  0.021996337    0.25051621           NA     AGL42
# 11 AT5G62940  0.000000000  0.32551922    0.1423515 0.5173541  0.026794172    0.08821128           NA    DOF5.6
# 6  AT4G14465 -0.029122906  0.02643625   -0.3113506 0.5150485 -0.027251150   -0.28038415           NA     AHL20
# 4  AT3G06740  0.091388846  0.01895477    0.3307985 0.5034917  0.072256535    0.27352397           NA    GATA15
# 2  AT1G64620  0.000000000  0.11792728    0.2198560 0.4519734  0.060525966    0.21949580           NA AT1G64620
# 9  AT5G16560  0.000000000  0.73158073    0.1300772 0.3951998  0.008438453    0.15668667           NA       KAN
# 12 AT5G65590  0.000000000  0.74044625   -0.1758766 0.3658925 -0.009455525   -0.30208884           NA AT5G65590
# 8  AT5G03790  0.000000000  0.19902620   -0.1720009 0.3526549 -0.013799422   -0.20096038           NA   ATHB-51

          #------------------------------------------------------------------------------
          # 12th model: 63 TFs with >= 95% motif match within 1000 (each) promoter.
          #------------------------------------------------------------------------------

tbl.atac <- data.frame(chrom=rep("Chr3", 2),
                       start=c(20431010,20436400),
                       end=c(20432010,20437400),
                       stringsAsFactors=FALSE)
tbl.tfs <- findCandidateTranscriptionFactorsByMotifInSequence(tp, tbl.atac, 95L)
candidate.tfs <- sort(unique(tbl.tfs$orf))
length(candidate.tfs)

recipe <- list(title="WBC19",
               type="noDNA.tfsSupplied",
               matrix=mtx3,
               candidateTFs=candidate.tfs,
               tfPool=getAllTranscriptionFactors(tp, "MotifDb"),
               tfPrefilterCorrelation=0.1,
               annotationDbFile=dbfile(org.At.tair.db),
               orderModelByColumn="rfScore",
               solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"),
               quiet=TRUE)

builder <- NoDnaModelBuilder(genome, targetGene, recipe, quiet=TRUE)
x12 <- build(builder)
lapply(x12, dim)
tbl.model.12 <- x12$model[1:10,]
tf.symbols <- unlist(lapply(tbl.model.12$gene, function(tf) getGeneNames(tp, tf)$symbol))
tbl.model.12$symbol <- tf.symbols

# tbl.model.12
#         gene    betaLasso  lassoPValue pearsonCoeff   rfScore   betaRidge spearmanCoeff bindingSites    symbol
# 17 AT4G17880  0.093382268 0.0008222258    0.4580876 1.2723674  0.04412493     0.4253061           NA AT4G17880
# 3  AT1G22985 -0.033356269 0.0099813571   -0.3858113 0.7122849 -0.04852291    -0.3998559           NA AT1G22985
# 13 AT3G15210  0.000000000 0.7570313477   -0.2229774 0.3615262 -0.01276041    -0.1080912           NA      ERF4
# 23 AT5G07680  0.000000000 0.0966796673   -0.2597433 0.3409373 -0.01130599    -0.1699400           NA   ANAC080
# 19 AT4G32800  0.000000000 0.0847828944    0.2511761 0.2765392  0.01980436     0.2693397           NA AT4G32800
# 15 AT4G14465 -0.018735267 0.0139181058   -0.3113506 0.2307609 -0.02534347    -0.2803842           NA     AHL20
# 14 AT4G00050  0.000000000 0.5831800806    0.1559773 0.2282816  0.01596675     0.2100840           NA     UNE10
# 26 AT5G43410  0.003556295 0.0252531286    0.2610044 0.2255470  0.01595780     0.3809364           NA AT5G43410
# 27 AT5G46760  0.000000000 0.8295434551    0.1201494 0.2081812  0.02031843     0.1872269           NA AT5G46760
# 29 AT5G62165  0.000000000 0.0632480856    0.2466379 0.1992148  0.01515128     0.2505162           NA     AGL42

          #------------------------------------------------------------------------------
          # 13th model: 52 TFs with >= 95% motif match within 1500 (each) promoter.
          #------------------------------------------------------------------------------

tbl.atac <- data.frame(chrom=rep("Chr3", 2),
                       start=c(20430000,20435500),
                       end=c(20431500,20437000),
                       stringsAsFactors=FALSE)
tbl.tfs <- findCandidateTranscriptionFactorsByMotifInSequence(tp, tbl.atac, 95L)
candidate.tfs <- sort(unique(tbl.tfs$orf))
length(candidate.tfs)

recipe <- list(title="WBC19",
               type="noDNA.tfsSupplied",
               matrix=mtx3,
               candidateTFs=candidate.tfs,
               tfPool=getAllTranscriptionFactors(tp, "MotifDb"),
               tfPrefilterCorrelation=0.1,
               annotationDbFile=dbfile(org.At.tair.db),
               orderModelByColumn="rfScore",
               solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"),
               quiet=TRUE)

builder <- NoDnaModelBuilder(genome, targetGene, recipe, quiet=TRUE)
x13 <- build(builder)
lapply(x13, dim)
tbl.model.13 <- x13$model[1:10,]
tf.symbols <- unlist(lapply(tbl.model.13$gene, function(tf) getGeneNames(tp, tf)$symbol))
tbl.model.13$symbol <- tf.symbols

# tbl.model.13
#         gene    betaLasso lassoPValue pearsonCoeff   rfScore   betaRidge spearmanCoeff bindingSites    symbol
# 2  AT1G22985 -0.072674574 0.005655387   -0.3858113 0.8585397 -0.05207582   -0.39985594           NA AT1G22985
# 11 AT3G15210  0.000000000 0.195240490   -0.2229774 0.4950218 -0.01760155   -0.10809124           NA      ERF4
# 20 AT5G07680  0.000000000 0.126299481   -0.2597433 0.4152436 -0.01300549   -0.16993998           NA   ANAC080
# 15 AT4G32800  0.000000000 0.060636290    0.2511761 0.3392273  0.02260782    0.26933974           NA AT4G32800
# 23 AT5G46760  0.000000000 0.402790027    0.1201494 0.2764723  0.02720795    0.18722689           NA AT5G46760
# 10 AT3G06740  0.004794159 0.044481922    0.3307985 0.2722499  0.04303362    0.27352397           NA    GATA15
# 7  AT1G69780  0.000000000 0.211205052   -0.1558149 0.2606508 -0.01934952   -0.24043217           NA    ATHB13
# 25 AT5G62940  0.000000000 0.322875673    0.1423515 0.2501957  0.01630500    0.08821128           NA    DOF5.6
# 12 AT3G47500  0.000000000 0.807130465    0.2092642 0.2484925  0.01951010    0.29901561           NA      CDF3
# 22 AT5G43410  0.000000000 0.052223010    0.2610044 0.2283790  0.01637408    0.38093637           NA AT5G43410

#no regulatory regions yet

          #------------------------------------------------------------------------------
          # 14th model: 42 TFs with >= 95% motif match within 1500 (each) promoter.
          #------------------------------------------------------------------------------

tbl.atac <- data.frame(chrom=rep("Chr3", 2),
                       start=c(20423050,20427384),
                       end=c(20424550,20428884),
                       stringsAsFactors=FALSE)
tbl.tfs <- findCandidateTranscriptionFactorsByMotifInSequence(tp, tbl.atac, 95L)
candidate.tfs <- sort(unique(tbl.tfs$orf))
length(candidate.tfs)

recipe <- list(title="WBC19",
               type="noDNA.tfsSupplied",
               matrix=mtx3,
               candidateTFs=candidate.tfs,
               tfPool=getAllTranscriptionFactors(tp, "MotifDb"),
               tfPrefilterCorrelation=0.1,
               annotationDbFile=dbfile(org.At.tair.db),
               orderModelByColumn="rfScore",
               solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"),
               quiet=TRUE)

builder <- NoDnaModelBuilder(genome, targetGene, recipe, quiet=TRUE)
x14 <- build(builder)
lapply(x14, dim)
tbl.model.14 <- x14$model[1:10,]
tf.symbols <- unlist(lapply(tbl.model.14$gene, function(tf) getGeneNames(tp, tf)$symbol))
tbl.model.14$symbol <- tf.symbols

# tbl.model.14
#         gene    betaLasso lassoPValue pearsonCoeff   rfScore    betaRidge spearmanCoeff bindingSites    symbol
# 10 AT4G14465 -0.032987913  0.02604191   -0.3113506 0.5811695 -0.032779964   -0.28038415           NA     AHL20
# 8  AT3G47500  0.000000000  0.10256329    0.2092642 0.5284426  0.038234129    0.29901561           NA      CDF3
# 3  AT1G69780  0.000000000  0.47126347   -0.1558149 0.4558065 -0.020464919   -0.24043217           NA    ATHB13
# 17 AT5G62940  0.000000000  0.44899545    0.1423515 0.4496875  0.019987754    0.08821128           NA    DOF5.6
# 7  AT3G06740  0.070921388  0.01895477    0.3307985 0.4338415  0.057434221    0.27352397           NA    GATA15
# 11 AT4G17460  0.009649258  0.04073222    0.3055506 0.4336755  0.024513398    0.24264106           NA      HAT1
# 2  AT1G64620  0.000000000  0.11742320    0.2198560 0.4100220  0.060567525    0.21949580           NA AT1G64620
# 12 AT4G37180 -0.025091221  0.05295391   -0.2804712 0.3833676 -0.044010154   -0.29142857           NA AT4G37180
# 15 AT5G16560  0.000000000  0.95486682    0.1300772 0.3681782  0.006399682    0.15668667           NA       KAN
# 5  AT2G46680  0.000000000  0.10509509   -0.2534576 0.3379815 -0.017638285   -0.26261705           NA    ATHB-7


          #---------------------------------------------------------------------------------------
          # 15th model: 20 TFs with >= 95% motif match within 2000 upstream, 500 downstream of TSS
          #---------------------------------------------------------------------------------------

tbl.atac <- data.frame(chrom=rep("Chr3", 2),
                       start=c(20432106,20434107),
                       end=c(20434106,20434607),
                       stringsAsFactors=FALSE)
tbl.tfs <- findCandidateTranscriptionFactorsByMotifInSequence(tp, tbl.atac, 95L)
candidate.tfs <- sort(unique(tbl.tfs$orf))
length(candidate.tfs)

recipe <- list(title="WBC19",
               type="noDNA.tfsSupplied",
               matrix=mtx3,
               candidateTFs=candidate.tfs,
               tfPool=getAllTranscriptionFactors(tp, "MotifDb"),
               tfPrefilterCorrelation=0.1,
               annotationDbFile=dbfile(org.At.tair.db),
               orderModelByColumn="rfScore",
               solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"),
               quiet=TRUE)

builder <- NoDnaModelBuilder(genome, targetGene, recipe, quiet=TRUE)
x15 <- build(builder)
lapply(x15, dim)
tbl.model.15 <- x15$model[1:10,]
tf.symbols <- unlist(lapply(tbl.model.15$gene, function(tf) getGeneNames(tp, tf)$symbol))
tbl.model.15$symbol <- tf.symbols

#         gene    betaLasso lassoPValue pearsonCoeff   rfScore   betaRidge spearmanCoeff bindingSites    symbol
# 6  AT3G54810 -0.039787610  0.04365073   -0.2940763 0.6846830 -0.06001076   -0.20384154           NA     GATA8
# 5  AT3G06740  0.124765638  0.01895477    0.3307985 0.6674693  0.08295289    0.27352397           NA    GATA15
# 10 AT5G62940  0.000000000  0.41461791    0.1423515 0.6670734  0.02197803    0.08821128           NA    DOF5.6
# 7  AT4G14465 -0.035092534  0.02604191   -0.3113506 0.6431163 -0.03020385   -0.28038415           NA     AHL20
# 2  AT1G69780  0.000000000  0.91517751   -0.1558149 0.6375511 -0.01767389   -0.24043217           NA    ATHB13
# 1  AT1G64620  0.003497219  0.12537262    0.2198560 0.6301531  0.06288908    0.21949580           NA AT1G64620
# 9  AT5G16560  0.000000000  0.54146029    0.1300772 0.5653065  0.01288556    0.15668667           NA       KAN
# 11 AT5G65590  0.000000000  0.53350580   -0.1758766 0.4449672 -0.01152127   -0.30208884           NA AT5G65590
# 3  AT1G80840  0.000000000  0.44366311   -0.1730834 0.3949913 -0.01562107   -0.14900360           NA    WRKY40
# 4  AT3G01470  0.000000000  0.85935301   -0.1515070 0.3905654 -0.01787426   -0.13373349           NA    ATHB-1

            #----------------------------------------------------------------------------------------
            # 16th model: 20 TFs with >= 95% motif match within 2000 upstream, 500 downstream of TSS
            # same region as #15 on atac.seq; Use bigger matrix
            #----------------------------------------------------------------------------------------

tbl.atac <- data.frame(chrom=rep("Chr3", 2),
                       start=c(20432112,20434113),
                       end=c(20434112,20434613),
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
x16 <- build(builder)
lapply(x16, dim)
tbl.model.16 <- x16$model[1:10,]
tf.symbols <- unlist(lapply(tbl.model.16$gene, function(tf) getGeneNames(tp, tf)$symbol))
tbl.model.16$symbol <- tf.symbols

#         gene    betaLasso   lassoPValue pearsonCoeff   rfScore   betaRidge spearmanCoeff bindingSites    symbol
# 7  AT4G14465 -0.631679182 1.919310e-162   -0.5628267 1392.8188 -0.45797827    -0.3313744           NA     AHL20
# 9  AT5G60200 -0.226439315  9.592281e-32   -0.3877858  663.6387 -0.24379551    -0.3551172           NA    DOF5.3
# 1  AT1G52150 -0.007069971  3.660479e-02   -0.2822424  510.0015 -0.10189183    -0.3354909           NA   ATHB-15
# 4  AT2G22430  0.534924218  3.009286e-43    0.3242114  484.3623  0.49642149     0.2671454           NA    ATHB-6
# 10 AT5G65590  0.196341104  5.449207e-39    0.3535122  453.4247  0.20074867     0.2939522           NA AT5G65590
# 6  AT3G01470  0.000000000  9.838639e-01   -0.2254802  318.4881 -0.10405340    -0.1848806           NA    ATHB-1
# 2  AT1G64620  0.032505573  1.688527e-02   -0.1025043  291.1077  0.11394430    -0.1231658           NA AT1G64620
# 5  AT2G37590  0.000000000  5.330965e-01   -0.2068212  260.5657 -0.08101563    -0.1980168           NA    DOF2.4
# 8  AT5G39660 -0.026584270  5.430603e-03   -0.2213106  253.8812 -0.12570706    -0.1426431           NA      CDF2
# 3  AT1G80840  0.000000000  3.945345e-01    0.1509693  251.0913  0.03387059     0.1306107           NA    WRKY40

#with mtx3:
#         gene   betaLasso lassoPValue pearsonCoeff   rfScore   betaRidge spearmanCoeff bindingSites    symbol
# 6  AT3G54810 -0.03167221  0.04365073   -0.2940763 0.8106728 -0.05785286   -0.20384154           NA     GATA8
# 2  AT1G69780  0.00000000  0.91517751   -0.1558149 0.6819254 -0.01714176   -0.24043217           NA    ATHB13
# 10 AT5G62940  0.00000000  0.41461791    0.1423515 0.6513227  0.02122323    0.08821128           NA    DOF5.6
# 7  AT4G14465 -0.02996968  0.02604191   -0.3113506 0.6449731 -0.02901716   -0.28038415           NA     AHL20
# 5  AT3G06740  0.10799807  0.01895477    0.3307985 0.6422672  0.07991416    0.27352397           NA    GATA15
# 9  AT5G16560  0.00000000  0.54146029    0.1300772 0.5825086  0.01236389    0.15668667           NA       KAN
# 1  AT1G64620  0.00000000  0.12537262    0.2198560 0.5129569  0.06036687    0.21949580           NA AT1G64620
# 11 AT5G65590  0.00000000  0.53350580   -0.1758766 0.4958118 -0.01116943   -0.30208884           NA AT5G65590
# 8  AT5G03790  0.00000000  0.18322862   -0.1720009 0.4236046 -0.01437446   -0.20096038           NA   ATHB-51
# 3  AT1G80840  0.00000000  0.44366311   -0.1730834 0.4227791 -0.01505340   -0.14900360           NA    WRKY40

          #----------------------------------------------------------------------------------------
          # 17th model: 28 TFs with >= 95% motif match within 2000 upstream, 2000 downstream of TSS
          # same region as #15 on atac.seq; Use bigger matrix
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
x17 <- build(builder)
lapply(x17, dim)
tbl.model.17 <- x17$model[1:10,]
tf.symbols <- unlist(lapply(tbl.model.17$gene, function(tf) getGeneNames(tp, tf)$symbol))
tbl.model.17$symbol <- tf.symbols

#         gene   betaLasso   lassoPValue pearsonCoeff   rfScore   betaRidge spearmanCoeff bindingSites    symbol
# 9  AT4G14465 -0.41358299 1.919310e-162   -0.5628267 1103.0008 -0.30944162    -0.3313744           NA     AHL20
# 13 AT5G62165 -0.21833735 8.382164e-104   -0.4991596  663.4349 -0.19320224    -0.3849297           NA     AGL42
# 2  AT1G53160  0.16486489  1.379180e-46    0.3917817  474.6536  0.16639273     0.3554472           NA      SPL4
# 10 AT5G18830 -0.30371523  7.029777e-31   -0.3776184  430.6054 -0.34049818    -0.3622195           NA      SPL7
# 12 AT5G60200 -0.13385030  5.107913e-22   -0.3877858  397.7115 -0.16566286    -0.3551172           NA    DOF5.3
# 6  AT2G22430  0.42007791  1.191427e-37    0.3242114  282.8667  0.38108896     0.2671454           NA    ATHB-6
# 1  AT1G52150  0.00000000  7.489972e-01   -0.2822424  281.1619 -0.04939149    -0.3354909           NA   ATHB-15
# 14 AT5G65590  0.12167810  4.984743e-33    0.3535122  253.8221  0.13395108     0.2939522           NA AT5G65590
# 8  AT3G01470  0.00000000  3.134900e-01   -0.2254802  193.3893 -0.01969008    -0.1848806           NA    ATHB-1
# 3  AT1G64620  0.01893602  1.611152e-02   -0.1025043  179.9023  0.09867434    -0.1231658           NA AT1G64620

          #----------------------------------------------------------------------------------------
          # 18th model: 104 TFs with >= 95% motif match within 5000 upstream, 5000 downstream of TSS
          # same region as #15 on atac.seq; Use bigger matrix
          #----------------------------------------------------------------------------------------

tbl.atac <- data.frame(chrom=rep("Chr3", 2),
                       start=c(20429112,20434113),
                       end=c(20434112,20439113),
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
x18 <- build(builder)
lapply(x18, dim)
tbl.model.18 <- x18$model[1:10,]
tf.symbols <- unlist(lapply(tbl.model.18$gene, function(tf) getGeneNames(tp, tf)$symbol))
tbl.model.18$symbol <- tf.symbols

#         gene   betaLasso   lassoPValue pearsonCoeff    rfScore   betaRidge spearmanCoeff bindingSites symbol
# 30 AT3G59060  0.20660505 7.586023e-257    0.6738150 1106.31242  0.09978381     0.5077358           NA   PIF5
# 33 AT4G14465 -0.14559575 7.720037e-103   -0.5628267  537.79490 -0.09335849    -0.3313744           NA  AHL20
# 32 AT4G00050  0.21111517 1.790158e-129    0.5765029  451.14938  0.14417470     0.4916015           NA  UNE10
# 38 AT4G36540  0.16916609  2.128292e-89    0.5333334  329.98386  0.12713276     0.4806764           NA   BEE2
# 23 AT3G01220 -0.08784450  1.524911e-46   -0.5233056  273.24071 -0.08700164    -0.3849471           NA ATHB20
# 51 AT5G62165 -0.03805766  3.058215e-24   -0.4991596  217.69778 -0.05817641    -0.3849297           NA  AGL42
# 9  AT1G53160  0.02257309  1.136179e-12    0.3917817   92.34605  0.04437246     0.3554472           NA   SPL4
# 50 AT5G60200 -0.02179131  5.127391e-12   -0.3877858   91.62244 -0.05630662    -0.3551172           NA DOF5.3
# 43 AT5G18830 -0.12205317  5.265642e-17   -0.3776184   82.37603 -0.15105459    -0.3622195           NA   SPL7
# 35 AT4G31920  0.00000000  9.739117e-01   -0.3963221   82.05983 -0.08120219    -0.3279602           NA  ARR10

      #----------------------------------------------------------------------------------------
      # 19th model: 98 TFs with >= 95% motif match within 5000 upstream, 5000 downstream of TSS
      # same region as #15 on atac.seq; Use bigger matrix
      #----------------------------------------------------------------------------------------
tbl.atac <- data.frame(chrom=rep("Chr3", 2),
                       start=c(20429112,20436288),
                       end=c(20434112,20441288),
                       stringsAsFactors=FALSE)
tbl.tfs <- findCandidateTranscriptionFactorsByMotifInSequence(tp, tbl.atac, 95L)
candidate.tfs <- sort(unique(tbl.tfs$orf))
length(candidate.tfs)

recipe <- list(title="WBC19",
               type="noDNA.tfsSupplied",
               matrix=mtx3,
               candidateTFs=candidate.tfs,
               tfPool=getAllTranscriptionFactors(tp, "MotifDb"),
               tfPrefilterCorrelation=0.1,
               annotationDbFile=dbfile(org.At.tair.db),
               orderModelByColumn="rfScore",
               solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"),
               quiet=TRUE)

builder <- NoDnaModelBuilder(genome, targetGene, recipe, quiet=TRUE)
x19 <- build(builder)
lapply(x19, dim)
tbl.model.19 <- x19$model[1:10,]
tf.symbols <- unlist(lapply(tbl.model.19$gene, function(tf) getGeneNames(tp, tf)$symbol))
tbl.model.19$symbol <- tf.symbols

#         gene   betaLasso   lassoPValue pearsonCoeff    rfScore   betaRidge spearmanCoeff bindingSites symbol
# 27 AT3G59060  0.20847974 7.586023e-257    0.6738150 1132.98732  0.10504491     0.5077358           NA   PIF5
# 30 AT4G14465 -0.14831565 7.720037e-103   -0.5628267  527.16775 -0.09513729    -0.3313744           NA  AHL20
# 29 AT4G00050  0.21213270 1.790158e-129    0.5765029  467.57087  0.14878185     0.4916015           NA  UNE10
# 21 AT3G01220 -0.08679869  1.524911e-46   -0.5233056  329.26474 -0.09229795    -0.3849471           NA ATHB20
# 34 AT4G36540  0.16911396  2.128292e-89    0.5333334  320.70080  0.12903040     0.4806764           NA   BEE2
# 47 AT5G62165 -0.03733050  3.058215e-24   -0.4991596  185.77783 -0.06039029    -0.3849297           NA  AGL42
# 46 AT5G60200 -0.02111709  5.127391e-12   -0.3877858   91.65985 -0.05759636    -0.3551172           NA DOF5.3
# 39 AT5G18830 -0.12003641  5.265642e-17   -0.3776184   86.14901 -0.15167134    -0.3622195           NA   SPL7
# 9  AT1G53160  0.02216566  1.136179e-12    0.3917817   84.55421  0.04472416     0.3554472           NA   SPL4
# 32 AT4G31920  0.00000000  9.427982e-01   -0.3963221   73.34116 -0.07899182    -0.3279602           NA  ARR10

# with smaller matrix mtx3:
#         gene   betaLasso  lassoPValue pearsonCoeff   rfScore    betaRidge spearmanCoeff bindingSites    symbol
# 29 AT4G17880  0.06445869 0.0008222258    0.4580876 0.9450456  0.034282571     0.4253061           NA AT4G17880
# 7  AT1G25550  0.00000000 0.0328008880   -0.3272428 0.5391818 -0.019705100    -0.2919088           NA AT1G25550
# 40 AT5G37020  0.10286850 0.0010999014    0.4500930 0.5314113  0.073007047     0.3855462           NA      ARF8
# 5  AT1G22985 -0.01027799 0.0142236209   -0.3858113 0.3790130 -0.038876945    -0.3998559           NA AT1G22985
# 37 AT5G07680  0.00000000 0.0518950825   -0.2597433 0.3328196 -0.010513774    -0.1699400           NA   ANAC080
# 22 AT3G15270  0.00000000 0.0593857696    0.3030271 0.2860886  0.013592940     0.2473469           NA      SPL5
# 21 AT3G15210  0.00000000 0.7331515557   -0.2229774 0.1963197 -0.007838267    -0.1080912           NA      ERF4
# 14 AT1G74930  0.00000000 0.6960000957   -0.2764086 0.1911578 -0.004877530    -0.2787515           NA    ERF018
# 41 AT5G43410  0.00000000 0.0281112451    0.2610044 0.1895305  0.016059188     0.3809364           NA AT5G43410
# 32 AT4G32800  0.00000000 0.0749800426    0.2511761 0.1830661  0.018230405     0.2693397           NA AT4G32800



#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
  runTests()

