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
# 10 AT5G65590  0.19838070  5.449207e-39    0.3535122  454.8653  0.19894077     0.2939522           NA AT5G65590 aka SCAP1
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
# 13 AT5G65590  0.14140722  6.260262e-34    0.3535122  297.9758  0.15024399     0.2939522           NA AT5G65590 aka SCAP1
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
# 16 AT5G65590  0.10574985  1.332595e-33    0.3535122  216.4876  0.12160481     0.2939522           NA AT5G65590 aka SCAP1
# 12 AT5G11510 -0.06323905  1.388600e-09   -0.2168312  177.7344 -0.09146526    -0.2150125           NA   MYB3R-4

#------------------------------------------------------------------------------------------------------------------------------
# Get high expression data
wbc19.vector <- mtx[c("AT3G55130", "AT3G59060"),]
plotWBC19 <- plot((mtx[c("AT3G55130"),]), (mtx[c("AT3G59060"),]))
new.wbc.vector <- t(wbc19.vector)
number1 <- new.wbc.vector[(new.wbc.vector[,1] > 2),]
number2 <- number1[(number1[,2] > 2),]
plot(number2)
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
# 13 AT5G65590  0.000000000  0.45605691   -0.1758766 0.3539472 -0.009495115   -0.30208884           NA AT5G65590 aka SCAP1
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
# 12 AT5G65590  0.000000000  0.74044625   -0.1758766 0.3658925 -0.009455525   -0.30208884           NA AT5G65590 aka SCAP1
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
# 11 AT5G65590  0.000000000  0.53350580   -0.1758766 0.4449672 -0.01152127   -0.30208884           NA AT5G65590 aka SCAP1
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
# 10 AT5G65590  0.196341104  5.449207e-39    0.3535122  453.4247  0.20074867     0.2939522           NA AT5G65590 aka SCAP1
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
# 11 AT5G65590  0.00000000  0.53350580   -0.1758766 0.4958118 -0.01116943   -0.30208884           NA AT5G65590 aka SCAP1
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
# 14 AT5G65590  0.12167810  4.984743e-33    0.3535122  253.8221  0.13395108     0.2939522           NA AT5G65590 aka SCAP1
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
      # same region as previous on atac.seq; Use smaller matrix
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
#Model transcription factors found and see what you find! Remember to use both the bigger and smaller matrix

          #----------------------------------------------------------------------------------------
          # 20th model: 48 TFs with >= 95% motif match within 2000 upstream, 500 downstream of TSS
          # use both matrices mtx and mtx3
          #----------------------------------------------------------------------------------------

targetGene <- canonicalizeName(tp, "AHL20")
tbl.atac <- data.frame(chrom=rep("Chr4", 2),
                       start=c(8318970,8320971),
                       end=c(8320970,8321471),
                       stringsAsFactors=FALSE)
tbl.tfs <- findCandidateTranscriptionFactorsByMotifInSequence(tp, tbl.atac, 95L)
candidate.tfs <- sort(unique(tbl.tfs$orf))
length(candidate.tfs)

recipe <- list(title="AHL20",
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
x20 <- build(builder)
lapply(x20, dim)
tbl.model.20 <- x20$model[1:10,]
tf.symbols <- unlist(lapply(tbl.model.20$gene, function(tf) getGeneNames(tp, tf)$symbol))
tbl.model.20$symbol <- tf.symbols

# with bigger mtx

#         gene   betaLasso   lassoPValue pearsonCoeff   rfScore   betaRidge spearmanCoeff bindingSites  symbol
# 24 AT5G65210  0.28749167 9.679734e-319    0.7271170 489.54621  0.17623964    0.59701574           NA    TGA1
# 13 AT3G01220  0.18996414 2.235051e-158    0.6135063 330.79385  0.14670230    0.44144524           NA  ATHB20
# 23 AT5G60200  0.17686271  2.253602e-80    0.4832782 141.30917  0.13456460    0.44778587           NA  DOF5.3
# 2  AT1G17950  0.07365367 9.087122e-100    0.5521620 128.18415  0.06900664    0.48992471           NA   MYB52
# 14 AT3G01470  0.16093068  1.016676e-59    0.5140175 127.37487  0.18588474    0.52632393           NA  ATHB-1
# 16 AT4G40060  0.13269051  9.650099e-54    0.5204446 117.28693  0.16293876    0.52720995           NA ATHB-16
# 4  AT1G52150  0.00000000  7.716481e-01    0.3620570  73.96752  0.03486273    0.41964467           NA ATHB-15
# 8  AT2G22430 -0.04970529  1.952835e-06   -0.1975364  43.26965 -0.13080218   -0.06137936           NA  ATHB-6
# 9  AT2G25000  0.03254845  4.784342e-07    0.1999678  41.89945  0.07458730    0.26489949           NA  WRKY60
# 21 AT5G49520  0.05100772  2.200487e-13    0.3075649  37.30326  0.08671889    0.25387406           NA  WRKY48

#with smaller mtx3
#         gene   betaLasso  lassoPValue pearsonCoeff  rfScore   betaRidge spearmanCoeff bindingSites  symbol
# 6  AT1G69780  0.52917750 9.251529e-08    0.6717750 9.866425  0.13695135     0.4410564           NA  ATHB13
# 16 AT4G40060  0.00000000 6.613187e-02    0.3992646 3.872053  0.15321376     0.2845138           NA ATHB-16
# 12 AT3G01470  0.00000000 9.751265e-01    0.4376918 3.548930  0.08064023     0.3113085           NA  ATHB-1
# 8  AT2G22430  0.00000000 1.042840e-01   -0.4260024 2.915782 -0.14902930    -0.4021609           NA  ATHB-6
# 9  AT2G25000  0.01184355 1.100339e-02    0.3901681 2.907999  0.11822496     0.3388715           NA  WRKY60
# 25 AT5G65210  0.00000000 2.394996e-02    0.4502085 2.695101  0.10478313     0.4290516           NA    TGA1
# 13 AT3G54810  0.09584137 4.146273e-03    0.4068371 2.419947  0.17859846     0.3986074           NA   GATA8
# 4  AT1G52150  0.00000000 1.860880e-01   -0.3088507 2.376620 -0.12262037    -0.3213926           NA ATHB-15
# 22 AT5G56270 -0.08530994 4.935621e-03   -0.3738075 2.367568 -0.17138477    -0.4034094           NA   WRKY2
# 24 AT5G62940  0.00000000 8.088673e-01   -0.1530649 2.075236 -0.02158921    -0.1006002           NA  DOF5.6

          #----------------------------------------------------------------------------------------
          # 21st model: 66 TFs with >= 95% motif match within 2000 upstream, 500 downstream of TSS
          # use both matrices mtx and mtx3
          #----------------------------------------------------------------------------------------

targetGene <- canonicalizeName(tp, "PIF5")
tbl.atac <- data.frame(chrom=rep("Chr3", 2),
                       start=c(21826187,21828188),
                       end=c(21828187,21828688),
                       stringsAsFactors=FALSE)
tbl.tfs <- findCandidateTranscriptionFactorsByMotifInSequence(tp, tbl.atac, 95L)
candidate.tfs <- sort(unique(tbl.tfs$orf))
length(candidate.tfs)

recipe <- list(title="PIF5",
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
x21 <- build(builder)
lapply(x21, dim)
tbl.model.21 <- x21$model[1:10,]
tf.symbols <- unlist(lapply(tbl.model.21$gene, function(tf) getGeneNames(tp, tf)$symbol))
tbl.model.21$symbol <- tf.symbols

#with bigger matrix
#         gene   betaLasso   lassoPValue pearsonCoeff   rfScore   betaRidge spearmanCoeff bindingSites    symbol
# 28 AT4G36540  0.31729995 7.689933e-146    0.5413038 1301.6056  0.22207458     0.4579084           NA      BEE2
# 36 AT5G62165 -0.20517737 9.475748e-170   -0.5731729 1169.2908 -0.15495633    -0.4259736           NA     AGL42
# 24 AT4G31920 -0.33374588  7.894880e-84   -0.4761674  558.2853 -0.31075981    -0.3834930           NA     ARR10
# 22 AT3G60490 -0.12682864  8.515071e-76   -0.4383758  403.6812 -0.09674127    -0.3911913           NA AT3G60490
# 35 AT5G60200 -0.16877910  5.074558e-47   -0.4158650  359.4931 -0.16318164    -0.3608002           NA    DOF5.3
# 31 AT5G25830 -0.21821714  6.054286e-93   -0.4750765  342.0631 -0.20263815    -0.3635412           NA    GATA12
# 1  AT1G17950 -0.08224473  1.652507e-97   -0.4788226  332.4711 -0.07033979    -0.3746657           NA     MYB52
# 6  AT1G53170  0.17270085  1.592222e-70    0.4460720  267.3134  0.13643533     0.3760144           NA      ERF8
# 21 AT3G57600  0.16206402  1.123689e-74    0.4555560  236.9006  0.13824439     0.3921747           NA AT3G57600
# 20 AT3G16770  0.16564353  1.974600e-35    0.3417913  179.2040  0.17261444     0.3205095           NA    RAP2-3

#with small matrix:
#         gene    betaLasso lassoPValue pearsonCoeff  rfScore   betaRidge spearmanCoeff bindingSites    symbol
# 1  AT1G06160  0.000000000 0.035329865    0.2210415 2.217605  0.03002289    0.03442977           NA    ERF094
# 10 AT2G22430  0.073019124 0.006538555    0.4177917 1.888240  0.07946312    0.38141657           NA    ATHB-6
# 30 AT5G62000  0.032641170 0.002867844    0.4235871 1.837604  0.03986967    0.44460984           NA      ARF2
# 3  AT1G28370  0.080778283 0.001006354    0.4512119 1.677484  0.03366724    0.46448980           NA     ERF11
# 23 AT5G08130 -0.169730460 0.004738418   -0.3754061 1.591577 -0.13639390   -0.27308523           NA      BIM1
# 32 AT5G65590  0.007382096 0.017227683    0.3482814 1.426294  0.03293105    0.32187275           NA AT5G65590 aka SCAP1
# 18 AT3G23240  0.000000000 0.457686819    0.3319695 1.276688  0.01306984    0.35001200           NA     ERF1B
# 9  AT1G80840  0.000000000 0.920894439    0.4096893 1.226684  0.02722907    0.38285714           NA    WRKY40
# 17 AT3G16770  0.005335604 0.022086968    0.2903235 1.150150  0.03497298    0.22698679           NA    RAP2-3
# 16 AT3G15210  0.000000000 0.359911361    0.3449356 1.053483  0.03197530    0.42770708           NA      ERF4

          #----------------------------------------------------------------------------------------
          # 22nd model: 26 TFs with >= 95% motif match within 2000 upstream, 500 downstream of TSS
          # use both matrices mtx and mtx3
          #----------------------------------------------------------------------------------------

targetGene <- canonicalizeName(tp, "UNE10")
tbl.atac <- data.frame(chrom=rep("Chr4", 2),
                       start=c(15862,17863),
                       end=c(17862,18363),
                       stringsAsFactors=FALSE)
tbl.tfs <- findCandidateTranscriptionFactorsByMotifInSequence(tp, tbl.atac, 95L)
candidate.tfs <- sort(unique(tbl.tfs$orf))
length(candidate.tfs)

recipe <- list(title="UNE10",
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
x22 <- build(builder)
lapply(x22, dim)
tbl.model.22 <- x22$model[1:10,]
tf.symbols <- unlist(lapply(tbl.model.22$gene, function(tf) getGeneNames(tp, tf)$symbol))
tbl.model.22$symbol <- tf.symbols

#with bigger matrix:
#         gene  betaLasso   lassoPValue pearsonCoeff  rfScore    betaRidge spearmanCoeff bindingSites    symbol
# 8  AT4G14465 -0.2352269 4.954686e-106   -0.4678254 376.5201 -0.194942574    -0.3154342           NA     AHL20
# 11 AT5G62165 -0.1605824 4.050781e-102   -0.4629587 339.9943 -0.142173176    -0.3505407           NA     AGL42
# 10 AT5G60200 -0.1010548  2.685066e-24   -0.3526809 175.9870 -0.122245187    -0.3010236           NA    DOF5.3
# 12 AT5G65590  0.1221619  8.165302e-37    0.3390436 171.3479  0.119279678     0.3121283           NA AT5G65590 aka SCAP1
# 2  AT1G52150  0.0000000  7.388956e-02   -0.2542754 155.6862 -0.055204319    -0.2906103           NA   ATHB-15
# 4  AT2G22430  0.2062901  8.318967e-26    0.2863719 125.5836  0.224416183     0.2595819           NA    ATHB-6
# 3  AT1G63480 -0.1759239  1.151318e-17   -0.1640331 114.9688 -0.167926041    -0.1373603           NA AT1G63480
# 6  AT2G46590  0.0000000  3.175261e-01   -0.1714005 112.2792  0.007433777    -0.1765432           NA      DAG2
# 9  AT4G40060  0.0000000  5.794193e-01   -0.1947899 103.9999 -0.049056939    -0.2129665           NA   ATHB-16
# 7  AT3G01470  0.0000000  1.467616e-01   -0.1471854 100.7755  0.026399529    -0.1529602           NA    ATHB-1

#with smaller matrix
#         gene  betaLasso  lassoPValue pearsonCoeff  rfScore   betaRidge spearmanCoeff bindingSites    symbol
# 10 AT4G40060 0.38676339 0.0003898992    0.4823011 3.436711  0.18352202     0.4773589           NA   ATHB-16
# 12 AT5G39660 0.19433705 0.0006938689    0.4612668 2.995625  0.12643245     0.3785354           NA      CDF2
# 1  AT1G47655 0.00000000 0.2856213499    0.2405810 2.256596  0.03568010     0.2492677           NA AT1G47655
# 9  AT4G14465 0.00000000 0.2353579376    0.2489781 1.789956  0.03119151     0.1388235           NA     AHL20
# 13 AT5G60200 0.00000000 0.0452494962    0.1770031 1.516223  0.10228559     0.1669628           NA    DOF5.3
# 7  AT3G01470 0.00000000 0.4349516377    0.1756762 1.469617  0.04019154     0.1840576           NA    ATHB-1
# 8  AT3G06740 0.02947161 0.0326759802    0.2637456 1.373348  0.13936260     0.2711710           NA    GATA15
# 4  AT1G63480 0.00000000 0.2213847405   -0.3395003 1.337998 -0.07550148    -0.2319808           NA AT1G63480
# 14 AT5G62165 0.00000000 0.0541638997   -0.3317078 1.224800 -0.03967390    -0.1822329           NA     AGL42
# 11 AT5G13180 0.00000000 0.1539337710    0.1325205 1.158471  0.06127619     0.1330612           NA   ANAC083

          #----------------------------------------------------------------------------------------
          # 23rd model: 29 TFs with >= 95% motif match within 2000 upstream, 500 downstream of TSS
          # use both matrices mtx and mtx3
          #----------------------------------------------------------------------------------------

targetGene <- canonicalizeName(tp, "BEE2")
tbl.atac <- data.frame(chrom=rep("Chr4", 2),
                       start=c(17241698,17243699),
                       end=c(17243698,17244199),
                       stringsAsFactors=FALSE)
tbl.tfs <- findCandidateTranscriptionFactorsByMotifInSequence(tp, tbl.atac, 95L)
candidate.tfs <- sort(unique(tbl.tfs$orf))
length(candidate.tfs)

recipe <- list(title="BEE2",
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
x23 <- build(builder)
lapply(x23, dim)
tbl.model.23 <- x23$model[1:10,]
tf.symbols <- unlist(lapply(tbl.model.23$gene, function(tf) getGeneNames(tp, tf)$symbol))
tbl.model.23$symbol <- tf.symbols
# with bigger matrix
#         gene   betaLasso  lassoPValue pearsonCoeff  rfScore   betaRidge spearmanCoeff bindingSites symbol
# 15 AT5G62165 -0.14654224 8.781013e-62   -0.3672063 386.6629 -0.11614451   -0.28922974           NA  AGL42
# 9  AT4G31920 -0.32184470 8.196087e-52   -0.3513091 343.8435 -0.28954908   -0.32828247           NA  ARR10
# 1  AT1G53160  0.14257230 3.025479e-64    0.3721841 309.8696  0.12612429    0.40741137           NA   SPL4
# 4  AT2G02540  0.16268244 1.387200e-67    0.3798380 261.2975  0.14506505    0.36150086           NA ATHB21
# 14 AT5G60200  0.00000000 5.376887e-02   -0.1971474 205.7040 -0.04499575   -0.25056895           NA DOF5.3
# 3  AT1G69780  0.19690706 1.589081e-17    0.1621895 201.0077  0.18118475    0.15222777           NA ATHB13
# 17 AT5G65310  0.27520837 4.586593e-20    0.1561768 173.0371  0.23328707    0.17642863           NA ATHB-5
# 10 AT4G32890 -0.13829996 7.193778e-26   -0.2485245 159.3232 -0.12355251   -0.25224104           NA  GATA9
# 11 AT5G16560 -0.04288807 2.020327e-07   -0.1283945 148.8646 -0.06416597   -0.08825239           NA    KAN
# 13 AT5G49520 -0.06936169 1.606297e-13   -0.2128161 148.5657 -0.08714175   -0.27128479           NA WRKY48

#with smaller matrix
#         gene   betaLasso  lassoPValue pearsonCoeff   rfScore   betaRidge spearmanCoeff bindingSites    symbol
# 4  AT1G69780 0.475263698 9.071561e-12    0.7900089 6.4627143  0.14877009     0.6877791           NA    ATHB13
# 12 AT4G31920 0.000000000 4.422050e-01   -0.2539478 1.7797840 -0.07590570    -0.2144058           NA     ARR10
# 16 AT5G60200 0.000000000 6.627797e-01   -0.3210403 1.5411595 -0.06501800    -0.2618487           NA    DOF5.3
# 14 AT5G16560 0.009221995 4.351465e-03    0.4371744 1.5307764  0.04171660     0.5045378           NA       KAN
# 6  AT2G22430 0.000000000 8.916344e-01   -0.5026041 1.4761931 -0.11028331    -0.4729412           NA    ATHB-6
# 3  AT1G63480 0.012537756 5.041249e-03    0.1390710 1.3466403  0.06012580     0.2322689           NA AT1G63480
# 5  AT1G80840 0.000000000 8.335613e-01   -0.4675158 1.1657589 -0.03904494    -0.3576951           NA    WRKY40
# 7  AT2G46830 0.000000000 8.288264e-01   -0.2435563 0.9598848 -0.01694524    -0.1338295           NA      CCA1
# 8  AT3G01470 0.000000000 6.096893e-01    0.3443222 0.7905707  0.04480747     0.3185114           NA    ATHB-1
# 10 AT3G15270 0.000000000 5.967193e-01   -0.3595525 0.7188257 -0.02229638    -0.1956783           NA      SPL5

          #----------------------------------------------------------------------------------------
          # 24th model: 47 TFs with >= 95% motif match within 2000 upstream, 500 downstream of TSS
          # use both matrices mtx and mtx3
          #----------------------------------------------------------------------------------------

targetGene <- canonicalizeName(tp, "AGL42")
tbl.atac <- data.frame(chrom=rep("Chr5", 2),
                       start=c(24963074,24965075),
                       end=c(24965074,24965574),
                       stringsAsFactors=FALSE)
tbl.tfs <- findCandidateTranscriptionFactorsByMotifInSequence(tp, tbl.atac, 95L)
candidate.tfs <- sort(unique(tbl.tfs$orf))
length(candidate.tfs)

recipe <- list(title="AGL42",
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
x24 <- build(builder)
lapply(x24, dim)
tbl.model.24 <- x24$model[1:10,]
tf.symbols <- unlist(lapply(tbl.model.24$gene, function(tf) getGeneNames(tp, tf)$symbol))
tbl.model.24$symbol <- tf.symbols
#with bigger mtx
#         gene   betaLasso   lassoPValue pearsonCoeff   rfScore   betaRidge spearmanCoeff bindingSites  symbol
# 14 AT4G14465  0.30093521 2.344620e-130    0.5125279 622.67856  0.18765823     0.4511495           NA   AHL20
# 16 AT4G31920  0.39061638  1.444580e-71    0.4203648 361.38089  0.30817815     0.3924070           NA   ARR10
# 23 AT5G60200  0.20384322  2.415366e-43    0.3851670 285.49324  0.16391718     0.3480110           NA  DOF5.3
# 17 AT4G36540 -0.11625581  3.054137e-41   -0.3672063 281.83242 -0.11514179    -0.2892297           NA    BEE2
# 1  AT1G52150  0.04116239  1.596194e-07    0.2863012 176.18715  0.08530244     0.3785580           NA ATHB-15
# 19 AT5G18830  0.00000000  9.281928e-01    0.2933309 156.37105  0.06125440     0.3515497           NA    SPL7
# 5  AT2G22430 -0.20208881  3.489263e-22   -0.2656049 129.99301 -0.21967644    -0.2478757           NA  ATHB-6
# 18 AT5G08130  0.21389193  1.763819e-10    0.1265168  99.08508  0.22683771     0.1324832           NA    BIM1
# 22 AT5G49520  0.00000000  2.478494e-02    0.1984116  92.59478  0.05805342     0.1869870           NA  WRKY48
# 9  AT2G46130  0.05209329  2.562067e-16    0.3187730  89.91458  0.06879770     0.2726951           NA  WRKY43

#with smaller mtx3
#         gene   betaLasso lassoPValue pearsonCoeff  rfScore   betaRidge spearmanCoeff bindingSites    symbol
# 10 AT3G54810 0.000000000 0.094217551   -0.2644196 6.867721 -0.14413012  -0.200000000           NA     GATA8
# 22 AT5G65590 0.000000000 0.420145989   -0.1392564 4.897552 -0.01964218  -0.325138055           NA AT5G65590 aka SCAP1
# 9  AT3G15270 0.102707048 0.006924263    0.3772227 4.677268  0.07139213   0.198559424           NA      SPL5
# 18 AT5G18830 0.000000000 0.124399548   -0.2049080 4.358970 -0.18306080  -0.091764706           NA      SPL7
# 11 AT4G14465 0.000000000 0.829154291   -0.2942996 4.109162 -0.03895698  -0.171668667           NA     AHL20
# 7  AT2G46130 0.000000000 0.196155136    0.1314209 3.775747  0.02131992  -0.006386555           NA    WRKY43
# 6  AT2G37590 0.006095311 0.039712174    0.2887062 3.765629  0.22335710   0.307274910           NA    DOF2.4
# 14 AT4G32890 0.053780624 0.021209033    0.3203169 3.368515  0.09258236   0.313229292           NA     GATA9
# 1  AT1G52150 0.000000000 0.064953297   -0.2188586 3.185514 -0.15943204  -0.124225690           NA   ATHB-15
# 4  AT1G69780 0.000000000 0.539715018   -0.2050077 2.765162 -0.01455309   0.056902761           NA    ATHB13

          #----------------------------------------------------------------------------------------
          # 25th model: 49 TFs with >= 95% motif match within 2000 upstream, 500 downstream of TSS
          # use both matrices mtx and mtx3
          #----------------------------------------------------------------------------------------

targetGene <- canonicalizeName(tp, "SPL7")
tbl.atac <- data.frame(chrom=rep("Chr5", 2),
                       start=c(6274114,6276115),
                       end=c(6276114,6276615),
                       stringsAsFactors=FALSE)
tbl.tfs <- findCandidateTranscriptionFactorsByMotifInSequence(tp, tbl.atac, 95L)
candidate.tfs <- sort(unique(tbl.tfs$orf))
length(candidate.tfs)

recipe <- list(title="SPL7",
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
x25 <- build(builder)
lapply(x25, dim)
tbl.model.25 <- x25$model[1:10,]
tf.symbols <- unlist(lapply(tbl.model.25$gene, function(tf) getGeneNames(tp, tf)$symbol))
tbl.model.25$symbol <- tf.symbols
# #with bigger mtx
#         gene  betaLasso  lassoPValue pearsonCoeff  rfScore    betaRidge spearmanCoeff bindingSites    symbol
# 3  AT1G52150 0.10143563 3.800640e-56    0.3712817 53.53663 0.0740524145     0.4772063           NA   ATHB-15
# 13 AT4G14465 0.04002867 1.704733e-92    0.4397313 44.82737 0.0275240169     0.4637010           NA     AHL20
# 11 AT3G01470 0.04394715 4.678820e-36    0.3598416 29.24611 0.0535552920     0.4344923           NA    ATHB-1
# 20 AT5G46760 0.18778745 2.933480e-58    0.3725163 28.04176 0.1431090218     0.3628578           NA AT5G46760
# 17 AT4G40060 0.03777721 5.064438e-48    0.3785698 23.80613 0.0440290409     0.4182756           NA   ATHB-16
# 23 AT5G65310 0.08498375 1.317547e-21    0.2923490 21.35004 0.0782008484     0.2661363           NA    ATHB-5
# 2  AT1G49560 0.03637750 2.537876e-08    0.1773350 21.25393 0.0529318604     0.1301963           NA AT1G49560
# 21 AT5G60200 0.00000000 9.504986e-01    0.3036750 17.63110 0.0142387163     0.3464370           NA    DOF5.3
# 14 AT4G34000 0.11099725 2.229469e-36    0.3064559 17.37278 0.0905053732     0.3412263           NA      ABF3
# 22 AT5G62165 0.00000000 4.002926e-01    0.2933309 16.76824 0.0005817714     0.3515497           NA     AGL42

#with smaller mtx3
#         gene   betaLasso lassoPValue pearsonCoeff   rfScore    betaRidge spearmanCoeff bindingSites    symbol
# 26 AT5G65590 0.000000000 0.034660547   -0.2762644 0.5059476 -0.018190685   -0.22487395           NA AT5G65590
# 15 AT4G34000 0.077155041 0.000558144    0.4708817 0.4976902  0.039810202    0.47716687           NA      ABF3
# 8  AT2G46590 0.005551830 0.023016935    0.2930937 0.3723662  0.022891704   -0.02722689           NA      DAG2
# 2  AT1G49560 0.000000000 0.229243435    0.2439457 0.2640093  0.026864711    0.11500600           NA AT1G49560
# 18 AT5G03790 0.000000000 0.037914342    0.2938346 0.2465082  0.016684373    0.30919568           NA   ATHB-51
# 14 AT4G17460 0.000000000 0.736322142   -0.3650563 0.2396181 -0.013813090   -0.27154862           NA      HAT1
# 19 AT5G08130 0.000000000 0.947897616    0.3043247 0.2312637  0.041961957    0.24436975           NA      BIM1
# 21 AT5G40350 0.000000000 0.035969700    0.2738329 0.2240069  0.013022035    0.22641056           NA     MYB24
# 10 AT3G01220 0.000000000 0.612019692   -0.1079931 0.2154041 -0.002190554   -0.08657863           NA    ATHB20
# 11 AT3G47500 0.005363208 0.022094394    0.3311737 0.1923897  0.025414752    0.21142857           NA      CDF3

          #----------------------------------------------------------------------------------------
          # 26th model: 42 TFs with >= 95% motif match within 2000 upstream, 500 downstream of TSS
          # use both matrices mtx and mtx3
          #----------------------------------------------------------------------------------------

targetGene <- canonicalizeName(tp, "ATHB20")
tbl.atac <- data.frame(chrom=rep("Chr3", 2),
                       start=c(71598,73599),
                       end=c(73598,74099),
                       stringsAsFactors=FALSE)
tbl.tfs <- findCandidateTranscriptionFactorsByMotifInSequence(tp, tbl.atac, 95L)
candidate.tfs <- sort(unique(tbl.tfs$orf))
length(candidate.tfs)

recipe <- list(title="ATHB20",
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
x26 <- build(builder)
lapply(x26, dim)
tbl.model.26 <- x26$model[1:10,]
tf.symbols <- unlist(lapply(tbl.model.26$gene, function(tf) getGeneNames(tp, tf)$symbol))
tbl.model.26$symbol <- tf.symbols

#with bigger mtx
#         gene   betaLasso   lassoPValue pearsonCoeff  rfScore   betaRidge spearmanCoeff bindingSites  symbol
# 10 AT4G31920  0.41772711 3.406814e-101    0.4580854 381.0669  0.34747429     0.4738247           NA   ARR10
# 14 AT5G25830  0.19961820  1.318223e-92    0.4445244 363.0023  0.17044102     0.4215165           NA  GATA12
# 8  AT4G16110  0.43145786  7.786821e-97    0.4506042 361.5957  0.35814863     0.4004716           NA    ARR2
# 3  AT1G69780  0.38328304  3.096823e-74    0.3917870 241.0333  0.29992065     0.3124582           NA  ATHB13
# 4  AT2G22430 -0.33991092  1.593083e-39   -0.3110957 123.9119 -0.27593860    -0.2621734           NA  ATHB-6
# 11 AT4G40060  0.16606026  1.371820e-26    0.2986074 121.7611  0.16935473     0.3205553           NA ATHB-16
# 20 AT5G60200  0.06673388  5.472103e-20    0.2569738 118.5619  0.06216835     0.2385958           NA  DOF5.3
# 1  AT1G52150  0.00000000  9.359605e-01    0.2043231 112.4625  0.02076009     0.2827836           NA ATHB-15
# 19 AT5G56270 -0.11981506  3.764489e-21   -0.2444904 105.9166 -0.15869759    -0.2278204           NA   WRKY2
# 12 AT5G12870  0.01968777  2.114678e-11    0.2610498 104.6406  0.04900232     0.2855275           NA   MYB46

#with smaller mtx3
#         gene  betaLasso lassoPValue pearsonCoeff   rfScore   betaRidge spearmanCoeff bindingSites    symbol
# 19 AT5G62000 -0.2862536 0.000191864   -0.5037641 13.112046 -0.11058991   -0.49848739           NA      ARF2
# 8  AT2G23320  0.0000000 0.454605160    0.1141251  8.021538  0.09733246    0.06160864           NA    WRKY15
# 4  AT1G64620  0.0000000 0.156291095    0.2507090  6.958213  0.17797144    0.17887155           NA AT1G64620
# 9  AT2G37590  0.1879416 0.017642741    0.3330993  6.609022  0.34492495    0.33618247           NA    DOF2.4
# 7  AT2G22430  0.0000000 0.670812321   -0.3491920  6.598530 -0.10096907   -0.38477791           NA    ATHB-6
# 18 AT5G60200  0.0000000 0.681120586   -0.2310342  4.441324 -0.10450002   -0.31130852           NA    DOF5.3
# 13 AT4G37180  0.0000000 0.955245566   -0.1257275  4.421769 -0.01297474   -0.11769508           NA AT4G37180
# 2  AT1G08010  0.0000000 0.711167545    0.2947692  3.809657  0.10327200    0.28384154           NA    GATA11
# 5  AT1G69780  0.0000000 0.200130094    0.3239425  3.625423  0.08096407    0.24504202           NA    ATHB13
# 3  AT1G52150  0.0000000 0.779012361   -0.2248019  3.554240 -0.10040755   -0.24504202           NA   ATHB-15

          #----------------------------------------------------------------------------------------
          # 27th model: 43 TFs with >= 95% motif match within 2000 upstream, 500 downstream of TSS
          # use both matrices mtx and mtx3
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
               matrix=mtx3,
               candidateTFs=candidate.tfs,
               tfPool=getAllTranscriptionFactors(tp, "MotifDb"),
               tfPrefilterCorrelation=0.1,
               annotationDbFile=dbfile(org.At.tair.db),
               orderModelByColumn="rfScore",
               solverNames=c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman"),
               quiet=TRUE)

builder <- NoDnaModelBuilder(genome, targetGene, recipe, quiet=TRUE)
x27 <- build(builder)
lapply(x27, dim)
tbl.model.27 <- x27$model[1:10,]
tf.symbols <- unlist(lapply(tbl.model.27$gene, function(tf) getGeneNames(tp, tf)$symbol))
tbl.model.27$symbol <- tf.symbols

#with bigger mtx
#         gene   betaLasso   lassoPValue pearsonCoeff   rfScore   betaRidge spearmanCoeff bindingSites    symbol
# 12 AT5G60200  0.21340936 2.382589e-140    0.5291238 180.25009  0.15742201     0.5104056           NA    DOF5.3
# 11 AT5G59780  0.12134591  9.156398e-89    0.4541901 134.24507  0.10896906     0.5001797           NA     MYB59
# 9  AT5G12870  0.10519445  1.720585e-64    0.4073818 113.38724  0.09952470     0.3741304           NA     MYB46
# 3  AT1G64620  0.17991695  1.974772e-73    0.4309093 104.51558  0.17044489     0.3959355           NA AT1G64620
# 6  AT2G37590  0.05159066  9.584809e-12    0.3314857  79.61398  0.09905095     0.3236278           NA    DOF2.4
# 1  AT1G19850  0.05712985  7.918163e-11    0.2524720  57.30713  0.07933547     0.2915014           NA      ARF5
# 7  AT3G01470  0.00000000  9.624542e-02    0.1885680  52.84702 -0.02167611     0.3273231           NA    ATHB-1
# 15 AT5G65310  0.02036068  2.148224e-03    0.1811693  49.68740  0.07261525     0.2085411           NA    ATHB-5
# 13 AT5G62165  0.00000000  9.062988e-01    0.2863012  42.93329  0.01439715     0.3785580           NA     AGL42
# 4  AT1G80840 -0.03598329  1.281630e-13   -0.1940822  39.17065 -0.04269959    -0.2079710           NA    WRKY40

#with smaller mtx3
# gene  betaLasso  lassoPValue pearsonCoeff   rfScore    betaRidge spearmanCoeff bindingSites symbol
# 20 AT5G62940 0.10801195 0.0006322369    0.4679888 0.8921942  0.065451480   0.135846339           NA DOF5.6
# 1  AT1G19850 0.00000000 0.0500041441    0.2097378 0.8527553  0.040026434   0.374597839           NA   ARF5
# 12 AT4G17460 0.00000000 0.6630970577   -0.1509035 0.7576329 -0.002103981   0.044321729           NA   HAT1
# 16 AT5G40350 0.04182587 0.0018847096    0.4204592 0.7557620  0.028697198   0.386506603           NA  MYB24
# 21 AT5G65310 0.06294694 0.0005086301    0.4738772 0.6652933  0.042665227   0.151500600           NA ATHB-5
# 4  AT2G22430 0.00000000 0.7015060035    0.2227988 0.4832946  0.016838179   0.122785114           NA ATHB-6
# 8  AT3G01470 0.00000000 0.8829875080   -0.1533313 0.3471583 -0.014339995   0.001776711           NA ATHB-1
# 5  AT2G42200 0.00000000 0.8771819918    0.1964781 0.3297144  0.012000219   0.330804322           NA   SPL9
# 10 AT3G15500 0.00000000 0.1752398645   -0.1873589 0.3038921 -0.014458599  -0.170900360           NA NAC055
# 2  AT1G53160 0.00000000 0.6867735006    0.1634598 0.2932131  0.016955481   0.004945978           NA   SPL4

          #----------------------------------------------------------------------------------------
          # 28th model: 72 TFs with >= 95% motif match within 2000 upstream, 500 downstream of TSS
          # use both matrices mtx and mtx3
          #----------------------------------------------------------------------------------------

targetGene <- canonicalizeName(tp, "ATHB6")
tbl.atac <- data.frame(chrom=rep("Chr2", 2),
                       start=c(9522469,9526470),
                       end=c(9526469,9526970),
                       stringsAsFactors=FALSE)
tbl.tfs <- findCandidateTranscriptionFactorsByMotifInSequence(tp, tbl.atac, 95L)
candidate.tfs <- sort(unique(tbl.tfs$orf))
length(candidate.tfs)

recipe <- list(title="ATHB6",
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
x28 <- build(builder)
lapply(x28, dim)
tbl.model.28 <- x28$model[1:10,]
tf.symbols <- unlist(lapply(tbl.model.28$gene, function(tf) getGeneNames(tp, tf)$symbol))
tbl.model.28$symbol <- tf.symbols

#with bigger matrix
# gene   betaLasso  lassoPValue pearsonCoeff  rfScore   betaRidge spearmanCoeff bindingSites    symbol
# 26 AT5G25830 -0.09483151 1.140345e-60   -0.3607472 43.16806 -0.06685793   -0.33905727           NA    GATA12
# 11 AT3G01470  0.22466705 1.211813e-22    0.1711951 42.22054  0.17356091    0.15181251           NA    ATHB-1
# 17 AT3G59060  0.01548082 7.926326e-40    0.3216937 34.07054  0.01428551    0.24764231           NA      PIF5
# 19 AT4G14465 -0.01061005 7.291409e-04   -0.1975364 28.95192 -0.01726298   -0.06137936           NA     AHL20
# 27 AT5G39660 -0.12605465 1.095966e-23   -0.2464525 24.99268 -0.10654957   -0.21454257           NA      CDF2
# 5  AT1G53170  0.05656884 4.104794e-15    0.2205160 24.07504  0.04516757    0.22032783           NA      ERF8
# 3  AT1G50640  0.03543838 9.952316e-13    0.1508054 21.46977  0.05594105    0.11563041           NA      ERF3
# 1  AT1G25550 -0.09614162 2.697151e-09   -0.1148429 21.07468 -0.10490239   -0.06847917           NA AT1G25550
# 30 AT5G62165 -0.01626665 2.727340e-15   -0.2656049 20.32351 -0.01924904   -0.24787569           NA     AGL42
# 24 AT4G37180  0.03308813 3.732311e-08    0.1568264 16.61269  0.04429686    0.16640427           NA AT4G37180

#with smaller matrix
#         gene   betaLasso  lassoPValue pearsonCoeff   rfScore   betaRidge spearmanCoeff bindingSites symbol
# 25 AT3G59060  0.00000000 6.013292e-03    0.4177917 1.4925109  0.04271919     0.3814166           NA   PIF5
# 2  AT1G06160  0.00000000 1.376525e-02   -0.3846738 1.4547389 -0.01850690    -0.5465066           NA ERF094
# 16 AT2G37590 -0.25457628 1.593588e-04   -0.5006159 0.6920075 -0.13581099    -0.3996639           NA DOF2.4
# 32 AT5G16560  0.00000000 4.887221e-01   -0.4302894 0.5912776 -0.01893960    -0.4443217           NA    KAN
# 9  AT1G53170  0.03921353 5.176469e-04    0.5012652 0.5717675  0.03725612     0.5334454           NA   ERF8
# 5  AT1G28370  0.02500702 5.105083e-05    0.5403776 0.5613470  0.01986683     0.4971429           NA  ERF11
# 27 AT4G31800  0.04091388 3.129935e-04    0.5047127 0.5010410  0.02207024     0.4797599           NA WRKY18
# 28 AT4G36540 -0.05821874 3.747159e-04   -0.5026041 0.4860311 -0.04020065    -0.4729412           NA   BEE2
# 34 AT5G25830  0.00000000 6.219192e-03   -0.3594084 0.4298235 -0.02886218    -0.4058103           NA GATA12
# 23 AT3G23240  0.00000000 1.194863e-01    0.4357799 0.3651550  0.02055299     0.4164706           NA  ERF1B

          #----------------------------------------------------------------------------------------
          # 29th model: 50 TFs with >= 95% motif match within 2000 upstream, 500 downstream of TSS
          # use both matrices mtx and mtx3
          #----------------------------------------------------------------------------------------
targetGene <- canonicalizeName(tp, "GATA8")
tbl.atac <- data.frame(chrom=rep("Chr3", 2),
                       start=c(20294956,20296957),
                       end=c(20296956,20297457),
                       stringsAsFactors=FALSE)
tbl.tfs <- findCandidateTranscriptionFactorsByMotifInSequence(tp, tbl.atac, 95L)
candidate.tfs <- sort(unique(tbl.tfs$orf))
length(candidate.tfs)

recipe <- list(title="GATA8",
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
x29 <- build(builder)
lapply(x29, dim)
tbl.model.29 <- x29$model[1:10,]
tf.symbols <- unlist(lapply(tbl.model.29$gene, function(tf) getGeneNames(tp, tf)$symbol))
tbl.model.29$symbol <- tf.symbols

#with bigger mtx
#         gene    betaLasso  lassoPValue pearsonCoeff  rfScore    betaRidge spearmanCoeff bindingSites    symbol
# 5  AT3G06740  0.192238208 5.337048e-36    0.2789892 91.33716  0.163749377    0.09585691           NA    GATA15
# 7  AT4G37260  0.156200032 1.473681e-28    0.2448333 66.03419  0.137607191    0.22819820           NA     MYB73
# 8  AT5G08130  0.174913830 2.326781e-16    0.1988544 60.12522  0.162337804    0.09287319           NA      BIM1
# 3  AT1G80840  0.053747950 2.394974e-15    0.1900307 52.04132  0.054437934    0.16466321           NA    WRKY40
# 4  AT2G23340 -0.082933899 7.309467e-14   -0.1798077 47.05292 -0.079265272   -0.13298754           NA     DEAR3
# 1  AT1G36060  0.009704538 3.783143e-11    0.1615724 46.93718  0.020403626    0.15544550           NA AT1G36060
# 9  AT5G16560  0.030179769 1.116185e-10    0.1482196 42.76163  0.033288525    0.14889430           NA       KAN
# 6  AT4G36540 -0.047619513 1.049123e-10   -0.1226409 41.97699 -0.040842893   -0.14314871           NA      BEE2
# 10 AT5G25810  0.000000000 4.827162e-01   -0.1014489 40.39362 -0.001452256   -0.04638677           NA      TINY
# 2  AT1G77200  0.008600916 2.797777e-04    0.1366130 39.17204  0.017254263    0.12835956           NA AT1G77200

#with smaller mtx3
#         gene betaLasso  lassoPValue pearsonCoeff   rfScore   betaRidge spearmanCoeff bindingSites    symbol
# 15 AT4G37260 0.1708379 8.031956e-07    0.6331905 1.7896975  0.08890945     0.5953902           NA     MYB73
# 14 AT4G36930 0.1438776 5.388445e-05    0.5527549 1.2571274  0.08813774     0.5863625           NA       SPT
# 20 AT5G62940 0.0000000 5.324920e-01   -0.2108329 0.9642275 -0.01741644    -0.1296038           NA    DOF5.6
# 17 AT5G16560 0.0000000 2.780080e-02    0.2876827 0.7651242  0.02414679     0.2155582           NA       KAN
# 10 AT3G16280 0.0000000 7.358323e-01   -0.2400175 0.6439938 -0.01895768    -0.2733733           NA AT3G16280
# 9  AT3G01470 0.0000000 9.614166e-01    0.2603260 0.5203657  0.01925202     0.2310204           NA    ATHB-1
# 8  AT2G40340 0.0000000 1.749215e-01    0.1445581 0.4774949  0.01651070     0.1782953           NA    DREB2C
# 12 AT4G16750 0.0000000 2.535999e-02    0.2803262 0.4495827  0.03932858     0.2551261           NA AT4G16750
# 7  AT2G36270 0.0000000 4.540776e-01    0.2403771 0.4448028  0.01199131     0.1084754           NA      ABI5
# 18 AT5G25810 0.0000000 1.631651e-01    0.2713764 0.4059476  0.01744686     0.1409364           NA      TINY

          #----------------------------------------------------------------------------------------
          # 30th model: 27 TFs with >= 95% motif match within 2000 upstream, 500 downstream of TSS
          # use both matrices mtx and mtx3
          #----------------------------------------------------------------------------------------
targetGene <- canonicalizeName(tp, "GATA15")
tbl.atac <- data.frame(chrom=rep("Chr3", 2),
                       start=c(2124658,2126659),
                       end=c(2126658,2127159),
                       stringsAsFactors=FALSE)
tbl.tfs <- findCandidateTranscriptionFactorsByMotifInSequence(tp, tbl.atac, 95L)
candidate.tfs <- sort(unique(tbl.tfs$orf))
length(candidate.tfs)

recipe <- list(title="GATA15",
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
x30 <- build(builder)
lapply(x30, dim)
tbl.model.30 <- x30$model[1:9,]
tf.symbols <- unlist(lapply(tbl.model.30$gene, function(tf) getGeneNames(tp, tf)$symbol))
tbl.model.30$symbol <- tf.symbols

#with bigger mtx
#       gene   betaLasso  lassoPValue pearsonCoeff   rfScore   betaRidge spearmanCoeff bindingSites symbol
# 1 AT1G42990 -0.33508042 1.351745e-54   -0.3428298 136.45933 -0.28399798   -0.33740759           NA BZIP60
# 3 AT4G16110 -0.11672350 3.434959e-07   -0.1189034 135.80445 -0.11646256   -0.11692962           NA   ARR2
# 2 AT2G22430 -0.10248342 1.859476e-08   -0.1354756 101.82278 -0.10607455   -0.09941991           NA ATHB-6
# 4 AT4G31920  0.07838795 4.844631e-05    0.1150265  93.92049  0.08577908    0.10381406           NA  ARR10
# 5 AT5G39660 -0.04991222 3.873510e-04   -0.1007136  93.64231 -0.06195625   -0.12276830           NA   CDF2

#with smaller mtx3
#       gene  betaLasso  lassoPValue pearsonCoeff   rfScore   betaRidge spearmanCoeff bindingSites    symbol
# 4 AT1G80840 -0.1237662 1.995854e-05   -0.5639998 1.7139994 -0.05265239    -0.5307690           NA    WRKY40
# 7 AT5G49520  0.0000000 8.609518e-01   -0.3637126 1.1599551 -0.02557242    -0.3389277           NA    WRKY48
# 5 AT2G22430  0.0000000 8.372288e-01   -0.2135799 0.9022090 -0.02643684    -0.2063435           NA    ATHB-6
# 1 AT1G42990  0.0000000 3.054324e-01   -0.3711281 0.8352476 -0.05579421    -0.3108358           NA    BZIP60
# 9 AT5G65590  0.0000000 2.991368e-01   -0.1733936 0.6920649 -0.01608158    -0.1782996           NA AT5G65590
# 8 AT5G60200  0.0000000 8.251771e-01   -0.1577859 0.5924579 -0.02278387    -0.1964513           NA    DOF5.3
# 2 AT1G49560  0.0000000 7.170586e-01   -0.1130704 0.5691463 -0.02138932    -0.1672069           NA AT1G49560
# 3 AT1G77920  0.0000000 3.590363e-01   -0.1727587 0.5415838 -0.02227929    -0.1709525           NA AT1G77920
# 6 AT5G39660  0.0000000 6.216311e-01   -0.1109419 0.4314407 -0.01836529    -0.1417081           NA      CDF2


          #----------------------------------------------------------------------------------------
          # 31st model: 29 TFs with >= 95% motif match within 2000 upstream, 500 downstream of TSS
          # use both matrices mtx and mtx3
          #----------------------------------------------------------------------------------------
targetGene <- canonicalizeName(tp, "DOF5.3")
tbl.atac <- data.frame(chrom=rep("Chr5", 2),
                       start=c(24239077,24241078),
                       end=c(24241077,24241578),
                       stringsAsFactors=FALSE)
tbl.tfs <- findCandidateTranscriptionFactorsByMotifInSequence(tp, tbl.atac, 95L)
candidate.tfs <- sort(unique(tbl.tfs$orf))
length(candidate.tfs)

recipe <- list(title="DOF5.3",
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
x31 <- build(builder)
lapply(x31, dim)
tbl.model.31 <- x31$model[1:10,]
tf.symbols <- unlist(lapply(tbl.model.31$gene, function(tf) getGeneNames(tp, tf)$symbol))
tbl.model.31$symbol <- tf.symbols
#with bigger mtx
#         gene    betaLasso   lassoPValue pearsonCoeff   rfScore   betaRidge spearmanCoeff bindingSites    symbol
# 8  AT4G14465  0.214155150 4.916909e-114    0.4832782 239.50481  0.15668213    0.44778587           NA     AHL20
# 6  AT2G37590  0.357741795  1.568545e-90    0.4384045 188.61084  0.31449635    0.42497559           NA    DOF2.4
# 2  AT1G64620  0.211534579  1.091095e-58    0.3815155 165.40192  0.19713607    0.41092741           NA AT1G64620
# 11 AT5G62165  0.070304535  7.690478e-39    0.3851670 130.60854  0.07760712    0.34801100           NA     AGL42
# 4  AT2G18380  0.080615872  2.501859e-21    0.2916731  94.26407  0.09837057    0.28940058           NA    GATA20
# 9  AT4G31920  0.000000000  6.329018e-02    0.1446562  78.76389 -0.04578754    0.20417106           NA     ARR10
# 12 AT5G62940  0.072017247  3.750608e-13    0.2035832  77.73078  0.08849514    0.17916316           NA    DOF5.6
# 10 AT5G25830  0.007394708  1.987848e-03    0.2516686  66.08657  0.04835927    0.26882677           NA    GATA12
# 5  AT2G22430  0.000000000  1.158494e-01   -0.1197237  64.66402  0.03705229   -0.09882664           NA    ATHB-6
# 3  AT1G80840 -0.001951971  1.292859e-02   -0.1299032  57.67901 -0.02536580   -0.13267635           NA    WRKY40

#with smaller mtx3
#         gene   betaLasso lassoPValue pearsonCoeff   rfScore   betaRidge spearmanCoeff bindingSites    symbol
# 4  AT1G80840 0.038016905  0.01034647    0.3594898 0.9927881  0.02919905    0.38986795           NA    WRKY40
# 14 AT5G62940 0.051637742  0.01879464    0.3284058 0.9061619  0.05611126    0.25570228           NA    DOF5.6
# 11 AT4G31920 0.000000000  0.64088629   -0.1142101 0.7730006 -0.02855464   -0.06228091           NA     ARR10
# 2  AT1G64620 0.000000000  0.09435306   -0.2271095 0.6201841 -0.05787497   -0.16600240           NA AT1G64620
# 12 AT5G16560 0.000000000  0.23323122   -0.2672106 0.5990362 -0.01935544   -0.36211285           NA       KAN
# 10 AT4G14465 0.000000000  0.87438075   -0.2027642 0.5986355 -0.01076157   -0.18943577           NA     AHL20
# 9  AT3G15270 0.000000000  0.15020663    0.2614142 0.5646575  0.01586613    0.24081633           NA      SPL5
# 13 AT5G62165 0.003856408  0.05513555    0.2602714 0.4449524  0.02247415    0.28172869           NA     AGL42
# 8  AT3G06740 0.000000000  0.69841228   -0.1577859 0.4156874 -0.02077967   -0.19645130           NA    GATA15
# 5  AT2G18380 0.000000000  0.35127576    0.1476492 0.4000534  0.01118307    0.02809124           NA    GATA20

          #----------------------------------------------------------------------------------------
          # 32nd model: 47 TFs with >= 95% motif match within 2000 upstream, 500 downstream of TSS
          # use both matrices mtx and mtx3
          #----------------------------------------------------------------------------------------
targetGene <- canonicalizeName(tp, "DOF1.8")
tbl.atac <- data.frame(chrom=rep("Chr1", 2),
                       start=c(24005287,24007288),
                       end=c(24007287,24007788),
                       stringsAsFactors=FALSE)
tbl.tfs <- findCandidateTranscriptionFactorsByMotifInSequence(tp, tbl.atac, 95L)
candidate.tfs <- sort(unique(tbl.tfs$orf))
length(candidate.tfs)

recipe <- list(title="DOF1.8",
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
x32 <- build(builder)
lapply(x32, dim)
tbl.model.32 <- x32$model[1:10,]
tf.symbols <- unlist(lapply(tbl.model.32$gene, function(tf) getGeneNames(tp, tf)$symbol))
tbl.model.32$symbol <- tf.symbols

#with bigger mtx
#         gene    betaLasso  lassoPValue pearsonCoeff   rfScore    betaRidge spearmanCoeff bindingSites symbol
# 11 AT5G60200  0.187105963 3.251037e-68    0.3815155 122.88080  0.143461584    0.41092741           NA DOF5.3
# 6  AT2G37590  0.195499776 6.374630e-33    0.3108857 108.78777  0.191921494    0.30756576           NA DOF2.4
# 9  AT4G37790  0.065707743 8.327646e-10    0.1851465  70.06470  0.073489682    0.21960667           NA  HAT22
# 12 AT5G62940  0.023270984 3.497178e-04    0.1522472  61.25121  0.041208829    0.15454709           NA DOF5.6
# 4  AT1G75080  0.000000000 5.510751e-01   -0.1180820  59.31977  0.004938775   -0.09912931           NA   BZR1
# 2  AT1G50640  0.020920599 4.737249e-03    0.1408610  57.08185  0.051081277    0.18047068           NA   ERF3
# 7  AT4G14465  0.009766064 3.472301e-06    0.2410530  55.09684  0.029934282    0.27975820           NA  AHL20
# 8  AT4G16110  0.000000000 6.983879e-01    0.1330114  50.03913  0.013573475    0.21742204           NA   ARR2
# 5  AT2G23340  0.002045172 5.886952e-02    0.1147439  46.57284  0.024416134    0.12424610           NA  DEAR3
# 3  AT1G53170 -0.032093620 1.722291e-05   -0.1726655  44.97815 -0.040581825   -0.17985127           NA   ERF8

#with smaller mtx3
#         gene    betaLasso lassoPValue pearsonCoeff   rfScore    betaRidge spearmanCoeff bindingSites    symbol
# 6  AT1G75080  0.000000000  0.07045807    0.2569291 0.8972444  0.065996112    0.26530612           NA      BZR1
# 18 AT5G65210  0.000000000  0.06127087    0.2385765 0.5332703  0.045748444    0.19788715           NA      TGA1
# 15 AT5G43410  0.008428567  0.02544334    0.3191712 0.5219050  0.022891551    0.28931573           NA AT5G43410
# 8  AT2G22430  0.000000000  0.05771801   -0.2717833 0.5065601 -0.038136257   -0.23697479           NA    ATHB-6
# 12 AT3G23230  0.007146285  0.02207873    0.3231495 0.4882881  0.015521526    0.32235294           NA AT3G23230
# 2  AT1G19210  0.000000000  0.19981887   -0.1720897 0.3969417 -0.014377444   -0.09147659           NA AT1G19210
# 17 AT5G60200  0.000000000  0.09561391   -0.2271095 0.3778582 -0.047719742   -0.16600240           NA    DOF5.3
# 1  AT1G06160  0.000000000  0.67509515    0.1302599 0.3752255  0.003923286    0.06055222           NA    ERF094
# 14 AT4G37790  0.000000000  0.09713631    0.2282825 0.3720663  0.066419448    0.26492197           NA     HAT22
# 13 AT4G16750 -0.016138806  0.02399706   -0.3204023 0.3581817 -0.026578976   -0.34588235           NA AT4G16750


          #----------------------------------------------------------------------------------------
          # 33rd model: 35 TFs with >= 95% motif match within 2000 upstream, 500 downstream of TSS
          # use both matrices mtx and mtx3
          #----------------------------------------------------------------------------------------
targetGene <- canonicalizeName(tp, "AT5G65590")
tbl.atac <- data.frame(chrom=rep("Chr5", 2),
                       start=c(26210083,26212084),
                       end=c(26212083,26212584),
                       stringsAsFactors=FALSE)
tbl.tfs <- findCandidateTranscriptionFactorsByMotifInSequence(tp, tbl.atac, 95L)
candidate.tfs <- sort(unique(tbl.tfs$orf))
length(candidate.tfs)

recipe <- list(title="AT5G65590",
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
x33 <- build(builder)
lapply(x33, dim)
tbl.model.33 <- x33$model[1:10,]
tf.symbols <- unlist(lapply(tbl.model.33$gene, function(tf) getGeneNames(tp, tf)$symbol))
tbl.model.33$symbol <- tf.symbols
#with bigger mtx
#         gene     betaLasso  lassoPValue pearsonCoeff  rfScore   betaRidge spearmanCoeff bindingSites  symbol
# 7  AT4G14465 -1.899293e-01 8.995185e-50   -0.3277134 309.2509 -0.13078474    -0.2360341           NA   AHL20
# 13 AT5G62165 -9.276661e-02 1.709055e-29   -0.2879464 252.7679 -0.08226675    -0.2327215           NA   AGL42
# 12 AT5G60200 -2.111557e-06 8.337027e-02   -0.1989609 204.1028 -0.04499716    -0.1928820           NA  DOF5.3
# 10 AT5G46350  1.105155e-01 2.792354e-10    0.1704134 189.7253  0.12265253     0.1650120           NA   WRKY8
# 1  AT1G52150  0.000000e+00 5.514104e-01   -0.1322944 182.0453  0.01505044    -0.1768898           NA ATHB-15
# 9  AT5G07680 -2.999106e-02 7.140412e-07   -0.2382704 178.5764 -0.04394351    -0.2189192           NA ANAC080
# 3  AT2G23340 -6.532837e-02 1.399493e-05   -0.1729476 177.4505 -0.09238407    -0.1817669           NA   DEAR3
# 5  AT2G46830  1.032759e-02 1.966855e-04    0.1492892 170.3037  0.02629840     0.1534758           NA    CCA1
# 8  AT4G40060  0.000000e+00 9.101853e-01   -0.1860643 168.0156 -0.03380036    -0.1864443           NA ATHB-16
# 2  AT1G80840  4.262706e-02 3.903982e-05    0.1508980 167.8314  0.06769133     0.1357928           NA  WRKY40

#with smaller mtx3
# gene betaLasso lassoPValue pearsonCoeff   rfScore   betaRidge spearmanCoeff bindingSites    symbol
# 7  AT4G16750 0.0000000  0.93452940    0.2044169 11.392812  0.06542909     0.2792317           NA AT4G16750
# 6  AT4G14465 0.0000000  0.20192799    0.1463848 10.096248  0.05083677     0.2896999           NA     AHL20
# 13 AT5G62165 0.0000000  0.56757660   -0.1392564  8.273597 -0.03793328    -0.3251381           NA     AGL42
# 8  AT4G17460 0.0000000  0.75047941   -0.1486456  8.234038 -0.02130059    -0.1741657           NA      HAT1
# 14 AT5G62940 0.0000000  0.37157142    0.1288301  7.951439  0.08740541     0.0657383           NA    DOF5.6
# 11 AT5G46350 0.1457086  0.02556695    0.3156192  7.527796  0.15276984     0.2780792           NA     WRKY8
# 10 AT5G07680 0.0000000  0.47114498    0.1410901  6.358087  0.03637233     0.1155822           NA   ANAC080
# 9  AT4G40060 0.0000000  0.15079891    0.1783725  5.394695  0.18413197     0.2516687           NA   ATHB-16
# 5  AT3G06740 0.0000000  0.33372932   -0.1733936  4.622696 -0.12740711    -0.1782996           NA    GATA15
# 2  AT1G64620 0.0000000  0.41067744   -0.1300946  3.475278 -0.11187940    -0.1661945           NA AT1G64620

          #----------------------------------------------------------------------------------------
          # 34th model: 28 TFs with >= 95% motif match within 2000 upstream, 500 downstream of TSS
          # use both matrices mtx and mtx3
          #----------------------------------------------------------------------------------------
targetGene <- canonicalizeName(tp, "MYC4")
tbl.atac <- data.frame(chrom=rep("Chr4", 2),
                       start=c(9931701,9933702),
                       end=c(9933701,9934202),
                       stringsAsFactors=FALSE)
tbl.tfs <- findCandidateTranscriptionFactorsByMotifInSequence(tp, tbl.atac, 95L)
candidate.tfs <- sort(unique(tbl.tfs$orf))
length(candidate.tfs)

recipe <- list(title="MYC4",
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
x34 <- build(builder)
lapply(x34, dim)
tbl.model.34 <- x34$model[1:10,]
tf.symbols <- unlist(lapply(tbl.model.34$gene, function(tf) getGeneNames(tp, tf)$symbol))
tbl.model.34$symbol <- tf.symbols

#with bigger mtx
#         gene  betaLasso  lassoPValue pearsonCoeff  rfScore  betaRidge spearmanCoeff bindingSites    symbol
# 6 AT5G62000 0.03293216 2.059567e-05    0.1116754 95.07190 0.04064757    0.05220811           NA      ARF2
# 3 AT3G47500 0.10271069 8.896262e-15    0.1744021 88.43036 0.09811551    0.15830419           NA      CDF3
# 2 AT3G06740 0.10785022 2.374480e-09    0.1340536 88.02321 0.10871368    0.11134859           NA    GATA15
# 5 AT5G39660 0.09237757 1.068600e-12    0.1690606 87.22162 0.09327233    0.12990382           NA      CDF2
# 4 AT4G14465 0.07248091 3.511876e-15    0.1774471 84.26109 0.06375901    0.20562914           NA     AHL20
# 7 AT5G62165 0.00709379 2.461137e-02    0.1026428 83.80349 0.01739119    0.12144408           NA     AGL42
# 1 AT1G49560 0.10077876 1.331365e-10    0.1442764 75.67232 0.09751514    0.13397775           NA AT1G49560

#with smaller mtx3
#         gene   betaLasso  lassoPValue pearsonCoeff  rfScore   betaRidge spearmanCoeff bindingSites    symbol
# 5  AT3G06740  0.40927216 0.0002450054    0.4965225 3.209359  0.21101236     0.5101203           NA    GATA15
# 6  AT3G47500  0.10391120 0.0138938349    0.3601382 2.830953  0.11856463     0.4494118           NA      CDF3
# 3  AT2G22430 -0.13576117 0.0131700679   -0.3625533 2.570702 -0.10784324    -0.2595438           NA    ATHB-6
# 8  AT5G16560  0.00000000 0.0756487998    0.3501211 2.274969  0.04497062     0.3777671           NA       KAN
# 9  AT5G39660  0.09777027 0.0079118978    0.3307048 2.254376  0.07885702     0.1379592           NA      CDF2
# 11 AT5G62940  0.00000000 0.1484288502    0.1361369 1.653449  0.05301587     0.1466026           NA    DOF5.6
# 2  AT1G64620  0.00000000 0.2869265118    0.2155694 1.575104  0.07660245     0.1059784           NA AT1G64620
# 4  AT2G37590  0.00000000 0.1679122457    0.2146886 1.373392  0.10137605     0.1597599           NA    DOF2.4
# 10 AT5G60200 -0.01546326 0.0570073327   -0.2469768 1.316233 -0.11839511    -0.1868427           NA    DOF5.3
# 1  AT1G49560  0.00000000 0.7869961954   -0.2202242 1.166599 -0.05085925    -0.2133493           NA AT1G49560

          #----------------------------------------------------------------------------------------
          # 35th model: 48 TFs with >= 95% motif match within 2000 upstream, 500 downstream of TSS
          # use both matrices mtx and mtx3
          #----------------------------------------------------------------------------------------
targetGene <- canonicalizeName(tp, "AT1G25550")
tbl.atac <- data.frame(chrom=rep("Chr1", 2),
                       start=c(8974643,8976644),
                       end=c(8976643,8977143),
                       stringsAsFactors=FALSE)
tbl.tfs <- findCandidateTranscriptionFactorsByMotifInSequence(tp, tbl.atac, 95L)
candidate.tfs <- sort(unique(tbl.tfs$orf))
length(candidate.tfs)

recipe <- list(title="AT1G25550",
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
x35 <- build(builder)
lapply(x35, dim)
tbl.model.35 <- x35$model[1:10,]
tf.symbols <- unlist(lapply(tbl.model.35$gene, function(tf) getGeneNames(tp, tf)$symbol))
tbl.model.35$symbol <- tf.symbols

#with bigger mtx
#         gene     betaLasso  lassoPValue pearsonCoeff  rfScore    betaRidge spearmanCoeff bindingSites    symbol
# 8  AT4G37180  0.1931458598 1.021758e-48    0.3242890 94.59979  0.163751975    0.30083819           NA AT4G37180
# 3  AT1G80840  0.1066981489 4.115582e-43    0.3056062 82.41714  0.088806587    0.30201622           NA    WRKY40
# 4  AT2G22430 -0.1848360668 3.291719e-12   -0.1148429 69.93445 -0.168051347   -0.06847917           NA    ATHB-6
# 7  AT4G36540  0.0522095676 4.995605e-09    0.1523035 60.67394  0.054368338    0.16142347           NA      BEE2
# 6  AT3G06740 -0.0842381150 5.731672e-09   -0.1691969 55.23117 -0.089748937   -0.15614525           NA    GATA15
# 9  AT4G38620  0.0151246394 2.285397e-05    0.1541935 45.27454  0.024560415    0.11727319           NA      MYB4
# 5  AT2G44840  0.0000000000 4.091548e-01    0.1630523 43.35259  0.002397413    0.14819748           NA     ERF13
# 10 AT5G49520  0.0161834476 1.974998e-03    0.1682330 42.05875  0.030072592    0.16698822           NA    WRKY48
# 1  AT1G08000 -0.0340783357 2.511792e-05   -0.1266205 38.88968 -0.028013797   -0.15010301           NA    GATA10
# 2  AT1G08010 -0.0009570583 1.000000e+00   -0.1266205 38.43543 -0.027596918   -0.15010301           NA    GATA11         

#with smaller mtx3
#         gene    betaLasso  lassoPValue pearsonCoeff   rfScore   betaRidge spearmanCoeff bindingSites symbol
# 9  AT3G06740 -0.299907577 6.439146e-05   -0.5342536 1.8577021 -0.15109740    -0.5212130           NA GATA15
# 4  AT1G80840  0.060545281 4.300083e-04    0.5025548 1.7573561  0.04383603     0.5300840           NA WRKY40
# 7  AT2G44840  0.000000000 8.838719e-01    0.3936641 1.1996919  0.02183806     0.4226170           NA  ERF13
# 21 AT5G16560  0.000000000 2.603671e-01   -0.2740614 1.1184633 -0.02773111    -0.3031453           NA    KAN
# 11 AT3G54810  0.009909955 2.564185e-02    0.2725065 0.9697879  0.08001210     0.2621369           NA  GATA8
# 13 AT4G00050  0.000000000 8.322380e-01   -0.1452447 0.9422956 -0.01770847    -0.1826170           NA  UNE10
# 12 AT3G59060  0.000000000 9.752550e-01    0.1444829 0.8469312  0.01177003     0.2466747           NA   PIF5
# 14 AT4G14465  0.000000000 1.443136e-01    0.1358985 0.8118220  0.02571154     0.1400720           NA  AHL20
# 5  AT2G01760  0.000000000 8.058925e-01    0.1600451 0.7753703  0.03482297     0.2307323           NA  ARR14
# 22 AT5G42630  0.013153562 1.765845e-02    0.2956409 0.7692770  0.02849112     0.2895078           NA    ATS

          #----------------------------------------------------------------------------------------
          # 36th model: 47 TFs with >= 95% motif match within 2000 upstream, 500 downstream of TSS
          # use both matrices mtx and mtx3
          #----------------------------------------------------------------------------------------
targetGene <- canonicalizeName(tp, "WRKY40")
tbl.atac <- data.frame(chrom=rep("Chr1", 2),
                       start=c(30381833,30383834),
                       end=c(30383833,30384334),
                       stringsAsFactors=FALSE)
tbl.tfs <- findCandidateTranscriptionFactorsByMotifInSequence(tp, tbl.atac, 95L)
candidate.tfs <- sort(unique(tbl.tfs$orf))
length(candidate.tfs)

recipe <- list(title="WRKY40",
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
x36 <- build(builder)
lapply(x36, dim)
tbl.model.36 <- x36$model[1:10,]
tf.symbols <- unlist(lapply(tbl.model.36$gene, function(tf) getGeneNames(tp, tf)$symbol))
tbl.model.36$symbol <- tf.symbols

#with bigger mtx

# gene   betaLasso   lassoPValue pearsonCoeff  rfScore   betaRidge spearmanCoeff bindingSites  symbol
# 6  AT4G31800  0.51525632 4.964320e-198    0.6100412 702.2947  0.38318834     0.5895222           NA  WRKY18
# 11 AT5G49520  0.20422611  4.091283e-24    0.3789303 312.8799  0.21953091     0.3569045           NA  WRKY48
# 5  AT4G14465 -0.16862619  5.477636e-25   -0.2222530 157.1492 -0.12068589    -0.1928905           NA   AHL20
# 8  AT4G40060 -0.08719258  2.176891e-11   -0.1802024 142.2125 -0.11461538    -0.1942995           NA ATHB-16
# 1  AT1G52150 -0.08624685  1.354613e-16   -0.1940822 138.6247 -0.10105451    -0.2079710           NA ATHB-15
# 7  AT4G31920 -0.02157285  2.252443e-03   -0.2124314 126.4829 -0.07640374    -0.2480439           NA   ARR10
# 4  AT3G54810  0.05375731  6.408844e-07    0.1900307 121.3865  0.11027880     0.1646632           NA   GATA8
# 12 AT5G60200 -0.04284976  4.494388e-07   -0.1299032 116.7443 -0.06072735    -0.1326763           NA  DOF5.3
# 3  AT3G01470  0.00000000  4.510046e-02   -0.2207950 112.6439 -0.08703087    -0.2035526           NA  ATHB-1
# 13 AT5G62165  0.00000000  2.372011e-01   -0.1474643 107.6476 -0.01547067    -0.1469365           NA   AGL42

#with smaller mtx3
#         gene   betaLasso  lassoPValue pearsonCoeff   rfScore   betaRidge spearmanCoeff bindingSites symbol
# 17 AT5G49520  0.30946554 2.050576e-08    0.6956844 16.702504  0.15112743     0.6376471           NA WRKY48
# 11 AT4G31800  0.14937448 7.223464e-06    0.6214864  8.312755  0.10574067     0.5374790           NA WRKY18
# 7  AT3G06740 -0.42878420 1.995821e-05   -0.5639998  5.991193 -0.31067054    -0.5307690           NA GATA15
# 9  AT3G50410  0.00000000 1.557519e-02    0.3614604  4.713995  0.12785945     0.4128211           NA   OBP1
# 15 AT5G16560  0.00000000 8.691337e-02   -0.4332895  3.753391 -0.05814861    -0.5171188           NA    KAN
# 2  AT1G69780 -0.02171626 3.770456e-03   -0.4954645  2.945440 -0.10626775    -0.4257863           NA ATHB13
# 3  AT2G22430  0.00000000 7.659822e-02    0.4055689  2.699469  0.12850483     0.3892917           NA ATHB-6
# 16 AT5G39660  0.00000000 9.408914e-01    0.1232301  2.434100  0.04390803     0.2660744           NA   CDF2
# 5  AT2G46830  0.00000000 5.002095e-01    0.1257317  1.782495  0.01610316     0.2545498           NA   CCA1
# 19 AT5G60200  0.00000000 1.777865e-01    0.3594898  1.689535  0.13931685     0.3898679           NA DOF5.3

          #----------------------------------------------------------------------------------------
          # 37th model: 63 TFs with >= 95% motif match within 2000 upstream, 500 downstream of TSS
          # use both matrices mtx and mtx3
          #----------------------------------------------------------------------------------------
targetGene <- canonicalizeName(tp, "HDG1")
tbl.atac <- data.frame(chrom=rep("Chr3", 2),
                       start=c(22628768,22630769),
                       end=c(22630768,22631269),
                       stringsAsFactors=FALSE)
tbl.tfs <- findCandidateTranscriptionFactorsByMotifInSequence(tp, tbl.atac, 95L)
candidate.tfs <- sort(unique(tbl.tfs$orf))
length(candidate.tfs)

recipe <- list(title="HDG1",
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
x37 <- build(builder)
lapply(x37, dim)
tbl.model.37 <- x37$model[1:10,]
tf.symbols <- unlist(lapply(tbl.model.37$gene, function(tf) getGeneNames(tp, tf)$symbol))
tbl.model.37$symbol <- tf.symbols

#with bigger mtx
#         gene    betaLasso   lassoPValue pearsonCoeff   rfScore   betaRidge spearmanCoeff bindingSites    symbol
# 15 AT4G14465 -0.154765390 2.596098e-115   -0.4856727 237.18354 -0.09040606    -0.4259938           NA     AHL20
# 28 AT5G62165 -0.143733663  8.555141e-99   -0.4655235 194.10680 -0.10233955    -0.4206713           NA     AGL42
# 4  AT1G53170  0.111083778  4.538903e-39    0.3885351 123.54931  0.10857248     0.3412866           NA      ERF8
# 27 AT5G60200 -0.003436034  3.476019e-03   -0.3182049  94.32478 -0.03456264    -0.3386358           NA    DOF5.3
# 14 AT3G16770  0.095642215  1.287627e-23    0.2962935  89.88464  0.09575802     0.3110551           NA    RAP2-3
# 23 AT5G42630  0.000000000  4.896472e-01   -0.2929949  85.81616 -0.02404363    -0.2472149           NA       ATS
# 2  AT1G43160 -0.038001593  1.805203e-35   -0.3635980  72.15085 -0.03897048    -0.3642711           NA    RAP2-6
# 19 AT5G16560  0.098398516  1.032615e-20    0.2124910  70.89221  0.10234046     0.1777708           NA       KAN
# 1  AT1G22985 -0.168688436  1.072297e-27   -0.2819869  65.52232 -0.13920117    -0.2649606           NA AT1G22985
# 26 AT5G59780  0.000000000  9.785459e-01   -0.3170906  61.38923 -0.02360199    -0.3422669           NA     MYB59

#with smaller mtx3
#         gene   betaLasso  lassoPValue pearsonCoeff   rfScore   betaRidge spearmanCoeff bindingSites symbol
# 2  AT1G19850 -0.15160166 0.0001122043   -0.5191075 2.7433970 -0.07520242    -0.6030732           NA   ARF5
# 13 AT3G06740  0.04224990 0.0058084445    0.3901419 1.2922685  0.07388667     0.4125909           NA GATA15
# 5  AT1G28370 -0.02436093 0.0025452370   -0.4350029 0.8931289 -0.02253084    -0.4381753           NA  ERF11
# 10 AT2G33810  0.00000000 0.2602239950   -0.2489882 0.6723020 -0.01894787    -0.3287875           NA   SPL3
# 20 AT4G14465  0.00000000 0.4723964614    0.3278602 0.6388994  0.01115103     0.3110204           NA  AHL20
# 25 AT5G16560  0.00530587 0.0085767546    0.3717312 0.5917463  0.02276504     0.2917167           NA    KAN
# 19 AT3G54810  0.09148248 0.0025801359    0.4166885 0.4226205  0.08816604     0.4079232           NA  GATA8
# 7  AT1G50640  0.00000000 0.0297311200   -0.1685020 0.3406988 -0.04450800    -0.1623529           NA   ERF3
# 15 AT3G15270  0.00000000 0.6326714126   -0.2212659 0.3184976 -0.01069675    -0.1128932           NA   SPL5
# 16 AT3G16770  0.00000000 0.0438706191    0.2069915 0.3161507  0.02201526     0.1376711           NA RAP2-3

          #----------------------------------------------------------------------------------------
          # 38th model: 39 TFs with >= 95% motif match within 2000 upstream, 500 downstream of TSS
          # use both matrices mtx and mtx3
          #----------------------------------------------------------------------------------------
targetGene <- canonicalizeName(tp, "OBP3")
tbl.atac <- data.frame(chrom=rep("Chr3", 2),
                       start=c(20525197,20527198),
                       end=c(20527197,20527698),
                       stringsAsFactors=FALSE)
tbl.tfs <- findCandidateTranscriptionFactorsByMotifInSequence(tp, tbl.atac, 95L)
candidate.tfs <- sort(unique(tbl.tfs$orf))
length(candidate.tfs)

recipe <- list(title="OBP3",
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
x38 <- build(builder)
lapply(x38, dim)
tbl.model.38 <- x38$model[1:10,]
tf.symbols <- unlist(lapply(tbl.model.38$gene, function(tf) getGeneNames(tp, tf)$symbol))
tbl.model.38$symbol <- tf.symbols

#with bigger mtx
#         gene   betaLasso   lassoPValue pearsonCoeff   rfScore   betaRidge spearmanCoeff bindingSites    symbol
# 2  AT1G52150  0.30032344 7.613967e-139    0.5288124 294.84890  0.21362574    0.52041400           NA   ATHB-15
# 20 AT5G60200  0.19500629 1.186395e-142    0.5328210 278.27108  0.15273043    0.49603968           NA    DOF5.3
# 13 AT4G14465  0.21858158 9.676539e-138    0.5271657 245.22635  0.12879812    0.49331909           NA     AHL20
# 8  AT2G03500  0.13145435  3.818057e-36    0.3335143  98.22692  0.11723378    0.34131853           NA AT2G03500
# 21 AT5G62165  0.06357558  1.869604e-32    0.4095408  85.46927  0.06649715    0.42466574           NA     AGL42
# 11 AT2G37590  0.06264759  1.633947e-13    0.3418685  70.43422  0.11443598    0.32495722           NA    DOF2.4
# 19 AT5G42630  0.00000000  6.830392e-01    0.2720615  69.88591  0.01116839    0.29793409           NA       ATS
# 17 AT5G16560 -0.02219210  2.446977e-05   -0.1551953  69.29002 -0.05366961   -0.08497416           NA       KAN
# 4  AT1G64620  0.03786746  1.173939e-09    0.3421884  67.12839  0.08976062    0.33598918           NA AT1G64620
# 12 AT3G01470  0.00000000  6.418641e-01    0.2890155  66.15996  0.04800213    0.34502592           NA    ATHB-1

#with smaller mtx3
#         gene   betaLasso  lassoPValue pearsonCoeff   rfScore   betaRidge spearmanCoeff bindingSites    symbol
# 4  AT1G70920  0.00000000 0.0549188851   -0.3220010 2.4308873 -0.02119084   -0.29114046           NA    ATHB18
# 16 AT5G42630  0.00000000 0.2775032985   -0.2703003 1.5977141 -0.01979709   -0.08235294           NA       ATS
# 6  AT2G01930  0.07720251 0.0096555402    0.3393510 1.5553011  0.06560244    0.24206483           NA      BPC1
# 3  AT1G63480 -0.08890034 0.0027604541   -0.4203422 1.4083491 -0.08789321   -0.33714286           NA AT1G63480
# 2  AT1G52150  0.00000000 0.5393362014    0.1994989 1.2331230  0.04341727   -0.06151261           NA   ATHB-15
# 11 AT4G17460  0.00000000 0.0360625159   -0.3557152 1.2322980 -0.04324252   -0.20364946           NA      HAT1
# 14 AT4G37790  0.00000000 0.0431609962   -0.3080342 1.0922844 -0.10320802   -0.21286915           NA     HAT22
# 20 AT5G65310  0.12564916 0.0005152703    0.4734606 0.9821801  0.06209034    0.14468187           NA    ATHB-5
# 19 AT5G62940  0.00000000 0.1403350982    0.3512926 0.9386974  0.06213477    0.12297719           NA    DOF5.6
# 7  AT2G03500  0.00000000 0.9089545184    0.1293064 0.7918994  0.01891004    0.17051621           NA AT2G03500
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
  runTests()

