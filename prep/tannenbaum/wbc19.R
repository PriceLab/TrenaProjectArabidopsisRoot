# All datasets generated in this study are available in GEO (GSE122772)
# (https:// www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122772).
#------------------------------------------------------------------------------------------------------------------------
library(igvR)
library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
library (RColorBrewer)
colors <- brewer.pal(8, "Dark2")
totalColorCount <- length(colors)
currentColorNumber <- 0
#------------------------------------------------------------------------------------------------------------------------
# library(FimoClient); FIMO_HOST <- "khaleesi"; FIMO_PORT <- 60029;
# fc <- FimoClient(FIMO_HOST, FIMO_PORT, quiet=FALSE)
#------------------------------------------------------------------------------------------------------------------------
library(trena)
library(MotifDb)

library(TrenaProjectArabidopsisRoot)
tp <- TrenaProjectArabidopsisRoot()
library(org.At.tair.db)
db <- org.At.tair.db
library(BSgenome.Athaliana.TAIR.TAIR9)
library(TxDb.Athaliana.BioMart.plantsmart28)
txdb <- TxDb.Athaliana.BioMart.plantsmart28

tbl.at <- select(db, keys=keys(db), columns=c("SYMBOL", "GENENAME"))
dim(tbl.at)  # 36537 3
if(!exists("igv")){
   igv <- igvR()
   setGenome(igv, "tair10")
   setBrowserWindowTitle(igv, "tair10")
   }

if(!exists("mtx")){
   mtx <- getExpressionMatrix(tp, "munoz.zinc.18993x42.canonicalRowNames")
   dim(mtx)
   }
displayBedTrack <- function(geneSymbol, shoulder=10000)
{
   geneID <- tbl.at[grep(geneSymbol, tbl.at$SYMBOL),"TAIR"][1]
   printf("%s: %s", geneSymbol, geneID)
   tbl.gene <- select(txdb, keys=geneID, columns=columns(txdb), keytype="GENEID")
   deleters <- which(is.na(tbl.gene$CDSSTART))
   if(length(deleters) > 0)
      tbl.gene <- tbl.gene[-deleters,]

   start.loc <- min(tbl.gene$CDSSTART, na.rm=TRUE) - shoulder
   end.loc <- max(tbl.gene$CDSEND, na.rm=TRUE) + shoulder
   chrom.loc <- sprintf("%s:%d-%d", tbl.gene$CDSCHROM[1], start.loc, end.loc)
   showGenomicRegion(igv, chrom.loc)

   files <- c("GSM3484757_Root_ATAC_1_peaks.narrowPeak",
              "GSM3484758_Root_ATAC_2_peaks.narrowPeak",
              "GSM3484759_Root_ATAC_3_peaks.narrowPeak")

   x <- getGenomicRegion(igv)

   for(i in 1:3){
      tbl.atac <- read.table(files[i], sep="\t", as.is=TRUE, nrow=-1)
      colnames(tbl.atac) <- c("chrom", "start", "end", "type", "score", "strand", "fc", "pScore", "qScore", "peakOffset")
      tbl.atac$chrom <- as.character(tbl.atac$chrom)
      print(dim(tbl.atac))
      x$chrom <- sub("chr", "", x$chrom)
      tbl.atac.gene <- subset(tbl.atac, chrom==x$chrom & start >= x$start & end <= x$end)[, c(1,2,3,5)]
      print(dim(tbl.atac.gene))
      track.name <- sprintf("ATAC.%d", i)
      #browser()
      currentColorNumber <<- (currentColorNumber %% totalColorCount) + 1
      color <- colors[currentColorNumber]
      track <- DataFrameQuantitativeTrack(track.name, tbl.atac.gene, color=color, autoscale=TRUE)
      displayTrack(igv, track)
      } # for i


} # displayBedTrack
#------------------------------------------------------------------------------------------------------------------------
capitalize <- function(chr.name)
{
   paste0(toupper(substr(chr.name, 1, 1)), substr(chr.name, 2, nchar(chr.name)))

} # capitalize
#------------------------------------------------------------------------------------------------------------------------
add.orfs.to.motif.table <- function(tbl.motifs)
{
   geneSymbols <- unlist(lapply(tbl.motifs$motifName, function(motifName) mcols(MotifDb[motifName])$geneSymbol))
   orfs <- unlist(lapply(geneSymbols, function(geneSymbol) canonicalizeName(tp, geneSymbol)))
   tbl.motifs$orf <- orfs

   tbl.motifs

} # add.orfs.to.motif.table
#------------------------------------------------------------------------------------------------------------------------
test_add.orfs.to.motif.table <- function()
{
   printf("--- test_add.orfs.to.motif.table")
     # create a small tbl.motifs, querying a small region just upstream of WBC19

   tbl.region <- data.frame(chrom="Chr3", start=20438630, end=20438660, stringsAsFactors=FALSE)
   pfms <- query(MotifDb, c("athaliana", "jaspar2018"))
   mm <- MotifMatcher("tair10", as.list(pfms), quiet=TRUE)

   tbl.motifs <- findMatchesByChromosomalRegion(mm, tbl.region, pwmMatchMinimumAsPercentage=98L)
   checkEquals(dim(tbl.motifs), c(2, 13))
   tbl.motifsWithOrfs <- add.orfs.to.motif.table(tbl.motifs)
   checkEquals(dim(tbl.motifsWithOrfs), c(2, 14))

   tbl.motifs <- findMatchesByChromosomalRegion(mm, tbl.region, pwmMatchMinimumAsPercentage=95L)
   checkEquals(dim(tbl.motifs), c(7, 13))
   tbl.motifsWithOrfs <- add.orfs.to.motif.table(tbl.motifs)
   checkEquals(dim(tbl.motifsWithOrfs), c(7, 14))

   mapped.to.orfs <- grep("^AT[1-5]", tbl.motifsWithOrfs$orf, value=TRUE)
   checkEquals(length(mapped.to.orfs), 6)
   map.failures <- setdiff(tbl.motifsWithOrfs$orf, mapped.to.orfs)
   checkEquals(length(map.failures), 1)

} # test_add.orfs.to.motif.table
#------------------------------------------------------------------------------------------------------------------------
matchMotifs <- function(display=FALSE)
{
   roi <- getGenomicRegion(igv)
   pfms <- query(MotifDb, c("athaliana", "jaspar2018"))
   length(pfms)
   mm <- MotifMatcher("tair10", as.list(pfms), quiet=TRUE)

   chrom <- capitalize(roi$chrom)
   start <- roi$start
   end <- roi$end
   tbl.regions <- data.frame(chrom=chrom, start=start, end=end, stringsAsFactors=FALSE)
     # printf("--- about to get TAIR10 sequence for %d bases", 1 + tbl.regions$end - tbl.regions$start)
     #   seq <- getSequence(mm, tbl.regions)
     # chrom    start      end                                                                                                   seq status
     #  Chr3 20438539 20438639 CTAAAAATAGGAGAAATAAGTGCAATTTCAAAGGAGACGTTAATGGCAAGTTTTTCTCTTCCAAAATTGACTAGTCGATGATGTTAATTGTTAATGTTGTC     wt

   tbl.motifs <- findMatchesByChromosomalRegion(mm, tbl.regions, pwmMatchMinimumAsPercentage=95L)
   dim(tbl.motifs)
   tbl.motifs <- add.orfs.to.motif.table(tbl.motifs)

   candidate.tfs <- unique(tbl.motifs$orf)
   candidate.tfs <- intersect(candidate.tfs, rownames(mtx))
   length(candidate.tfs)

   lapply(candidate.tfs, function(tf) cor(mtx["AT3G55130", ], mtx[tf,]))
   tbl.regions <-
   sequence <- as.list(as.character(getSeq(BSgenome.Athaliana.TAIR.TAIR9, chrom, start, end)))
   #requestMatch(fc, sequence, pvalThreshold=0.00001)
   #if(display){


} # matchMotifs
#------------------------------------------------------------------------------------------------------------------------
findCorrelatedExpression <- function(target="ATWBC19")
{
   x <- apply(mtx.munoz, 1, function(row) cor(mtx.munoz[target,], row))
}
#------------------------------------------------------------------------------------------------------------------------
# displayBedTrack("WBC19")
