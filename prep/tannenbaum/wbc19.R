# All datasets generated in this study are available in GEO (GSE122772)
# (https:// www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122772).
#------------------------------------------------------------------------------------------------------------------------
library(igvR)
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

library (RColorBrewer)
colors <- brewer.pal(8, "Dark2")
totalColorCount <- length(colors)
currentColorNumber <- 0
#------------------------------------------------------------------------------------------------------------------------
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
# displayBedTrack("WBC19")
