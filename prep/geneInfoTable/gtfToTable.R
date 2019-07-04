# post hoc reconstruction of the commands needed to create a legit geneInfoTable
tbl <- import("Arabidopsis_thaliana.TAIR10.44.gtf.gz")
tbl.at <- as.data.frame(tbl)
tbl.at <- subset(tbl.at, type=="gene")
colnames(tbl.at)[11] <- "geneSymbol"
colnames(tbl.at)[1] <- "chr"
tbl.at$tss <- tbl.at$start
revStrand <- which(tbl.at$strand == "-")
tbl.at$tss[revStrand] <- tbl.at$end[revStrand]
tbl.at$chr <- as.character(tbl.at$chr)
tbl.at$type <- as.character(tbl.at$type)
tbl.at$source <- as.character(tbl.at$source)
tbl.at$strand[tbl.at$strand == 2] <- -1


save(tbl.at, file="../../inst/extdata/geneInfoTable_tair10.RData")
