library(TrenaProjectArabidopsisRoot)
library(RCyjs)
rcy <- RCyjs(title="WBC19", quiet=TRUE)
tbl <- read.table("highAndLowWBC19regulators.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)[1:6,]
tbl
tbl.big <- read.table("highAndLowWBC19regulators.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)[,1:3]
# build expression table from genes in each model
# tmp <- c("AT3G55130", tbl.model[1,1], tbl.model[3,1], tbl.model[2,1], tbl.model[4,1], tbl.model[6,1], tbl.model[5,1])
# top3 <- highAndLowMtx[tmp, 1:3]
# low3 <- highAndLowMtx[tmp, 51:53]
# testNetworkHL <- cbind(top3, low3)
load("testNetworkHL.RData")
mtx <- testNetworkHL

wbc19ModelHL <- load("WBC19ModelHighLow.1.RData")
rownames(mtx)[2:7] <- c(tbl.model[1,9], tbl.model[3,9], tbl.model[2,9], tbl.model[4,9], tbl.model[6,9], tbl.model[5,9])
rownames(mtx)[1] <- "WBC19"


nodes <- unique(c(tbl.big$Regulators, tbl.big$Target.Gene))
nodes
g <- new("graphNEL", edgemode = "directed")
edgeDataDefaults(g, attr="edgeType") <- "undefined"
nodeDataDefaults(g, attr="expression") <- 0

g <- addNode(nodes, g)
g <- addEdge(tbl.big$Regulators, tbl.big$Target.Gene, g)

edgeData(g, tbl.big$Regulators, tbl.big$Target, "edgeType") <- tbl.big$Regulation

setGraph(rcy, g)
layout(rcy, "cola")

#for(layout in getLayoutStrategies(rcy)){
#   print(layout)
#   layout(rcy, layout)
#   Sys.sleep(2)
#   }

loadStyleFile(rcy, "style.js")
fit(rcy, 20)
layout(rcy, "breadthfirst")

tpar <- TrenaProjectArabidopsisRoot()
get.all.gene.names.for.expression <- getGeneNames(tp, unlist(tbl.big$Regulators))$orf
mtx.expression.for.model.p1 <- highAndLowMtx[get.all.gene.names.for.expression, 1:3]
mtx.expression.for.model.p2 <- highAndLowMtx[get.all.gene.names.for.expression, 51:53]
mtx.expression.for.model <- cbind(mtx.expression.for.model.p1, mtx.expression.for.model.p2)
rownames(mtx.expression.for.model) <- unlist(lapply(rownames(mtx.expression.for.model), function(orf) getGeneNames(tpar, orf)$symbol))

# rownames(mtx2) <- unlist(lapply(rownames(mtx2), function(orf) getGeneNames(tpar, orf)$symbol))
# rownames(mtx2)[6] <- "SCAP1" # symbol for "AT5G65590", however, the symbol does not show in the 

setNodeAttributes(rcy, "expression", rownames(mtx.expression.for.model), as.numeric(mtx.expression.for.model[, 1]))
# sample with highest expression
Sys.sleep(4)
setNodeAttributes(rcy, "expression", rownames(mtx.expression.for.model), as.numeric(mtx.expression.for.model[, 2]))
# sample with second highest expression
Sys.sleep(4)
setNodeAttributes(rcy, "expression", rownames(mtx.expression.for.model), as.numeric(mtx.expression.for.model[, 3]))
# sample with third highest expression
Sys.sleep(4)
setNodeAttributes(rcy, "expression", rownames(mtx.expression.for.model), as.numeric(mtx.expression.for.model[, 4]))
# sample with lowest expression
Sys.sleep(4)
setNodeAttributes(rcy, "expression", rownames(mtx.expression.for.model), as.numeric(mtx.expression.for.model[, 5]))
# sample with second lowest expression
Sys.sleep(4)
setNodeAttributes(rcy, "expression", rownames(mtx.expression.for.model), as.numeric(mtx.expression.for.model[, 6]))
# sample with third high expression




# saveLayout(rcy, "preferredLayout")
# restoreLayout(rcy, "preferredLayout")



