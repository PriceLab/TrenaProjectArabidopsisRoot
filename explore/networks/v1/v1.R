# biocGet("RCyjs")
#------------------------------------------------------------------------------------------------------------------------
#       gene   betaLasso lassoPValue pearsonCoeff  rfScore   betaRidge spearmanCoeff     symbol
#  AT3G54810 -0.03315773  0.08019973   -0.2644196 8.106659 -0.16569641   -0.20000000      GATA8
#  AT2G37590  0.19950526  0.04200000    0.2887062 6.680566  0.24873045    0.30727491     DOF2.4
#  AT5G65590  0.00000000  0.24946301   -0.1392564 6.477805 -0.03827349   -0.32513806  AT5G65590
#  AT4G14465 -0.06861681  0.03802656   -0.2942996 6.385959 -0.06734547   -0.17166867      AHL20
#  AT5G18830  0.00000000  0.12425412   -0.2049080 6.275481 -0.23378909   -0.09176471       SPL7
#  AT1G69780  0.00000000  0.95502237   -0.2050077 4.572243 -0.04328460    0.05690276     ATHB13
#------------------------------------------------------------------------------------------------------------------------
library(TrenaProjectArabidopsisRoot)
library(RCyjs)
rcy <- RCyjs(title="WBC19", quiet=TRUE)
tbl <- read.table("wbc19.regulators.root.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)[1:6,]
tbl
tbl.big <- read.table("wbc19.regulators.root.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)[,1:3]
tbl.big$Regulators[3] <- "SCAP1"
load("testNetwork.mtx.RData")
load("testNetwork1.mtx.RData")
mtx2 <- testNetwork1.mtx
mtx <- testNetwork.mtx
wbc19Model <- load("WBC19Model.small.RData")
rownames(tbl.model.17)

rownames(mtx2)[2:7] <- tbl.model.17$symbol
rownames(mtx2)[1] <- "WBC19"
rownames(mtx2)[6] <- "SCAP1"

nodes <- unique(c(tbl.big$Regulators, tbl.big$Target.Gene))
nodes
g <- new("graphNEL", edgemode = "directed")
edgeDataDefaults(g, attr="edgeType") <- "undefined"
nodeDataDefaults(g, attr="expression") <- 0

g <- addNode(nodes, g)
g <- addEdge(tbl.big$Regulators, tbl.big$Target.Gene, g)

edgeData(g, tbl.big$Regulators, tbl.big$Target, "edgeType") <- tbl.big$Regulation

setGraph(rcy, g)
layout(rcy, "breadthfirst")
layout(rcy, "cola")

#for(layout in getLayoutStrategies(rcy)){
#   print(layout)
#   layout(rcy, layout)
#   Sys.sleep(2)
#   }

loadStyleFile(rcy, "style.js")
fit(rcy, 20)

tpar <- TrenaProjectArabidopsisRoot()
# rownames(mtx2) <- unlist(lapply(rownames(mtx2), function(orf) getGeneNames(tpar, orf)$symbol))
# rownames(mtx2)[6] <- "SCAP1" # symbol for "AT5G65590", however, the symbol does not show in the 

setNodeAttributes(rcy, "expression", rownames(mtx.expression.for.model), as.numeric(mtx.expression.for.model[, 1]))
Sys.sleep(3)
setNodeAttributes(rcy, "expression", rownames(mtx.expression.for.model), as.numeric(mtx.expression.for.model[, 2]))
Sys.sleep(3)
setNodeAttributes(rcy, "expression", rownames(mtx.expression.for.model), as.numeric(mtx.expression.for.model[, 3]))


get.all.gene.names.for.expression <- getGeneNames(tp, unlist(tbl.big$Regulators))$orf
mtx.expression.for.model <- mtx3[get.all.gene.names.for.expression, 1:3]
rownames(mtx.expression.for.model) <- unlist(lapply(rownames(mtx.expression.for.model), function(orf) getGeneNames(tpar, orf)$symbol))

# saveLayout(rcy, "preferredLayout")
# restoreLayout(rcy, "preferredLayout")



