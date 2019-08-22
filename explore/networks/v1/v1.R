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
tbl <- read.table("wbc19.regulators.root.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)[1:4,]
tbl


nodes <- unique(c(tbl$Regulators, tbl$Target.Gene))
nodes
g <- new("graphNEL", edgemode = "directed")
edgeDataDefaults(g, attr="edgeType") <- "undefined"
nodeDataDefaults(g, attr="expression") <- 0

g <- addNode(nodes, g)
g <- addEdge(tbl$Regulators, tbl$Target.Gene, g)

edgeData(g, tbl$Regulators, tbl$Target, "edgeType") <- tbl$Regulation

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

rownames(mtx) <- unlist(lapply(rownames(mtx), function(orf) getGeneNames(tpar, orf)$symbol))
rownames(mtx)[4] <- "AT5G65590"

setNodeAttributes(rcy, "expression", rownames(mtx), as.numeric(mtx[, 1]))
setNodeAttributes(rcy, "expression", rownames(mtx), as.numeric(mtx[, 2]))
setNodeAttributes(rcy, "expression", rownames(mtx), as.numeric(mtx[, 3]))






