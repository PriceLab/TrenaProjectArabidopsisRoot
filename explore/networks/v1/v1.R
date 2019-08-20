# biocGet("RCyjs")
library(RCyjs)
rcy <- RCyjs(title="WBC19", quiet=TRUE)
tbl <- read.table("wbc19.regulators.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)

nodes <- unique(c(tbl$Regulators, tbl$Target.Gene))
nodes
g <- new("graphNEL", edgemode = "directed")
edgeDataDefaults(g, attr="edgeType") <- "undefined"

g <- addNode(nodes, g)
g <- addEdge(tbl$Regulators, tbl$Target.Gene, g)

edgeData(g, tbl$Regulators, tbl$Target, "edgeType") <- tbl$Regulation

setGraph(rcy, g)
getLayoutNames(rcy)
layout(rcy, "breadthfirst")
layout(rcy, "cola")
for(layout in getLayoutStrategies(rcy)){
   print(layout)
   layout(rcy, layout)
   Sys.sleep(2)
   }
loadStyleFile(rcy, "style.js")
fit(rcy, 20)


