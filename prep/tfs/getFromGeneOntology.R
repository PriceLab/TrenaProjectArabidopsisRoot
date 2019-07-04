library(GO.db)
library(org.At.tair.db)

go.term <- "GO:0003700"

 # DNA-binding transcription factor activity
 # https://www.ebi.ac.uk/QuickGO/term/GO:0003700

columns(GO.db)  #  "DEFINITION" "GOID"       "ONTOLOGY"   "TERM"
keytypes(GO.db) # "DEFINITION" "GOID"       "ONTOLOGY"   "TERM"
tbl.go <- select(GO.db, keys=go.term, keytype="GOID", columns=columns(GO.db))
as.data.frame(t(select(GO.db, keys=go.term, keytype="GOID", columns=columns(GO.db))))

tbl.go <- select(org.At.tair.db, keys=go.term, columns=c("TAIR", "SYMBOL"), keytype="GOALL")
tfs <- sort(unique(tbl.go$TAIR))
length(tfs)   # 1670
write(tfs, file="../../inst/extdata/geneSets/GO-0003700-DNA-binding-transcription-factor-activity.txt")
