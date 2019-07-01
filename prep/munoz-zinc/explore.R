print(load("mtx.zinc.21201x42.RData"))
dim(mtx)
  # rownames are genes, colnames are samples.
  # is our gene of interest in the matrix?
grep("WBC19", rownames(mtx), value=TRUE, ignore.case=TRUE)

  # extract the vector of expression values, then explore them
wbc19.vector <- mtx["ATWBC19",]
fivenum(wbc19.vector)
plot(wbc19.vector)   # what do you see?  what possibilities does it suggest?

   # now look at another gene of interest, FRD3
grep("FRD3", rownames(mtx), value=TRUE, ignore.case=TRUE)

  # extract the vector of expression values, then explore them
frd3.vector <- mtx["FRD3",]
fivenum(frd3.vector)
plot(frd3.vector)   # different pattern from wbc19?

  # trena uses two techniques for identifying regulatory TF candidates
  #   1) tf is a transcription factor, and has plausible binding site in a regulatory region
  #   2) tf and targetGene are strongly correlated, or strongly anti-correlated in expression
  # the first step in this is to identify arabidopsis transcription factor genes
  # geneontology is a good approach.  read up!  https://en.wikipedia.org/wiki/Gene_ontology
  # and then look at this:  https://www.ebi.ac.uk/QuickGO/term/GO:0003700

