print(load("mtx.zinc.21201x42.RData")) #-loading the file listed within the string and is printing out the file name, which is "mtx"
dim(mtx) #-is displaying the dimensions of the matrix. There are 21201 rows and 42 columns
  # rownames are genes, colnames are samples.
  # is our gene of interest in the matrix?
grep("WBC19", rownames(mtx), value=TRUE, ignore.case=TRUE) #- grep returns the WBC19 vector. this is done by telling R to search in the rownames of the matrix. "value=TRUE" allows you to get the vector, rather than the position of where it occured (rows, columns). "Ignore.case=TRUE" allows you to ignore the case of the vector so it won't be case sensitive
  # extract the vector of expression values, then explore them
wbc19.vector <- mtx["ATWBC19",] #sets the vector for wbc19 by setting values within the matrix and align with the row name "ATWBC19" from the previous line
fivenum(wbc19.vector) # returns the minumum, lower-hinge, median, upper-hinge, and maximum values (there are 42 values) for the wbc19 vector.
plot(wbc19.vector)   # what do you see?  what possibilities does it suggest?
#^ there are no wbc19 samples between 7-9 of the matrix
#^ the first 24 samples range around 6 while the last 18 range around 10
   # now look at another gene of interest, FRD3
grep("FRD3", rownames(mtx), value=TRUE, ignore.case=TRUE)
#grep returns the FRD3 vector from the rownames in the matrix."value=TRUE" allows you to get the vector, rather than the position of where it occured (rows, columns). "Ignore.case=TRUE" allows you to ignore the case of the vector so it won't be case sensitive
  # extract the vector of expression values, then explore them
frd3.vector <- mtx["FRD3",] #sets the vector for frd3 by setting values within the matrix and align with the row name "FRD3" from the previous line
fivenum(frd3.vector) # returns the minumum, lower-hinge, median, upper-hinge, and maximum values (there are 42 values) for the frd3 vector.
plot(frd3.vector)   # different pattern from wbc19?
#^the plot for frd3 is more linear and indicates a relationship
#^ as the sample positions/index increase, the matrix values decrease
#^ all samples range between 6.9-9.4
  # trena uses two techniques for identifying regulatory TF candidates
  #   1) tf is a transcription factor, and has plausible binding site in a regulatory region
  #   2) tf and targetGene are strongly correlated, or strongly anti-correlated in expression
  # the first step in this is to identify arabidopsis transcription factor genes
  # geneontology is a good approach.  read up!  https://en.wikipedia.org/wiki/Gene_ontology
  # and then look at this:  https://www.ebi.ac.uk/QuickGO/term/GO:0003700

