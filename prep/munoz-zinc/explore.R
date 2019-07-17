print(load("mtx.zinc.21201x42.RData")) 
#-loading the file listed within the string and is printing out the file name, which is "mtx"
dim(mtx) 
#-is displaying the dimensions of the matrix. There are 21201 rows and 42 columns
  # rownames are genes, colnames are samples.
  # is our gene of interest in the matrix?
grep("WBC19", rownames(mtx), value=TRUE, ignore.case=TRUE) 
  #- grep returns the WBC19 vector. this is done by telling R to search in the rownames of the matrix. "value=TRUE" allows you to get the vector, rather than the position of where it occured (rows, columns). "Ignore.case=TRUE" allows you to ignore the case of the vector so it won't be case sensitive
  # extract the vector of expression values, then explore them
wbc19.vector <- mtx["ATWBC19",] 
#sets the vector for wbc19 by setting values within the matrix and align with the row name "ATWBC19" from the previous line
fivenum(wbc19.vector) 
# returns the minumum, lower-hinge, median, upper-hinge, and maximum values (there are 42 values) for the wbc19 vector.
plot(wbc19.vector)   
# what do you see?  what possibilities does it suggest?
#^ there are no wbc19 samples between 7-9 of the matrix
#^ the first 24 samples range around 6 while the last 18 range around 10
   # now look at another gene of interest, FRD3
grep("FRD3", rownames(mtx), value=TRUE, ignore.case=TRUE)
#grep returns the FRD3 vector from the rownames in the matrix."value=TRUE" allows you to get the vector, rather than the position of where it occured (rows, columns). "Ignore.case=TRUE" allows you to ignore the case of the vector so it won't be case sensitive
  # extract the vector of expression values, then explore them
frd3.vector <- mtx["FRD3",] 
#sets the vector for frd3 by setting values within the matrix and align with the row name "FRD3" from the previous line
fivenum(frd3.vector) 
# returns the minumum, lower-hinge, median, upper-hinge, and maximum values (there are 42 values) for the frd3 vector.
plot(frd3.vector)   #
different pattern from wbc19?
#^the plot for frd3 is more linear and indicates a relationship
#^ as the sample positions/index increase, the matrix values decrease
#^ all samples range between 6.9-9.4
  # trena uses two techniques for identifying regulatory TF candidates
  #   1) tf is a transcription factor, and has plausible binding site in a regulatory region
  #   2) tf and targetGene are strongly correlated, or strongly anti-correlated in expression
  # the first step in this is to identify arabidopsis transcription factor genes
  # geneontology is a good approach.  read up!  https://en.wikipedia.org/wiki/Gene_ontology
# and then look at this:  https://www.ebi.ac.uk/QuickGO/term/GO:0003700

# plot WBC19 expression - and this time use covariate data to color the points
tbl.cov <- read.table("experimentalVariables.tsv", sep="\t", header=TRUE, as.is=TRUE)

# plot all samples with the same color, same solid dot
plot(wbc19.vector, main="WBC19 expression", col="red", pch=16)

# recall the covariate column names in preparation for coloring by zinc treatment: -, +, ++
head(tbl.cov)
zincTreatmentColor <- rep("black", nrow(tbl.cov))
zinc.plus <- which(tbl.cov$zinc == "+")
zinc.plus.plus <- which(tbl.cov$zinc == "++")
zincTreatmentColor[zinc.plus] <- "pink"
zincTreatmentColor[zinc.plus.plus] <- "red"
plot(wbc19.vector, main="WBC19 expression", pch=16, col=zincTreatmentColor)

# now color by Root/Shoot tissue differences
tissueColor <- rep("black", nrow(tbl.cov))
tissue.root <- which(tbl.cov$tissue == "Root")
tissue.shoot <- which(tbl.cov$tissue == "Shoot")
tissueColor[tissue.root] <- "brown"
tissueColor[tissue.shoot] <- "#2E8925"
plot(wbc19.vector, main="WBC19 expression", pch=16, col=tissueColor)

# what pch (point character) values are available?
max <- 25
plot(1:max, pch=1:max)

# set shape to either  16 (solid dot) or 7 (an x in a rectangle)
genetics.shape <- rep(16, nrow(tbl.cov))
mutants <- which(tbl.cov$genetics == "Mut")
genetics.shape[mutants] <- 7

plot(wbc19.vector, main="WBC19 expression", col=tissueColor, pch=genetics.shape)

