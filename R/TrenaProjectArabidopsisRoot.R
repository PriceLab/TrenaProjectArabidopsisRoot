#----------------------------------------------------------------------------------------------------
#' @import methods
#' @import TrenaProject
#' @importMethodsFrom TrenaProject getAllTranscriptionFactors

#' @importFrom AnnotationDbi select
#' @import org.Hs.eg.db
#'
#' @title TrenaProjectArabidopsisRoot-class
#'
#' @name TrenaProjectArabidopsisRoot-class
#' @rdname TrenaProjectArabidopsisRoot-class
#' @aliases TrenaProjectArabidopsisRoot
#' @exportClass TrenaProjectArabidopsisRoot
#'

.TrenaProjectArabidopsisRoot <- setClass("TrenaProjectArabidopsisRoot",
                                  contains="TrenaProject")

#----------------------------------------------------------------------------------------------------
#' Define an object of class TrenaProjectArabidopsisRoot
#'
#' @description
#' Expression, variant and covariate data for the genes of interest (perhaps unbounded) for pre-term birth studies
#'
#' @rdname TrenaProjectArabidopsisRoot-class
#'
#' @export
#'
#' @return An object of the TrenaProjectArabidopsisRoot class
#'

TrenaProjectArabidopsisRoot <- function(quiet=TRUE)
{
   genomeName <- "tair10"

   directory <- system.file(package="TrenaProjectArabidopsisRoot", "extdata", "geneSets")
   geneSet.files <- list.files(directory)
   geneSets <- list()
   for(file in geneSet.files){
      full.path <- file.path(directory, file)
      genes <- scan(full.path, sep="\t", what=character(0), quiet=TRUE)
      geneSet.name <- sub(".txt", "", file)
      geneSets[[geneSet.name]] <- genes
      }

   geneInfoTable.path <- system.file(package="TrenaProjectArabidopsisRoot", "extdata", "geneInfoTable_tair10.RData")
   stopifnot(file.exists(geneInfoTable.path))

   footprintDatabaseNames <- NA_character_;
   footprintDatabaseHost <- NA_character_;
   dataDirectory <- system.file(package="TrenaProjectArabidopsisRoot", "extdata")

   stopifnot(file.exists(dataDirectory))

   obj <- .TrenaProjectArabidopsisRoot(TrenaProject("TrenaProjectArabidopsisRoot",
                                                    supportedGenes=geneSets[["rootMetalTransport"]],
                                                    genomeName=genomeName,
                                                    geneInfoTable.path=geneInfoTable.path,
                                                    footprintDatabaseHost=footprintDatabaseHost,
                                                    footprintDatabaseNames=footprintDatabaseNames,
                                                    packageDataDirectory=dataDirectory,
                                                    quiet=quiet
                                                    ))
   tbl.orfsAndSyms <- get(load(system.file(package="TrenaProjectArabidopsisRoot", "extdata", "misc", "geneIdMap.RData")))
   obj@state$tbl.geneNames <- tbl.orfsAndSyms
   obj@state$tbl.aliases <- get(load(system.file(package="TrenaProjectArabidopsisRoot", "extdata", "misc", "tbl.aliases.RData")))

   obj

} # TrenaProjectArabidopsisRoot, the constructor
#----------------------------------------------------------------------------------------------------
#' get all genes annotated by GO to
#'
#' @rdname getAllTranscriptionFactors
#' @aliases getAllTranscriptionFactors
#'
#' @param obj An object of class TrenaProject
#'
#' @export

setMethod('getAllTranscriptionFactors', 'TrenaProjectArabidopsisRoot',

   function(obj, source) {
      source <- tolower(source)
      stopifnot(source %in% c("geneontology", "gene ontology", "motifdb"))
      if(grepl("ontology", source)){
         full.path <- system.file(package="TrenaProjectArabidopsisRoot", "extdata", "geneSets",
                                 "GO-0003700-DNA-binding-transcription-factor-activity.txt")
         return(scan(full.path, sep="\t", what=character(0), quiet=TRUE))
         }
      else{ # must be MotifDb
         tfs.mdb <-mcols(query(MotifDb, c("athali", "jaspar2018")))$geneSymbol
         tfs.mdb2 <- lapply(tfs.mdb, function(tf) getGeneNames(tp, tf)$orf)
         return(sort(unique(unlist(tfs.mdb2))))
         }
      })

#------------------------------------------------------------------------------------------------------------------------
setGeneric('getGeneNames', signature='obj', function(obj, name) standardGeneric ('getGeneNames'))
setGeneric('canonicalizeName', signature='obj', function(obj, name) standardGeneric ('canonicalizeName'))
setGeneric('findCandidateTranscriptionFactors', signature='obj', function(obj, tbl.regions, pwmMatchMinimumAsPercentage)
              standardGeneric ('findCandidateTranscriptionFactors'))

#------------------------------------------------------------------------------------------------------------------------
#' get orf and geneSymbol for orf or geneSymbol
#'
#' @rdname getGeneNames
#' @aliases getGeneNames
#'
#' @param obj An object of class TrenaProject
#' @param name An character string, either an orf or a gene symbol
#'
#' @export

setMethod('getGeneNames', 'TrenaProjectArabidopsisRoot',

    function(obj, name) {
       name <- tolower(name)
       tbl <- obj@state$tbl.geneNames
       symbol.match <- match(name, tolower(tbl$geneSymbol), nomatch=FALSE)
       orf.match    <- match(name, tolower(tbl$orf), nomatch=FALSE)
       if(symbol.match)
          return(list(symbol=tbl$geneSymbol[symbol.match], orf=tbl$orf[symbol.match]))
       if(orf.match)
          return(list(symbol=tbl$geneSymbol[orf.match], orf=tbl$orf[orf.match]))
       if(grepl("^at", name))
          return(list(symbol=toupper(name), orf=toupper(name)))
       return(NULL)
       })

#------------------------------------------------------------------------------------------------------------------------
#' map every possible gene name and/or symbol to a standard AT orf name
#'
#' @rdname canonicalizeName
#' @aliases canonicalizeName
#'
#' @param obj An object of class TrenaProject
#' @param name A character string: gene symbol, orf name, odd alias
#'
#' @export

setMethod('canonicalizeName', 'TrenaProjectArabidopsisRoot',

      function(obj, name) {
        tbl.aliases <- obj@state$tbl.aliases
        tbl.orfSym <- obj@state$tbl.geneNames
        #browser()
        #xyz <- "in canoncializeName"
        if(!name %in% tbl.orfSym$geneSymbol){
           if(grepl("^AT[A-Z]", "ATWBC19"))
              name <- sub("^AT", "", name)
           }
        name <- strsplit(name, " ")[[1]][1]
        if(grepl("^AT[1-5]", name))
           return(name)
        if(name %in% tbl.orfSym$geneSymbol){
           index <- match(name, tbl.orfSym$geneSymbol)
           return(tbl.orfSym$orf[index])
           }
        if(name %in% tbl.aliases$symbol){
           index <- match(name, tbl.orfSym$geneSymbol)
           return(tbl.orfSym$orf[index])
           }
        index <- grep(name, tbl.aliases$tokens)
        if(length(index) > 0)
           return(tbl.aliases$locus_name[index])
        return(name)
        })

#------------------------------------------------------------------------------------------------------------------------
#' identify transcription factors (using orf names) with motif match above threshold in the specified regions
#'
#' @description
#' uses the "Bioconductor motif matcher" aka Biostrings::matchPWM as wrapped by trena's MotifMatcher class
#'
#' @rdname findCandidateTranscriptionFactors
#' @aliases findCandidateTranscriptionFactors
#'
#' @param obj An object of class TrenaProject
#' @param tbl.regions A data.frame with chrom, start and end columsn
#' @param pwmMatchMinimumAsPercentage a long integer (eg, 95L) expression minimum match threshold percentage
#'
#' @export

setMethod('findCandidateTranscriptionFactors', 'TrenaProjectArabidopsisRoot',

    function(obj, tbl.regions, pwmMatchMinimumAsPercentage){

       capitalize <- function(chr.name){
          paste0(toupper(substr(chr.name, 1, 1)), substr(chr.name, 2, nchar(chr.name)))
          }
       add.orfs.to.motif.table <- function(tbl.motifs){
         geneSymbols <- unlist(lapply(tbl.motifs$motifName, function(motifName) mcols(MotifDb[motifName])$geneSymbol))
         orfs <- unlist(lapply(geneSymbols, function(geneSymbol) canonicalizeName(tp, geneSymbol)))
         tbl.motifs$orf <- orfs
         tbl.motifs
         }

      pfms <- query(MotifDb, c("athaliana", "jaspar2018"))
      mm <- MotifMatcher("tair10", as.list(pfms), quiet=TRUE)
      tbl.motifs <- findMatchesByChromosomalRegion(mm, tbl.regions, pwmMatchMinimumAsPercentage=95L)
      if(nrow(tbl.motifs) == 0)
         return(vector(mode="character", length=0))
      tbl.motifs <- add.orfs.to.motif.table(tbl.motifs)
      unique(tbl.motifs$orf)
      })

#------------------------------------------------------------------------------------------------------------------------

