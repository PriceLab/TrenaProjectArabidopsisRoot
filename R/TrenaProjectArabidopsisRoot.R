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
   tbl.names <- get(load(system.file(package="TrenaProjectArabidopsisRoot", "extdata", "misc", "geneIdMap.RData")))
   obj@state$tbl.geneNames <- tbl.names

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
       return(NULL)
       })

#------------------------------------------------------------------------------------------------------------------------

