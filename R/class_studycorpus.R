#' Title
#'
#' @slot corpus_dir character.
#' @slot study_type character.
#' @slot studies list.
#'
#' @return StudyCorpus object
#' @export
#' @include class_gwas_sumstats.R
#'
StudyCorpus <- setClass(
  Class = "StudyCorpus",
  contains = "GWASsumstats",
  slots = list(
    corpus_dir = "character",
    study_type = "character",
    studies = "list"
  ),
  prototype = list(
    corpus_dir = character(),
    study_type = character(),
    studies = list()
  )
)

setMethod(
  f = "initialize",
  signature = "StudyCorpus",
  definition = function(.Object, corpus_dir, study_type, ...) {

    print("init StudyCorpus")
    print(names(as.list(match.call()[-1])))

    # peel off the StudyCorpus specific slots
    .Object@corpus_dir <- corpus_dir
    .Object@study_type <- study_type

    # send the rest down to GWASsumstats
    .Object <- callNextMethod(.Object, ...)

    # create the list of studies from the directory
    .Object <- create_study_list(.Object, study_type)

    validObject(.Object)
    return(.Object)
  }
)

# setValidity(
#   Class = "StudyCorpus",
#   method = function(object) {
#
#     valid_study_types <- c('GWASsumstats')
#     stopifnot("`dir` must be a valid directory" = dir.exists(object@corpus_dir))
#     stopifnot("`study_type` must be one of `c('GWASsumstats')" = object@study_type %in% valid_study_types)
#
#   }
# )

setGeneric("create_study_list", function(object, study_type) standardGeneric("create_study_list"))
setMethod(
  f = "create_study_list",
  signature = c("StudyCorpus", "character"),
  definition = function(object, study_type) {

    object@studies <- list()

    study_directories <- list.dirs(object@corpus_dir, recursive=FALSE)

    for(study_dir in study_directories) {

      if(study_type == "GWASsumstats") {

        study <- GWASsumstats(dir = study_dir,
                              pre_qc_dir = object@pre_qc_dir,
                              post_qc_dir = object@post_qc_dir,
                              file_structure = object@file_structure,
                              ref_path = object@ref_path)

        object@studies[[basename(study_dir)]] <- study

      } else if(study_type == "AnotherType") {

        ##

      } else {

        rlog::log_error("create_study_list(StudyCorpus) error: shouldn't be here, validity has failed.")
        stop("create_study_list() `study_type` error")

      }

    }

    validObject(object)
    return(object)
  }
)
