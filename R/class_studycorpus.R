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

    call <- paste0(names(as.list(match.call()[-1])),collapse=', ')
    rlog::log_trace(glue::glue("Init StudyCorpus: call({call})"))

    # peel off the StudyCorpus specific slots
    .Object@corpus_dir <- corpus_dir
    .Object@study_type <- study_type

    # send the rest down to GWASsumstats
    .Object <- callNextMethod(.Object, ...)

    # create the list of studies from the directory
    rlog::log_trace(glue::glue("Creating StudyCorpus study list..."))
    .Object <- create_study_list(.Object, study_type)
    rlog::log_trace(glue::glue("Finished creating StudyCorpus study list"))

    validObject(.Object)
    rlog::log_trace(glue::glue("Init StudyCorpus: complete"))
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


setMethod(
  f = "show",
  signature = "StudyCorpus",
  definition = function(object) {

    file_path_list <- function(df_list_el) {
      out <- vector("list", length(df_list_el))
      names(out) <- names(df_list_el)
      for(i in seq_along(df_list_el)) {
        if(is.list(df_list_el[[i]])) {
          out[[i]] <- file_path_list(df_list_el[[i]])
        }
      }
      return(results)
    }

    l <- file_path_list(object@data_files)
    print(Hmisc::list.tree(l, size=F))

  }
)

