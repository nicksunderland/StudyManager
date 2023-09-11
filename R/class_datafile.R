#' Title
#'
#' @slot path character.
#' @slot data data.table.
#' @slot mapping ColMap
#'
#' @return a DataFile object
#' @export
#' @importFrom methods new validObject
#' @importFrom rlang :=
#' @importFrom checkmate test_path_for_output
#' @importFrom data.table data.table fread fwrite set is.data.table
#'
DataFile <- setClass(
  Class = "DataFile",
  slots = list(
    path = "character",
    data = "data.table",
    active_cols = "character",
    mapping = "ColMap"
  ),
  prototype = list(
    path = NA_character_,
    data = data.table::data.table(NULL),
    active_cols = character(),
    mapping = StudyManager::ColumnMapping
  )
)


setMethod(
  f = "initialize",
  signature = "DataFile",
  definition = function(.Object, path=NA_character_,
                                 data=data.table::data.table(NULL),
                                 active_cols=character(),
                                 mapping=StudyManager::ColumnMapping) {

    if(ncol(data) > 0) {

      header <- colnames(data)

    } else if(!is.na(path)) {

      # get the current file column names
      h <- data.table::fread(path, nrows=0)
      header <- colnames(h)

    } else {

      stop("Initializse DataFile path or data input error")

    }

    .Object@data <- data
    .Object@path <- path
    .Object@mapping <- mapping
    .Object@active_cols <- parse_col_names(.Object@mapping, header)

    return(.Object)
  }
)


setGeneric("free", function(x) standardGeneric("free"))
setMethod("free", "DataFile", function(x) {
  data.table::set(x@data, j=names(x@data), value=NULL)
  validObject(x)
  x
})

setGeneric("write_file", function(x, file_path, append=FALSE, sep="\t", ...) standardGeneric("write_file"))
setMethod(
  f = "write_file",
  signature = c("DataFile", "character"),
  definition = function(x, file_path, append=FALSE, sep="\t", ...) {

    if(length(data) > 0) {

      data.table::fwrite(x@data, file=file_path, append=append, sep=sep, ...)

    } else {

      warning("Data empty - do you need to call `extract` first?")

    }
  }
)


setGeneric("free", function(x) standardGeneric("free"))
setMethod("free", "DataFile", function(x) {
  data.table::set(x@data, j=names(x@data), value=NULL)
  validObject(x)
  x
})


setGeneric("extract", function(x) standardGeneric("extract"))
setMethod(
  f = "extract",
  signature = "DataFile",
  definition = function(x) {

    if(ncol(x@data) > 0) {

      warning("Already extracted data, if you want to re-extract first call free() to delete the current data.table")

    } else {

      # get missing columns
      missing_cols <- names(x@active_cols[is.na(x@active_cols)])
      non_missing_cols <- x@active_cols[!is.na(x@active_cols)]

      if(length(missing_cols) > 0 & is.null(col_fill(x@mapping))) {

        stop("Error importing data, missing columns with no fill provided")

      }

      x@data <- data.table::fread(x@path, select=unname(non_missing_cols))

      # new names
      data.table::setnames(x@data, non_missing_cols, names(non_missing_cols))

      if(length(missing_cols) > 0) {

        x@data[, (missing_cols) := col_fill(x@mapping)]

      }

      col_type_check_funcs <- col_type_funcs(x@mapping)

      for(i in seq_along(col_type_check_funcs)) {
          col_name <- names(col_type_check_funcs)[i]
          x@data[, (col_name) := col_type_check_funcs[[i]](get(col_name))]
      }

    }
    validObject(x)
    x
  }
)


#' data
#'
#' @param x object
#'
#' @return a data.table
#' @export
#'
setGeneric("file_data", function(x) standardGeneric("file_data"))
#' Title
#'
#' @param x obj
#'
#' @return object
#' @export
#'
setMethod("file_data", "DataFile", function(x) x@data)


setGeneric("set_data<-", function(x, value) standardGeneric("set_data<-"))
setMethod("set_data<-", "DataFile", function(x, value) {
  x@data <- value
  validObject(x)
  x
})


setGeneric("file_name", function(x) standardGeneric("file_name"))
setMethod("file_name", "DataFile", function(x) basename(x@path))


setGeneric("path", function(x) standardGeneric("path"))
setMethod("path", "DataFile", function(x) x@path)


setGeneric("path<-", function(x, value) standardGeneric("path<-"))
setMethod("path<-", "DataFile", function(x, value) {
  x@path <- value
  validObject(x)
  x
})


setGeneric("mapping", function(x) standardGeneric("mapping"))
setMethod("mapping", "DataFile", function(x) x@mapping)


setGeneric("mapping<-", function(x, value) standardGeneric("mapping<-"))
setMethod("mapping<-", "DataFile", function(x, value) {
  x@mapping <- value
  validObject(x)
  x
})

setGeneric("active_cols", function(x) standardGeneric("active_cols"))
setMethod("active_cols", "DataFile", function(x) x@active_cols)


setValidity(
  Class = "DataFile",
  method = function(object) {

    # path
    if(is.na(object@path)) {

      stopifnot("If `path` is NA/not given, `data` must be a valid data.table" = data.table::is.data.table(object@data))

    }

    # path and no data
    if(!is.na(object@path) & ncol(object@data) == 0) {

      stopifnot("`path` must be valid file path" = file.exists(object@path))

    }

    # path and data
    if(!is.na(object@path) & ncol(object@data) > 0) {

      stopifnot("`path` must be a valid and writable file path if `data` present" = checkmate::test_path_for_output(object@path, overwrite = TRUE))

    }

    # mapping object
    stopifnot("Mapping must be a valid ColMap object" = inherits(object@mapping, "ColMap"))
    stopifnot("No active mapping columns, adjust the `ColMap@active` slot" = length(object@active_cols) > 0)

  }
)





