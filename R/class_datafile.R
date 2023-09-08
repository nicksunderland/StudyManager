#' Title
#'
#' @slot path character.
#' @slot data data.table.
#' @slot col_names character.
#' @slot col_types list.
#' @slot col_fill ANY.
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
    col_names = "character",
    col_types = "list",
    col_fill = "ANY"
  ),
  prototype = list(
    path = NA_character_,
    data = data.table::data.table(NULL),
    col_names = character(),
    col_types = list(),
    col_fill = NA
  )
)


setGeneric("free", function(x) standardGeneric("free"))
setMethod("free", "DataFile", function(x) {
  data.table::set(x@data, j=names(x@data), value=NULL)
  validObject(x)
  x
})

setGeneric("write_file", function(x, file_path, append, sep, ...) standardGeneric("write_file"))
setMethod("write_file", "DataFile", function(x, file_path, append=FALSE, sep="\t", ...) {
  data.table::fwrite(x@data, file=file_path, append=append, sep=sep, ...)
})


setGeneric("free", function(x) standardGeneric("free"))
setMethod("free", "DataFile", function(x) {
  data.table::set(x@data, j=names(x@data), value=NULL)
  validObject(x)
  x
})


setGeneric("extract", function(x) standardGeneric("extract"))
setMethod("extract", "DataFile", function(x) {

  if(ncol(x@data) == 0) {

    # see if missing columns
    if(length(x@col_names) == 0) {

      x@data <- data.table::fread(x@path)

    } else {

      h <- data.table::fread(x@path, nrows=0)
      header <- colnames(h)
      missing_cols <- setdiff(x@col_names, header)
      non_missing_cols <- x@col_names[!x@col_names %in% missing_cols]

      if(length(missing_cols) > 0 & length(x@col_fill) == 0) {

        stop("Error importing data, missing columns with no fill provided")

      }

      x@data <- data.table::fread(x@path, select=unname(non_missing_cols))

      if(length(missing_cols) > 0) {

        x@data[, (missing_cols) := rep(x@col_fill, length(missing_cols))]

      }

      if(length(x@col_types) > 0 & all(sapply(x@col_types, is.character))) {

        type.convert <- list(
          "character" = as.character,
          "numeric" = as.numeric,
          "logical" = as.logical,
          "integer" = as.integer,
          "complex" = as.complex,
          "raw" = as.raw
        )

        for(i in seq_along(x@col_types)) {
          col_name <- names(x@col_types)[i]
          x@data[, (col_name) := type.convert[[ x@col_types[[col_name]] ]](get(col_name))]
        }

      } else if(length(x@col_types) > 0 & all(sapply(x@col_types, is.function))) {

        for(i in seq_along(x@col_types)) {
          col_name <- names(x@col_types)[i]
          x@data[, (col_name) := x@col_types[[i]](get(col_name))]
        }

      }

      # new names provided
      if(length(x@col_names) > 0 & rlang::is_named(x@col_names)) {

        data.table::setnames(x@data, x@col_names, names(x@col_names))

      }

    } # end else col_names provided

  } # end if data doesn't exist

  validObject(x)
  x
})


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


setGeneric("path", function(x) standardGeneric("path"))
setMethod("path", "DataFile", function(x) x@path)

setGeneric("path<-", function(x, value) standardGeneric("path<-"))
setMethod("path<-", "DataFile", function(x, value) {
  x@path <- value
  validObject(x)
  x
})


setGeneric("col_names", function(x) standardGeneric("col_names"))
setMethod("col_names", "DataFile", function(x) x@col_names)

setGeneric("col_names<-", function(x, value) standardGeneric("col_names<-"))
setMethod("col_names<-", "DataFile", function(x, value) {
  x@col_names <- value
  validObject(x)
  x
})


setGeneric("col_types", function(x) standardGeneric("col_types"))
setMethod("col_types", "DataFile", function(x) x@col_types)

setGeneric("col_types<-", function(x, value) standardGeneric("col_types<-"))
setMethod("col_types<-", "DataFile", function(x, value) {
  x@col_types <- value
  validObject(x)
  x
})


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

    # col full
    if(length(object@col_fill) > 0) {

      stopifnot("`col_fill` must be a length==1 atomic type" = is.atomic(object@col_fill) & length(object@col_fill)==1)

    }

    # col_names
    if(length(object@col_names) > 0) {

      stopifnot("`col_names` must be a character vector (names optional)" = is.character(object@col_names))

    }

    # col_types
    if(length(object@col_types) > 0) {

      # col_types must have names
      stopifnot("`col_types` must be a named list" = rlang::is_named(object@col_types))

      # names must be in the col_names
      stopifnot("names(`col_types`) and `col_names` must be the same" = setequal(names(object@col_types), unname(object@col_names)))

      # col_types must be valid atomic character, or function
      if(all(sapply(object@col_types, is.character))) {

        stopifnot("If `col_types` are characters, they must an atomic type: 'character', 'numeric', 'integer', 'logical', 'raw', 'complex'" = all(object@col_types %in% c('character', 'numeric', 'integer', 'logical', 'raw', 'complex')))

      } else if(all(sapply(object@col_types, is.function))) {

        #TODO: testing functions work

      } else {

        stop("`col_types` should all be a list of characters or a list of functions to apply to each column")

      }

    }

    # both col_types and col_names
    if(length(object@col_types) > 0 | length(object@col_names) > 0) {

      if(ncol(object@data) == 0) {

        tryCatch(
          expr = {
            header <- data.table::fread(object@path, nrows=0)
            cols <- colnames(header)
          },
          error = function(e){
            stop(e)
          }
        )
      } else {

        cols <- colnames(object@data)

      }

      # col names must be in the data file
      if(length(object@col_fill) == 0) {

        stopifnot("All `col_names` must be valid column names, or a `col_fill` provided" = all(object@col_names %in% cols) | all(names(object@col_names) %in% cols))

      } else {

        stopifnot("At least one of `col_names` must be a valid column name" = any(object@col_names %in% cols))

      }

    }
  }

)


#
# d <- DataFile(
#   path = "/Users/nicholassunderland/Downloads/hermes_progression/bioshift_triumph/bioshift_triumph.females.allcause.gz"
#   # col_names = c("PVAL"="P", "POSITION"="POS", "Foo"="Foo"),
#   # col_types = list("POS"=function(x) paste(x,"test"), "P"=as.character, "Foo"=as.character),
#   # col_fill = NA_real_
# )
#
# d <- extract(d)
# foo<-data(d)
# head(data(d))



