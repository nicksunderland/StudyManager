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
    mapping = "ColMap"
  ),
  prototype = list(
    path = NA_character_,
    data = data.table::data.table(NULL),
    mapping = StudyManager:::ColumnMapping
  )
)

setValidity(
  Class = "DataFile",
  method = function(object) {

    stopifnot("`path`(s) must be writable file path(s)"
              = all(sapply(object@path, checkmate::test_path_for_output, overwrite=TRUE)))
    stopifnot("`data` must be a data.table" = data.table::is.data.table(object@data))
    stopifnot("`mapping` must be a valid ColMap object" = inherits(object@mapping, "ColMap"))

  }
)

# see generic set in 'ColMap' class
# setGeneric("col_names", function(x) standardGeneric("col_names"))
setMethod("col_names", "DataFile", function(x) names(col_map(x@mapping, only.active=TRUE)) )

# see generic set in 'ColMap' class
# setGeneric("col_map", function(x, ...) standardGeneric("col_map"))
setMethod(
  f = "col_map",
  signature = "DataFile",
  definition = function(x) {

    # get file header
    header <- c()
    for(p in x@path) {

      h <- data.table::fread(p, nrows=0)
      cols <- colnames(h)

      if(any(duplicated(cols))) {

        rlog::log_fatal(glue::glue("Duplicate column names [{paste0(cols[duplicated(cols)],collapse=',')}] found in {basename(x@path)}"))
        stop()
      }

      header <- unique(c(header, cols))

    }
    raw_cols = header

    # use the 'ColMap' col_names signature method
    active_col_map <- col_map(x@mapping,
                              input_col_names=raw_cols,
                              only.active=TRUE,
                              ignore.case=TRUE)

    return(active_col_map)
  }
)

setGeneric("path", function(x) standardGeneric("path"))
setMethod("path", "DataFile", function(x) x@path)

setGeneric("path<-", function(x, value) standardGeneric("path<-"))
setMethod("path<-", "DataFile", function(x, value) {

  x@path <- value

  if(data_exists(x)) {

    warning("DataFile path changed with existing data extracted -> existing data cleared, you should re-extract")
    x <- free(x)

  }

  validObject(x)
  x
})

setGeneric("file_name", function(x) standardGeneric("file_name"))
setMethod("file_name", "DataFile", function(x) basename(x@path))

setGeneric("mapping", function(x) standardGeneric("mapping"))
setMethod("mapping", "DataFile", function(x) x@mapping)
setMethod("mapping", "list", function(x) lapply(x, mapping))

setGeneric("mapping<-", function(x, value) standardGeneric("mapping<-"))
setMethod(
  f = "mapping<-",
  signature = "list",
  definition = function(x, value) {

    stopifnot("mapping value must be length==1 or length==length(x)" = length(value)==1 | length(value)==length(x))

    value <- rep(as.list(value), length.out=length(x))

    for(i in seq_along(x)) {

      mapping(x[[i]]) <- value[[i]]

    }
    return(x)
  }
)
setMethod("mapping<-", "DataFile", function(x, value) {

  x@mapping <- value

  if(data_exists(x)) {

    warning("DataFile mapping changed with existing data extracted -> existing data cleared, you should re-extract")
    x <- free(x)

  }

  validObject(x)
  x
})

setGeneric("data_exists", function(x) standardGeneric("data_exists"))
setMethod(
  f = "data_exists",
  signature = "DataFile",
  definition = function(x){

    if(ncol(x@data) > 0) {
      return(TRUE)
    } else {
      return(FALSE)
    }

  }
)

#' @title get_data
#' @param x a DataFile object
#' @return a data.table
#' @rdname get_data
#' @export
setGeneric("get_data", function(x, cols=NULL) standardGeneric("get_data"))
#' @rdname get_data
setMethod(
  f = "get_data",
  signature = "DataFile",
  definition = function(x, cols=NULL){

    if(data_exists(x)) {

      # cols to return, default all
      active_cols <- col_names(x)

      # Insert col_fill columns for any missing entries in mapping vector
      current_cols <- colnames(x@data)
      missing_cols <- setdiff(active_cols, current_cols)

      for(col in missing_cols) {

        x@data[, (col) := col_fill(x@mapping)]

      }

      # optional col select
      if(!is.null(cols)) {

        stopifnot("`cols` must be a character vector of column names" = all(sapply(cols, is.character)))

        new_dt <- x@data[, ..cols]

        return(new_dt)

      } else {

        return(x@data)

      }

    } else {

      warning("Data empty - do you need to call `extract()` first?")
      return(NULL)

    }
  }
)

setGeneric("set_data<-", function(x, value) standardGeneric("set_data<-"))
setMethod("set_data<-", "DataFile", function(x, value) {
  x@data <- value
  validObject(x)
  x
})


setGeneric("write_file", function(x, file_path, overwrite=FALSE, append=FALSE, sep="\t", ...) standardGeneric("write_file"))
setMethod(
  f = "write_file",
  signature = c("DataFile", "character"),
  definition = function(x, file_path, overwrite=FALSE, append=FALSE, sep="\t", ...) {

    if(!data_exists(x)) {

      rlog::log_warn("Trying to write from empty DataFile - do you need to call `extract()` first?")
      return(NULL)

    }

    if(!checkmate::test_path_for_output(file_path, overwrite=overwrite)) {

      rlog::log_error("If file_path exists you must confirm overwrite == TRUE ")
      stop()

    }

    data.table::fwrite(x@data, file=file_path, append=append, sep=sep, ...)
  }
)

setGeneric("free", function(x) standardGeneric("free"))
setMethod("free", "list", function(x) lapply(x, free))
setMethod("free", "DataFile", function(x) {
  if(data_exists(x)) {
    x@data[, names(x@data) := NULL]
  }
  validObject(x)
  x
})










setGeneric("extract", function(x, ...) standardGeneric("extract"))
setMethod(
  f = "extract",
  signature = "DataFile",
  definition = function(x) {

    if(length(x@path) > 1) {

      rlog::log_warn(glue::glue("Multiple file paths set, calling extract on all and merging: {paste0(basename(x@path),collapse=', ')}"))

      dfile_list <- lapply(x@path, function(p, m) DataFile(path=p, mapping=m), m=x@mapping)

      x <- extract(dfile_list) # different signature of 'list'

      return(x)

    } else {

      rlog::log_info(glue::glue("Extract from DataFile: ../{basename(dirname(x@path))}/{basename(x@path)}"))

      # active columns
      active_col_map <- col_map(x)
      rlog::log_debug(glue::glue("Active mapped columns: {paste0(names(active_col_map),'=',active_col_map,collapse=',')}"))

      # non-missing cols in the datafile
      is_missing <- is.na(active_col_map)
      missing_cols <- active_col_map[is_missing]
      non_missing_cols <- active_col_map[!is_missing]
      active_col_map[is_missing] <- names(active_col_map)[is_missing]

      # column types
      col_type <- col_types( x@mapping )[!is_missing]
      raw_name_and_type <- stats::setNames(col_type, active_col_map[names(col_type)])

      # potential custom functions
      funcs <- col_type_funcs( x@mapping )[!is_missing]

      # read the data (only store non-missing, as when we request data we'll fill it in later)
      rlog::log_info(glue::glue("Reading file..."))
      x@data <- data.table::fread(x@path, select=raw_name_and_type)
      rlog::log_debug(glue::glue("Finished reading file: {nrow(x@data)} rows by {ncol(x@data)} cols"))

      # add missing
      if(length(missing_cols) > 0) {

        rlog::log_debug(glue::glue("Missing mapping columns: {paste0(names(missing_cols),'=',missing_cols, collapse=',')}. These will be filled with {col_fill(x@mapping)}"))

        for(std_name in names(missing_cols)) {

          x@data[, (std_name) := col_fill(x@mapping)]

        }

      }

      # set standard names
      data.table::setnames(x@data, active_col_map, names(active_col_map))

      # enforce column types
      rlog::log_debug(glue::glue("Ensuring column types:  {paste0(names(col_type),'=',col_type,collapse=',')}"))
      x@data[ , names(funcs) := Map(function(f, x) f(x), funcs, .SD), .SDcols=names(funcs)]
    }

    validObject(x)
    return(x)
  }
)

setMethod(
  f = "extract",
  signature = "list",
  definition = function(x, merge_col="JOIN") {

    stopifnot(length(x) > 0)

    names = names(x)
    if(is.null(names)) {
        names = as.character(seq_along(x))
    }

    rlog::log_info(glue::glue("Extract from multiple ({length(x)}) DataFiles. Concatenating: {paste0(names, collapse=', ')}"))

    # for each DataFile
    for(i in seq_along(x)) {

      x[[i]] <- extract( x[[i]] )

    }

    # ?previous nested merges
    num_merged <- sum(grepl(merge_col, colnames(x[[1]]@data)))
    if(num_merged > 0) {

      merge_col <- paste0(merge_col, "_", as.character(num_merged))

    }

    # update the mapping to have the merge column
    new_map <- add_col(x = x[[1]]@mapping,
                       col_name = merge_col,
                       aliases = list(merge_col),
                       col_type = "character",
                       func = as.character,
                       active = TRUE,
                       overwrite = TRUE)

    # create a merged datafile and overwrite x
    x <- DataFile(path = c(sapply(x, function(y) y@path)),
                  data = data.table::rbindlist(lapply(x, function(x) x@data), idcol=merge_col),
                  mapping = new_map)

    rlog::log_debug(glue::glue("Data size: {nrow(x@data)} rows by {ncol(x@data)} cols"))

    # return
    validObject(x)
    return(x)
  }
)












