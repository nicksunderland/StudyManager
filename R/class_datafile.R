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
    mapping = StudyManager::ColumnMapping
  )
)

setValidity(
  Class = "DataFile",
  method = function(object) {

    path_missing = all(is.na(object@path))
    data_missing = ncol(object@data) == 0

    if(!path_missing) {
      stopifnot("`path` must be valid file path" = all(file.exists(object@path)))
    }

    if(!path_missing & !data_missing) {
        stopifnot("`path` must be a valid and writable file path if `data` present" =
                    all(sapply(object@path, checkmate::test_path_for_output, overwrite=TRUE)))
    }

    stopifnot("`data` must be a data.table" = data.table::is.data.table(object@data))
    stopifnot("Mapping must be a valid ColMap object" = inherits(object@mapping, "ColMap"))
  }
)

setGeneric("get_raw_col_names", function(x) standardGeneric("get_raw_col_names"))
setMethod(
  f = "get_raw_col_names",
  signature = "DataFile",
  definition = function(x){

    if(!all(is.na(x@path))) {

      header <- c()
      for(p in x@path) {

        h <- data.table::fread(p, nrows=0)
        header <- unique(c(header, colnames(h)))

      }
      return( header )

    # data was directly inserted to the data slot, without pathname
    } else if(ncol(x@data) > 0) {

        return( colnames(x@data) )

    # object validity (data and/or path) should stop us getting here
    } else {

      rlog::log_error("ERROR: get_raw_col_names() ... shouldn't be here")
      stop()

    }
  }
)

setGeneric("path", function(x) standardGeneric("path"))
setMethod("path", "DataFile", function(x) x@path)

setGeneric("path<-", function(x, value) standardGeneric("path<-"))
setMethod("path<-", "DataFile", function(x, value) {
  x@path <- value
  validObject(x)
  x
})

setGeneric("file_name", function(x) standardGeneric("file_name"))
setMethod("file_name", "DataFile", function(x) basename(x@path))

setGeneric("mapping", function(x) standardGeneric("mapping"))
setMethod("mapping", "DataFile", function(x) x@mapping)

setGeneric("mapping<-", function(x, value) standardGeneric("mapping<-"))
setMethod("mapping<-", "DataFile", function(x, value) {
  x@mapping <- value
  warning("DataFile mapping changed, existing data cleared, you should re-extract")
  x <- free(x)
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
      warning("Data empty - do you need to call `extract()` first?")
      return(FALSE)
    }

  }
)

setGeneric("get_data", function(x) standardGeneric("get_data"))
setMethod(
  f = "get_data",
  signature = "DataFile",
  definition = function(x){

    raw_cols <- get_raw_col_names(x)
    active_col_map <- col_names(x@mapping, input_col_names=raw_cols)

    if(data_exists(x)) {

      # rename columns using the mapping vector
      data.table::setnames(x@data, unname(active_col_map), names(active_col_map), skip_absent=TRUE)

      # Insert col_fill columns for any missing entries in mapping vector
      missing_cols <- setdiff(names(active_col_map), names(x@data))
      for(col in missing_cols) {

          x@data[, (col) := col_fill(x@mapping)]

      }

      # drop the others
      unwanted_cols <- setdiff(names(x@data), names(active_col_map))
      if(length(unwanted_cols) > 0) {

        x@data[, (unwanted_cols) := NULL]

      }

      return(x@data)

    } else {

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


setGeneric("write_file", function(x, file_path, append=FALSE, sep="\t", ...) standardGeneric("write_file"))
setMethod(
  f = "write_file",
  signature = c("DataFile", "character"),
  definition = function(x, file_path, append=FALSE, sep="\t", ...) {

    if(data_exists(x)) {

      data.table::fwrite(x@data, file=file_path, append=append, sep=sep, ...)

    }
  }
)

setGeneric("free", function(x) standardGeneric("free"))
setMethod("free", "DataFile", function(x) {
  if(data_exists(x)) {
    x@data[, names(x@data) := NULL]
  }
  validObject(x)
  x
})

setGeneric("extract", function(data_file_list, ...) standardGeneric("extract"))

setMethod(
  f = "extract",
  signature = "list",
  definition = function(data_file_list, merge_col=NULL) {

    stopifnot(length(data_file_list) > 0)

    names = names(data_file_list)
    if(is.null(names)) {
        names = as.character(seq_along(data_file_list))
    }

    rlog::log_info(glue::glue("Extract from multiple ({length(data_file_list)}) DataFiles. Concatenating: {paste0(names, collapse=', ')}"))

    # for each DataFile
    for(i in seq_along(data_file_list)) {

      data_file_list[i] <- extract(data_file_list[i])

    }

    dt_list <- lapply(data_file_list, function(y) get_data(extract(y)))
    names(dt_list) <- names

    if(!is.null(merge_col)) {

      mapping(data_file_list[[1]]) <- add_col(x = x[[i]]@mapping,
                                              col_name = merge_col,
                                              aliases = list(merge_col),
                                              col_type = "character",
                                              func = as.character,
                                              active = TRUE,
                                              overwrite = TRUE)

      set_data(data_file_list[[1]]) <- rbind(unlist(dt_list), idcol=merge_col)

    } else {

      set_data(data_file_list[[1]]) <- rbind(unlist(dt_list))

    }
    rlog::log_debug(glue::glue("Data size: {nrow(data_file_list[[1]]@data)} rows by {ncol(data_file_list[[i]]@data)} cols"))

    data_file_list[[1]]@path <- sapply(data_file_list, function(y) y@path)

    # return
    validObject(data_file_list[[1]])
    return(data_file_list[[1]])
  }
)


setMethod(
  f = "extract",
  signature = "DataFile",
  definition = function(x) {

    rlog::log_info(glue::glue("Extract from DataFile: {basename(x@path)}"))

    # active columns from the mapping
    raw_cols <- get_raw_col_names(x)
    active_col_map <- col_names(x@mapping, input_col_names=raw_cols)
    rlog::log_debug(glue::glue("Mapping columns: {paste0(names(active_col_map),'=',active_col_map,collapse=',')}"))

    # non-missing cols in the datafile
    missing_cols <- is.na(active_col_map)
    non_missing_cols <- active_col_map[!missing_cols]

    if(any(missing_cols)) {

      rlog::log_debug(glue::glue("Missing mapping columns: {paste0(names(active_col_map[missing_cols]),'=',active_col_map[missing_cols], collapse=',')}. These will be added when calling `get_data()`}"))
    }

    # read the data (only store non-missing, as when we request data we'll fill it in later)
    x@data <- data.table::fread(x@path, select=unname(non_missing_cols))
    rlog::log_debug(glue::glue("Data file {basename(x@path)} size: {nrow(x@data)} rows by {ncol(x@data)} cols"))
    rlog::log_debug(glue::glue("Data file {basename(x@path)} cols: {paste0(colnames(x@data), collapse=',')}"))

    # renames
    data.table::setnames(x@data, active_col_map, names(active_col_map), skip_absent=TRUE)
    rlog::log_debug(glue::glue("Renamed extracted columns to: {paste0(names(x@data),collapse=',')}"))

    # get the column type function
    col_type_check_funcs <- col_type_funcs(mapping(x))
    col_type_str <- col_types(mapping(x))
    rlog::log_debug(glue::glue("Ensuring column types: {paste0(names(col_type_str),'=',col_type_str,collapse=',')}"))

    # enforce function on each column
    for(j in seq_along(col_type_check_funcs)) {

      col_name <- names(col_type_check_funcs)[j]
      x@data[, (col_name) := col_type_check_funcs[[j]](get(col_name))]

    }

    validObject(x)
    return(x)
  }
)









