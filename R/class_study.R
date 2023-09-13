#' Title
#'
#' @slot dir character.
#' @slot pmid integer.
#' @slot file_structure list.
#' @slot data_files list.
#'
#' @return object
#' @export
#'
Study <- setClass(
  Class = "Study",
  contains = "VIRTUAL",
  slots = list(
    dir = "character",
    pre_qc_dir = "character",
    post_qc_dir = "character",
    pmid = "integer",
    file_structure = "list",
    mapping = "ColMap",
    data_files = "list",
    qc_data_files = "list"
  ),
  prototype = list(
    dir = character(),
    pre_qc_dir = "pre_qc",
    post_qc_dir = "post_qc",
    pmid = integer(),
    file_structure = list(),
    mapping = StudyManager::ColumnMapping,
    data_files = list(),
    qc_data_files = list()
  )
)

setMethod(
  f = "initialize",
  signature = "Study",
  definition = function(.Object, ...) {
    .Object <- callNextMethod(.Object, ...)
    .Object@data_files <- create_data_files(.Object)
    .Object@qc_data_files <- copy_file_structure(.Object@data_files)
    return(.Object)
  }
)

setGeneric("pre_qc_data_dir", function(x) standardGeneric("pre_qc_data_dir"))
setMethod("pre_qc_data_dir", "Study", function(x) file.path(x@dir, x@pre_qc_dir))

setGeneric("post_qc_data_dir", function(x) standardGeneric("post_qc_data_dir"))
setMethod("post_qc_data_dir", "Study", function(x) file.path(x@dir, x@post_qc_dir))


setClassUnion("listORcharacterORmissing", c("list", "character", "missing"))
setGeneric("create_data_files", function(object, x=NULL) standardGeneric("create_data_files"))
setMethod(
  f = "create_data_files",
  signature = c("Study", "listORcharacterORmissing"),
  definition = function(object, x=NULL) {

    if(is.null(x)) {
      file_struct_element <- object@file_structure
    } else {
      file_struct_element <- x
    }

    if (is.list(file_struct_element)) {

      data_file_list <- list()

      # recurse down the lists
      for (el_name in names(file_struct_element)) {

        data_file_list[[el_name]] <- create_data_files(object, file_struct_element[[el_name]])

      }

      # the end result
      return(data_file_list)

    # the base case - reached a non-list element
    } else {

      # get the full file path by regex
      file_name_regex = file_struct_element
      file_paths = list.files(pre_qc_data_dir(object),  full.names=TRUE)

      file_path = file_paths[grepl(file_name_regex, file_paths, perl=TRUE)]

      # only allow one file per regex (each file should be named), or no file found
      if(length(file_path) == 0) {

        warning(paste0("file regex `", file_name_regex, "` did not find any files"))
        return(NULL)

      } else if(length(file_path) > 1) {

        stop(paste0("file regex must yeild one file (or no file) but found multiple using `", file_name_regex, "`: \n",
                    paste0(file_path, collapse="\n")))

      } else {

        data_file <- DataFile(path=file_path, mapping=object@mapping)
        return(data_file)

      }
    }
  }
)

# recursion via different signatures
setGeneric("get_data_keys", function(input, prefix=NULL) standardGeneric("get_data_keys"))
# if 'Study' call the function on the study data_file object
setMethod("get_data_keys", "Study", function(input, prefix = NULL) get_data_keys(input@data_files))
# if list, i.e. the data_file object, recursively return into that list
setMethod(
  f = "get_data_keys",
  signature = "list",
  definition = function(input, prefix = NULL) {

    result <- list()

    for (name in names(input)) {

      new_prefix <- c(prefix, name)
      result <- c(result, get_data_keys(input[[name]], new_prefix))

    }
    return(result)
  }
)
# the base case, essentially 'if(DataFile)', return the prefix (i.e. the name of the DataFile)
setMethod("get_data_keys", "DataFile",  function(input, prefix = NULL) list(prefix))

setGeneric("keys_valid", function(object, keys, ...) standardGeneric("keys_valid"))
setMethod(
  f = "keys_valid",
  signature = c("Study", "list"),
  definition = function(object, keys){
    test_all = all(sapply(keys, keys_valid, object=object))
    return(test_all)
  }
)
setMethod(
  f = "keys_valid",
  signature = c("Study", "character"),
  definition = function(object, keys, ...){

    tryCatch(
      expr = {
        # ... more keys
        test <- object@data_files[[ c(keys, ...) ]]
        if(is.null(test)) {
          stop()
        } else {
          return(TRUE)
        }
      },
      error = function(e) {
        rlog::log_error(glue::glue("invalid keys: c({paste0('\"', c(keys, ...), '\"', collapse=', ')}). Possible keys: {paste(get_data_keys(object), collapse=', ')}"))
        return(FALSE)
      }
    )
  }
)

setGeneric("copy_file_structure", function(file_structure) standardGeneric("copy_file_structure"))
setMethod(
  f = "copy_file_structure",
  signature = c("list"),
  definition = function(file_structure) {

    result <- vector("list", length(file_structure))
    names(result) <- names(file_structure)

    for (i in seq_along(file_structure)) {
      if(is.list(file_structure[[i]])) {
        result[[i]] <- copy_file_structure(file_structure[[i]])
      }
    }
    return(result)
  }
)

#
# setGeneric("get_data", function(object, ...) standardGeneric("get_data"))
# setMethod(
#   f = "get_data",
#   signature = c("Study"),
#   definition = function(object, ...) {
#
#     if(!keys_valid(object@data_files, ...)) {
#       stop("All arguments in '...' must be valid ordered keys into the `data_files` structure")
#     }
#
#     object@data_files[[ c(...) ]] <- extract( object@data_files[[ c(...) ]] )
#
#     return( file_data(object@data_files[[ c(...) ]] ))
#
#   }
# )


setGeneric("get_data_file", function(object, ...) standardGeneric("get_data_file"))
setMethod(
  f = "get_data_file",
  signature = c("Study"),
  definition = function(object, ...) {

    if(!keys_valid(object@data_files, ...)) {
      stop("All arguments in '...' must be valid ordered keys into the `data_files` structure")
    }

    return( object@data_files[[ c(...) ]] )

  }
)

setGeneric("get_qc_data_file", function(object, ...) standardGeneric("get_qc_data_file"))
setMethod(
  f = "get_qc_data_file",
  signature = c("Study"),
  definition = function(object, ...) {

    if(!keys_valid(object@qc_data_files, ...)) {
      stop("All arguments in '...' must be valid ordered keys into the `qc_data_files` structure")
    }

    return( object@qc_data_files[[ c(...) ]] )

  }
)





