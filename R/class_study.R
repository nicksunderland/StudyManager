#' Study virtual class
#'
#' @slot dir character.
#' @slot file_structure list.
#' @slot data_files list.
#' @slot pre_qc_dir .
#' @slot post_qc_dir .
#' @slot mapping .
#' @slot qc_data_files .
#'
#' @return object
#' @importFrom methods .valueClassTest slotNames
#' @export
#'
Study <- setClass(
  Class = "Study",
  contains = "VIRTUAL",
  slots = list(
    dir = "character",
    pre_qc_dir = "character",
    post_qc_dir = "character",
    file_structure = "list",
    mapping = "ColMap",
    data_files = "list",
    qc_data_files = "list"
  ),
  prototype = list(
    dir = character(),
    pre_qc_dir = "pre_qc",
    post_qc_dir = "post_qc",
    file_structure = list(),
    mapping = NULL,
    data_files = list(),
    qc_data_files = list()
  )
)

setMethod(
  f = "initialize",
  signature = "Study",
  definition = function(.Object, mapping=StudyManager::base_column_mapping, ...) {

    call <- paste0(names(as.list(match.call()[-1])),collapse=', ')
    rlog::log_trace(glue::glue("Init Study: call({call})"))

    # cant seem to have default data as a slot, so set here
    .Object@mapping <- mapping

    # send the remaining to the default init method to set slots
    .Object <- callNextMethod(.Object, ...)

    # form the rest of the slots, if we are a Study or parent of a study (e.g. GWASsumstats)
    # this is a hack to test if we are in fact a StudyCorpus, which inherits from
    # the parents GWASsumstats etc. To test we just see if we have the slot 'corpus_dir'
    if(!"corpus_dir" %in% slotNames(.Object)) {
      .Object@data_files <- create_data_files(.Object)
      .Object@qc_data_files <- create_qc_data_files(.Object)
    }

    validObject(.Object)
    rlog::log_trace(glue::glue("Init Study: complete"))
    return(.Object)
  }
)

#' pre_qc_data_dir
#'
#' @param x .
#'
#' @return .
#' @export
#'
setGeneric("pre_qc_data_dir", function(x) standardGeneric("pre_qc_data_dir"))
#' @rdname pre_qc_data_dir
setMethod("pre_qc_data_dir", "Study", function(x) file.path(x@dir, x@pre_qc_dir))

#' post_qc_data_dir
#'
#' @param x .
#'
#' @return .
#' @export
#'
setGeneric("post_qc_data_dir", function(x) standardGeneric("post_qc_data_dir"))
#' @rdname post_qc_data_dir
setMethod("post_qc_data_dir", "Study", function(x) file.path(x@dir, x@post_qc_dir))

setClassUnion("listORcharacterORmissing", c("list", "character", "missing"))
#' create_data_files
#'
#' @param object .
#' @param x .
#' @param append_regex .
#'
#' @return .
#' @export
#'
setGeneric("create_data_files", function(object, x=NULL, append_regex="") standardGeneric("create_data_files"))
#' @rdname create_data_files
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

        rlog::log_warn(glue::glue("Dir: {basename(object@dir)}"))
        rlog::log_warn(glue::glue("file regex `{file_name_regex}` did not find any files"))
        return(NULL)

      } else if(length(file_path) > 1) {

        rlog::log_warn(glue::glue("Dir: {basename(object@dir)}"))
        rlog::log_error(glue::glue("file regex must yeild one file (or no file) but found multiple using `{file_name_regex}`: \n{paste0(file_path, collapse='\n')}"))
        stop("create_data_files() multi file regex error")

      } else {

        data_file <- DataFile(path=file_path, mapping=object@mapping)
        return(data_file)

      }
    }
  }
)

# recursion via different signatures; always returns a list of keys (character vectors)
#' get_data_keys
#'
#' @param input .
#' @param ... .
#' @param prefix description
#'
#' @return .
#' @export
#'
setGeneric("get_data_keys", valueClass="list", function(input, ...) standardGeneric("get_data_keys"))
# if 'Study' call the function on the study data_files object (a list)\
#' @rdname get_data_keys
setMethod("get_data_keys", "Study", function(input, prefix = NULL) get_data_keys(input@data_files))
# if list, i.e. the data_file object, recursively return into that list
#' @rdname get_data_keys
setMethod("get_data_keys", "list", function(input, prefix = NULL) {

    result <- list()

    for (name in names(input)) {

      new_prefix <- c(prefix, name) # add the name of the list, i.e. build up a the names of each level
      result <- c(result, get_data_keys(input[[name]], new_prefix))

    }
    return(result)
  })
# the base case, essentially 'if(DataFile)', return the prefix (i.e. the name of the DataFile)
#' @rdname get_data_keys
setMethod("get_data_keys", "DataFile",  function(input, prefix = NULL) list(prefix))
# a different case where actual list level names are passed as characters or character vectors
#' @rdname get_data_keys
setMethod("get_data_keys", "character", function(input, ...) {

    # must be character
    stopifnot("All inputs must be characters or character vectors: level1, level2, ..., levelN" = is.character(c(input, ...)))

    # in order: input is level 1, then subsequent characters level2, level3, etc
    levels_list <- c(list(input), list(...))

    # get all combinations of the upper levels with the lower levels
    keys_df <- expand.grid(levels_list)

    # convert this to a list of keys (character vectors)
    keys_list <- unname(as.list(as.data.frame(t(keys_df))))

    # return
    return(keys_list)
  })

#' keys_valid
#'
#' @param object .
#' @param keys .
#' @param ... .
#'
#' @return .
#' @export
#'
setGeneric("keys_valid", function(object, keys, ...) standardGeneric("keys_valid"))
#' @rdname keys_valid
setMethod("keys_valid", c("Study", "list"), function(object, keys){
    stopifnot("An empty list is not a valid key" = length(keys) > 0)
    test_all = all(sapply(keys, keys_valid, object=object))
    return(test_all)
  })
#' @rdname keys_valid
setMethod("keys_valid", c("Study", "character"), function(object, keys, ...){

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
        levels_list <- c(list(keys), list(...))
        keys_df <- expand.grid(levels_list)
        keys_list <- unname(as.list(as.data.frame(t(keys_df))))

        # return
        return(keys_list)


        rlog::log_error(glue::glue("invalid keys: c({paste0('\"', c(keys, ...), '\"', collapse=', ')}). Possible keys: {paste(get_data_keys(object), collapse=', ')}"))
        return(FALSE)
      }
    )
  })

#' copy_file_structure
#'
#' @param file_structure .
#'
#' @return .
#' @export
#'
setGeneric("copy_file_structure", function(file_structure) standardGeneric("copy_file_structure"))
#' @rdname copy_file_structure
setMethod("copy_file_structure", "list", function(file_structure) {

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


#' create_qc_data_files
#'
#' @param object .
#' @param ... .
#' @param input_list .
#' @param parent_names .
#'
#' @return .
#' @export
#'
setGeneric("create_qc_data_files", function(object, ...) standardGeneric("create_qc_data_files"))
#' @rdname create_qc_data_files
setMethod("create_qc_data_files", "Study", function(object, input_list=NULL, parent_names=character(0)) {

  if(is.null(input_list)) {

    create_qc_data_files(object, object@file_structure)

  } else {

    for (key in names(input_list)) {

      value <- input_list[[key]]
      current_names <- c(parent_names, key)

      if (is.list(value)) {

        input_list[[key]] <- create_qc_data_files(object, value, current_names)

      } else {

        file_name_regex <- paste0(c(basename(object@dir), current_names, "post_qc.gz$"), collapse="_")
        file_path <- list.files(post_qc_data_dir(object),  full.names=TRUE)
        file_path = file_path[grepl(file_name_regex, file_path, perl=TRUE)]

        # only allow one file per regex (each file should be named), or no file found
        if(length(file_path) == 0) {

          rlog::log_warn(glue::glue("Dir: {basename(object@dir)}"))
          rlog::log_warn(glue::glue("file regex `{file_name_regex}` did not find any post_qc files"))
          input_list[[key]] <- NULL

        } else if(length(file_path) > 1) {

          rlog::log_warn(glue::glue("Dir: {basename(object@dir)}"))
          rlog::log_error(glue::glue("More than one post-qc file found multiple with `{file_name_regex}`: \n{paste0(file_path, collapse='\n')}"))
          stop("create_qc_data_files() multi file regex error")

        } else {
          input_list[[key]] <- DataFile(path=file_path, mapping=object@mapping)
        }
      }
    }
    return(input_list)
  }
}
)


# setValidity(
#   Class = "Study",
#   method = function(object) {
#
#     stopifnot("`dir` must be a valid directory path" = dir.exists(object@dir))
#
#   }
# )
