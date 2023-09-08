#' Title
#'
#' @slot dir character.
#' @slot pmid integer.
#' @slot file_structure list.
#' @slot col_names character.
#' @slot col_types list.
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
    pmid = "integer",
    file_structure = "list",
    col_names = "character",
    col_types = "list",
    data_files = "list"
  ),
  prototype = list(
    dir = character(),
    pmid = integer(),
    file_structure = list(),
    col_names = character(),
    col_types = list(),
    data_files = list()
  )
)


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
      file_paths = list.files(object@dir, full.names=TRUE)
      file_path = file_paths[grepl(file_name_regex, file_paths, perl=TRUE)]

      # only allow one file per regex (each file should be named), or no file found
      if(length(file_path) == 0) {

        warning(paste0("file regex `", file_name_regex, "` did not find any files"))
        return(NULL)

      } else if(length(file_path) > 1) {

        stop(paste0("file regex must yeild one file (or no file) but found multiple using `", file_name_regex, "`: \n",
                    paste0(file_path, collapse="\n")))

      } else {

        data_file <- DataFile(path=file_path, col_names=object@col_names, col_types=object@col_types)
        return(data_file)

      }
    }
  }
)


setGeneric("get_data", function(object, ...) standardGeneric("get_data"))
setMethod(
  f = "get_data",
  signature = c("Study"),
  definition = function(object, ...) {

    if (!all(sapply(list(...), is.character))) {
      stop("All arguments in '...' must be keys of class 'character'")
    }

    object@data_files[[ c(...) ]] <- extract( object@data_files[[ c(...) ]] )

    return( file_data(object@data_files[[ c(...) ]] ))

  }
)





