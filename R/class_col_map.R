#' ColMap
#'
#' @slot active logical.
#' @slot aliases list.
#' @slot col_types character.
#' @slot funcs list.
#' @slot col_fill ANY.
#'
#' @return a ColMap object
#' @export
#'
ColMap <- setClass(
  Class = "ColMap",
  slots = list(
    active = "logical",
    aliases = "list",
    col_types = "character",
    funcs = "list",
    col_fill = "ANY"
  ),
  prototype = list(
    active = logical(),
    aliases = list(),
    col_types = character(),
    funcs = list(),
    col_fill = NA_character_
  )
)

setValidity(
  Class = "ColMap",
  method = function(object) {

    stopifnot("Column aliases list must be named with the standardised column name" = rlang::is_named(object@aliases))
    stopifnot("Column type vector must be named with the standardised column name" = rlang::is_named(object@col_types))
    stopifnot("Column active vector must be named with the standardised column name" = rlang::is_named(object@active))
    stopifnot("Column function list must be named with the standardised column name" = rlang::is_named(object@funcs))
    stopifnot("Column fill must be either NULL or an atomic type length 1" = length(object@col_fill)==1 & is.atomic(object@col_fill))
    stopifnot("All vector slot names must be the same" = all(setequal(names(object@aliases), names(object@col_types)),
                                                             setequal(names(object@aliases), names(object@funcs)),
                                                             setequal(names(object@aliases), names(object@active))))
    stopifnot("Mapping cannot contain duplicate column names" = !any(duplicated(names(object@aliases))))


    # duplicate aliases check
    aliases <- unlist(object@aliases)
    if(any(duplicated(aliases))) {
      i <- which(duplicated(aliases))
      print("Duplicate alias:")
      print(aliases[i])
      stop("Mapping cannot contain duplicate aliases")
    }
  }
)

#' col_fill
#'
#' @param x a Colmap object
#'
#' @return the fill value
#' @export
#'
setGeneric("col_fill", function(x) standardGeneric("col_fill"))
#' @rdname col_fill
setMethod("col_fill", "ColMap", function(x) x@col_fill )


#' all_cols
#'
#' @param x a Colmap object
#'
#' @return the names of the active columns
#' @export
#'
setGeneric("all_cols", function(x) standardGeneric("all_cols"))
#' @rdname all_cols
setMethod("all_cols", "ColMap", function(x) names(x@active))


#' col_types
#'
#' @param x a Colmap object
#'
#' @return the column types of the active columns
#' @export
#'
setGeneric("col_types", function(x) standardGeneric("col_types"))
#' @rdname col_types
setMethod("col_types", "ColMap", function(x) x@col_types[x@active] )

#' col_type_funcs
#'
#' @param x a Colmap object
#'
#' @return a list of functions to apply to the columns on extracting
#' @export
#'
setGeneric("col_type_funcs", function(x) standardGeneric("col_type_funcs"))
#' @rdname col_type_funcs
setMethod("col_type_funcs", "ColMap", function(x) x@funcs[x@active] )

#' col_aliases
#'
#' @param x a Colmap object
#'
#' @return a list of alias vectors
#' @export
#'
setGeneric("col_aliases", function(x) standardGeneric("col_aliases"))
#' @rdname col_aliases
setMethod("col_aliases", "ColMap", function(x) x@aliases[x@active] )

#' col_active
#'
#' @param x a Colmap object
#'
#' @return logical vectors of possible columns (nameS) and whether active or not
#' @export
#'
setGeneric("col_active", function(x) standardGeneric("col_active"))
#' @rdname col_active
setMethod("col_active", "ColMap", function(x) x@active[x@active] )

#' set_active
#'
#' @param x a Colmap object
#' @param col_names the names to activate
#' @param rest.off whether to turn the rest of the columns off
#'
#' @return a Colmap object
#' @export
#'
setGeneric("set_active", function(x, col_names, rest.off=TRUE) standardGeneric("set_active"))
#' @rdname set_active
setMethod(
  f = "set_active",
  signature = c("list", "character"),
  definition = function(x, col_names, rest.off=TRUE) {
    stopifnot("x must be a list of ColMap objects" = all(sapply(x, inherits, what="ColMap")))
    x <- lapply(x, set_active, col_names=col_names, rest.off=rest.off)
    return(x)
  }
)
#' @rdname set_active
setMethod(
  f = "set_active",
  signature = c("ColMap", "character"),
  definition = function(x, col_names, rest.off=TRUE) {

    if(!all(col_names %in% names(x@aliases))) {
      stop(paste0("All `col_names` must be ColMap columns. Not found: ",
                  paste0(col_names[!col_names %in% names(x@aliases)], collapse=", ")))
    }

    for(i in seq_along(x@active)) {

      if(names(x@active[i]) %in% col_names) {

        x@active[[i]] <- TRUE

      } else if(rest.off) {

        x@active[[i]] <- FALSE

      } else {

        next

      }

    }

    validObject(x)
    x
  }
)

#' col_names
#'
#' @param x a Colmap object
#'
#' @return names of the active columns
#' @export
#'
setGeneric("col_names", function(x) standardGeneric("col_names"))
#' @rdname col_names
setMethod("col_names", "ColMap", function(x) names(col_map(x, only.active=TRUE)) )

#' col_map
#'
#' @param x a Colmap object
#' @param input_col_names NULL,
#' @param only.active .
#' @param ignore.case .
#' @param ... .
#'
#' @return the col name map
#' @export
#'
setGeneric("col_map", function(x, ...) standardGeneric("col_map"))
#' @rdname col_map
setMethod(
  f = "col_map",
  signature = "ColMap",
  definition = function(x, input_col_names=NULL, only.active=TRUE, ignore.case=TRUE) {

    stopifnot("`input_col_names` must be a character vector" = is.character(input_col_names) | is.null(input_col_names))
    stopifnot("`ignore.case` must be logical" = is.logical(ignore.case) & length(ignore.case)==1)

    # cols to get
    if(only.active) {
      std_columns <- names(x@active)[x@active]
    } else {
      std_columns <- names(x@active)
    }

    # generate the output named vector; std_name = old_name OR na
    output_cols <- rep(NA_character_, length(std_columns))
    names(output_cols) <- std_columns

    # for each standard column
    for(i in seq_along(std_columns)) {

      # see if there is a matching input col
      for(j in seq_along(input_col_names)) {

        # does the input column exist in the aliases
        alias_present <- grepl(pattern = paste0("^", input_col_names[[j]], "$"),
                               x = x@aliases[[ std_columns[[i]] ]],
                               ignore.case = ignore.case)

        # add input name to std_name if is an alias
        if(any(alias_present)) {

          output_cols[[ std_columns[[i]] ]] <- input_col_names[[j]]
          break

        }

      } # end cycle input_col_names

    } # end cycle std_columns

    # if requesting only active warn if not found
    if(only.active & any(is.na(output_cols)) & !is.null(input_col_names)) {

      rlog::log_warn(glue::glue("Some active columns where not found in the input_col_names: {paste0(names(output_cols[is.na(output_cols)]), collapse=', ')}"))

    }

    # return
    return(output_cols)
  }
)

#' add_col
#'
#' @param x ColMap object
#' @param col_name .
#' @param col_type .
#' @param func .
#' @param aliases .
#' @param active .
#' @param overwrite .
#'
#' @return ColMap object
#' @export
#'
setGeneric("add_col", function(x, col_name, col_type, func, aliases=list(col_name), active=TRUE, overwrite=FALSE) standardGeneric("add_col"))
#' @rdname add_col
setMethod(
  f = "add_col",
  signature = c("ColMap", "character", "character", "function", "list"),
  definition = function(x, col_name, col_type, func, aliases=list(col_name), active=TRUE, overwrite=FALSE) {

    stopifnot("`col_name` must be a character" = is.character(col_name) & length(col_name)==1)
    stopifnot("`aliases` must be a list of characters" = is.list(aliases) & all(sapply(aliases, is.character)))
    stopifnot("`col_type` must be an atomic type, one of c('numeric', 'character', 'integer', 'logical)" = col_type %in% c('numeric', 'character', 'integer', 'logical'))
    stopifnot("`overwrite` must be a logical" = is.logical(overwrite) & length(overwrite)==1)
    stopifnot("`active` must be a logical" = is.logical(active) & length(active)==1)

    if(!col_name %in% names(x@aliases)) {

      names(aliases) <- col_name
      x@aliases = c(x@aliases, aliases)
      names(col_type) <- col_name
      x@col_types = c(x@col_types, col_type)
      func <- list(func)
      names(func) <- col_name
      x@funcs = c(x@funcs, func)
      names(active) <- col_name
      x@active <- c(x@active, active)

    } else if(col_name %in% names(x@aliases) & overwrite) {

      x@aliases[[col_name]] = aliases
      x@col_types[[col_name]] = col_type
      x@funcs[[col_name]] = func
      x@active[[col_name]] = active

    } else {

      stop("`col_name` already exists, set `overwrite=TRUE` to proceed anyway")

    }

    validObject(x)
    return(x)
  }
)

#' remove_col
#'
#' @param x ColMap object
#' @param col_name .
#'
#' @return a Colmap object
#' @export
#'
setGeneric("remove_col", function(x, col_name) standardGeneric("remove_col"))
#' @rdname remove_col
setMethod(
  f = "remove_col",
  signature = c("ColMap", "character"),
  definition = function(x, col_name) {

    if(!all(col_name %in% names(x@aliases))) {

      warning("All `col_name(s)` do not exist in ColMap")

    } else {

      name_filter <- names(x@aliases)[!names(x@aliases) %in% col_name]
      a <- x@aliases[name_filter]
      t <- x@col_types[name_filter]
      f <- x@funcs[name_filter]
    }

    x <- ColMap(active = x@active[name_filter],
                aliases = x@aliases[name_filter],
                col_types = x@col_types[name_filter],
                funcs = x@funcs[name_filter],
                col_fill = x@col_fill)

    validObject(x)
    return(x)
  }
)

#' add_alias
#'
#' @param x .
#' @param col_name .
#' @param alias .
#'
#' @return a Datafile obj
#' @export
#'
setGeneric("add_alias", function(x, col_name, alias) standardGeneric("add_alias"))
#' @rdname add_alias
setMethod(
  f = "add_alias",
  signature = c("ColMap", "character", "character"),
  definition = function(x, col_name, alias) {

    if(!col_name %in% names(x@aliases)) {

      stop(paste0("Column name `", col_name, "` not found in ColMap names: ", paste0(names(x@aliases), collapse=", ")))

    } else {

      x@aliases[col_name] = c(x@aliases[[col_name]], alias)

    }

    validObject(x)
    return(x)
  }
)


#' remove_alias
#'
#' @param x .
#' @param col_name .
#' @param alias .
#'
#' @return a Datafile obj
#' @export
#'
setGeneric("remove_alias", function(x, col_name, alias) standardGeneric("remove_alias"))
#' @rdname remove_alias
setMethod(
  f = "remove_alias",
  signature = c("ColMap", "character", "character"),
  definition = function(x, col_name, alias) {

    if(!col_name %in% names(x@aliases)) {

      stop(paste0("Column name `", col_name, "` not found in ColMap names: ", paste0(names(x@aliases), collapse=", ")))

    } else {

      new_aliases <- x@aliases[[col_name]][!grepl(paste0("^",alias,"$"), x@aliases[[col_name]])]
      x@aliases[[col_name]] <- new_aliases

    }

    validObject(x)
    return(x)
  }
)
