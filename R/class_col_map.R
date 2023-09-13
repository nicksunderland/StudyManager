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
    stopifnot("Mapping cannot contain duplicate aliases" = !any(duplicated(unname(unlist(object@aliases)))))

  }
)

setGeneric("col_fill", function(x) standardGeneric("col_fill"))
setMethod("col_fill", "ColMap", function(x) x@col_fill )

setGeneric("all_cols", function(x) standardGeneric("all_cols"))
setMethod("all_cols", "ColMap", function(x) names(x@active))

setGeneric("col_types", function(x) standardGeneric("col_types"))
setMethod("col_types", "ColMap", function(x) x@col_types[x@active] )

setGeneric("col_type_funcs", function(x) standardGeneric("col_type_funcs"))
setMethod("col_type_funcs", "ColMap", function(x) x@funcs[x@active] )

setGeneric("col_aliases", function(x) standardGeneric("col_aliases"))
setMethod("col_aliases", "ColMap", function(x) x@aliases[x@active] )

setGeneric("col_active", function(x) standardGeneric("col_active"))
setMethod("col_active", "ColMap", function(x) x@active[x@active] )

setGeneric("set_active", function(x, col_names, rest.off=TRUE) standardGeneric("set_active"))
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

setGeneric("col_names", function(x, ...) standardGeneric("col_names"))
setMethod(
  f = "col_names",
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

      warning(paste0("Some active columns where not found in the input_col_names: ",
                     paste0(names(output_cols[is.na(output_cols)]), collapse=", ")))

    }

    # return
    return(output_cols)
  }
)

setGeneric("add_col", function(x, col_name, aliases, col_type, func, active=TRUE, overwrite=FALSE) standardGeneric("add_col"))
setMethod(
  f = "add_col",
  signature = c("ColMap", "character", "list", "character", "function", "logical", "logical"),
  definition = function(x, col_name, aliases, col_type, func, active=TRUE, overwrite=FALSE) {

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

setGeneric("remove_col", function(x, col_name) standardGeneric("remove_col"))
setMethod(
  f = "remove_col",
  signature = c("ColMap", "character"),
  definition = function(x, col_name) {

    if(!col_name %in% names(x@aliases)) {

      warning("`col_name` does not exist in ColMap")

    } else {

      name_filter <- names(x@aliases != col_name)

      x@aliases <- x@aliases[name_filter]
      x@col_type <- x@col_type[name_filter]
      x@funcs <- x@funcs[name_filter]
    }

    validObject(x)
    return(x)
  }
)

setGeneric("add_alias", function(x, col_name, alias) standardGeneric("add_alias"))
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

