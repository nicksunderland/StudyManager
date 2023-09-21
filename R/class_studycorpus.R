#' Title
#'
#' @slot corpus_dir character.
#' @slot study_type character.
#' @slot studies list.
#' @slot results .
#' @slot meta_quantitative .
#' @slot meta_indel_alleles .
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
    studies = "list",
    results = "list",
    meta_quantitative = "logical",
    meta_indel_alleles = "logical"
  ),
  prototype = list(
    corpus_dir = character(),
    study_type = character(),
    studies = list(),
    results = list(),
    meta_quantitative = TRUE,
    meta_indel_alleles = TRUE
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

#' Title
#'
#' @param x ..
#' @param i .
#'
#' @return .
#' @export
#'
setMethod("[", "StudyCorpus", function(x, i) x@studies[i])

#' Title
#'
#' @param x ..
#' @param i ..
#'
#' @return .
#' @export
#'
setMethod("[[", "StudyCorpus", function(x, i) x@studies[[i]])

#' studies
#'
#' @param object .
#' @param ... .
#' @param index .
#'
#' @return .
#' @export
#'
setGeneric("studies", function(object, ...) standardGeneric("studies"))
#' @rdname studies
setMethod("studies", "StudyCorpus", function(object, index=NULL) {

  if(is.null(index)) {

    return(object@studies)

  } else {

    tryCatch(
      expr = {
        study <- object@studies[index]
        if(is.null(study) | length(study)==0) stop()
        return(study)

      },
      error = function(e) {

        rlog::log_error(glue::glue("Index error, object@studies[[`{index}`]] is not defined"))

      }
    )
  }
})

#' studies-set
#'
#' @param x .
#' @param value .
#' @param ... . .
#' @param index .
#'
#' @return .
#' @export
#'
setGeneric("studies<-", function(x, value, ...) standardGeneric("studies<-"))
#' @rdname studies-set
setMethod("studies<-", "StudyCorpus", function(x, value, index=NULL) {

  if(is.null(index)) {

    x@studies <- value
    validObject(x)
    return(x)

  } else {

    tryCatch(
      expr = {

        x@studies[[index]] <- value
        validObject(x)
        return(x)

      },
      error = function(e) {

        rlog::log_error(glue::glue("Index error, object@studies[[`{index}`]] is not defined"))

      }
    )
  }
})

#' length
#'
#' @param x ..
#'
#' @return .
#' @export
#'
setMethod("length", "StudyCorpus", function(x) length(x@studies))


#' create_study_list
#'
#' @param object .
#' @param study_type .
#'
#' @return .
#' @export
#'
setGeneric("create_study_list", function(object, study_type) standardGeneric("create_study_list"))
#' @rdname create_study_list
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
                              ref_path = object@ref_path,
                              qc_freq_threshold = object@qc_freq_threshold,
                              qc_low_events_n = object@qc_low_events_n,
                              qc_freq_low_events_thresh = object@qc_freq_low_events_thresh,
                              qc_freq_ambig_thresh = object@qc_freq_ambig_thresh,
                              qc_qual_threshold = object@qc_qual_threshold,
                              qc_call_rate_threshold = object@qc_call_rate_threshold,
                              qc_freq_diff_threshold = object@qc_freq_diff_threshold)

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




#' run_gwama
#'
#' @param corpus .
#' @param output_dir .
#' @param ... .
#'
#' @return .
#' @export
#'
setGeneric("run_gwama", function(corpus, output_dir, ...) standardGeneric("run_gwama"))
#' @rdname run_gwama
setMethod(
  f = "run_gwama",
  signature = c("StudyCorpus", "character"),
  definition = function(corpus, output_dir, ...) {

    stopifnot(dir.exists(output_dir))

    # clean up previous runs
    old_run_files <- list.files(output_dir, pattern="_gwama.in", full.names=TRUE)
    unlink(old_run_files)

    # do for each study
    for(i in seq_along(corpus@studies)) {

      # extract the study to work on
      study <- corpus[[i]]

      # get the keys into the data_files structure
      if(length(c(...)) == 0) {

        data_file_keys <- get_data_keys(study)

      } else {

        data_file_keys <- get_data_keys(...)

      }
      rlog::log_debug(glue::glue("DataFile keys to process for meta-analysis: {basename(study@dir)}, {paste0(lapply(data_file_keys, paste0, collapse=','), collapse=' | ')}"))
      stopifnot(keys_valid(study, data_file_keys))

      # for each key, process, each key will need a separate meta-analysis
      for(keys in data_file_keys) {

        # No data, skip
        if(length(study@qc_data_files[[ keys ]]) == 0) {

          rlog::log_warn(glue::glue("No post-QC DataFile found for {basename(study@dir)}{paste0('[',keys,']',collapse='')}: Pre-QC: {length(study@data_files[[ keys ]])} files, Post-QC: {length(study@qc_data_files[[ keys ]])} files --> skipping..."))
          next

        } else {

          rlog::log_info(glue::glue("Preparing meta-analaysis data for {basename(study@dir)}{paste0('[',keys,']',collapse='')}"))

        }

        # define the gwama.in file open an append connect to it
        gwama.in <- file.path(output_dir, glue::glue("{paste0(keys,collapse='_')}_gwama.in"))
        if(!file.exists(gwama.in)) {

          gwama.in.connection <- file(gwama.in, "w")

        } else {

          gwama.in.connection <- file(gwama.in, "a")

        }

        # set the mapping with the required columns for meta-analysis analysis
        req_cols <- c("cptid","N","EFFECT_ALLELE","OTHER_ALLELE","FRQ","BETA","SE")
        map <- StudyManager::base_column_mapping
        map <- set_active(map, req_cols)

        # the post-QC data; potentially combining multiple files depending on keys
        rlog::log_debug(glue::glue("Extracting post-qc plot data: {basename(study@dir)} qc_data_files{paste0('[',keys,']',collapse='')}"))
        free( study@qc_data_files[[ keys ]] )
        mapping(study@qc_data_files[[ keys ]]) <- map
        qc_data_file <- extract( study@qc_data_files[[ keys ]], merge_col="CHR_FCT")

        # write new (maybe coombined file) out; then free the memory
        tmp_data_path <- tempfile(glue::glue("{basename(study@dir)}_{paste0(keys,collapse='_')}_"))
        write_file(qc_data_file, tmp_data_path)
        free( qc_data_file )

        # add the file to the gwama input
        writeLines(tmp_data_path, gwama.in.connection)

        # close connection
        close(gwama.in.connection)

      } # end for each key
    } # end for each study

    # GWAMA runs
    gwama.in_file_list = list.files(output_dir, pattern="_gwama.in", full.names=TRUE)

    # run a meta-analysis for each of the keys (outcomes usually)
    for(gwama.in_file in gwama.in_file_list) {

      # which key set
      key_str <- sub("_gwama.in", "", basename(gwama.in_file))

      # separate output directory; overwrite previous runs
      output_sub_dir <- file.path(output_dir, key_str)
      if(dir.exists(output_sub_dir)) unlink(output_sub_dir)
      dir.create(output_sub_dir, showWarnings=FALSE)
      gwama.out_file <- file.path(output_sub_dir,key_str)

      # GWAMA executable


      # the GWAMA command
      gwama_cmd <- paste("/Users/xx20081/Downloads/GWAMA_v2/GWAMA",  # need to change this
                         "--filelist", gwama.in_file,
                         "--output", gwama.out_file,
                         ifelse(corpus@meta_quantitative, "--quantitative", ""),
                         ifelse(corpus@meta_indel_alleles, "--indel_alleles", ""),
                         "--name_marker cptid",
                         "--name_n N",
                         "--name_ea EFFECT_ALLELE",
                         "--name_nea OTHER_ALLELE",
                         "--name_eaf FRQ",
                         "--name_beta BETA",
                         "--name_se SE"
      )

      # run the meta-analysis
      system(gwama_cmd)

      # GWAMA mapping
      gwama_map <- StudyManager::base_column_mapping
      gwama_map <- remove_col(gwama_map, c("EFFECT_ALLELE","OTHER_ALLELE"))
      gwama_map <- add_col(gwama_map, "EFFECT_ALLELE", "character", as.character, list("reference_allele"), active=TRUE)
      gwama_map <- add_col(gwama_map, "OTHER_ALLELE", "character", as.character, list("other_allele"), active=TRUE)
      gwama_map <- set_active(gwama_map, c("SNP","EFFECT_ALLELE","OTHER_ALLELE","FRQ","BETA","SE","BETA_95L","BETA_95U","Z","P","LOG10_P","Q_STATISTIC","Q_P_VALUE","HETISQT","N","N_CAS","DIRECTION"))

      # create the results file
      corpus@results[[key_str]] <- DataFile(path = gwama.out_file,
                                            mapping = gwama_map)

    }

    validObject(corpus)
    return(corpus)

  }
)


cluster <- function(on, cluster=NULL, cores=NA_real_) {

  stopifnot(is.logical(on) & (is.numeric(cores) | is.na(cores)))

  if(on) {

    cores_detected <- parallel::detectCores()
    if(!is.na(cores)) {

      if(cores <= cores_detected) {

        cores <- as.integer(cores)

      } else {

        rlog::log_warn(glue::glue("Specified number of cores [{cores}] is more than the detected cores [{cores_detected}]. Using {cores_detected-1}..."))
        cores <- cores_detected - 1

      }

    } else {

      cores <- cores_detected - 1

    }

    cl <- parallel::makeForkCluster(cores)
    doParallel::registerDoParallel(cl)
    return(cl)

  } else {

    parallel::stopCluster(cluster)
    return(NULL)

  }
}

#' @rdname run_qc
setMethod(
  f = "run_qc",
  signature = "StudyCorpus",
  definition = function(object, ..., index=NULL, parallel_cores=NA_integer_) {

    rlog::log_info("Running run_qc(`StudyCorpus`)...")

    stopifnot("`parallel_cores` must be an integer" = (is.numeric(parallel_cores) | is.na(parallel_cores)))

    if(!is.na(parallel_cores)) {

      log_path <- file.path(object@corpus_dir, '__log_run_qc_plots.txt')
      writeLines(c(""), log_path)
      cl <- cluster(on=TRUE)
      `%do_or_dopar%` <- foreach::`%dopar%`

    } else {

      `%do_or_dopar%` <- foreach::`%do%`

    }

    # if index NULL do all of them
    if(is.null(index)) {
      index <- 1:length(object)
    }

    # run the corpus
    study_names <- names(studies(object))
    studies(object) <- foreach::foreach(study    = studies(object),
                                        name     = names(studies(object)),
                                        i        = 1:length(object),
                                        .inorder = TRUE,
                                        .combine = 'c') %do_or_dopar% {

      # if parallel, turn everything on, logging file will be put in corpus main directory
      if(!is.na(parallel_cores)) {

       sink(log_path, append=TRUE)

      }

      if(i %in% index) {

        return( list(run_qc(study, ...)) ) # return as list, then 'c' to combine ensures that loop of n=1 still returns a list

      } else {

        return( list(study) ) # unprocessed

      }


    }
    names(studies(object)) <- study_names

    # if was parallel, turn everything off
    if(!is.na(parallel_cores)) {

      cluster(on=FALSE, cluster=cl)
      sink()

    }

    # return
    validObject(object)
    return(object)
  }
)


#' @rdname run_qc_plots
setMethod(
  f = "run_qc_plots",
  signature = c("StudyCorpus", "character"),
  definition = function(object, output_dir, ..., index=NULL, parallel_cores=NA_integer_) {

    rlog::log_info("Running run_qc_plots(`StudyCorpus`)...")

    stopifnot("`parallel_cores` must be an integer" = (is.numeric(parallel_cores) | is.na(parallel_cores)))

    if(!is.na(parallel_cores)) {

      log_path <- file.path(object@corpus_dir, '__log_run_qc_plots.txt')
      writeLines(c(""), log_path)
      cl <- cluster(on=TRUE)
      `%do_or_dopar%` <- foreach::`%dopar%`

    } else {

      `%do_or_dopar%` <- foreach::`%do%`

    }

    foreach::foreach(study = studies(object, index)) %do_or_dopar% {

      if(!is.na(parallel_cores)) {

        sink(log_path, append=TRUE)

      }

      run_qc_plots(study, output_dir, ...)

    }

    if(!is.na(parallel_cores)) {

      cluster(on=FALSE, cluster=cl)
      sink()

    }

    # return
    validObject(object)
    return(object)
  }
)

