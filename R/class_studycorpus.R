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

#' results
#'
#' @param object .
#' @param ... .
#' @param index .
#'
#' @return .
#' @export
#'
setGeneric("results", function(object, ...) standardGeneric("results"))
#' @rdname results
setMethod("results", "StudyCorpus", function(object, index=NULL) {

  if(is.null(index)) {

    return(object@results)

  } else {

    tryCatch(
      expr = {
        result <- object@results[index]
        if(is.null(result) | length(result)==0) stop()
        return(result)

      },
      error = function(e) {

        rlog::log_error(glue::glue("Index error, object@results[`{index}`] is not defined"))

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

#' create_results_list
#'
#' @param object .
#' @param study_type .
#' @param results_dir description
#'
#' @return .
#' @export
#'
setGeneric("create_results_list", function(corpus, study_type, results_dir) standardGeneric("create_results_list"))
#' @rdname create_results_list
setMethod(
  f = "create_results_list",
  signature = c("StudyCorpus", "character", "character"),
  definition = function(corpus, study_type, results_dir) {

    corpus@results <- list()

    results_directories <- list.dirs(results_dir, recursive=FALSE)

    for(dir in results_directories) {

      if(study_type == "GWASsumstats") {

        result_file <- list.files(dir, "\\.out", full.names = TRUE)
        result_file <- result_file[!grepl("(\\.err\\.out|\\.gc\\.out|\\.log\\.out)$", result_file)]

        if(length(result_file)==1) {
          key <- basename(sub("\\.out$", "", result_file))
          corpus@results[[key]] <- DataFile(path = result_file)
        } else {
          next
        }
      }
    }

    validObject(corpus)
    return(corpus)
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
  definition = function(corpus, output_dir, ..., index=NULL, parallel_cores=NA_integer_) {

    stopifnot(dir.exists(output_dir))
    stopifnot("`parallel_cores` must be an integer" = (is.numeric(parallel_cores) | is.na(parallel_cores)))

    # clean up previous runs
    old_run_files <- list.files(output_dir, pattern="_gwama.in", full.names=TRUE)
    unlink(old_run_files)

    # parallel?
    if(!is.na(parallel_cores)) {

      log_path <- file.path(corpus@corpus_dir, '__log_run_qc_plots.txt')
      writeLines(c(""), log_path)
      cl <- cluster(on=TRUE)
      `%do_or_dopar%` <- foreach::`%dopar%`

    } else {

      `%do_or_dopar%` <- foreach::`%do%`

    }

    # if index NULL do all of them
    if(is.null(index)) {
      index <- 1:length(corpus)
    }

    # define the studies
    studies <- corpus[index]

    # get all the possible keys
    possible_keys <- c()
    for(s in studies) {

      # get the keys into the data_files structure
      if(length(c(...)) == 0) {

        data_file_keys <- get_data_keys(s)

      } else {

        data_file_keys <- get_data_keys(...)

      }
      rlog::log_debug(glue::glue("DataFile keys to process for meta-analysis: {basename(s@dir)}, {paste0(lapply(data_file_keys, paste0, collapse=','), collapse=' | ')}"))
      stopifnot(keys_valid(s, data_file_keys))

      possible_keys <- unique(c(possible_keys, data_file_keys))
    }

    # run for each key / outcome
    for(key in possible_keys) {

      # for each study
      for(s in studies) {

        # check valid key for this study
        if(!keys_valid(s, key) | length(s@qc_data_files[[key]]) == 0) {

          rlog::log_warn(glue::glue("No post-QC DataFile found for {basename(s@dir)}{paste0('[',key,']',collapse='')} --> skipping..."))
          next

        }

        # get the data
        rlog::log_info(glue::glue("Preparing meta-analaysis data for {basename(s@dir)}{paste0('[',key,']',collapse='')}"))

        # define the gwama.in file open an append connect to it
        gwama.in <- file.path(output_dir, glue::glue("{paste0(key,collapse='_')}_gwama.in"))
        if(!file.exists(gwama.in)) {

          gwama.in.connection <- file(gwama.in, "w")

        } else {

          gwama.in.connection <- file(gwama.in, "a")

        }

        # set the mapping with the required columns for meta-analysis analysis
        req_cols <- c("cptid","EFFECT_ALLELE","OTHER_ALLELE","BETA","SE","N","FRQ","STRAND")
        map <- StudyManager::base_column_mapping
        map <- set_active(map, req_cols, rest.off=TRUE)

        # the post-QC data; potentially combining multiple files depending on keys
        rlog::log_debug(glue::glue("Extracting post-qc plot data: {basename(s@dir)} qc_data_files{paste0('[',key,']',collapse='')}"))
        free( s@qc_data_files[[ key ]] )
        mapping(s@qc_data_files[[ key ]]) <- map
        qc_data_file <- extract( s@qc_data_files[[ key ]])


        # write new (maybe combined file) out; then free the memory
        tmp_data_path <- tempfile(glue::glue("{basename(s@dir)}_{paste0(key,collapse='_')}_"))
        write_file(qc_data_file, tmp_data_path, na=".", quote=FALSE) # GWAMA doesn't like empty fields, so set empty to "."?; also doesnt like quoted string headers
        free( qc_data_file )

        # add the file to the gwama input
        writeLines(tmp_data_path, gwama.in.connection)

        # close connection
        close(gwama.in.connection)

      } # end for each study
    } # end for each key

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

      # the GWAMA command
      gwama_exe <- system.file("bin", "GWAMA_v2", "GWAMA", package="StudyManager")
      gwama_cmd <- paste(gwama_exe,
                         "--filelist", gwama.in_file,
                         "--output", gwama.out_file,
                         ifelse(corpus@meta_quantitative, "--quantitative", ""),
                         ifelse(corpus@meta_indel_alleles, "--indel_alleles", ""),
                         "--name_marker cptid",
                         "--name_ea EFFECT_ALLELE",
                         "--name_nea OTHER_ALLELE",
                         "--name_beta BETA",
                         "--name_se SE",
                         "--name_eaf FRQ",
                         "--name_n N",
                         "--name_strand STRAND",
                         "" #--random
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
      corpus@results[[key_str]] <- DataFile(path = paste0(gwama.out_file, ".out"),
                                            mapping = gwama_map)

      # delete the tmp files
      gwama.in_files <- readLines(gwama.in_file)
      unlink(gwama.in_files)

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

      # generate the plots as side effects
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







#' run_meta_plots
#'
#' @param object .
#' @param output_dir .
#' @param ... .
#' @param index .
#' @param parallel_cores .
#'
#' @return True
#' @export
#'
setGeneric("run_meta_plots", function(object, output_dir, index=NULL, parallel_cores=NA_integer_) standardGeneric("run_meta_plots"))
#' @rdname run_meta_plots
setMethod(
  f = "run_meta_plots",
  signature = c("StudyCorpus", "character"),
  definition = function(object, output_dir, index=NULL, parallel_cores=NA_integer_) {

    rlog::log_info("Running run_meta_plots(`StudyCorpus`)...")

    stopifnot("`parallel_cores` must be an integer" = (is.numeric(parallel_cores) | is.na(parallel_cores)))

    if(!is.na(parallel_cores)) {

      log_path <- file.path(object@corpus_dir, '__log_run_qc_plots.txt')
      writeLines(c(""), log_path)
      cl <- cluster(on=TRUE)
      `%do_or_dopar%` <- foreach::`%dopar%`

    } else {

      `%do_or_dopar%` <- foreach::`%do%`

    }

    outcome_names = names(results(object, index))
    foreach::foreach(result = results(object, index), name = outcome_names) %do_or_dopar% {

      if(!is.na(parallel_cores)) {

        sink(log_path, append=TRUE)

      }

      # do the processing of the result
      stopifnot(dir.exists(output_dir))

      # set the mapping with the required columns for analysis
      gwama_map <- StudyManager::base_column_mapping
      gwama_map <- remove_col(gwama_map, c("EFFECT_ALLELE","OTHER_ALLELE"))
      gwama_map <- add_col(gwama_map, "EFFECT_ALLELE", "character", as.character, list("reference_allele"), active=TRUE)
      gwama_map <- add_col(gwama_map, "OTHER_ALLELE", "character", as.character, list("other_allele"), active=TRUE)
      gwama_map <- set_active(gwama_map, c("SNP","P", "EFFECT_ALLELE", "OTHER_ALLELE", "BETA"))

      # set the mapping into the results DataFile and extract
      rlog::log_debug(glue::glue("Extracting meta-analysis results data for DataFile: {result@path}"))
      free(result)
      mapping(result) <- gwama_map
      result <- extract(result)

      # get the data.table
      dt <- get_data(result)
      dt[, CHR :=  sub(":.*", "", SNP)]
      dt[, BP :=  as.numeric(sub(".+:([0-9]+)(?::.*)?", "\\1", SNP))]


      # PLOTTING
      rlog::log_debug(glue::glue("Plotting - output dir: {output_dir}"))
      output_filename <- sub(".out", "", basename(result@path))
      output_path <- file.path(output_dir, output_filename)


      # GET TOP HITS
      # wont need the is.finite once I sort the QC procedure.....
      top_hits <- dt[P<5e-8 & is.finite(BETA), ]
      data.table::fwrite(top_hits, paste0(output_path, "_meta_tophits.tsv"), sep="\t")


      # CREATE Manhattan
      create_manhattan(dt = dt[, c("SNP","CHR","BP","P")],
                       file_path = paste0(output_path, "_meta_manhattan.png"),
                       highlight = top_hits[,SNP],
                       annotate = top_hits[,SNP],
                       title = paste0("Manhattan Plot - ", output_filename))

      # CREATE FOREST
      for(hit_idx in seq_along(nrow(top_hits))) {

        study_contrib <- data.table::data.table(cptid=character(),
                                                CHR=character(),
                                                BP=integer(),
                                                P=numeric(),
                                                BETA=numeric(),
                                                SE=numeric(),
                                                EFFECT_ALLELE=character(),
                                                OTHER_ALLELE=character(),
                                                study=character())
        hit <- top_hits[hit_idx, ]

        key = c(sub("(.*?_.*?)_.*", "\\1", name), sub("(?:.*?_.*?)_(.*)", "\\1", name))

        for(study in studies(corpus)) {

          empty_row <- list(study = basename(study@dir),
                            cptid = NA_character_,
                            CHR = hit$CHR,
                            BP = hit$BP,
                            BETA = 0,
                            SE = 0,
                            P = NA_real_,
                            EFFECT_ALLELE = hit$EFFECT_ALLELE,
                            OTHER_ALLELE = hit$OTHER_ALLELE)

          if(keys_valid(study, key)) {
            dat <- get_data(extract( study@qc_data_files[[ key ]] ))
            dat <- dat |> dplyr::select(cptid, CHR, BP, P, BETA, SE, EFFECT_ALLELE, OTHER_ALLELE)
            dat <- dat[CHR==hit$CHR & BP==hit$BP, ]
            dat[ , study := .(basename(study@dir)) ]

            if(nrow(dat)==0) {
              dat <- data.table(t(empty_row))
            }
            study_contrib <- rbind(study_contrib, dat)
          } else {
            study_contrib <- rbind(study_contrib, empty_row, fill=T)
          }

        }

        study_contrib[ , c("upper", "lower") := list(as.numeric(BETA)-as.numeric(SE), as.numeric(BETA)+as.numeric(SE))]

        create_forestplot(study_contrib,
                          file_path = paste0(output_path, "_", gsub(":","_",top_hits[hit_idx, SNP]), "_meta_forestplot.png"),
                          title = paste0("Forest Plot - ", top_hits[hit_idx, SNP], " - ", output_filename))

      }

      # CREATE ZOOM PLOTS
      for(hit_idx in seq_along(nrow(top_hits))) {
        file_path <- paste0(output_path, "_", gsub(":","_",top_hits[hit_idx, SNP]), "_meta_locuszoom.png")
        create_locuszoom(dt,
                         snp = top_hits[hit_idx, SNP],
                         file_path = file_path,
                         title = paste0("LocusZoom - ", top_hits[hit_idx, SNP], " - ", output_filename))
      }

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


#' run_filter_summary_plots
#'
#' @param object .
#' @param output_dir .

#'
#' @return True
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr summarise across starts_with rename_with rename group_by
#' @export
#'
setGeneric("run_filter_summary_plots", function(object, output_dir, parallel_cores=NA_integer_) standardGeneric("run_filter_summary_plots"))
#' @rdname run_filter_summary_plots
setMethod(
  f = "run_filter_summary_plots",
  signature = c("StudyCorpus", "character"),
  definition = function(object, output_dir, parallel_cores=NA_integer_) {

    rlog::log_info("Running run_filter_summary_plots(`StudyCorpus`)...")

    stopifnot("`parallel_cores` must be an integer" = (is.numeric(parallel_cores) | is.na(parallel_cores)))

    if(!is.na(parallel_cores)) {

      log_path <- file.path(object@corpus_dir, '__log_run_filter_summary_plots.txt')
      writeLines(c(""), log_path)
      cl <- cluster(on=TRUE)
      `%do_or_dopar%` <- foreach::`%dopar%`

    } else {

      `%do_or_dopar%` <- foreach::`%do%`

    }

    # loop over studies
    summary_dt <- foreach::foreach(study = studies(object), .combine="rbind") %do_or_dopar% {
      if(!is.na(parallel_cores)) {

        sink(log_path, append=TRUE)

      }
      dt <- qc_filter_summary(study)
    }


    if(!is.na(parallel_cores)) {

      cluster(on=FALSE, cluster=cl)
      sink()

    }

    # PLOTTING
    rlog::log_debug(glue::glue("Plotting - output dir: {output_dir}"))

    # heatmap - filter overview per study
    summary_dt |>
      dplyr::group_by(study) |>
      dplyr::summarise(n = sum(n, na.rm=T),
                       dplyr::across(dplyr::starts_with("FILTER_"), ~ sum(.x))) |>
      dplyr::mutate(dplyr::across(dplyr::starts_with("FILTER_"), ~ .x / n)) |>
      dplyr::select(-n) |>
      dplyr::rename_with(~ sub("_n$","",.), dplyr::starts_with("FILTER_")) |>
      tidyr::pivot_longer(starts_with("FILTER_"), names_to="FILTER", values_to="FAIL_RATE") |>
      data.table::as.data.table() |>
      create_heatmap(x="study", y="FILTER", z="FAIL_RATE",
                     x_lab="Study", y_lab="Filter", fill_lab="Fail rate %", fill_lims=c(0,0.3),
                     file_path=file.path(output_dir, "filter_summary_studies.png"),
                     title = "Filter summary - overview")

    # heatmap - filter by outcome per study
    summary_dt |>
      dplyr::rename(outcome=key_1) |>
      dplyr::group_by(study, outcome) |>
      dplyr::summarise(n = sum(n, na.rm=T),
                       dplyr::across(dplyr::starts_with("FILTER_"), ~ sum(.x))) |>
      dplyr::mutate(dplyr::across(dplyr::starts_with("FILTER_"), ~ .x / n)) |>
      dplyr::select(-n) |>
      dplyr::rename_with(~ sub("_n$","",.), dplyr::starts_with("FILTER_")) |>
      tidyr::pivot_longer(starts_with("FILTER_"), names_to="FILTER", values_to="FAIL_RATE") |>
      data.table::as.data.table() |>
      create_heatmap(x="study", y="FILTER", z="FAIL_RATE", fy="outcome",
                     x_lab="Study", y_lab="Filter", fill_lab="Fail rate %",fill_lims=c(0,0.3),
                     file_path=file.path(output_dir, "filter_summary_studies_outcome.png"),
                     title = "Filter summary - outcome")

    # heatmap - filter by chromosome type per study
    summary_dt |>
      dplyr::rename(chrom_type=key_2) |>
      dplyr::group_by(study, chrom_type) |>
      dplyr::summarise(n = sum(n, na.rm=T),
                       dplyr::across(dplyr::starts_with("FILTER_"), ~ sum(.x))) |>
      dplyr::mutate(dplyr::across(dplyr::starts_with("FILTER_"), ~ .x / n)) |>
      dplyr::select(-n) |>
      dplyr::rename_with(~ sub("_n$","",.), dplyr::starts_with("FILTER_")) |>
      tidyr::pivot_longer(starts_with("FILTER_"), names_to="FILTER", values_to="FAIL_RATE") |>
      data.table::as.data.table() |>
      create_heatmap(x="study", y="FILTER", z="FAIL_RATE", fy="chrom_type",
                     x_lab="Study", y_lab="Filter", fill_lab="Fail rate %",fill_lims=c(0,0.3),
                     file_path=file.path(output_dir, "filter_summary_studies_chrom_type.png"),
                     title = "Filter summary - autosomes vs. X-chr")

    # heatmap - filter by chromosome number per study
    summary_dt |>
      dplyr::mutate(CHR = as.character(as.integer(CHR))) |>
      dplyr::mutate(CHR = dplyr::coalesce(CHR, key_2)) |>
      dplyr::mutate(CHR = factor(CHR, levels=c(as.character(1:25), unique(summary_dt$key_2)))) |>
      dplyr::group_by(study, CHR) |>
      dplyr::summarise(n = sum(n, na.rm=T),
                       dplyr::across(dplyr::starts_with("FILTER_"), ~ sum(.x))) |>
      dplyr::mutate(dplyr::across(dplyr::starts_with("FILTER_"), ~ .x / n)) |>
      dplyr::select(-n) |>
      dplyr::rename_with(~ sub("_n$","",.), dplyr::starts_with("FILTER_")) |>
      tidyr::pivot_longer(starts_with("FILTER_"), names_to="FILTER", values_to="FAIL_RATE") |>
      data.table::as.data.table() |>
      create_heatmap(x="CHR", y="FILTER", z="FAIL_RATE", fy="study",
                     x_lab="Chromosome", y_lab="Filter", fill_lab="Fail rate %",fill_lims=c(0,0.3),
                     file_path=file.path(output_dir, "filter_summary_studies_chromosome.png"),
                     h = 4800,
                     title = "Filter summary - chromosome")

    # return
    validObject(object)
    return(object)
  }
)
