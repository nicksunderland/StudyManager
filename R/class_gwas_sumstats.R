#' Title
#'
#' @slot parameter numeric.
#'
#' @return object
#' @export
#' @importFrom methods callNextMethod
#' @include class_study.R
#'
GWASsumstats <- setClass(
  Class = "GWASsumstats",
  contains = "Study",
  slots = list(
    qc_freq_threshold = "numeric",
    qc_low_events_n = "numeric",
    qc_freq_low_events_thresh = "numeric",
    qc_freq_ambig_thresh = "numeric",
    qc_qual_threshold = "numeric",
    qc_call_rate_threshold = "numeric",
    qc_freq_diff_threshold = "numeric",
    ref_path = "character",
    ref_data_file = "DataFile"
  ),
  prototype = list(
    qc_freq_threshold = 0.01,
    qc_low_events_n = 200,
    qc_freq_low_events_thresh = 0.05,
    qc_freq_ambig_thresh = 0.40,
    qc_qual_threshold = 0.5,
    qc_call_rate_threshold = 0.98,
    qc_freq_diff_threshold = 0.2,
    ref_path = character(),
    ref_data_file = NULL
  )
)

setMethod(
  f = "initialize",
  signature = "GWASsumstats",
  definition = function(.Object, ref_path, ref_data_file=NULL, ...) {

    print("init GWASsumstats")
    print(names(as.list(match.call()[-1])))

    # peel off GWASsumstats specific slots
    .Object@ref_path <- ref_path
    if(is.null(.Object@ref_data_file)) {
      .Object@ref_data_file <- DataFile(path=.Object@ref_path)
    }

    # send the rest down to Study
    .Object <- callNextMethod(.Object, ...)

    # .Object picked up a mapping from 'Study' initilialiser, so use to set ref map
    ref_cols <- c("cptid","SNP","CHR","BP","MARKER_TYPE","A0","OTHER_ALLELE","EUR_FRQ")
    mapping(.Object@ref_data_file) <- set_active(.Object@mapping, ref_cols)

    # required gwas cols, now set the mapping to these
    gwas_cols <- c("SNP","SE","P","N_CAS","N","INFO","FRQ","CHR","BP","BETA","EFFECT_ALLELE","OTHER_ALLELE","STRAND")
    .Object@mapping <- set_active(.Object@mapping, gwas_cols)  # see mapping slot of 'Study'

    validObject(.Object)
    return(.Object)
  }
)


setGeneric("run_qc", function(object, ...) standardGeneric("run_qc"))
setMethod(
  f = "run_qc",
  signature = "StudyCorpus",
  definition = function(object, ...) {

    `%dopar%` <- foreach::`%dopar%`
    cores <- parallel::detectCores()
    cl <- parallel::makeForkCluster(cores)
    doParallel::registerDoParallel(cl)

    foreach::foreach(study = object@studies) %dopar% {

      run_qc(study, ...)

    }
  }
)
setMethod(
  f = "run_qc",
  signature = "GWASsumstats",
  definition = function(object, ...) {

    # extract the ref data, enforcing types and doing checks etc.
    rlog::log_info("Reading reference data into GWASsumstats object")
    object@ref_data_file <- extract( object@ref_data_file )

    # write out
    tmp_ref_path <- tempfile()
    rlog::log_info(glue::glue("Writing parsed reference data to tmp file: {basename(tmp_ref_path)}"))
    write_file(object@ref_data_file, file_path=tmp_ref_path)

    # free the memory / old data
    object@ref_data_file <- free( object@ref_data_file )

    # update the the ref_data_file with this new tmp file path as the source
    path(object@ref_data_file) <- tmp_ref_path

    # create an output directory if it doesn't exist; or clean up (delete) if it does
    output_dir <- file.path(object@dir, object@post_qc_dir)
    if(!dir.exists(output_dir)) {

      rlog::log_info(glue::glue("Output directory not found, creating one at: {output_dir}"))
      dir.create(output_dir, showWarnings=TRUE, recursive=TRUE)

    } else {

      output_file_regex <- "post_qc(.[0-9]+)?(.rep|.ecf|.out|.notinref.txt|.gz|)"
      old_run_files <- list.files(output_dir, pattern=output_file_regex, full.names=TRUE)

      if(length(old_run_files) > 0) {

        rlog::log_info(glue::glue("Output directory found containing {length(old_run_files)} files matching pattern '{output_file_regex}' --> DELETING"))
        unlink(old_run_files)

      }
    }

    # get the keys into the data_files structure
    if(length(c(...)) == 0) {

      data_file_keys <- get_data_keys(object)

    } else {

      data_file_keys <- list(c(...))

    }
    stopifnot(keys_valid(object, data_file_keys))

    # process each data_file
    for(keys in data_file_keys) {

      args = c(list(object), as.list(keys))

      object <- do.call("easy_qc", args)

    }

    # return
    validObject(object)
    return(object)
  }
)



# TODO: clean up function
# delete temp files and reset
# clean up reference file, reset reference
# unlink(tmp_ref_path)
# path(object@ref_data_file) <- ref_path



# don't export this method, just for internal class use
#' Title
#'
#' @param object GWASsumstats oobject
#' @param ... data_file keys
#'
#' @return GWASsumstats oobject
#' @import rlog
#' @export
#'
setGeneric("easy_qc", function(object, ...) standardGeneric("easy_qc"))

#' @rdname easy_qc
setMethod(
  f = "easy_qc",
  signature = c("GWASsumstats"),
  definition = function(object, ...) {

    # ensure valid keys into the data_file structure
    keys <- c(...)
    rlog::log_info(glue::glue("Running EASYQC for DataFile{paste0('[',keys,']',collapse='')}"))
    stopifnot(keys_valid(object, keys))

    # the DataFile data
    rlog::log_info("Extracting file data")
    object@data_files[[ keys ]] <- extract( object@data_files[[ keys ]] )

    # get the column mapping details for EasyQC
    map <- mapping(object@data_files[[ keys ]]) # from each DataFile map as could be different
    input_types <- col_types(map)
    input_cols <- names(input_types)

    # create a written temp file from the DataFile (appropriate col names enforced when we extracted)
    input_path <- tempfile()
    rlog::log_debug(glue::glue("Writting standardised data file to: {basename(input_path)}"))
    fwrite(x = get_data( object@data_files[[ keys ]] ),
           file = input_path,
           sep = "\t")

    # free the data
    rlog::log_debug("Freeing memory / DataFile data")
    object@data_files[[ keys ]] <- free( object@data_files[[ keys ]] )

    # create output path and check
    output_filename <- paste0(c(basename(object@dir), keys, "post_qc"), collapse="_")
    output_path <- file.path(object@dir, object@post_qc_dir, output_filename)
    rlog::log_debug(glue::glue("EASYQC output filename set to: ../{output_filename}"))
    stopifnot("`output_path` should be a writable file path" = checkmate::check_path_for_output(output_path, overwrite=TRUE))

    # get the written temp file from the reference DataFile
    tmp_ref_path <- path(object@ref_data_file)
    ref_map <- mapping(object@ref_data_file)
    ref_input_types <- col_types(ref_map)
    ref_input_cols <- names(ref_input_types)

    # the ecf config contents
    ecf_str <- glue::glue(
      "
      DEFINE
          --pathOut { dirname(output_path) }
          --strSeparator TAB
          --acolIn { paste0(input_cols, collapse=';') }
          --acolInClasses { paste0(input_types, collapse=';') }
          --acolNewName { paste0(input_cols, collapse=';') }

      EASYIN
          --fileIn { input_path }
          --fileInShortName { basename(output_path) }

      START EASYQC

      FILTER
          --rcdFilter INFO >= { object@qc_qual_threshold }
          --strFilterName good_qual_score

      # FILTER
          --rcdFilter ((FRQ >= { object@qc_freq_threshold }) & (FRQ <= { 1-object@qc_freq_threshold }) & (N_CAS > { object@qc_low_events_n } | is.na(N_CAS))) | ((FRQ >= { object@qc_freq_low_events_thresh }) & (FRQ <= { 1-object@qc_freq_low_events_thresh }) & (N_CAS <= { object@qc_low_events_n } | is.na(N_CAS)))
          --strFilterName good_freq

      ADDCOL
          --rcdAddCol N / max(N)
          --colOut CALL_RATE
          --blnOverwrite 1

      FILTER
          --rcdFilter CALL_RATE >= { object@qc_call_rate_threshold }
          --strFilterName good_call_rate

      ADDCOL
          --rcdAddCol ((EFFECT_ALLELE == 'A') & (OTHER_ALLELE == 'T')) | ((EFFECT_ALLELE == 'T') & (OTHER_ALLELE == 'A')) | ((EFFECT_ALLELE == 'C') & (OTHER_ALLELE == 'G')) | ((EFFECT_ALLELE == 'G') & (OTHER_ALLELE == 'C'))
          --colOut AMBIGUOUS
          --blnOverwrite 1

      FILTER
          --rcdFilter (! AMBIGUOUS) | (AMBIGUOUS & ((FRQ <= { object@qc_freq_ambig_thresh }) | (FRQ >= { 1-object@qc_freq_ambig_thresh })))
          --strFilterName no_common_ambiguous

      HARMONIZEALLELES
          --colInA1 EFFECT_ALLELE
          --colInA2 OTHER_ALLELE

      CREATECPTID
          --colInMarker SNP
          --colInA1 EFFECT_ALLELE
          --colInA2 OTHER_ALLELE
          --colInChr CHR
          --colInPos BP
          --blnUseInMarker 0

      ADDCOL
          --rcdAddCol '+'
          --colOut STRAND
          --blnOverwrite 0

      MERGE
          --colInMarker cptid
          --fileRef { tmp_ref_path }
          --acolIn { paste0(ref_input_cols, collapse=';')  }
          --acolInClasses { paste0(ref_input_types, collapse=';') }
          --colRefMarker cptid
          --strRefSuffix _REF
          --blnInAll 0
          --blnRefAll 0
          --blnWriteNotInRef 1

      ADJUSTALLELES
          --colRefA1 OTHER_ALLELE_REF
          --colRefA2 A0_REF
          --colInStrand STRAND
          --colInA1 EFFECT_ALLELE
          --colInA2 OTHER_ALLELE
          --colInFreq FRQ
          --colInBeta BETA
          --blnMetalUseStrand 1
          --blnWriteMismatch 1
          --blnRemoveMismatch 1
          --blnWriteInvalid 1
          --blnRemoveInvalid 1

      FILTER
          --rcdFilter abs(FRQ-EUR_FRQ_REF) <= { object@qc_freq_diff_threshold }
          --strFilterName eaf_difference_within_tol

      WRITE
          --strMode gz
          --strSep TAB

      STOP EASYQC
      "
    )

    rlog::log_info("Starting EASYQC")
    rlog::log_info(glue::glue("ECF file:\n{ecf_str}"))

    # write the .ecf file
    ecf_path <- paste0(output_path, ".ecf")
    writeLines(ecf_str, ecf_path)

    # run EasyQC
    easy_qc_success <- EasyX::EasyX(ecf_path)

    # clean up
    unlink(input_path)

    # get the QC data n.b. no file is written if empty data.table so check for existance
    if(!easy_qc_success) {

      rlog::log_error("EasyQC failed")

    } else {

      # save DataFile
      new_post_qc_map <- data.table::copy(object@mapping)
      new_post_qc_map <- set_active(new_post_qc_map, unique(c(col_names(new_post_qc_map),
                                                              "CALL_RATE","AMBIGUOUS","STRAND","cptid")))

      easyqc_output_file <- paste0(output_path, ".gz")

      if(file.exists(easyqc_output_file)) {

        rlog::log_debug(glue::glue("EasyQC created file: {basename(easyqc_output_file)}"))
        rlog::log_debug(glue::glue("Creating qc_data_file with columns: {paste0(col_names(new_post_qc_map), collapse=',')}"))

        object@qc_data_files[[ keys ]] <- DataFile(path=paste0(output_path, ".gz"),
                                                   mapping=new_post_qc_map)

      } else {

        rlog::log_debug(glue::glue("EasyQC failed to created file, likely no SNPs found"))
        rlog::log_debug(glue::glue("Creating qc_data_file with columns: {paste0(col_names(new_post_qc_map), collapse=',')}"))
        rlog::log_debug(glue::glue("Writing the data.table to: {basename(easyqc_output_file)}"))

        empty_dt <- data.table::data.table(matrix(NA, 0, length(col_names(new_post_qc_map))))
        data.table::setnames(empty_dt, col_names(new_post_qc_map))
        object@qc_data_files[[ keys ]] <- DataFile(path=easyqc_output_file,
                                                   data=empty_dt,
                                                   mapping=new_post_qc_map)
        write_file(object@qc_data_files[[ keys ]], file_path=easyqc_output_file, overwrite=TRUE)
      }
    }

    validObject(object)
    rlog::log_info("Finished EASYQC")
    return(object)
  }
)


setGeneric("run_qc_plots", function(object, output_dir, ...) standardGeneric("run_qc_plots"))
setMethod(
  f = "run_qc_plots",
  signature = c("StudyCorpus", "character"),
  definition = function(object, output_dir, ...) {

    `%dopar%` <- foreach::`%dopar%`
    cores <- parallel::detectCores()
    cl <- parallel::makeForkCluster(cores)
    doParallel::registerDoParallel(cl)

    foreach::foreach(study = object@studies, output_dir=rep(output_dir, length(object@studies))) %dopar% {

      run_qc_plots(study, output_dir)

    }
  }
)
setMethod(
  f = "run_qc_plots",
  signature = c("GWASsumstats", "character"),
  definition = function(object, output_dir, ...) {

    print("herre333")
    print(object@qc_freq_threshold)

    stopifnot(dir.exists(output_dir))

    # get the keys into the data_files structure
    if(length(c(...)) == 0) {

      data_file_keys <- get_data_keys(object)

    } else {

      data_file_keys <- list(c(...))

    }
    stopifnot(keys_valid(object, data_file_keys))

    # set the reference mapping with the required columns for analysis
    ref_map <- mapping(object@ref_data_file)
    ref_map <- set_active(ref_map, c("cptid","EUR_FRQ"))
    object@ref_data_file <- free( object@ref_data_file )
    mapping(object@ref_data_file) <- ref_map

    # # extract the required ref data
    # object@ref_data_file <- extract( object@ref_data_file )
    #
    # # set key for faster processing
    # data.table::setkey(get_data(object@ref_data_file), "cptid")
    #
    # # for each data source
    # for(keys in data_file_keys) {
    #
    #   # set the mapping with the required columns for analysis
    #   req_cols <- c("cptid","SNP","CHR","BP","P","BETA","INFO","FRQ")
    #   maps <- mapping(object@data_files[[ keys ]])
    #   maps <- lapply(maps, set_active, col_names=req_cols)
    #
    #   # the pre-QC data, chromosomes merged
    #   rlog::log_debug(glue::glue("Extracting pre-qc plot data data_files{paste0('[',keys,']',collapse='')}"))
    #   free( object@data_files[[ keys ]] )
    #   mapping(object@data_files[[ keys ]]) <- maps
    #   data_file <- extract( object@data_files[[ keys ]], merge_col="CHR_FCT")
    #
    #   # the pre-QC data, chromosomes merged
    #   rlog::log_debug(glue::glue("Extracting post-qc plot data qc_data_files{paste0('[',keys,']',collapse='')}"))
    #   free( object@qc_data_files[[ keys ]] )
    #   mapping(object@qc_data_files[[ keys ]]) <- maps
    #   qc_data_file <- extract( object@qc_data_files[[ keys ]], merge_col="CHR_FCT" )
    #
    #   # which SNPs were filtered out
    #   filt_snps <- get_data(data_file)[!get_data(qc_data_file), on="SNP"][["SNP"]]
    #
    #   # bind the pre and post QC data
    #   rlog::log_info("Creating plotting data...")
    #   dt <- rbind(
    #       get_data(data_file)   [, "QC_STATUS" := "PRE-QC" ],
    #       get_data(qc_data_file)[, "QC_STATUS" := "POST-QC"]
    #   )[, QC_STATUS := factor(QC_STATUS, levels=c("PRE-QC", "POST-QC"))]
    #   data.table::setkey(dt, "cptid")
    #
    #   # Join the reference data and take the European frequency column
    #   dt[get_data(object@ref_data_file), EUR_FRQ := EUR_FRQ]
    #
    #   diff_levels <- c("No reference data", "EAF outlier", "EAF within tolerance")
    #   dt[, FRQ_DIFF_FCT := data.table::fcase(is.na(EUR_FRQ), factor("No reference data", levels=diff_levels),
    #                                          abs(FRQ-EUR_FRQ) > object@qc_freq_diff_threshold, factor("EAF outlier", levels=diff_levels),
    #                                          abs(FRQ-EUR_FRQ) <=object@qc_freq_diff_threshold, factor("EAF within tolerance", levels=diff_levels),
    #                                          TRUE, factor(NA_character_)) ]
    #
    #   # PLOTTING
    #   rlog::log_debug(glue::glue("Plotting - output dir: {output_dir}"))
    #
    #   output_filename <- paste0(c(basename(object@dir), keys), collapse="_")
    #   output_path <- file.path(output_dir, output_filename)
    #
    #   # Manhattan
    #   create_manhattan(dt = dt[, c("SNP","CHR","BP","P")],
    #                    file_path = paste0(output_path, "_qc_manhattan.png"),
    #                    highlight = filt_snps,
    #                    title = paste0(c("Manhattan Plot -", basename(object@dir), keys), collapse=" "))
    #
    #   # QQ-plot
    #   create_qq(dt = dt[, c("P","QC_STATUS")],
    #             file_path = paste0(output_path, "_qc_qq_plot.png"),
    #             title = paste0(c("QQ plot -", basename(object@dir), keys), collapse=" "))
    #
    #   # P value histogram
    #   create_pvalhist(dt = dt[, c("P","QC_STATUS")],
    #                   file_path = paste0(output_path, "_qc_pval_hist.png"),
    #                   title = paste0(c("P value histogram -", basename(object@dir), keys), collapse=" "))
    #
    #   # Effect size histogram
    #   create_eshist(dt = dt[, c("BETA","QC_STATUS")],
    #                 file_path = paste0(output_path, "_qc_es_hist.png"),
    #                 title = paste0(c("Effect size histogram -", basename(object@dir), keys), collapse=" "))
    #
    #   # Imputation quality histogram
    #   create_impqualhist(dt = dt[, c("INFO","QC_STATUS")],
    #                      file_path = paste0(output_path, "_qc_impqual_hist.png"),
    #                      title = paste0(c("Imputation quality histogram -", basename(object@dir), keys), collapse=" "))
    #
    #   # allele frequency plots
    #   create_eafplot(dt = dt[, c("FRQ","FRQ_DIFF_FCT","EUR_FRQ","QC_STATUS","CHR_FCT")],
    #                  threshold = object@qc_freq_diff_threshold,
    #                  file_path = paste0(output_path, "_qc_eaf_plot.png"),
    #                  title = paste0(c("Allele frequency -", basename(object@dir), keys), collapse=" "))
    #
    #   return(TRUE)
    #
    #  }

  }
)

setGeneric("create_impqualhist", function(dt, file_path=getwd(), ...) standardGeneric("create_impqualhist"))
setMethod(
  f = "create_impqualhist",
  signature = c("data.table"),
  definition = function(dt, file_path=getwd(), ...) {

    # logging info
    rlog::log_info("Effect size histogram...")
    rlog::log_info(glue::glue("Writing image to: {file_path}"))

    # choose a pallete
    palette <- c(wesanderson::wes_palette("GrandBudapest1"), wesanderson::wes_palette("GrandBudapest2"))
    #c("#F1BB7B", "#FD6467", "#5B1A18", "#D67236", "#E6A0C4", "#C6CDF7", "#D8A499", "#7294D4")

    # plot
    plot <- ggplot2::ggplot(dt, ggplot2::aes(x = INFO)) +
      ggplot2::geom_histogram(color = palette[7],
                              fill = palette[5],
                              alpha = 0.9,
                              binwidth = 0.01) +
      ggplot2::labs(...,
           x = "INFO",
           y = "Count") +
      ggplot2::theme_minimal() +
      ggplot2::facet_grid(cols = ggplot2::vars(QC_STATUS))

    # save plot
    ggplot2::ggsave(plot=plot, filename=file_path, device="png", bg="white", height=15, width=30, units="cm")

    # return
    invisible(plot)
  }
)


setGeneric("create_eafplot", function(dt, file_path=getwd(), threshold=0.2, ...) standardGeneric("create_eafplot"))
setMethod(
  f = "create_eafplot",
  signature = c("data.table"),
  definition = function(dt, file_path=getwd(), threshold=0.2, ...) {

    # logging info
    rlog::log_info("Allele frequency plot...")
    rlog::log_info(glue::glue("Writing image to: {file_path}"))

    # choose a pallete
    palette <- c(wesanderson::wes_palette("GrandBudapest1"), wesanderson::wes_palette("GrandBudapest2"))
    #c("#F1BB7B", "#FD6467", "#5B1A18", "#D67236", "#E6A0C4", "#C6CDF7", "#D8A499", "#7294D4")

    colour_labels <- setNames(c(palette[3],palette[2],palette[8]),
                              c("No reference data", "EAF outlier", "EAF within tolerance"))

    # plot
    plot <- ggplot2::ggplot(dt, ggplot2::aes(x = EUR_FRQ, y = FRQ, color = FRQ_DIFF_FCT)) +
      ggplot2::geom_point(size = 0.1) +
      ggplot2::geom_abline(slope=1, intercept= threshold, color=palette[3], linetype = "dotted") +
      ggplot2::geom_abline(slope=1, intercept=-threshold, color=palette[3], linetype = "dotted") +
      ggplot2::lims(x=c(-0.08,1), y=c(0,1)) +
      ggplot2::labs(...,
           x = "Ref EAF",
           y = "Study EAF",
           color = NULL) +
      ggplot2::scale_color_manual(values=colour_labels) +
      ggplot2::theme_minimal() +
      ggplot2::theme(aspect.ratio = 1) +
      ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 5))) +
      ggplot2::facet_grid(rows = ggplot2::vars(CHR_FCT), cols = ggplot2::vars(QC_STATUS))

    # save plot
    ggplot2::ggsave(plot=plot, filename=file_path, device="png", bg="white", height=30, width=25, units="cm")

    # return
    invisible(plot)
})


setGeneric("create_eshist", function(dt, file_path=getwd(), ...) standardGeneric("create_eshist"))
setMethod(
  f = "create_eshist",
  signature = c("data.table"),
  definition = function(dt, file_path=getwd(), ...) {

    # logging info
    rlog::log_info("Effect size histogram...")
    rlog::log_info(glue::glue("Writing image to: {file_path}"))

    # choose a pallete
    palette <- c(wesanderson::wes_palette("GrandBudapest1"), wesanderson::wes_palette("GrandBudapest2"))
    #c("#F1BB7B", "#FD6467", "#5B1A18", "#D67236", "#E6A0C4", "#C6CDF7", "#D8A499", "#7294D4")

    # plot
    plot <- histogram_w_outlier_bins(data = dt,
                             col_name = "BETA",
                             col_group = "QC_STATUS",
                             bin_cutoffs = seq(-2.5, 2.5, 0.05),
                             plot_frequency = TRUE,
                             x_axis_title = "BETA",
                             y_axis_title = "Frequency",
                             x_tick_marks = seq(-2.5, 2.5, 0.5),
                             x_tick_mark_labels = seq(-2.5, 2.5, 0.5),
                             outlier_bin_fill_color = palette[2],
                             non_outlier_bin_fill_color = palette[6],
                             border_color = "grey") +
      ggplot2::labs(...) +
      ggplot2::theme_minimal() +
      ggplot2::facet_grid(cols = ggplot2::vars(QC_STATUS))

    # save plot
    ggplot2::ggsave(plot=plot, filename=file_path, device="png", bg="white", height=15, width=30, units="cm")

    # return
    invisible(plot)
})


setGeneric("create_pvalhist", function(dt, file_path=getwd(), ...) standardGeneric("create_pvalhist"))
setMethod(
  f = "create_pvalhist",
  signature = c("data.table"),
  definition = function(dt, file_path=getwd(), ...) {

    # logging info
    rlog::log_info("P value histogram plot...")
    rlog::log_info(glue::glue("Writing image to: {file_path}"))

    # choose a pallete
    palette <- c(wesanderson::wes_palette("GrandBudapest1"), wesanderson::wes_palette("GrandBudapest2"))
    #c("#F1BB7B", "#FD6467", "#5B1A18", "#D67236", "#E6A0C4", "#C6CDF7", "#D8A499", "#7294D4")

    plot <- ggplot2::ggplot(dt, ggplot2::aes(x = P)) +
      ggplot2::geom_histogram(binwidth = 0.01,
                              color=palette[6],
                              fill=palette[8]) +
      ggplot2::labs(...,
           x = "P value",
           y = "Count") +
      ggplot2::theme_minimal() +
      ggplot2::facet_grid(cols = ggplot2::vars(QC_STATUS))

    # save plot
    ggplot2::ggsave(plot=plot, filename=file_path, device="png", bg="white", height=15, width=30, units="cm")

    # return
    invisible(plot)
  }
)



setGeneric("create_qq", function(dt, file_path=getwd(), ...) standardGeneric("create_qq"))
setMethod(
  f = "create_qq",
  signature = c("data.table"),
  definition = function(dt, file_path=getwd(), ...) {

    # ensure there is some data
    stopifnot(ncol(dt) > 0)
    stopifnot(nrow(dt) > 0)

    # choose a pallete
    palette <- c(wesanderson::wes_palette("GrandBudapest1"), wesanderson::wes_palette("GrandBudapest2"))
    #c("#F1BB7B", "#FD6467", "#5B1A18", "#D67236", "#E6A0C4", "#C6CDF7", "#D8A499", "#7294D4")

    # create dir if needed
    dir.create(dirname(file_path), recursive=TRUE, showWarnings=FALSE)

    # logging info
    rlog::log_info("QQ plot...")
    rlog::log_info(glue::glue("Writing image to: {file_path}"))

    # generate PP-plot but with ggplot instead of qqman
    plot <- dt |>
      dplyr::filter(!is.na(P),
             !is.nan(P),
             !is.null(P),
             is.finite(P),
             P < 1.0,
             P > 0.0) |>
      dplyr::group_by(QC_STATUS) |>
      dplyr::mutate(observed = -log10(sort(P, decreasing=FALSE)),
                    expected = -log10(stats::ppoints(dplyr::n()))) |>
      dplyr::ungroup() |>
      ggplot2::ggplot(ggplot2::aes(x = expected, y = observed)) +
      ggplot2::geom_point(size = 0.5, color=palette[8]) +
      ggplot2::geom_abline(slope=1, intercept=0, color="darkred", linetype = "dotted") +
      ggplot2::labs(...,
           x = "Expected (-log10 P)",
           y = "Observed (-log10 P)",
           color = NULL) +
      ggplot2::theme_minimal() +
      ggplot2::theme(aspect.ratio = 1) +
      ggplot2::facet_grid(cols = ggplot2::vars(QC_STATUS))

    # save plot
    ggplot2::ggsave(plot=plot, filename=file_path, device="png", bg="white", height=15, width=30, units="cm")

    # return
    invisible(plot)
  }
)

setGeneric("create_manhattan", function(dt, file_path=getwd(), highlight=NULL, ...) standardGeneric("create_manhattan"))
setMethod(
  f = "create_manhattan",
  signature = c("data.table"),
  definition = function(dt, file_path=getwd(), highlight=NULL, ...) {

    # ensure there is some data
    stopifnot(ncol(dt) > 0)
    stopifnot(nrow(dt) > 0)

    # choose a pallete
    palette <- c(wesanderson::wes_palette("GrandBudapest1"), wesanderson::wes_palette("GrandBudapest2"))
    #c("#F1BB7B", "#FD6467", "#5B1A18", "#D67236", "#E6A0C4", "#C6CDF7", "#D8A499", "#7294D4")

    # create dir if needed
    dir.create(dirname(file_path), recursive=TRUE, showWarnings=FALSE)

    # logging info
    rlog::log_info("Manhattan plot...")
    rlog::log_info(glue::glue("Writing image to: {file_path}"))

    dt <- dt |> dplyr::mutate(CHR = factor(CHR, levels=c(as.character(1:23),"X")))

    # Prepare the dataset
    d <- dt |>

      # Compute chromosome size
      dplyr::group_by(CHR) |>
      dplyr::summarise(chr_len=as.numeric(max(BP))) |>

      # Calculate cumulative position of each chromosome
      dplyr::mutate(tot=cumsum(chr_len)-chr_len) |>
      dplyr::select(-chr_len) |>

      # Add this info to the initial dataset
      dplyr::left_join(dt, ., by=c("CHR"="CHR")) |>

      # Add a cumulative position of each SNP
      dplyr::arrange(CHR, BP) |>
      dplyr::mutate( BPcum=BP+tot ) |>

      # Add highlight and annotation information
      dplyr::mutate( is_highlight=ifelse(SNP %in% highlight, "yes", "no")) |>
      dplyr::mutate( is_annotate=ifelse(-log10(P)>4, "yes", "no"))

    # Prepare X axis
    axisdf <- d |>
      dplyr::group_by(CHR) |>
      dplyr::summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

    # Make the plot
    plot <- ggplot2::ggplot(d, ggplot2::aes(x=BPcum, y=-log10(P))) +

      # Show all points
      ggplot2::geom_point( ggplot2::aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
      ggplot2::scale_color_manual(values = rep(c(palette[5], palette[6]), 22 )) +

      # custom X axis:
      ggplot2::scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
      ggplot2::scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis

      # Add highlighted points
      ggplot2::geom_point(data=subset(d, is_highlight=="yes"), color="orange", size=2) +

      # Add label using ggrepel to avoid overlapping
      ggrepel::geom_label_repel( data=subset(d, is_annotate=="yes"), ggplot2::aes(label=SNP), size=2) +

      # Custom the theme:
      ggplot2::theme_minimal() +
      #ggplot2::theme_bw() +
      ggplot2::labs(...,
                    x = "Chromosome",
                    y = "-log10(P)") +
      ggplot2::theme(
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
      )

    # save plot
    ggplot2::ggsave(plot=plot, filename=file_path, device="png", bg="white", height=15, width=35, units="cm")

    # return
    invisible(plot)
  }
)


