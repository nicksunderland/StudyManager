#' Title
#'
#' @slot parameter numeric.
#'
#' @return object
#' @export
#' @importFrom methods callNextMethod
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
    reference_path = "character",
    reference_data_file = "DataFile",
    mapping = "ColMap"
  ),
  prototype = list(
    qc_freq_threshold = 0.01,
    qc_low_events_n = 200,
    qc_freq_low_events_thresh = 0.05,
    qc_freq_ambig_thresh = 0.40,
    qc_qual_threshold = 0.5,
    qc_call_rate_threshold = 0.98,
    qc_freq_diff_threshold = 0.2,
    reference_path = character(),
    reference_data_file = NULL,
    mapping = StudyManager::ColumnMapping
  )
)


setMethod(
  f = "initialize",
  signature = "GWASsumstats",
  definition = function(.Object, reference_path,
                                 reference_data_file = NULL,
                                 mapping = StudyManager::ColumnMapping, ...) {

    # required ref cols
    ref_map <- mapping
    ref_cols <- c("cptid","SNP","CHR","BP","MARKER_TYPE","A0","OTHER_ALLELE","EUR_FRQ")
    ref_map <- turn_on(ref_map, ref_cols)

    # create the ref data file
    .Object@reference_path <- reference_path
    .Object@reference_data_file <- DataFile(path=.Object@reference_path, mapping=ref_map)

    # required gwas cols , set mapping back to these
    col_map <- mapping
    gwas_sumstats_cols <- c("SNP","SE","P","N_CAS","N","INFO","FRQ","CHR","BP","BETA","EFFECT_ALLELE","OTHER_ALLELE","STRAND")
    .Object@mapping <- turn_on(col_map, gwas_sumstats_cols)

    # initialise the rest
    .Object <- callNextMethod(.Object, ...)
    validObject(.Object)
    return(.Object)
  }
)


setValidity(
  Class = "GWASsumstats",
  method = function(object) {

# stopifnot("Minimum required reference columns: c('cptid','SNP','CHR','BP','MARKER_TYPE','A0','OTHER_ALLELE','EUR_FRQ')" =
#             all(c('cptid','SNP','CHR','BP','MARKER_TYPE','A0','OTHER_ALLELE','EUR_FRQ') %in% names(active_cols(object@reference_data_file))))

  }
)


setGeneric("run_qc", function(object, ...) standardGeneric("run_qc"))
setMethod(
  f = "run_qc",
  signature = c("GWASsumstats"),
  definition = function(object, ...) {

    # extract the ref data
    rlog::log_info("Reading reference data")
    object@reference_data_file <- extract( object@reference_data_file )
    tmp_ref_path <- tempfile()
    write_file(object@reference_data_file, tmp_ref_path)
    object@reference_data_file <- free( object@reference_data_file )
    ref_map <- mapping(object@reference_data_file)
    object@reference_data_file <- DataFile(path=tmp_ref_path, mapping=ref_map)

    # create an output directory if it doesn't exist; clean up
    output_dir <- file.path(object@dir, object@post_qc_dir)
    if(!dir.exists(output_dir)) {
      dir.create(output_dir, showWarnings=TRUE, recursive=TRUE)
    } else {
      old_run_files <- list.files(output_dir, pattern="post_qc(.[0-9]+)?(.rep|.ecf|.out|.notinref.txt|.gz)", full.names=TRUE)
      unlink(old_run_files)
    }

    # single data file use
    if(length(c(...) > 0)) {

      object <- easy_qc(object, ...)

    # recursively process all data_files
    } else {

      key_list <- get_data_keys(object@data_files)
      for(keys in key_list) {
        args = c(list(object), as.list(keys))
        object <- do.call("easy_qc", args)
      }

    }

    validObject(object)
    return(object)
  }
)


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

    # print(c(...))
    # checks failing ?????

    # stopifnot("Internal code error, no `data_file` keys passed to `easy_qc()`" = length(c(...)) > 0)
    # if(!keys_valid(object@data_files, ...)) {
    #   stop("All arguments in '...' must be valid ordered keys into the `data_files` structure")
    # }

    # the DataFile data
    rlog::log_info("Extracting file data")
    object@data_files[[ c(...) ]] <- extract( object@data_files[[ c(...) ]] )

    # create a written temp file from the DataFile (appropriate col names enforced etc.)
    mapping <- mapping(object@data_files[[ c(...) ]])
    input_path <- tempfile()
    input_types <- col_types(mapping)
    input_cols <- names(input_types)

    rlog::log_info("Writting standardised file")
    fwrite(x = file_data( object@data_files[[ c(...) ]] ),
           file = input_path,
           sep = "\t")

    # free the data
    rlog::log_info("Freeing memory / DataFile data")
    object@data_files[[ c(...) ]] <- free( object@data_files[[ c(...) ]] )

    # create output path
    output_path <- file.path(object@dir, object@post_qc_dir, paste0(c(basename(object@dir), ..., "post_qc"), collapse="_"))

    # make sure the output_path is useable
    stopifnot("`output_path` should be a writable file path" = checkmate::check_path_for_output(output_path, overwrite=TRUE))

    # get the written temp file from the reference DataFile
    tmp_ref_path <- path(object@reference_data_file)
    ref_mapping <- mapping(object@reference_data_file)
    ref_input_types <- col_types(ref_mapping)
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

    # write the file
    ecf_path <- paste0(output_path, ".ecf")
    writeLines(ecf_str, ecf_path)

    # run EasyQC
    easy_qc_success <- EasyQC::EasyQC(ecf_path)

    # get the QC data n.b. no file is written if empty data.table so check forexistance
    if(!easy_qc_success) {

      rlog::log_error("EasyQC failed")

    } else {

      # save DataFile
      post_qc_map <- object@mapping
      map_names <- col_names(object@mapping)
      easy_qc_names <- c("CALL_RATE","AMBIGUOUS","STRAND","cptid")
      post_qc_map <- turn_on(post_qc_map, unique(c(map_names, easy_qc_names)))
      rlog::log_info(glue::glue("EasyQC output DataFiles mapping set to: {paste0(col_names(post_qc_map), collapse=',')}"))

      if(file.exists(paste0(output_path, ".gz"))) {

        object@qc_data_files[[ c(...) ]] <- DataFile(path=paste0(output_path, ".gz"),
                                                     mapping=post_qc_map)

      } else {

        empty_dt <- data.table::data.table(matrix(NA, 0, length(col_names(post_qc_map))))
        data.table::setnames(empty_dt, col_names(post_qc_map))
        object@qc_data_files[[ c(...) ]] <- DataFile(path=paste0(output_path, ".gz"),
                                                     data=empty_dt,
                                                     mapping=post_qc_map)
        write_file(object@qc_data_files[[ c(...) ]], paste0(output_path, ".gz"))

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
  signature = c("GWASsumstats", "character"),
  definition = function(object, output_dir, ...) {

    stopifnot(dir.exists(output_dir))

    # specific data files used
    if(length(c(...) > 0)) {

      key_list <- list(c(...))
      rlog::log_debug(glue::glue("Keys provdied: {paste0(unlist(key_list), collapse=',')}"))

      # recursively process all data_files individually
    } else {

      key_list <- get_data_keys(object@data_files)
      rlog::log_debug(glue::glue("No keys provdied, using all: {paste0(unlist(key_list), collapse=',')}"))

    }

    # for each data source
    for(keys in key_list) {

      # set up arguments for plotting
      rlog::log_debug(glue::glue("Extracting pre-qc plot data with keys: '{paste0(unlist(keys), collapse=',')}'"))
      data_file <- extract( get_data_file(object, unlist(key_list) ), merge_extract="CHR_FCT" )
      rlog::log_debug(glue::glue("Extracting post-qc plot data with keys: {paste0(unlist(keys), collapse=',')}"))
      qc_data_file <- extract( get_qc_data_file(object, unlist(key_list) ), merge_extract="CHR_FCT" )
      pre_snps  <- file_data(data_file)[, get("SNP")]
      post_snps <- file_data(qc_data_file)[, get("SNP")]
      filtered_snps_logi <- !pre_snps %in% post_snps
      filtered_snps <- pre_snps[filtered_snps_logi]
      file_path <- file.path(output_dir, paste0(paste0(c(basename(object@dir), unlist(keys)), collapse="_")))

      # subsets
      req_cols <- c("cptid","P","BETA","INFO","FRQ","CHR_FCT")
      req_ref_cols <- c("cptid","EUR_FRQ")
      pre_dt <- file_data(turn_on(data_file, req_cols))
      dt <- rbind(
        cbind(,
              list(QC_STATUS = rep("PRE-QC",  nrow(file_data(data_file))))),
        cbind(file_data(qc_data_file)[, req_cols,, with=FALSE],
              list(QC_STATUS = rep("POST-QC", nrow(file_data(qc_data_file)))))
      )

      #########

      object@reference_data_file <- extract( object@reference_data_file )
      print("ref:")
      warning(paste0(colnames( file_data(object@reference_data_file) ), collapse=","))
      print("dt:")
      warning(paste0(colnames( dt ), collapse=","))
      dt[file_data(object@reference_data_file), on='cptid']
      print(head(dt))
        # A[B, on = 'a', bb := i.b].   https://stackoverflow.com/questions/34598139/left-join-using-data-table

      # compute extra columns
      rlog::log_info("Creating plotting data...")
      dt <- dplyr::bind_rows(file_data(data_file)    |> dplyr::select(req_cols) |> dplyr::mutate(QC_STATUS="PRE-QC"),
                             file_data(qc_data_file) |> dplyr::select(req_cols) |> dplyr::mutate(QC_STATUS="POST-QC")) |>
        dplyr::mutate(FRQ_DIFF_FCT = dplyr::case_when(is.na(EUR_FRQ)
                                                        ~ factor("No reference data"),
                                                      abs(FRQ-EUR_FRQ)> object@qc_freq_diff_threshold
                                                        ~ factor("EAF outlier"),
                                                      abs(FRQ-EUR_FRQ)<=object@qc_freq_diff_threshold
                                                        ~ factor("EAF within tolerance"),
                                                      TRUE
                                                        ~ factor(NA_character_)))
      rlog::log_debug(glue::glue("Plotting data.table columns: {paste0(dt, collapse=', ')}"))

      # Manhattan
      create_manhattan(dt = file_data(data_file)[, c("SNP","CHR","BP","P")],
                       file_path = paste0(file_path, "_qc_manhattan.png"),
                       highlight = filtered_snps,
                       main = paste0(c("Manhattan Plot -", unlist(keys)), collapse=" "))

      # QQ-plot
      create_qq(dt = dt[, c("P","QC_STATUS")],
                file_path = paste0(file_path, "_qc_qq_plot.png"),
                title = paste0(c("QQ plot -", unlist(keys)), collapse=" "))

      # P value histogram
      create_pvalhist(dt = dt[, c("P","QC_STATUS")],
                      file_path = paste0(file_path, "_qc_pval_hist.png"),
                      title = paste0(c("P value histogram -", unlist(keys)), collapse=" "))

      # Effect size histogram
      create_eshist(dt = dt[, c("BETA","QC_STATUS")],
                    file_path = paste0(file_path, "_qc_es_hist.png"),
                    title = paste0(c("Effect size histogram -", unlist(keys)), collapse=" "))

      # Imputation quality histogram
      create_impqualhist(dt = dt[, c("INFO","QC_STATUS")],
                         file_path = paste0(file_path, "_qc_impqual_hist.png"),
                         title = paste0(c("Imputation quality histogram -", unlist(keys)), collapse=" "))

      # allele frequency plots
      create_eafplot(dt = dt[, c("FRQ","FRQ_DIFF_FCT","QC_STATUS","CHR_FCT")],
                     file_path = paste0(file_path, "_qc_eaf_plot.png"),
                     title = paste0(c("Allele frequency -", unlist(keys)), collapse=" "))




    }

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

    # start PNG
    png(file=file_path, width=800, height=350)

    # plot
    ggplot2::ggplot(dt, aes(x = QUAL_SCORE)) +
      geom_histogram(color = palette[7],
                     fill = palette[5],
                     alpha = 0.9,
                     binwidth = 0.01) +
      labs(...,
           x = "QUAL_SCORE",
           y = "Count") +
      theme_minimal() +
      facet_grid(cols = vars(QC_STATUS))

    # stop PND
    dev.off()
  }
)


setGeneric("create_eafplot", function(dt, file_path=getwd(), ...) standardGeneric("create_eafplot"))
setMethod(
  f = "create_eafplot",
  signature = c("data.table"),
  definition = function(dt, file_path=getwd(), ...) {

    # logging info
    rlog::log_info("Allele frequency plot...")
    rlog::log_info(glue::glue("Writing image to: {file_path}"))

    # choose a pallete
    palette <- c(wesanderson::wes_palette("GrandBudapest1"), wesanderson::wes_palette("GrandBudapest2"))
    #c("#F1BB7B", "#FD6467", "#5B1A18", "#D67236", "#E6A0C4", "#C6CDF7", "#D8A499", "#7294D4")

    colour_labels <- setNames(c(palette[3],palette[2],palette[8]),
                              c("No reference data", "EAF outlier", "EAF within tolerance"))

    # start PNG
    png(file=file_path, width=800, height=1050)

    # plot
    ggplot2::ggplot(dt, aes(x = EUR_FRQ, y = FRQ, color = FRQ_DIFF_FCT)) +
      geom_point(size = 0.1) +
      geom_abline(slope=1, intercept= freq_diff_threshold, color=palette[3], linetype = "dotted") +
      geom_abline(slope=1, intercept=-freq_diff_threshold, color=palette[3], linetype = "dotted") +
      lims(x=c(-0.08,1), y=c(0,1)) +
      labs(...,
           x = "Ref EAF",
           y = "Study EAF",
           color = NULL) +
      scale_color_manual(values=colour_labels) +
      theme_minimal() +
      guides(color = guide_legend(override.aes = list(size = 5))) +
      facet_grid(rows = vars(CHR_FCT), cols = vars(QC_STATUS))

dev.off()
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

    # start PNG
    png(file=file_path, width=800, height=350)

    # plot
    histogram_w_outlier_bins(data = dt,
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
      labs(...) +
      theme_minimal() +
      facet_grid(cols = vars(QC_STATUS))

    # stop PND
    dev.off()
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

    # start PNG
    png(file=file_path, width=800, height=350)

    ggplot(dt, aes(x = P)) +
      geom_histogram(binwidth = 0.01,
                     color=palette[6],
                     fill=palette[8]) +
      labs(...,
           x = "P value",
           y = "Count") +
      theme_minimal() +
      facet_grid(cols = vars(QC_STATUS))

    # stop PND
    dev.off()
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

    # start PNG
    png(file=file_path, width=800, height=350)

    # generate PP-plot but with ggplot instead of qqman
    plot <- dt |>
      filter(!is.na(P),
             !is.nan(P),
             !is.null(P),
             is.finite(P),
             P < 1.0,
             P > 0.0) |>
      group_by(QC_STATUS) |>
      mutate(observed = -log10(sort(P, decreasing=FALSE)),
             expected = -log10(ppoints(n()))) |>
      ungroup() |>
      ggplot(aes(x = expected, y = observed)) +
      geom_point(size = 0.5, color=palette[8]) +
      geom_abline(slope=1, intercept=0, color="darkred", linetype = "dotted") +
      labs(...,
           x = "Expected (-log10 P)",
           y = "Observed (-log10 P)",
           color = NULL) +
      theme_minimal() +
      facet_grid(cols = vars(QC_STATUS))

    plot

    # stop PND
    dev.off()
  }
)

setGeneric("create_manhattan", function(dt, file_path=getwd(), ...) standardGeneric("create_manhattan"))
setMethod(
  f = "create_manhattan",
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
    rlog::log_info("Manhattan plot...")
    rlog::log_info(glue::glue("Writing image to: {file_path}"))

    # start PNG
    png(file=file_path, width=800, height=450)

    # plot
    qqman::manhattan(x = dt |> dplyr::mutate(CHR = dplyr::if_else(CHR=="X", 23, as.numeric(CHR))),
                     chr = "CHR",
                     bp = "BP",
                     p = "P",
                     snp = "SNP",
                     ylim = c(0, 10), cex = 0.6, cex.axis = 0.9,
                     col = palette[5:6],
                     suggestiveline = -log10(1e-5),
                     genomewideline = -log10(5e-8),
                     ...)

    # stop PND
    dev.off()
  }
)


