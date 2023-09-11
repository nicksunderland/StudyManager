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
    ref_cols <- c("CPTID","SNP","CHR","BP","MARKER_TYPE","A0","OTHER_ALLELE","EUR_FRQ")
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

stopifnot("Minimum required reference columns: c('CPTID','SNP','CHR','BP','MARKER_TYPE','A0','OTHER_ALLELE','EUR_FRQ')" =
            all(c('CPTID','SNP','CHR','BP','MARKER_TYPE','A0','OTHER_ALLELE','EUR_FRQ') %in% names(active_cols(object@reference_data_file))))

  }
)

setGeneric("run_qc", function(object, ...) standardGeneric("run_qc"))
setMethod(
  f = "run_qc",
  signature = c("GWASsumstats"),
  definition = function(object, ...) {

    # extract the ref data
    object@reference_data_file <- extract( object@reference_data_file )
    tmp_ref_path <- tempfile()
    print(tmp_ref_path)
    write_file(object@reference_data_file, tmp_ref_path)
    object@reference_data_file <- free( object@reference_data_file )
    ref_map <- mapping(object@reference_data_file)
    object@reference_data_file <- DataFile(path=tmp_ref_path, mapping=ref_map)

    # single data file use
    if(length(c(...) > 0)) {

      object <- easy_qc(object, ...)

    # recursively process all data_files
    } else {

      key_list <- get_keys(object@data_files)
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

    # log
    rlog::log_info("Starting EASYQC")

    # ensure valid keys into the data_file structure

    # print(c(...))
    # checks failing ?????

    # stopifnot("Internal code error, no `data_file` keys passed to `easy_qc()`" = length(c(...)) > 0)
    # if(!keys_valid(object@data_files, ...)) {
    #   stop("All arguments in '...' must be valid ordered keys into the `data_files` structure")
    # }

    # the DataFile data
    object@data_files[[ c(...) ]] <- extract( object@data_files[[ c(...) ]] )

    # create a written temp file from the DataFile (appropriate col names enforced etc.)
    mapping <- mapping(object@data_files[[ c(...) ]])
    input_path <- tempfile()
    input_types <- col_types(mapping)
    input_cols <- names(input_types)
    fwrite(x = file_data( object@data_files[[ c(...) ]] ),
           file = input_path,
           sep = "\t")

    # free the data
    object@data_files[[ c(...) ]] <- free( object@data_files[[ c(...) ]] )

    # create an output directory if it doesn't exist
    output_dir <- file.path(object@dir, object@post_qc_dir)
    output_path <- file.path(output_dir, paste0(c(...), collapse="_"))
    if(!dir.exists(output_dir)) {
      dir.create(output_dir, showWarnings=TRUE, recursive=TRUE)
    } else {
      old_run_files <- list.files(output_dir, pattern=paste0(basename(output_path),"(.[0-9]+)?(.rep|.ecf|.out|.notinref.txt|.gz)"), full.names=TRUE)
      unlink(old_run_files)
    }

    # make sure the output_path is useable
    stopifnot("`output_path` should be a writable file path" = checkmate::check_path_for_output(output_path, overwrite=TRUE))

    # get the written temp file from the reference DataFile
    tmp_ref_path <- path(object@reference_data_file)
    print(tmp_ref_path)
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
          --colRefMarker CPTID
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

    rlog::log_info(glue::glue("ECF file\n:{ecf_str}"))

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
      ref_map_names <- col_names(mapping(object@reference_data_file))
      easy_qc_names <- c("CALL_RATE", "AMBIGUOUS", "STRAND", "CPTID")
      out_map <- turn_on(post_qc_map, unique(c(map_names, ref_map_names, easy_qc_names)))

      if(file.exists(paste0(output_path, ".gz"))) {

        object@qc_data_files[[ c(...) ]] <- DataFile(path=paste0(output_path, ".gz"),
                                                     mapping=post_qc_map)

      } else {

        empty_dt <- data.table::data.table(matrix(NA, 0, length(col_names(post_qc_map))))
        data.table::setnames(empty_dt, col_names(post_qc_map))
        object@qc_data_files[[ c(...) ]] <- DataFile(data=empty_dt,
                                                     mapping=post_qc_map)
        write_file(object@qc_data_files[[ c(...) ]], paste0(output_path, ".gz"))

      }
    }

    validObject(object)
    rlog::log_info("Finished EASYQC")
    return(object)
  }
)

