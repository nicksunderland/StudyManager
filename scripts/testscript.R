load_all()

# d <- DataFile(path="/Users/nicholassunderland/Downloads/hermes_progression/bioshift_triumph/pre_qc/bioshift_triumph.females.allcause.gz")
# d<- extract(d)
# head(d@data)
# d<- free(d)
# head(d@data)
# #
# #
# #
# ref_data_file <- DataFile(
#   path = "/Users/nicholassunderland/Downloads/genome_reference/ref_justX_1000GP_Phase3_maf_biallelic.gz",
#   mapping = turn_on(StudyManager::ColumnMapping, c("CPTID","SNP","CHR","BP","MARKER_TYPE","A*","OTHER_ALLELE","EUR_FRQ"))
# )
# ref_data_file <-extract(ref_data_file)

s <- GWASsumstats(dir = "/Users/nicholassunderland/Downloads/hermes_progression/bioshift_triumph",
                  pre_qc_dir = "pre_qc",
                  post_qc_dir = "post_qc",
                  file_structure = list(
                    "allcause_death" = list(
                      "autosomes" = "(?i)^(?!.*(?:fe)?male).*allcause.*",
                      "xchr_male" = "(?i)^(?=.*allcause)(?=.*male)(?!.*female).*",
                      "xchr_female" = "(?i)^(?=.*allcause)(?=.*female).*"
                    ),
                    "composite_1" = list(
                      "autosomes" = ".*allcause_death_post_qc.tsv"
                    )
                  ),
                  reference_path = "/Users/nicholassunderland/Downloads/genome_reference/ref_justX_1000GP_Phase3_maf_biallelic.gz"
)


#### EasyQC
s <- run_qc(s)
# s <- run_qc(s, "allcause_death", "xchr_female")


#### DATAFILE

# d <- DataFile(path="/Users/nicholassunderland/Downloads/hermes_progression/bioshift_triumph/pre_qc/bioshift_triumph.females.allcause.gz",
#               col_names = c("P","POS","CHR"))
# d <- DataFile(
#   path = "/Users/nicholassunderland/Downloads/hermes_progression/bioshift_triumph/bioshift_triumph.females.allcause.gz"
#   # col_names = c("PVAL"="P", "POSITION"="POS", "Foo"="Foo"),
#   # col_types = list("POS"=function(x) paste(x,"test"), "P"=as.character, "Foo"=as.character),
#   # col_fill = NA_real_
# )
#
# d <- extract(d)
# foo<-data(d)
# head(data(d))


#### Get data
# foo2 <- get_data(s, "allcause_death", "xchr_female")


#
# create_named_nested_list <- function(names_vector) {
#   if (length(names_vector) == 1) {
#     element <- "NA"
#     return(list(element))
#   } else {
#     nested_list <- create_named_nested_list(names_vector[-1])
#     names(nested_list) <- names_vector[-1][1]
#     final_list <- list(nested_list)
#     names(final_list) <- names_vector[1]
#     return(final_list)
#   }
# }
#
# # Example character vector
# names_vector <- c("level1", "level2", "level3", "level4")
#
# # Create the named nested list
# nested_list <- create_named_nested_list(names_vector)
#
# # Print the structure of the nested list
# str(nested_list)




