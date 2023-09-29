load_all()

s <- GWASsumstats(dir = "/Users/xx20081/Downloads/hermes_progression/solid",
                  ref_path = "/Users/xx20081/Downloads/ref_1000GP_Phase3/ref_1000GP_Phase3_legend_cptid.gz",
                  pre_qc_dir = "pre_qc",
                  post_qc_dir = "post_qc",
                  file_structure = list(
                    "allcause_death" = list(
                      "autosomes" = "(?i)^(?!.*(?:fe)?male).*allcause.*",
                      "xchr_male" = "(?i)^(?=.*allcause)(?=.*male)(?!.*female).*",
                      "xchr_female" = "(?i)^(?=.*allcause)(?=.*female).*"
                    ),
                    "composite_1" = list(
                      "autosomes" = "(?i)^(?!.*(?:fe)?male).*comp(?:osite)?1.*",
                      "xchr_male" = "(?i)^(?=.*comp(?:osite)?1)(?=.*male)(?!.*female).*",
                      "xchr_female" = "(?i)^(?=.*comp(?:osite)?1)(?=.*female).*"
                    ),
                    "composite_2" = list(
                      "autosomes" = "(?i)^(?!.*(?:fe)?male).*comp(?:osite)?2.*",
                      "xchr_male" = "(?i)^(?=.*comp(?:osite)?2)(?=.*male)(?!.*female).*",
                      "xchr_female" = "(?i)^(?=.*comp(?:osite)?2)(?=.*female).*"
                    )
                  )


)

#### EasyQC
s <- run_qc(s, "allcause_death", "autosomes")

#### QC plots
load_all()
run_qc_plots(s, "/Users/xx20081/Downloads/figures", "allcause_death", "autosomes") #, "composite_1", "composite_2"))


pre <- data.table::fread("/Users/xx20081/Downloads/hermes_progression/solid/pre_qc/allcause_chrN_solid.txt.gz")
post <- data.table::fread("/Users/xx20081/Downloads/hermes_progression/solid/post_qc/solid_allcause_death_autosomes_post_qc.gz")

check <- sum(post$ORI_OTHER_ALLELE==post$OTHER_ALLELE)










corpus <- StudyCorpus(corpus_dir ="/Users/xx20081/Downloads/hermes_progression", #/Users/xx20081/Documents/local_data/hermes_progression", #
                      study_type = "GWASsumstats",
                      ref_path = "/Users/xx20081/Downloads/ref_1000GP_Phase3/ref_1000GP_Phase3_legend_cptid.gz",
                      mapping = StudyManager::base_column_mapping,
                      file_structure = list(
                        "allcause_death" = list(
                          "autosomes" = "(?i)^(?!.*(?:fe)?male).*allcause.*",
                          "xchr_male" = "(?i)^(?=.*allcause)(?=.*male)(?!.*female).*",
                          "xchr_female" = "(?i)^(?=.*allcause)(?=.*female).*"
                        ),
                        "composite_1" = list(
                          "autosomes" = "(?i)^(?!.*(?:fe)?male).*comp(?:osite)?1.*",
                          "xchr_male" = "(?i)^(?=.*comp(?:osite)?1)(?=.*male)(?!.*female).*",
                          "xchr_female" = "(?i)^(?=.*comp(?:osite)?1)(?=.*female).*"
                        ),
                        "composite_2" = list(
                          "autosomes" = "(?i)^(?!.*(?:fe)?male).*comp(?:osite)?2.*",
                          "xchr_male" = "(?i)^(?=.*comp(?:osite)?2)(?=.*male)(?!.*female).*",
                          "xchr_female" = "(?i)^(?=.*comp(?:osite)?2)(?=.*female).*"
                        )
                      )
)





#
# summary_dt <- run_filter_summary_plots(corpus, "/Users/xx20081/Downloads/figures", parallel_cores=12)
#
#
# # jobs, rerun all of the QC as changed filters
# # create the summary figures
# # create the zoomed fi
#
#
#
# corpus <- run_qc(corpus, index=2) # "allcause_death", "xchr_male",
#
# corpus <- run_qc_plots(corpus, "/Users/xx20081/Downloads/figures", c("allcause_death"), index=1)#"composite_1", "composite_2"
#
#
#
# load_all()

# exclude
excluded <- "ephesus"
included_idx <- which(!names(studies(corpus)) %in% excluded)

corpus <- run_gwama(corpus, "/Users/xx20081/Downloads/meta_analysis_output", index=included_idx) #, parallel_cores=4)

#corpus <- create_results_list(corpus, "GWASsumstats", "/Users/xx20081/Downloads/meta_analysis_output")

corpus <- run_meta_plots(corpus, "/Users/xx20081/Downloads/figures")










#studies(corpus, "bioshift_triumph") <- run_qc( studies(corpus, "bioshift_triumph") )




x = data.table::fread("/Users/xx20081/Downloads/hermes_progression/solid/pre_qc/allcause_chrN_solid.txt.gz")













foo_pre<-data.table::fread("/Users/xx20081/Downloads/hermes_progression/cathgen/pre_qc/composite1_chrN_cathgen_maf5.txt.gz")
View(head(foo_pre))
foo<- data.table::fread("/Users/xx20081/Downloads/hermes_progression/cathgen/post_qc/cathgen_composite_1_autosomes_post_qc.notinref.txt")
View(head(foo))












#
# corpus <- run_qc(corpus, parallel_cores=4)
#
# run_qc_plots(corpus, "/Users/xx20081/git/thesis/vignettes/figures", c("allcause_death", "composite_1", "composite_2") ) #, parallel_cores=4)
#
#




load_all()











# d1 <- DataFile(
#   path="/Users/nicholassunderland/Downloads/hermes_progression/bioshift_triumph/pre_qc/bioshift_triumph.females.allcause.gz"
# )
#
# d1@mapping <- set_active(d1@mapping, c("SIGNED_SUMSTAT"), rest.off = F)
#
# d1 <- extract( d1 )
#
# head( d1@data )
#
# foo <- get_data( d1 )



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


#
# load_all()
#
#
# d <- DataFile(
#   path="/Users/nicholassunderland/Downloads/hermes_progression/bioshift_triumph/pre_qc/bioshift_triumph.females.allcause.gz"
# )
#
# d2 <- DataFile(
#   path="/Users/nicholassunderland/Downloads/hermes_progression/bioshift_triumph/pre_qc/bioshift_triumph.males.allcause.gz"
# )
#
# df_list=list("F"=d,"M"=d2)
#
# d_single <- extract(d)
#
# d_list <- extract(df_list)
#
# head(get_data(d_single))
# head(get_data(d_list))
#
# d_single <- free(d_single)
# d_list <- free(d_list)
#
# head(get_data(d_single))
# head(get_data(d_list))

# n=100000
# s@data_files$allcause_death$xchr_male <- extract( s@data_files$allcause_death$xchr_male )
# dt = head(s@data_files$allcause_death$xchr_male@data, n=n) |>
#   mutate(CHR_FCT = "autosomes",
#          QC_STATUS = c(rep("PRE-QC",n/2),rep("POST-QC",n/2)),
#          FRQ_DIFF_FCT = sample(c("Testing","3ds","fds"),n,replace=TRUE),
#          INFO = runif(n),
#          EUR_FRQ = runif(n))

