load_all()

s <- GWASsumstats(dir = "/Users/nicholassunderland/Downloads/hermes_progression/bioshift_triumph/pre_qc",
                  file_structure = list(
                    "allcause_death" = list(
                      "autosomes" = "(?i)^(?!.*(?:fe)?male).*allcause.*",
                      "xchr_male" = "(?i)^(?=.*allcause)(?=\\.*male)(?!.*female).*",
                      "xchr_female" = "(?i)^(?=.*allcause)(?=.*female).*"
                    ),
                    "composite_1" = list(
                      "autosomes" = ".*allcause_death_post_qc.tsv"
                    )
                  ),
                  col_names = c("PVAL"="P", "POSITION"="POS", "Foo"="Foo"),
                  col_types = list("POS"=function(x) paste(x,"test"), "P"=as.character, "Foo"=as.character)
)

foo2 <- get_data(s, "allcause_death", "xchr_female")
