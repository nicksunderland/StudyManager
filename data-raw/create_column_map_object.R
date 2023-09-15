
a = list()
t = character()
f = list()
l = logical()
map = MungeSumstats:::sumstatsColHeaders
map <- map |>
  rbind(c("QUAL_SCORE", "INFO"),
        c("STRAND", "STRAND"),
        c("N_EVENTS", "N_CAS"),
        c("N_SAMPLE", "N"),
        c("TYPE", "MARKER_TYPE"),
        c("EUR", "EUR_FRQ"),
        c("cptid", "cptid"),
        c("CPTID", "cptid"),
        c("CALL_RATE","CALL_RATE"),
        c("AMBIGUOUS", "AMBIGUOUS")
        ) |>
  dplyr::mutate(Corrected = dplyr::if_else(Corrected=="A1", "OTHER_ALLELE", Corrected),
                Corrected = dplyr::if_else(Corrected=="A2", "EFFECT_ALLELE", Corrected),
                Corrected = dplyr::if_else(Corrected=="A*", "A0", Corrected)) |>
  dplyr::filter(!(Uncorrected=="MAJORALLELE" & Corrected=="EFFECT_ALLELE"))

for(name in unique(map$Corrected)) {

  uncorrected <- map[map$Corrected==name,"Uncorrected"]
  uncorrected <- c(uncorrected, name)
  aliases = list( unique(uncorrected) )
  names(aliases)<- name
  a<-c(a,aliases)
  stopifnot("Mapping cannot contain duplicate aliases" = !any(duplicated(unname(unlist(aliases)))))

  types = c("Z"="numeric",
            "SNP"="character",
            "SIGNED_SUMSTAT"="numeric",
            "SE"="numeric",
            "P"="numeric",
            "OR"="numeric",
            "NSTUDY"="integer",
            "N_CON"="integer",
            "N_CAS"="integer",
            "N"="integer",
            "LOG_ODDS"="numeric",
            "INFO"="numeric",
            "HETPVAL"="numeric",
            "HETISQT"="numeric",
            "HETDF"="numeric",
            "HETCHISQ"="numeric",
            "FRQSE"="numeric",
            "FRQMIN"="numeric",
            "FRQMAX"="numeric",
            "FRQ"="numeric",
            "DIRECTION"="character",
            "STRAND" = "character",
            "CHR"="character",
            "BP"="integer",
            "BETA"="numeric",
            "AC"="integer",
            "OTHER_ALLELE"="character",
            "EFFECT_ALLELE"="character",
            "A0"="character",
            "MARKER_TYPE"="character",
            "EUR_FRQ"="numeric",
            "cptid"="character",
            "AMBIGUOUS"="logical",
            "CALL_RATE"="numeric")
  t<-c(t,types[name])

  func = list("Z"=as.numeric,
              "SNP"=as.character,
              "SIGNED_SUMSTAT"=as.numeric,
              "SE"=as.numeric,
              "P"=as.numeric,
              "OR"=as.numeric,
              "NSTUDY"=as.integer,
              "N_CON"=as.integer,
              "N_CAS"=as.integer,
              "N"=as.integer,
              "LOG_ODDS"=as.numeric,
              "INFO"=as.numeric,
              "HETPVAL"=as.numeric,
              "HETISQT"=as.numeric,
              "HETDF"=as.numeric,
              "HETCHISQ"=as.numeric,
              "FRQSE"=as.numeric,
              "FRQMIN"=as.numeric,
              "FRQMAX"=as.numeric,
              "FRQ"=as.numeric,
              "DIRECTION"=as.character,
              "STRAND" = as.character,
              "CHR"=as.character,
              "BP"=as.integer,
              "BETA"=as.numeric,
              "AC"=as.integer,
              "EFFECT_ALLELE"=as.character,
              "OTHER_ALLELE"=as.character,
              "A0"=as.character,
              "MARKER_TYPE"=as.character,
              "EUR_FRQ"=as.numeric,
              "cptid"=as.character,
              "AMBIGUOUS"=as.logical,
              "CALL_RATE"=as.numeric)
  f<-c(f,func[name])

  if(name %in%  c("cptid","SNP","SE","P","N_CAS","N","INFO","FRQ","STRAND","CHR","BP","BETA","EFFECT_ALLELE","OTHER_ALLELE")) {
    logi <- TRUE
  } else {
    logi <- FALSE
  }
  names(logi) <- name
  l <- c(l,logi)

}

ColumnMapping = ColMap(active = l, aliases = a, col_types = t, funcs = f)

usethis::use_data(ColumnMapping,internal = T, overwrite=TRUE)

