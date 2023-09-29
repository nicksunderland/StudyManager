
a = list()
t = character()
f = list()
l = logical()
map = MungeSumstats:::sumstatsColHeaders
map <- map |>
  # binding in other mappings that I've found
       #c(UNCORRECTED, CORRECTED)
  rbind(c("QUAL_SCORE", "INFO"),
        c("STRAND", "STRAND"),
        c("N_EVENTS", "N_CAS"),
        c("N_SAMPLE", "N"),
        c("TYPE", "MARKER_TYPE"),
        c("EUR", "EUR_FRQ"),
        c("EUR_FRQ_REF", "EUR_FRQ"),
        c("cptid", "cptid"),
        c("CPTID", "cptid"),
        c("CALL_RATE","CALL_RATE"),
        c("AMBIGUOUS", "AMBIGUOUS"),
        c("rs_number", "SNP"),
        c("beta_95L", "BETA_95L"),
        c("beta_95U", "BETA_95U"),
        c("_-log10_p-value", "LOG10_P"),
        c("q_statistic", "Q_STATISTIC"),
        c("q_p-value", "Q_P_VALUE"),
        c("n_studies", "NSTUDY"),
        c("effects", "DIRECTION"),
        c("i2", "HETISQT"),
        c("IMPUTED", "IMPUTED"),
        c("Imputed", "IMPUTED"),
        c("imputed", "IMPUTED"),
        c("ORI_OTHER_ALLELE", "ORI_OTHER_ALLELE"),
        c("ORI_OTHER_ALLELE_REF", "ORI_OTHER_ALLELE_REF"),
        c("ORI_EFFECT_ALLELE", "ORI_EFFECT_ALLELE"),


        c("OTHER_ALLELE_REF", "OTHER_ALLELE_REF"),
        c("A0_REF", "A0"),
        c("ORI_A0", "ORI_A0"),
        c("ORI_A0_REF", "ORI_A0"),
        c("REF_RSID", "REF_RSID")) |>

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
            "LOG10_P"="numeric",
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
            "BETA_95L"="numeric",
            "BETA_95U"="numeric",
            "Q_STATISTIC"="numeric",
            "Q_P_VALUE"="numeric",
            "AC"="integer",
            "OTHER_ALLELE"="character",
            "EFFECT_ALLELE"="character",
            "A0"="character",
            "MARKER_TYPE"="character",
            "EUR_FRQ"="numeric",
            "cptid"="character",
            "AMBIGUOUS"="logical",
            "CALL_RATE"="numeric",
            "IMPUTED"="logical",
            "OTHER_ALLELE_REF"="character",
            "ORI_OTHER_ALLELE"="character",
            "ORI_OTHER_ALLELE_REF"="character",
            "ORI_EFFECT_ALLELE"="character",
            "ORI_A0"="character",
            "REF_RSID"="character")


  t<-c(t,types[name])

  func = list("Z"=as.numeric,
              "SNP"=as.character,
              "SIGNED_SUMSTAT"=as.numeric,
              "SE"=as.numeric,
              "P"=as.numeric,
              "LOG10_P"=as.numeric,
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
              "BETA_95L"=as.numeric,
              "BETA_95U"=as.numeric,
              "Q_STATISTIC"=as.numeric,
              "Q_P_VALUE"=as.numeric,
              "AC"=as.integer,
              "EFFECT_ALLELE"=as.character,
              "OTHER_ALLELE"=as.character,
              "A0"=as.character,
              "MARKER_TYPE"=as.character,
              "EUR_FRQ"=as.numeric,
              "cptid"=as.character,
              "AMBIGUOUS"=as.logical,
              "CALL_RATE"=as.numeric,
              "IMPUTED"=as.logical,
              "OTHER_ALLELE_REF"=as.character,
              "ORI_OTHER_ALLELE"=as.character,
              "ORI_EFFECT_ALLELE"=as.character,
              "ORI_OTHER_ALLELE_REF"=as.character,
              "ORI_A0"=as.character,
              "REF_RSID"=as.character)
  f<-c(f,func[name])

  if(name %in%  c("cptid","SNP","SE","P","N_CAS","N","INFO","FRQ","STRAND","CHR","BP","BETA","EFFECT_ALLELE","OTHER_ALLELE","IMPUTED")) {
    logi <- TRUE
  } else {
    logi <- FALSE
  }
  names(logi) <- name
  l <- c(l,logi)

}

base_column_mapping = ColMap(active = l, aliases = a, col_types = t, funcs = f)

usethis::use_data(base_column_mapping, internal = F, overwrite=TRUE)

