Study <- setClass(
  Class = "Study",
  contains = "VIRTUAL",
  slots = list(
    pmid = "integer"
  ),
  prototype = list(
    pmid = NA_integer_
  )
)

GWAS <- setClass(
  Class = "GWAS",
  contains = "Study",
  slots = list(
  ),
  prototype = list(
  )
)

s <- GWAS()
