StudyCorpus <- setClass(
  Class = "StudyCorpus",
  slots = list(
    studies = "list"
  ),
  prototype = list(
    studies = list()
  )
)

corpus <- StudyCorpus()
