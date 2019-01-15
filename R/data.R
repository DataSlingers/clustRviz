#' Log transformed word count of presidential speeches
#'
#' A dataset of the top 75 most variable log-transformed word counts for
#' each US president aggregated over several speeches
#' (Inaugural, State of the Union, etc.).
#' Stop words have been removed and words have been stemmed.
#'
#' @format A data.frame with 44 rows (one for each president) and 75 columns (log transformed word counts)
#' @details Grover Cleveland was elected president twice (1892 and 1884). For our purposes his speeches are combined.
#' @source \url{http://www.presidency.ucsb.edu}
"presidential_speech"

#' Word Count Data from Four English-Language Authors
#' 
#' This data set (\eqn{n=841, p = 69}) consists of counts of common words
#' appearing in texts written by four popular English-language authors 
#' (Jane Austen, Jack London, William Shakespeare, and John Milton).
#' The row names are the authors (true cluster labels) and the column
#' names are the words (slightly processed). 
"authors"

#' Log-Transformed Level III RPKM Gene Expression Levels for 438 Breast-Cancer Patients
#' 
#' This data set (\eqn{n = 438, p = 353}) contains log-transformed Level III RPKM gene 
#' expression levels for 438 breast-cancer patients collected by the Cancer Genome Atlas
#' Network. The Luminal A and Luminal B subtypes have been combined. The row names give the
#' clinically diagnosed subtype (true cluster labels) and the column names are the gene IDs.
#' 
#' @references
#' The Cancer Genome Atlas Network. "Comprehensive Molecular Portraits of Human Breast Tumours"
#' Nature 490, p.61-70. 2012. \doi{10.1038/nature11412}
"tcga_breast"
