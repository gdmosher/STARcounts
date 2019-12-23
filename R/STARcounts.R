#' STAR counts
#'
#' This function loads STAR Aligner count files *ReadsPerGene.out.tab$ from results_path, accumulates the
#' data into a data frame, then writes the results to table_name.
#' Some code from biostars 241602.
#' Some code from stackoverflow 13762224.
#' annotation.gtf must already exist before the mapping step, ./gffread -T annotation.gff -o annotation.gtf
#' _SampleName_.ReadsPerGene.out.tab must already exist (created during mapping step)
#' @param results_path Path to the STAR Aligner count files
#' @param table_name Name of the output file including the .xls extension
#' @param dbg Optional debug level defaults to 0
#'
#' @return
#' @export
#'
#' @examples
#' write.table_STARcounts(results_path="./results", table_name="STARcountDFeByg.xls")
#' countDF <- read.delim("results/STARcountDFeByg.xls", row.names=1, check.names=FALSE)
#' head(countDF)
write.table_STARcounts <- function(results_path="./results", table_name="STARcountDFeByg.xls", dbg=0) {
  details <- file.info(list.files( path = results_path, pattern = "*ReadsPerGene.out.tab$", full.names = TRUE, recursive = TRUE ))
  details <-  details[with(details, order(as.POSIXct(ctime))), ]
  ff <- rownames(details) # from stackoverflow 13762224, filenames sorted by created time
  if (length(ff) == 0) {
    stop("No STAR count files found in results_path.")
  }
  counts.files <- lapply( ff, read.table, skip = 4 )
  ## column 1: gene ID
  ## column 2: counts for unstranded RNA-seq
  ## column 3: counts for the 1st read strand aligned with RNA (htseq-count option -s yes)
  ## column 4: counts for the 2nd read strand aligned with RNA (htseq-count option -s reverse)
  counts <- as.data.frame( sapply( counts.files, function(x) x[ , 'V2' ] ) ) # V2 is column 2: counts for unstranded RNA-seq
  ## prepare colnames for xls file
  ff <- gsub( "[.]ReadsPerGene[.]out[.]tab", "", ff ) # get rid of all of the filename except the SampleName
  ff <- gsub( "[.]/results/", "", ff )
  .last <- function(x) { return( x[length(x)] ) } # or use tail(x,1) since ns slower is ok
  xxx <- lapply( ff, function(x) { .last(strsplit(x, split="\\/")[[1]]) } ) # keep all that follows last slash, the SampleName
  ff <- unlist(xxx)
  colnames(counts) <- ff
  row.names(counts) <- counts.files[[1]]$V1 # V1 are the gene ID's, see STAR Manual --quantMode
  write.table(counts, paste0(results_path, "/", table_name), col.names=NA, quote=FALSE, sep="\t")
}


## Usage:
# write.table_STARcounts(results_path="./results", table_name="STARcountDFeByg.xls")

