
load_input_file <- function (file)
{
  skipl <- length(grep("^#", readLines(file)))
  
  result <- read.delim(file, quote = "",  skip=skipl, stringsAsFactors=FALSE)
  #result <- join(subset(cebs_df, select=-CAS), subset(mapping, select=c(CAS, Tox21.ID, Chemical.Name,StructureID)), by="Tox21.ID",  match = "first")
  #result$readout <- paste(readout, substr(result$Library, nchar(result$Library), nchar(result$Library) ), sep=".")
  #result$pathway <- pathway
  #result$uniqueID <- paste(result$Tox21.ID, "@", result$Cmpd_Library, sep="")
  #result <- result[order(result$uniqueID, result$readout),]
  #replicate <- names(table(result$uniqueID)[table(result$uniqueID) > 3])
  #times <- as.vector(t(table(result[result$uniqueID %in% replicate, c("uniqueID", "readout")])))
  #result[result$uniqueID %in% replicate,]$readout <- paste(result[result$uniqueID %in% replicate,]$readout,  "_", unlist(sapply(times, seq)), sep="")
  if (is.null(result$parent)) 
  {
    result$parent <- ''
  } else
  {
    result$parent[is.na(result$parent)] <- ''
  }
  if (is.null(result$mask))
  {
    result[, "mask"] <- ''
  } else
  {
    result[is.na(result$mask), "mask"] <- ''
  }
  return(result)
  
}