
load_input_file <- function (file)
{
  skipl <- length(grep("^#", readLines(file)))
  
  result <- read.delim(file, quote = "",  skip=skipl, stringsAsFactors=FALSE)
  
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
  
  result <- result[,which(colSums(is.na(result)) != nrow(result))]
  
  return(result)
  
}