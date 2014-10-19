
load_input_file <- function (file)
{
  skipl <- length(grep("^#", readLines(file)))
  
  result <- read.delim(file, quote = "",  skip=skipl, stringsAsFactors=FALSE)
  result <- check_data_input(result)
#   if (is.null(result$parent)) 
#   {
#     result$parent <- ''
#   } else
#   {
#     result$parent[is.na(result$parent)] <- ''
#   }
#   if (is.null(result$mask))
#   {
#     result[, "mask"] <- ''
#   } else
#   {
#     result[is.na(result$mask), "mask"] <- ''
#   }
#   
#   result <- result[,which(colSums(is.na(result)) != nrow(result))]
#   
#   return(result)
  
}

check_data_input <- function (df)
{
  if (is.null(df$parent)) 
  {
    df$parent <- ''
  } else
  {
    df$parent[is.na(df$parent)] <- ''
  }
  if (is.null(df$mask))
  {
    df[, "mask"] <- ''
  } else
  {
    df[is.na(df$mask), "mask"] <- ''
  }
  
  df <- df[,which(colSums(is.na(df)) != nrow(df))]
  
  return(df)
}