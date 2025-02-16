T2dRowName <- function(db_path,tableName){
  library(DBI)
  library(RSQLite)
  con <- dbConnect(SQLite(), dbname = db_path)
  query <- paste0("SELECT * FROM ", tableName)
  data <- dbGetQuery(con, query)
  res <- data
  return(rownames(res));
}

T2dResCol <- function(db_path,tableName,colInx){
  library(DBI)
  library(RSQLite)
  con <- dbConnect(SQLite(), dbname = db_path)
  query <- paste0("SELECT * FROM ", tableName)
  data <- dbGetQuery(con, query)
  res <- data[[colInx]];
  hit.inx <- is.na(res) | res == ""; # note, must use | for element-wise operation
  res[hit.inx] <- "N/A";
  return(res);
}

