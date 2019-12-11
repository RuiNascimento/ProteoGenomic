library(httr)
UniprotConvertion <- function(query, from, to, format = 'tab') {
  if (length(query)>1) {query <- paste0(query, collapse = " ")}
  b2 <- "https://www.uniprot.org/uploadlists/"
  response <- POST(b2, body = list('from' = from, 'to' = to, 'format' = format,
                                   'query' = query))
  read.delim(text = content(response, "text", encoding = "UTF-8"))
}

UniprotConvertion(c("296083299", "225431846", "147773334", "225431854", "100241610"), "P_GI", "ID")


getUniprotGoodies <- function(query, columns = c("id","entry name","reviewed","protein names","genes","organism","length"))
{
  ## Set columns with some defaults
  ## query and columns start as a character vectors
  qstring <- paste(query, collapse="+or+")
  cstring <- paste(columns, collapse=",")
  uri <- 'http://www.uniprot.org/uniprot/?query='
  fullUri <- paste0(uri,qstring,'&format=tab&columns=',cstring)
  dat <- read.delim(URLencode(fullUri), stringsAsFactors=FALSE)
  ## now remove things that were not in the specific original query...
  dat <- dat[dat[,1] %in% query,]
  dat
}

getUniprotGoodies(c("A5AQ75","F6H6U5"), columns = c("id","entry name","reviewed","protein names","genes","organism","length"))
