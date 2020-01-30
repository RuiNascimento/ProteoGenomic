library(httr)
library(data.table)
library(stringr)


UniprotConvertion <- function(query, from, to, format = 'tab') {
  ## Convert identifiers with the help of UNIPROT
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
  ## https://www.uniprot.org/help/uniprotkb_column_names to check the possible column names
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



getGI <- function(s) {
  
  ## Helper function to extract GI from fasta headers in MAxQuant results
  ex <- str_extract_all(string = s, pattern = "gi\\|[:digit:]+\\|")[[1]]
  str_extract_all(string = ex, pattern = "[:digit:]+",simplify = T)[,1]
}

get1GI <- function(s) {
  
  ## Helper function to extract GI from fasta headers in MaxQuant results
  ex <- str_extract_all(string = s, pattern = "gi\\|[:digit:]+\\|")[[1]][[1]]
  str_extract(string = ex, pattern = "[:digit:]+")
}

s <- "gi|359491770|ref|XP_003634320.1|;gi|225439342|ref|XP_002270170.1|;gi|225437695|ref|XP_002279878.1|;gi|147853527|emb|CAN80663.1|;gi|147836392|emb|CAN75420.1|;gi|147834215|emb|CAN70884.1|;gi|147767872|emb|CAN71283.1|;gi|302143958|emb|CBI23063.3|;gi|297740140|emb|CBI30322.3|;gi|225455266|ref|XP_002273568.1|;gi|147769165|emb|CAN60770.1|;gi|147802377|emb|CAN77120.1|;gi|297736166|emb|CBI24204.3|;gi|225465030|ref|XP_002265864.1|;gi|297735960|emb|CBI23934.3|;gi|302143596|emb|CBI22349.3|;gi|731432763|ref|XP_010644407.1|;gi|225465609|ref|XP_002266370.1|;gi|731432761|ref|XP_010644406.1|;gi|731432759|ref|XP_010644405.1|;gi|731432757|ref|XP_010644404.1|;gi|147805226|emb|CAN64480.1|;gi|225465615|ref|XP_002267017.1|;gi|147834511|emb|CAN71997.1|;gi|302143601|emb|CBI22354.3|;gi|731417327|ref|XP_002267206.2|;gi|297744041|emb|CBI37011.3|;gi|76559884|tpe|CAI56329.1|;gi|297742776|emb|CBI35456.3|"

s2 <- gsub("gi", "\ngi", s)


res.6hpi$GI <- lapply(res.6hpi$protein, get1GI)

cat(s2)
library(data.table)
sdt <- fread(text = s2)

sdt[, c("V1", "V5"):=NULL]

dcast(sdt, (V2~V3))

getGI(s)

query <- getGI(s)
converted <- UniprotConvertion(query = query, from = "P_GI", "ACC")
head(converted)
getUniprotGoodies(converted$To)

spl <- str_split(s, ";")
spl
q <- str_split_fixed(spl[[1]], "\\|",n = 5)



#q <- str_split_fixed(t[[1]], "\\|",n = 3)
head(q)
df <- data.frame("gi" = q[,2], "ref" = q[,3], "id" = q[,4])

head(df)

