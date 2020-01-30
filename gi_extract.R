library(org.Vvinifera.eg.db)
library(AnnotationDbi)
library(stringr)
library(topGO)
api_key <- readLines("~/.ncbi_api")

head(data_results)

head(data_results$name)

XP <- data_results$name[1]
NO_XP <- data_results$name[3]

# Criar uma coluna com os GI, de preferencia o referenciado ao XP

extract_GI <- function(x){
 # Get 1 GI, priority to GI with XP annotation, otherwise get the first GI 
  if (str_detect(x, "XP")) {
    gi <- str_extract(string = x, pattern = "(?<=gi\\|)[:number:]+(?=\\|ref\\|XP)")
  } else {
      gi <- str_extract(string = x, pattern = "(?<=gi\\|)[:number:]+")
    }
  return(gi)
  }

data_results$GI <- sapply(data_results$name, extract_GI, simplify = TRUE)

# Obter GI do gene correspondente

protein2gene <- function(x, api_key = api_key){
  gene_gi <- entrez_link(dbfrom="protein", db="gene", id = x, api_key = api_key)$links$protein_gene[1]
  if (is.null(gene_gi)){gene_gi <- ""}
  return(gene_gi[[1]])
}

data_results$GENE_GI <- sapply(data_results$GI, protein2gene, api_key = api_key, simplify = TRUE, USE.NAMES = FALSE)


# Para obter os GO ID

df <- select(org.Vvinifera.eg.db, keys = data_results$GENE_GI, column = c("go_idall", "ontologyall"), keytype = "gene_id")

BP <- df[df$ontologyall == "BP",]

select(GO.db, keys = BP$go_idall, columns = "TERM", keytype = "GOID")


library(topGO)

# 6h

geneList6h <- data_results$Ri6h_vs_Rm6h_p.val
names(geneList6h) <- data_results$GENE_GI
head(geneList6h)


GOdata6h <- new("topGOdata",
              ontology = "BP",
              allGenes = geneList6h,
              geneSelectionFun = function(x) (x < 0.05),
              annot = annFUN.org , mapping = "org.Vvinifera.eg.db",
              ID = "entrez")

resultKS <- runTest(GOdata6h, algorithm = "weight01", statistic = "fisher")

tab <- GenTable(GOdata6h, raw.p.value = resultKS, topNodes = length(resultKS@score), numChar = 120)
# Convert p.value to numeric
tab$raw.p.value <- as.numeric(tab$raw.p.value)

head(tab, 15)

par(cex = 0.6)
showSigOfNodes(GOdata6h, score(resultKS), firstSigNodes = 5, useInfo = "def")
par(cex = 1)

# 12h

geneList12h <- data_results$Ri12h_vs_Rm12h_p.val
names(geneList12h) <- data_results$GENE_GI
head(geneList12h)

GOdata12h <- new("topGOdata",
                ontology = "BP",
                allGenes = geneList12h,
                geneSelectionFun = function(x) (x < 0.05),
                annot = annFUN.org , mapping = "org.Vvinifera.eg.db",
                ID = "entrez")

resultKS12h <- runTest(GOdata12h, algorithm = "weight01", statistic = "fisher")

tab12h <- GenTable(GOdata12h, raw.p.value = resultKS12h, topNodes = length(resultKS@score), numChar = 120)
# Convert p.value to numeric
tab12h$raw.p.value <- as.numeric(tab12h$raw.p.value)

head(tab12h, 15)

par(cex = 0.6)
showSigOfNodes(GOdata12h, score(resultKS12h), firstSigNodes = 5, useInfo = "def")
par(cex = 1)

# 24h

geneList24h <- data_results$Ri24h_vs_Rm24h_p.val
names(geneList24h) <- data_results$GENE_GI
head(geneList24h)

GOdata24h <- new("topGOdata",
                 ontology = "BP",
                 allGenes = geneList24h,
                 geneSelectionFun = function(x) (x < 0.05),
                 annot = annFUN.org , mapping = "org.Vvinifera.eg.db",
                 ID = "entrez")

resultKS24h <- runTest(GOdata24h, algorithm = "weight01", statistic = "fisher")

tab24h <- GenTable(GOdata24h, raw.p.value = resultKS24h, topNodes = length(resultKS@score), numChar = 120)
# Convert p.value to numeric
tab24h$raw.p.value <- as.numeric(tab24h$raw.p.value)

head(tab24h, 15)

par(cex = 0.6)
showSigOfNodes(GOdata24h, score(resultKS24h), firstSigNodes = 5, useInfo = "def")
par(cex = 1)
