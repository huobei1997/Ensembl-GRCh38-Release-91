library(AnnotationDbi)
library(org.Hs.eg.db)
library(readxl)
library(writexl)
library(dplyr)
library(clusterProfiler)
library(tidyr)
library(biomaRt)
library(AnnotationHub)

##########################01-read gtf_id #58302####################################
ID1 <- read.table("ENSMBL_hg38_release_91.txt",header = T)

#############02-AnnotationDbi-mapIds #35879##############################
ID1$gene_id = mapIds(org.Hs.eg.db, 
                 keys = ID1$gene_name, 
                 column = "ENTREZID", 
                 keytype = "ENSEMBL")

ID2<-ID1%>%dplyr::filter(!is.na(ID1$gene_id))
ID3<-ID1%>%dplyr::filter(is.na(ID1$gene_id))#22423

####################03-ClusterProfiler-bitr #0##############################
g <- setdiff(ID1$gene_name, ID2$gene_name)
valid_symbols <- keys(org.Hs.eg.db, keytype = "ENSEMBL")
g <- intersect(valid_symbols,g)

##################################04-syngo#707##########################
#https://www.syngoportal.org/convert

syngo<-read_xlsx("idmap.xlsx")
syngo<-syngo[,c(1,2,4)]
ID4<-syngo%>%dplyr::filter(!is.na(syngo$entrezgene))
ID5<-syngo%>%dplyr::filter(is.na(syngo$entrezgene))

##########################05-AnnotationDbi-select#0########################
mart <- useMart("ensembl")
dataset <- useDataset("hsapiens_gene_ensembl", mart)
gene_ids <- getBM(attributes = c("ensembl_gene_id","entrezgene_id"),
                  mart = dataset)
gene_ids<-gene_ids%>%dplyr::filter(!is.na(gene_ids$entrezgene_id))

g<-intersect(gene_ids$ensembl_gene_id,ID5$query)

#####################06-EnsDb.Hsapiens.v91#793##############################
#https://support.bioconductor.org/p/9150642/

ah <- AnnotationHub()
query(ah, "EnsDb.Hsapiens.v91")
edb <- ah[["AH60773"]]
prots <- genes(edb)

ID6<- prots %>% data.frame()
ID6$gene_id_2 <- as.numeric(as.character(ID6$entrezid))
ID6<-ID6%>%dplyr::filter(!is.na(ID6$gene_id_2))
ID6<-ID6[,c(6,14)]
colnames(ID6)<-c("query","gene_id")

ID7<-left_join(ID5,ID6,by="query")
ID8<-ID7%>%dplyr::filter(!is.na(ID7$gene_id))
ID9<-ID7%>%dplyr::filter(is.na(ID7$gene_id))

################################07-ensembl-biomart#975##############################
mart<-read.table("mart_export.txt",sep="\t",quote = "",fill = T,header = T)
mart<- mart[!duplicated(mart$Gene.stable.ID),]
mart<-mart[,c(1,7)]
mart<-mart%>%dplyr::filter(!is.na(mart$NCBI.gene..formerly.Entrezgene..ID))
colnames(mart)<-c("query","gene_id")

ID10<-left_join(ID9,mart,by="query")
ID11<-ID10%>%dplyr::filter(!is.na(ID10$gene_id.y))
ID12<-ID10%>%dplyr::filter(is.na(ID10$gene_id.y))

################################08-add type and description##############################
id<-rbind(ID2,ID4,ID8,ID11,ID12)
id<-id[,c(1,3)]
colnames(id)<-c("gene_id","id")
##
ah <- AnnotationHub()
edb <- ah[["AH60773"]]
prots <- genes(edb)
ID<- prots %>% data.frame()
ID<-ID[,c(6,8,10,12)]
##
res<-left_join(id,ID,by="gene_id")
colnames(res)<-c("ENSEMBL","ENTREZID","gene_biotype","description","SYMBOL")
writexl::write_xlsx(res,"ID.xlsx")




