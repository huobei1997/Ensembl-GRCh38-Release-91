library(dplyr)
library(tidyr)
library(AnnotationHub)

##########################01-extract gtf ENSEMBL and SYMBOL #58302###################
ID1 <- read.table("id.txt", header = FALSE, sep = "\t", quote = "", fill = TRUE)
ID1 <- ID1[!duplicated(ID1$V1), c(1, 4)]
colnames(ID1) <- c("gene_name", "gene_symbol")
ID1[] <- lapply(ID1, function(x) sub(".*: ", "", x))

##########################02-extract EnsDb.Hsapiens.v91##############################
#https://support.bioconductor.org/p/9150642/

ah <- AnnotationHub()
query(ah, "EnsDb.Hsapiens.v91")
edb <- ah[["AH60773"]]
prots <- genes(edb)

ID2<- prots %>% data.frame()
ID2<-ID2[,c(6,8,10,13)]
colnames(ID2)<-c("gene_name","gene_biotype","description","gene_id")
ID2<-left_join(ID1,ID2,by="gene_name")

#########################03-check data##############################################
ID3<-ID2%>%dplyr::filter(!is.na(ID2$gene_id)) #25282
ID3$gene_id_2 <- as.numeric(as.character(ID3$gene_id))
na_values <- ID3[is.na(ID3$gene_id_2), ] #288

ID3$gene_id<-ID3$gene_id_2
ID3<-ID3[,c(1:5)]
ID3<- ID3[!is.na(ID3$gene_id),]

##https://www.syngoportal.org/convert和https://www.genecards.org/

ID4<-readxl::read_xlsx("ID_ENTREZID.xlsx")
ID4<-rbind(ID3,ID4)
ID5<-readxl::read_xlsx("ID_no_ENTREZID.xlsx")
ID<-rbind(ID4,ID5)           
writexl::write_xlsx(ID,"ID.xlsx")








