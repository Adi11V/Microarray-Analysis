download_dir <- "C:/Users/adity/Documents/RScripts/GSE25724/"
if (!dir.exists(download_dir)) {
  dir.create(download_dir, recursive = TRUE)}
getGEOSuppFiles("GSE25724", baseDir = download_dir)
untar("GSE25724_RAW.tar", exdir = "C:/Users/adity/Documents/RScripts/GSE25724")

#READING THE CEL DATA AND ATTACHING THE EXPERIMENTAL DATA (PDATA)

cel_data <- ReadAffy(celfile.path = "C:/Users/adity/Documents/RScripts/GSE25724")
pData(cel_data)

proc_data <- getGEO('GSE25724') #For the phenotype data
proc_data
length(proc_data)
class(proc_data[[1]])
proc_data <- proc_data[[1]]
varLabels(proc_data)
proc_data$supplementary_file
pd <- pData(proc_data)
pd['cel_file'] <- str_split(pd$supplementary_file,"/") %>% map_chr(tail,1)

list.files("C:/Users/adity/Documents/RScripts/GSE25724", full.names = TRUE)
cel_files <- list.files("C:/Users/adity/Documents/RScripts/GSE25724", pattern = "*.CEL.gz", full.names = TRUE)

celdata <- read.celfiles(filenames = cel_files, phenoData = phenoData(proc_data))
pData(celdata)
varLabels(celdata)
varLabels(celdata)
phenotype_data <- pData(celdata)
exprs_celdata <- exprs(celdata)

disease_state <- phenotype_data[["disease state:ch1"]]

#PCA Pre normalization

log_celdata <- log2(exprs(celdata) + 1)
pca_celdata <- prcomp(t(log_celdata), scale. =TRUE)
pca_celdata_df <- data.frame(PC1 = pca_celdata$x[, 1], PC2 = pca_celdata$x[, 2], Sample = colnames(celdata), DiseaseState = disease_state)
pca_pre_plot <- ggplot(pca_celdata_df, aes(x = PC1, y = PC2, color = DiseaseState)) +
  geom_point(size = 3) +  
  ggtitle("PCA of Pre-Normalized Data Colored by Disease") +
  xlab(paste("PC1 -", round(summary(pca_celdata)$importance[2, 1] * 100, 1), "% variance")) +
  ylab(paste("PC2 -", round(summary(pca_celdata)$importance[2, 2] * 100, 1), "% variance")) +
  scale_color_manual(values = c("type 2 diabetes" = "red", "non-diabetic" = "blue")) +
  theme_minimal()
pca_pre_plot
summary(pca_celdata)

 boxplot((exprs_celdata), main = "Raw Data", outline = FALSE)

#RMA Normalization of the data
rma_norm_data <- rma(celdata)
phenotype_rma <- pData(rma_norm_data)
exprs_rma <- exprs(rma_norm_data)

boxplot_rma <- boxplot((exprs_rma), main = "Normalized Data", outline = FALSE)

#PCA of rma normalized data
pca_norm <- prcomp(t(exprs_rma), scale. = TRUE)
disease_state_n <- phenotype_rma[['disease state:ch1']]
summary(pca_norm)
pca_norm_df <- data.frame(PC1 = pca_norm$x[, 1], PC2 = pca_norm$x[, 2], Sample = colnames(rma_norm_data), DiseaseState = disease_state_n)
pca_rma_plot <- ggplot(pca_norm_df, aes(x = PC1, y = PC2, color = DiseaseState)) +
  geom_point(size = 3) +  
  ggtitle("PCA of Normalized Data Colored by Sample") +
  xlab(paste("PC1 -", round(summary(pca_norm)$importance[2, 1] * 100, 1), "% variance")) +
  ylab(paste("PC2 -", round(summary(pca_norm)$importance[2, 2] * 100, 1), "% variance")) +
  scale_color_manual(values = c("type 2 diabetes" = "red", "non-diabetic" = "blue")) +
  theme_minimal()
pca_rma_plot

library("biomaRt")
mart <- biomaRt::useEnsembl(biomart =  "ensembl",
                            dataset= "hsapiens_gene_ensembl",
)
attributes <- listAttributes(mart) 

probeIDs <- rownames(exprs_rma)
Probeset_gene_id <- biomaRt::getBM(attributes = c("hgnc_symbol", "ensembl_gene_id",
                                                  "affy_hg_u133a_2"),
                                   filters = "affy_hg_u133a_2",
                                   values = probeIDs,
                                   mart = mart)

match_rows <- match(rownames(exprs_rma), Probeset_gene_id$affy_hg_u133a_2)     
matched_rows <- Probeset_gene_id[match_rows, ]

exprs_anno <- exprs_rma
rownames(exprs_anno)<- matched_rows$hgnc_symbol
nrow(exprs_anno)

exprs_anno <- avereps(x=exprs_anno, ID=rownames(exprs_anno))
sum(is.na(rownames(exprs_anno)))
sum(is.na(exprs_anno))

exprs_anno <- exprs_anno[!is.na(rownames(exprs_anno )), ]
exprs_anno <- na.omit(exprs_anno)
exprs_anno <- exprs_anno[rownames(exprs_anno)!= "",]


#DEG identification

library("limma")
design <- model.matrix(~ 0 + disease_state_n)
colnames(design) <- levels(as.factor(disease_state_n))  
design
colnames(design) <- make.names(colnames(design))
contrast_matrix <- makeContrasts(T2DvsControl = type.2.diabetes -  non.diabetic,
                                 levels = design)

fit1 <- lmFit(exprs_anno, design)
fit2 <- contrasts.fit(fit1, contrast_matrix)
fit2 <- eBayes(fit2)

top_genes <- topTable(fit2, adjust = "fdr", number = Inf)
head(top_genes)

summary(decideTests(fit2, lfc=1, p.value=0.05))

deg_significant <- topTable(fit2, number = Inf, p.value = 0.05, lfc = 1)

top_genes1 <- topTable(fit2, adjust = "none", number = Inf)

volcanoplot(fit2, coeff=2) #rma data
interesting_genes1 <- topTable(fit2,number=Inf,p.value = 0.05,lfc > 1)
points(interesting_genes1[['logFC']],-log10(interesting_genes1[['P.Value']]),col='red')





##END##
##FOR REF##




#gcrma normalization
if (!requireNamespace("gcrma", quietly = TRUE)) {
  install.packages("gcrma", repos = "http://bioconductor.org/packages/release/bioc", type = "source")
}
library(gcrma)
gcrma_data <- gcrma(celdataA) 
exprs_gcrma <- exprs(gcrma_data)
phenotype_gcrma <- pData(gcrma_data)

#pca of gcrma data

pca_gcrma<- prcomp(t(exprs_gcrma))
disease_state_gc <- phenotype_gcrma[['disease state:ch1']]
summary(pca_gcrma)
pca_gcrma_df <- data.frame(PC1 = pca_gcrma$x[, 1], PC2 = pca_gcrma$x[, 2], Sample = colnames(gcrma_data), DiseaseState = disease_state_gc)
pca_gcrma_plot <- ggplot(pca_gcrma_df, aes(x = PC1, y = PC2, color = DiseaseState)) +
  geom_point(size = 3) +  
  ggtitle("PCA of GCRMA Normalized Data Colored by Sample") +
  xlab(paste("PC1 -", round(summary(pca_gcrma)$importance[2, 1] * 100, 1), "% variance")) +
  ylab(paste("PC2 -", round(summary(pca_gcrma)$importance[2, 2] * 100, 1), "% variance")) +
  theme_minimal()
pca_gcrma_plot

#DEG
library(limma)
design <- model.matrix(~ 0 + disease_state_n)  #designing the matrix
colnames(design) <- levels(as.factor(disease_state_n))  
design
colnames(design) <- make.names(colnames(design))
contrast_matrix <- makeContrasts(T2DVsNormal = type.2.diabetes - non.diabetic,
                                 levels = design)
#Fitting limma
fit <- lmFit(exprs_rma, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)
#Upregulated
ps2 <- topTable(fit2,number=Inf,p.value = 0.05,lfc=1)
ps2_up <- rownames(ps2[ps2$logFC > 0,])
df <- AnnotationDbi::select(hgu133a.db,ps2_up,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")
mapping_summary <- table(df$PROBEID)
head(mapping_summary)
up_unique <- df %>%
  distinct(PROBEID, .keep_all = TRUE) #for duplicates removal
up_filtered <- up_unique %>%
  filter(!is.na(SYMBOL) & !is.na(GENENAME)) #removing non annotated data
head(up_filtered)
#downregulated
ps2_dn <- rownames(ps2[ps2$logFC < 0,])
df_dn <- AnnotationDbi::select(hgu133a.db,ps2_dn,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")
mapping_summary <- table(df_dn$PROBEID)
head(mapping_summary)
dn_unique <- df_dn %>%
  distinct(PROBEID, .keep_all = TRUE)
dn_filtered <- dn_unique %>%
  filter(!is.na(SYMBOL) & !is.na(GENENAME))
head(dn_filtered)
#DEG for Gcrma data
design1 <- model.matrix(~ 0 + disease_state_gc)  
colnames(design1) <- levels(as.factor(disease_state_gc))  
design1
contrast_matrix1 <- makeContrasts(T2DVsNormal = type.2.diabetes - non.diabetic,
                                 levels = design1)
fit_g <- lmFit(exprs_gcrma, design1)
colnames(design1) <- make.names(colnames(design1))

fit_g2 <- contrasts.fit(fit_g, contrast_matrix1)
fit_g2 <- eBayes(fit_g2)
summary(decideTests(fit_g2, lfc=1, p.value=0.05)) #summary of degs gcrma
summary(decideTests(fit2, lfc=1, p.value=0.05)) #summary of degs rma

volcanoplot(fit_g2, coeff=2) #gcrma 
interesting_genes <- topTable(fit_g2,number=Inf,p.value = 0.05,lfc=1)
volcanoplot(fit_g2, coef=2, main=sprintf("%d features pass our cutoffs",nrow(interesting_genes)))
points(interesting_genes[['logFC']],-log10(interesting_genes[['P.Value']]),col='red')
ps_g <- topTable(fit_g2,number=Inf,p.value = 0.05,lfc=1.5)

volcanoplot(fit2, coeff=2) #rma data
interesting_genes1 <- topTable(fit2,number=Inf,p.value = 0.05,lfc > 1)
points(interesting_genes1[['logFC']],-log10(interesting_genes1[['P.Value']]),col='red')
