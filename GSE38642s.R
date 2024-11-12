cel_data <- ReadAffy(celfile.path = "C:/Users/adity/Documents/RScripts/GSE38642")
pData(cel_data)
proc_data <- getGEO('GSE38642')
proc_data <- proc_data[[1]]
proc_data$supplementary_file
pd <- pData(proc_data)
pd['cel_file'] <- str_split(pd$supplementary_file,"/") %>% map_chr(tail,1)
celdata <- read.celfiles(paste0('new/',pd$cel_file),phenoData=phenoData(proc_data))
print(pd$cel_file)

list.files("C:/Users/adity/Documents/RScripts/new", full.names = TRUE)
cel_files <- list.files("C:/Users/adity/Documents/RScripts/new", pattern = "*.CEL.gz", full.names = TRUE)

celdata <- read.celfiles(filenames = cel_files, phenoData = phenoData(proc_data))
pData(celdata)
varLabels(celdata)
varLabels(celdata)
phenotype_data <- pData(celdata)
exprs_celdata <- exprs(celdata)
disease_state <- phenotype_data[["status:ch1"]]
log_celdata <- log2(exprs(celdata) + 1)
pca_celdata <- prcomp(t(log_celdata), scale. =TRUE)
pca_celdata_df <- data.frame(PC1 = pca_celdata$x[, 1], PC2 = pca_celdata$x[, 2], Sample = colnames(celdata), DiseaseState = disease_state)
pca_pre_plot <- ggplot(pca_celdata_df, aes(x = PC1, y = PC2, color = DiseaseState)) +
  geom_point(size = 3) +  
  ggtitle("PCA of Pre-Normalized Data Colored by Sample") +
  xlab(paste("PC1 -", round(summary(pca_celdata)$importance[2, 1] * 100, 1), "% variance")) +
  ylab(paste("PC2 -", round(summary(pca_celdata)$importance[2, 2] * 100, 1), "% variance")) +
  theme_minimal()+
  scale_color_viridis_d()
pca_pre_plot
summary(pca_celdata)


#RMA Normalization of the data
rma_norm_data <- rma(celdata)
phenotype_rma <- pData(rma_norm_data)
exprs_rma <- exprs(rma_norm_data)
#PCA of rma normalized data
pca_norm <- prcomp(t(exprs_rma), scale. = TRUE)
disease_state_n <- phenotype_rma[['status:ch1']]
pca_norm_df <- data.frame(PC1 = pca_norm$x[, 1], PC2 = pca_norm$x[, 2], Sample = colnames(rma_norm_data), DiseaseState = disease_state_n)
pca_rma_plot <- ggplot(pca_norm_df, aes(x = PC1, y = PC2, color = DiseaseState)) +
  geom_point(size = 3) +  
  ggtitle("PCA of Normalized Data Colored by Sample") +
  xlab(paste("PC1 -", round(summary(pca_norm)$importance[2, 1] * 100, 1), "% variance")) +
  ylab(paste("PC2 -", round(summary(pca_norm)$importance[2, 2] * 100, 1), "% variance")) +
  theme_minimal()
pca_rma_plot

#annotation of the genes
probeIDs <- rownames(exprs_rma)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")

library("biomaRt")
mart <- biomaRt::useEnsembl(biomart =  "ensembl",
                       dataset= "hsapiens_gene_ensembl",
                       )

Probeset_gene_id <- biomaRt::getBM(attributes = c("hgnc_symbol", "ensembl_gene_id",
                                                 "affy_hugene_1_0_st_v1"),
                                  filters = "affy_hugene_1_0_st_v1",
                                  values = probeIDs,
                                  mart = mart)

match_rows <- match(rownames(exprs_rma), Probeset_gene_id$affy_hugene_1_0_st_v1)     
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
contrast_matrix <- makeContrasts(ControlVsT2D = non.diabetic.donors - T2D.donors,
                                 levels = design)

fit1 <- lmFit(exprs_anno, design)
fit2 <- contrasts.fit(fit1, contrast_matrix)
fit2 <- eBayes(fit2)

top_genes <- topTable(fit2, adjust = "fdr", number = Inf)
head(top_genes)

summary(decideTests(fit2, lfc=0.5, p.value=0.05))

deg_significant <- topTable(fit2, number = Inf, p.value = 0.05, lfc = 0.5)

top_genes1 <- topTable(fit2, adjust = "none", number = Inf)

