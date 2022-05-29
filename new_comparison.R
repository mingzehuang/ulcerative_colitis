library(Rtsne)
library(latentcor)
source("~/ulcerative_colitis/ulcerative_colitis/analysis.r")
train.Fib.seur <- readRDS("~/ulcerative_colitis/ulcerative_colitis/train.Fib.seur.rds")
train.Fib.updated = UpdateSeuratObject(object = train.Fib.seur)

###################################################################################

# Calculate tSNE
print('Calculating tSNE')
tsne.rot = Rtsne(train.Fib.updated@reductions[['pca']]@cell.embeddings %*% diag(1/((train.Fib.updated@reductions[['pca']]@stdev)**2)), do.fast=TRUE, max_iter=1000, perplexity=25, verbose=T)$Y
rownames(tsne.rot) = rownames(train.Fib.updated@reductions[['pca']]@cell.embeddings)
train.Fib.updated_new_tsne = CreateDimReducObject(
  embeddings = tsne.rot,
  assay = 'RNA',
  key = 'tSNE_'
)

# This plot does not agree.
plot(train.Fib.updated_new_tsne@cell.embeddings, col=as.numeric(train.Fib.updated@meta.data[["Cluster"]]), main = "t-sne for fib with original")
plot(train.Fib.updated@reductions[["tsne"]]@cell.embeddings, col=as.numeric(train.Fib.updated@meta.data[["Cluster"]]), main = "t-sne for fib in original t-sne")


###################################################################################

# # Batch correction
#   batch.use = NULL
#   print('Batch correcting with ComBat')
#   
#   # Get batch vector
#   if(is.null(batch.use)){
#     batch.use = train.Fib.updated@active.ident
#   }
#   batch.use = batch.use[names(train.Fib.updated@active.ident)]
#   
#   # Batch correct with ComBat
#   bc.data = ComBat(dat=as.matrix(train.Fib.updated@assays[['RNA']]@data), batch.use, par.prior=TRUE, prior.plots=FALSE)
#   pc.data = t(scale(t(bc.data), center=F))


  
pc.data = as.matrix(train.Fib.updated@assays[['RNA']]@data)

# Fast PCA

print('Calculating PCA')

pc.data = pc.data[intersect(rownames(train.Fib.updated@reductions[['pca']]@feature.loadings), rownames(pc.data)),]
fwrite(as.data.table(pc.data), file='pc.data.txt', sep='\t')

# pc.data = train.Fib.updated@assays[['RNA']]@data[intersect(rownames(train.Fib.updated@reductions[['pca']]@feature.loadings), rownames(train.Fib.updated@assays[['RNA']]@data)),]
pca.obj = rpca(t(pc.data), center=TRUE, scale=TRUE, retx=TRUE, k=30)

train.Fib.updated_new_pca =  CreateDimReducObject(
  embeddings = pca.obj$x,
  loadings = pca.obj$rotation,
  assay = 'RNA',
  stdev = pca.obj$sdev,
  key = 'PC_',
  misc = list(total.variance = sum(pca.obj$sdev))
)
plot(train.Fib.updated_new_pca@cell.embeddings[,1], (train.Fib.updated@reductions[['pca']]@cell.embeddings %*% diag(1/((train.Fib.updated@reductions[['pca']]@stdev)**2)))[,1])
# Calculate tSNE
print('Calculating tSNE')
tsne.rot = Rtsne(train.Fib.updated_new_pca@cell.embeddings, do.fast=TRUE, max_iter=1000, perplexity=25, verbose=T)$Y
rownames(tsne.rot) = rownames(train.Fib.updated_new_pca@cell.embeddings)
train.Fib.updated_new_tsne = CreateDimReducObject(
  embeddings = tsne.rot,
  assay = 'RNA',
  key = 'tSNE_'
)
plot(train.Fib.updated_new_tsne@cell.embeddings, col=as.numeric(train.Fib.updated@meta.data[["Cluster"]]), main = "From data to t-sne")
plot(train.Fib.updated@reductions[["tsne"]]@cell.embeddings, col=as.numeric(train.Fib.updated@meta.data[["Cluster"]]), main = "pre-stored t-sne")


pc.data_latentcor = read.table("pc_latentcor.txt", sep=" ", header=FALSE)
rownames(pc.data_latentcor) = colnames(pc.data_latentcor) = colnames(pc.data)
pca.obj = rpca(pc.data_latentcor, center=TRUE, scale=TRUE, retx=TRUE, k=30)
train.Fib.updated_new_pca =  CreateDimReducObject(
  embeddings = pca.obj$x,
  loadings = pca.obj$rotation,
  assay = 'RNA',
  stdev = pca.obj$sdev,
  key = 'PC_',
  misc = list(total.variance = sum(pca.obj$sdev))
)
print('Calculating tSNE')
tsne.rot = Rtsne(train.Fib.updated_new_pca@cell.embeddings, do.fast=TRUE, max_iter=1000, perplexity=25, verbose=T)$Y
rownames(tsne.rot) = rownames(train.Fib.updated_new_pca@cell.embeddings)
train.Fib.updated_new_tsne = CreateDimReducObject(
  embeddings = tsne.rot,
  assay = 'RNA',
  key = 'tSNE_'
)
plot(train.Fib.updated_new_tsne@cell.embeddings, col=as.numeric(train.Fib.updated@meta.data[["Cluster"]]), main = "From data to t-sne (latentcor before pca)")

print('Calculating tSNE')
tsne.rot = Rtsne(pc.data_latentcor, do.fast=TRUE, max_iter=1000, perplexity=25, verbose=T)$Y
rownames(tsne.rot) = rownames(train.Fib.updated_new_pca@cell.embeddings)
train.Fib.updated_new_tsne = CreateDimReducObject(
  embeddings = tsne.rot,
  assay = 'RNA',
  key = 'tSNE_'
)
plot(train.Fib.updated_new_tsne@cell.embeddings, col=as.numeric(train.Fib.updated@meta.data[["Cluster"]]), main = "From data to t-sne")
