
  source("~/ulcerative_colitis/ulcerative_colitis/analysis.r")
  train.Fib.seur <- readRDS("~/ulcerative_colitis/ulcerative_colitis/train.Fib.seur.rds")
  train.Fib.updated = UpdateSeuratObject(object = train.Fib.seur)
  
  # Fast PCA
  print('Calculating PCA')
  
  pc.data = train.Fib.updated@assays[['RNA']]@data[intersect(rownames(train.Fib.updated@reductions[['pca']]@feature.loadings), rownames(train.Fib.updated@assays[['RNA']]@data)),]
  # pc.data = latentcor(pc.data, types = rep("tru", ncol(pc.data)), use.nearPD = FALSE)$R
  pca.obj = rpca(t(pc.data), center=TRUE, scale=TRUE, retx=TRUE, k=30)

  train.Fib.updated_new_pca =  CreateDimReducObject(
    embeddings = pca.obj$x %*% diag(pca.obj$sdev**2),
    loadings = pca.obj$rotation,
    assay = 'RNA',
    stdev = pca.obj$sdev,
    key = 'PC_',
    misc = list(total.variance = sum(pca.obj$sdev))
  )

  # Calculate tSNE
  print('Calculating tSNE')
  tsne.rot = Rtsne(train.Fib.updated_new_pca@cell.embeddings[,1:30], do.fast=TRUE, max_iter=1000, perplexity=25, verbose=T)$Y
  rownames(tsne.rot) = rownames(train.Fib.updated_new_pca@cell.embeddings)
  train.Fib.updated_new_tsne = CreateDimReducObject(
    embeddings = tsne.rot,
    assay = 'RNA',
    key = 'tSNE_'
  )
  plot(train.Fib.updated_new_tsne@cell.embeddings, col=as.numeric(train.Fib.updated@meta.data[["Cluster"]]), main = "t-sne for fib with original")
  plot(train.Fib.updated@reductions[["tsne"]]@cell.embeddings, col=as.numeric(train.Fib.updated@meta.data[["Cluster"]]), main = "t-sne for fib in original t-sne")
  
  train.Imm.seur <- readRDS("~/ulcerative_colitis/ulcerative_colitis/train.Imm.seur.rds")
  train.Imm.updated = UpdateSeuratObject(object = train.Imm.seur)
  
  # Fast PCA
  print('Calculating PCA')
  
  pc.data = train.Imm.updated@assays[['RNA']]@data[intersect(rownames(train.Imm.updated@reductions[['pca']]@feature.loadings), rownames(train.Imm.updated@assays[['RNA']]@data)),]
  # pc.data = latentcor(pc.data, types = rep("tru", ncol(pc.data)), use.nearPD = FALSE)$R
  pca.obj = rpca(t(pc.data), center=TRUE, scale=TRUE, retx=TRUE, k=30)
  
  
  train.Imm.updated_new_pca =  CreateDimReducObject(
    embeddings = pca.obj$x %*% diag(pca.obj$sdev**2),
    loadings = pca.obj$rotation,
    assay = 'RNA',
    stdev = pca.obj$sdev,
    key = 'PC_',
    misc = list(total.variance = sum(pca.obj$sdev))
  )
  
  # Calculate tSNE
  print('Calculating tSNE')
  tsne.rot = Rtsne(train.Imm.updated_new_pca@cell.embeddings[,1:30], do.fast=TRUE, max_iter=1000, perplexity=25, verbose=T)$Y
  rownames(tsne.rot) = rownames(train.Imm.updated_new_pca@cell.embeddings)
  train.Imm.updated_new_tsne = CreateDimReducObject(
    embeddings = tsne.rot,
    assay = 'RNA',
    key = 'tSNE_'
  )
  plot(train.Imm.updated_new_tsne@cell.embeddings, col=as.numeric(train.Imm.updated@meta.data[["Cluster"]]), main = "t-sne for imm with original")
  plot(train.Imm.updated@reductions[["tsne"]]@cell.embeddings, col=as.numeric(train.Imm.updated@meta.data[["Cluster"]]), main = "t-sne for imm in original t-sne")

  train.Epi.seur <- readRDS("~/ulcerative_colitis/ulcerative_colitis/train.Epi.seur.rds")
  train.Epi.updated = UpdateSeuratObject(object = train.Epi.seur)
  
  # Fast PCA
  print('Calculating PCA')
  
  pc.data = train.Epi.updated@assays[['RNA']]@data[intersect(rownames(train.Epi.updated@reductions[['pca']]@feature.loadings), rownames(train.Epi.updated@assays[['RNA']]@data)),]
  # pc.data = latentcor(pc.data, types = rep("tru", ncol(pc.data)), use.nearPD = FALSE)$R
  pca.obj = rpca(t(pc.data), center=TRUE, scale=TRUE, retx=TRUE, k=30)
  
  
  train.Epi.updated_new_pca =  CreateDimReducObject(
    embeddings = pca.obj$x %*% diag(pca.obj$sdev**2),
    loadings = pca.obj$rotation,
    assay = 'RNA',
    stdev = pca.obj$sdev,
    key = 'PC_',
    misc = list(total.variance = sum(pca.obj$sdev))
  )
  
  # Calculate tSNE
  print('Calculating tSNE')
  tsne.rot = Rtsne(train.Epi.updated_new_pca@cell.embeddings[,1:30], do.fast=TRUE, max_iter=1000, perplexity=25, verbose=T)$Y
  rownames(tsne.rot) = rownames(train.Epi.updated_new_pca@cell.embeddings)
  train.Epi.updated_new_tsne = CreateDimReducObject(
    embeddings = tsne.rot,
    assay = 'RNA',
    key = 'tSNE_'
  )
  plot(train.Epi.updated_new_tsne@cell.embeddings, col=as.numeric(train.Epi.updated@meta.data[["Cluster"]]), main = "t-sne for Epi with original")
  plot(train.Epi.updated@reductions[["tsne"]]@cell.embeddings, col=as.numeric(train.Epi.updated@meta.data[["Cluster"]]), main = "t-sne for Epi in original t-sne")
  
  
  # 
  # # Cluster cells with Phenograph
  # print('Clustering with Phenograph')
  # clusters = run_phenograph(seur@reductions[['pca']]@cell.embeddings[,1:num_pcs], k=k)
  # seur@meta.data[,colnames(clusters)] = clusters
