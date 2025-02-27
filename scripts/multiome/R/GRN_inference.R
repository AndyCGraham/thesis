#' Calculates Gene-Regulatory Networks (GRNs) for a set of transcription factors (TFs), using extreme gradient boosting regression.
#' 
#' @param TFs Character vector of names of TFs (present in rownames of chromvarMat), to infer GRNs for.
#' @param rnaMatrix Normalised gene count matrix. 
#' @param peakMatrix Normalised peak count matrix..
#' @param motifMat Motif x Peak matrix denoting which peaks contain which binding motifs.
#' @param chromvarMat Motif x Cell matrix, containing chromvar motif activity score in each cell.
#' @param annotation Genome annotation.
#' @param genesets List of genesets for enrichment analysis.
#' @param enrichment Boolean specifying whether to conduct gene-set over enrichment analysis of GRN target genes (default: TRUE).
#' @param distance Distance from gene promotors to find regulatory peaks, in base pairs (default: 200000).
#' @return List containing GRN stats for every TF, with accompanying enrichment analysis results if enrichment=TRUE.
#' @export
#' 
#' @importFrom GenomicRanges  GRanges
#' @importFrom dplyr  inner_join
#' @importFrom clusterProfiler  enricher 
#' 
FindTFregulons = function(TFs, rnaMatrix, peakMatrix, motifMat, boundRegions=NULL, chromvarMat,chromvar.obj, annotation, genesets, 
                          enrichment = T, distance = 150000, ...){
  
  #Get gene locations relative to peaks
  gene.coords <- CollapseToLongestTranscript(
    ranges = annotation
  )
  #Remove peaks not near a gene
  Peaks = rownames(peakMatrix)
  Peaks= gsub(":", "-", Peaks) #Convert alternate peakname variations
  Peaks = GRanges(seqnames=sapply(strsplit(Peaks, "-"), `[`, 1),
                                 ranges=IRanges(start = as.numeric(sapply(strsplit(Peaks, "-"), `[`, 2)), 
                                                end = as.numeric(sapply(strsplit(Peaks, "-"), `[`, 3))))
  
  gene.coords = gene.coords[gene.coords$gene_name %in% rownames(rnaMatrix),]
  peak_distance_matrix <- DistanceToTSS(
    peaks = Peaks,
    genes = gene.coords,
    distance = distance
  )
  
  #Only want peaks in peak matrix which are within the window of at least one gene
  peakMatrix = peakMatrix[match(rownames(peak_distance_matrix),rownames(peakMatrix)), ]
  peakMatrix = peakMatrix[Matrix::rowSums(peak_distance_matrix) > 0,]
  
  motifMat = motifMat[, which(colnames(motifMat) %in% TFs), drop=FALSE ]
  
  boundRegions = lapply(boundRegions, \(TF) TF[TF %in% rownames(peakMatrix)])
  
  if(is.null(boundRegions)){
    message("Finding Peaks Correlating to TF Activities")
    TFtoPEAK = lapply(setNames(TFs, TFs), \(TF) TFtoPEAK(TF, peakMatrix, motifMat, chromvarMat[rownames(chromvarMat) == TF,], ))
  } else {
    message("Recomputing TF Activities From Bound Regions")
    for(TF in TFs){
      motifMat[!rownames(motifMat) %in% boundRegions[[TF]], TF, drop=FALSE] = FALSE
    }
    chromvarMat = chromVAR::computeDeviations(object = chromvar.obj, annotations = motifMat) #Calculate average accessibility of bound peaks
    chromvarMat = SummarizedExperiment::assays(chromvarMat)[[2]]
    chromvarMat[is.na(chromvarMat)] = 0
    
    TFtoPEAK = lapply(setNames(TFs, TFs), function(TF){
      corr = map_dbl(boundRegions[[TF]], \(region) cor(peakMatrix[rownames(peakMatrix) == region,], chromvarMat[TF,]) )
      data.frame(peak = boundRegions[[TF]], TFtoPeakIMP = 100, TFtoPeakR = corr, impXr=ifelse(corr < 0, -100, 100))
      })
  }
  
  message("Finding Genes Correlating to TF Activities")
  TFtoGENE = lapply(setNames(TFs, TFs), \(TF) TFtoGENE(TF, rnaMatrix, chromvarMat[TF,], 
                                                       distance=distance, gene.coords,  TFtoPEAK[[TF]]$peak))
  message("Inferring Links between Peaks and Genes")
  PEAKtoGENE = PEAKtoGENE(rnaMatrix, peakMatrix, distance=distance, gene.coords,
                          unique(unlist(sapply(TFtoGENE, \(TF) TF$gene))),
                          unique(unlist(sapply(TFtoPEAK, \(TF) TF$peak))))
  
  
  GRNs = lapply(setNames(names(TFtoGENE),names(TFtoGENE)), function(TF){
      
    regulatoryRegions = rownames(motifMat[, which(colnames(motifMat) %in% TF), drop=FALSE ])
      
      if(all(nrow(TFtoGENE[[TF]])>5, nrow(TFtoPEAK[[TF]])>5, sum(PEAKtoGENE$peak %in% TFtoPEAK[[TF]]$peak)>5)){
      
      GRN = inner_join(TFtoGENE[[TF]], PEAKtoGENE[PEAKtoGENE$peak %in% regulatoryRegions, 
                                                  colnames(PEAKtoGENE) %in%c("gene","peak", "PeaktoGeneIMP", "PeaktoGeneR", "PeaktoGeneIMPxR")],
                       by = "gene")
      
      GRN = inner_join(GRN, TFtoPEAK[[TF]],
                       by = "peak")
      
      GRN$TF = rep(TF)
      
      GRN_pos = GRN[GRN$TFtoGeneR > 0 & GRN$PeaktoGeneR > 0 & GRN$TFtoPeakR > 0,]
      GRN_neg = GRN[GRN$TFtoGeneR < 0 & GRN$PeaktoGeneR > 0 & GRN$TFtoPeakR < 0,]
      
      TERM2GENE = data.frame(GSname = rep(names(genesets), times = sapply(genesets, length)), Genes = unlist(genesets))
      
      if(all(enrichment, length(unique(GRN_pos$gene)) > 5)){
        
        Result = enricher(unique(GRN_pos$gene),
                          universe = rownames(rnaMatrix), minGSSize = 5, pAdjustMethod = "BH",
                          maxGSSize = 1000, qvalueCutoff = 0.1, TERM2GENE=TERM2GENE)
        GRN_pos = list(links = GRN_pos, GSEA = Result)
      }
      
      if(all(enrichment, length(unique(GRN_neg$gene)) > 5)){
        
        Result = enricher(unique(GRN_neg$gene),
                          universe = rownames(rnaMatrix), minGSSize = 5, pAdjustMethod = "BH",
                          maxGSSize = 1000, qvalueCutoff = 0.1, TERM2GENE=TERM2GENE)
        GRN_neg = list(links = GRN_neg, GSEA = Result)
      }
      
      
      return(list(pos = GRN_pos, neg = GRN_neg))
      
    } else {
      return(NULL)
    }
  })
  
  return(list(GRNs = GRNs, TFactivities = chromvarMat))
}
  

#' Calculates peaks correlating to the activity score of a transcription factor (TFs), using extreme gradient boosting 
#' regression.
#' 
#' @param TF TF name.. 
#' @param peakMatrix Normalised peak count matrix..
#' @param motifMat Motif x Peak matrix denoting which peaks contain which binding motifs.
#' @param TFactivity Numeric vector, containing chromvar motif activity score of TF in each cell.
#' @param distance Distance from gene promotors to find regulatory peaks, in base pairs (default: 200000).
#' @return Dataframe containing TF activity to peak stats for TF.
#' @export
#' 
#' @importFrom caret  varImp  train
#' @importFrom purrr  map_dbl
#' 
TFtoPEAK = function(TF, peakMatrix, motifMat, TFactivity, distance = distance, ...){

  motifMat = motifMat[, which(colnames(motifMat) %in% TF), drop=FALSE ]
  peakList = rownames(motifMat[motifMat[,1], 1, drop=FALSE])
  peakMatrix = peakMatrix[rownames(peakMatrix) %in% peakList,]
  
  #Find peaks significantly positively correlated to TF activity
  #importance = varImp(train(x=as.matrix(t(peakMatrix)),y=TFactivity,ntree=500, method = "rf", nodesize=100, sampsize=1000), scale=T)$importance[,1] 
  
  params=expand.grid(nrounds=c(500), max_depth = c(250),eta =c(0.01),gamma=c(0),
              colsample_bytree=c(0.9),min_child_weight=c(15),subsample=c(0.9))  
  importance = varImp(train(x=as.matrix(t(peakMatrix)),y=TFactivity, method = "xgbTree",
                            tuneGrid=params, trControl = caret::trainControl(
                              method = "none",
                              verboseIter = FALSE, # no training log
                              allowParallel = TRUE # FALSE for reproducible results 
                            )), scale=T)$importance
  importance = importance[match(rownames(peakMatrix), rownames(importance)),]
  
  #peakList = peakList[which(importance > quantile(importance, 0.9))]
  #importance = importance[which(importance > quantile(importance, 0.9))]
  corr = map_dbl(rownames(peakMatrix), \(gene) cor(peakMatrix[rownames(peakMatrix) == gene,],TFactivity) )
  CorrtoTFexpression = F
  
  out = data.frame(peak = rownames(peakMatrix), TFtoPeakIMP = importance,
                   TFtoPeakR = corr, impXr=ifelse(corr < 0, -importance, importance))
  return(out)
}

#' Calculates genes correlating to the activity score of a transcription factor (TFs), using extreme gradient boosting 
#' regression.
#' 
#' @param TF TF name.. 
#' @param rnaMatrix Normalised gene count matrix..
#' @param TFactivity Numeric vector, containing chromvar motif activity score of TF in each cell.
#' @param distance Distance from gene promotors to find regulatory peaks, in base pairs (default: 200000).
#' @param gene.coords Dataframe containing genomic coordinates of genes in rnaMatrix.
#' @param peakList Character vector of peaks potentially regulated by TF.
#' @return Dataframe containing TF activity to gene stats for TF.
#' @export
#' 
#' @importFrom caret  varImp train
#' @importFrom GenomicRanges  GRanges
#' @importFrom purrr  map_dbl
#' 
TFtoGENE = function(TF, rnaMatrix, distance = distance, TFactivity, gene.coords, peakList, ...){
  
  #Get gene locations relative to peaks
  Peaks = peakList
  Peaks= gsub(":", "-", Peaks) #Convert alternate peakname variations
  Peaks = GRanges(seqnames=sapply(strsplit(Peaks, "-"), `[`, 1),
                                 ranges=IRanges(start = as.numeric(sapply(strsplit(Peaks, "-"), `[`, 2)), 
                                                end = as.numeric(sapply(strsplit(Peaks, "-"), `[`, 3))))
  
  gene.coords = gene.coords[gene.coords$gene_name %in% rownames(rnaMatrix),]
  peak_distance_matrix <- DistanceToTSS(
    peaks = Peaks,
    genes = gene.coords,
    distance = distance
  )
  
  #Only want genes in RNA matrix which are in annot and next to at least one peak harbouring TF Motif
  rnaMatrix = rnaMatrix[match(colnames(peak_distance_matrix),rownames(rnaMatrix)), ]
  rnaMatrix = rnaMatrix[Matrix::colSums(peak_distance_matrix) > 0,]
  
  #Set Up xgboost params
  params = expand.grid(nrounds=c(500), max_depth = c(min(250, max(nrow(rnaMatrix)-5, 5))),eta =c(0.01),gamma=c(0),
                       colsample_bytree=c(0.9),min_child_weight=c(15),subsample=c(0.9))    

  #Find genes significantly positively correlated to TF activity
  TFtoGeneImportance = varImp(train(x=t(as.matrix(rnaMatrix)),y=TFactivity, method = "xgbTree",
                                    tuneGrid=params, trControl = caret::trainControl(
                                      method = "none",
                                      verboseIter = F, # no training log
                                      allowParallel = TRUE # FALSE for reproducible results 
                                    )), scale=T)$importance
  TFtoGeneImportance = TFtoGeneImportance[match(rownames(rnaMatrix), rownames(TFtoGeneImportance)),]
  
  selected = which(TFtoGeneImportance > quantile(TFtoGeneImportance, 0.75))
  rnaMatrix = rnaMatrix[selected,]
  
  TFtoGeneImportance = TFtoGeneImportance[selected]
  names(TFtoGeneImportance) = rownames(rnaMatrix)
  
  #Iterate over each gene within specified distance from peak
  corr = purrr::map_dbl(rownames(rnaMatrix), \(gene) cor(rnaMatrix[rownames(rnaMatrix) == gene,],TFactivity) )
  
  TFtoGene = data.frame(gene=rownames(rnaMatrix), TFtoGeneIMP = TFtoGeneImportance, TFtoGeneR = corr,
                        TFtoGeneIMPxR= ifelse(corr < 0, -TFtoGeneImportance, TFtoGeneImportance))
  return(TFtoGene)
}

#' Infers peaks correlating to the aexpression of nearby genes.
#' 
#' @param TF TF name.. 
#' @param rnaMatrix Normalised gene count matrix..
#' @param peakMatrix Normalised peak count matrix..
#' @param distance Distance from gene promotors to find regulatory peaks, in base pairs (default: 200000).
#' @param gene.coords Dataframe containing genomic coordinates of genes in rnaMatrix.
#' @param peakList Character vector of peaks potentially regulated by any TF.
#' @return Dataframe containing peak to gene stats.
#' @export
#' 
#' @importFrom caret  varImp train
#' @importFrom GenomicRanges  GRanges
#' @importFrom purrr  map_dfr  map_dbl
#' 
PEAKtoGENE = function(rnaMatrix, peakMatrix, distance = distance, gene.coords, geneList, peakList, ...){
  
  #Get gene locations relative to peaks
  Peaks = rownames(peakMatrix)
  Peaks= gsub(":", "-", Peaks) #Convert alternate peakname variations
  Peaks = GRanges(seqnames=sapply(strsplit(Peaks, "-"), `[`, 1),
                                 ranges=IRanges(start = as.numeric(sapply(strsplit(Peaks, "-"), `[`, 2)), 
                                                end = as.numeric(sapply(strsplit(Peaks, "-"), `[`, 3))))
  
  rnaMatrix = rnaMatrix[rownames(rnaMatrix) %in% geneList,]
  gene.coords = gene.coords[gene.coords$gene_name %in% rownames(rnaMatrix),]
  peak_distance_matrix <- DistanceToTSS(
    peaks = Peaks,
    genes = gene.coords,
    distance = distance
  )
  
  #Only want genes in RNA matrix which are in annot and next to at least one peak 
  rnaMatrix = rnaMatrix[match(colnames(peak_distance_matrix),rownames(rnaMatrix)), ]
  rnaMatrix = rnaMatrix[Matrix::colSums(peak_distance_matrix) > 0,]
  peakMatrix = peakMatrix[match(rownames(peak_distance_matrix),rownames(peakMatrix)), ]
  
  params = expand.grid(nrounds=c(200), max_depth = min(20, max(nrow(rnaMatrix)-5, 5)),eta =c(0.01),gamma=c(0),
                       colsample_bytree=c(0.9),min_child_weight=c(10),subsample=c(0.9))  
  
  #Now iterate over each gene 
  PEAKtoGene = purrr::map_dfr(rownames(rnaMatrix), function(gene) {
    
    #Find nearby peaks which are important for predicting that gene's expression
    nearPeaks=peakMatrix[peak_distance_matrix[,colnames(peak_distance_matrix) == gene] == 1,]
    
    #If there's only one nearby peak, just do a linear regression
    if(is.vector(nearPeaks) == 1){
      corr = cor(rnaMatrix[rownames(rnaMatrix) == gene,],nearPeaks) 
      if(abs(corr) > 0.01){
        out = data.frame(gene=gene, peak = rownames(peakMatrix)[peak_distance_matrix[,colnames(peak_distance_matrix) == gene] == 1], 
                         PeaktoGeneIMP = 100, PeaktoGeneR = corr, PeaktoGeneIMPxR=ifelse(corr < 0, -0.1, 0.1))
        return(out)
      } else {
        return(NULL)
      }
      
    }
    
    #Otherwise use XGBoost to identify which peaks are the most imporatant predictors
    OneGeneRes = varImp(train(x=t(as.matrix(nearPeaks)),
                              y=as.numeric(rnaMatrix[rownames(rnaMatrix)==gene,]), method = "xgbTree",
                              tuneGrid=params, trControl = caret::trainControl(
                                method = "none",
                                verboseIter = F, # no training log
                                allowParallel = TRUE # FALSE for reproducible results 
                              )), scale=T)$importance
    OneGeneRes = OneGeneRes[match(rownames(nearPeaks), rownames(OneGeneRes)),]
    
    OneGeneResMotif = OneGeneRes[which(OneGeneRes > Rfast::nth(OneGeneRes, min(5, nrow(OneGeneRes)), descending = T) 
                                       & rownames(nearPeaks) %in% peakList )]
    
    if(length(OneGeneResMotif) < 1 ){
      return(NULL)
    }
    
    regulatoryRegions = rownames(nearPeaks)[which(OneGeneRes > Rfast::nth(OneGeneRes, min(5, length(OneGeneRes)), descending = T) 
                                                  & rownames(nearPeaks) %in% peakList )]
    corr = purrr::map_dbl(regulatoryRegions, 
                          \(peak) cor(rnaMatrix[rownames(rnaMatrix) == gene,],nearPeaks[rownames(nearPeaks) == peak,]) )
    
    #Iterate over each gene within specified distance from peak
    out = data.frame(gene=gene, peak = regulatoryRegions, PeaktoGeneIMP = OneGeneResMotif,
                     PeaktoGeneR = corr, PeaktoGeneIMPxR=ifelse(corr < 0, -OneGeneResMotif, OneGeneResMotif))
    return(out)
    
  } )
  
  return(PEAKtoGene)
  
} 


#' From signac
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @import data.table
#' @export
CollapseToLongestTranscript <- function(ranges) {
  range.df <- as.data.frame(x = ranges, row.names = NULL)
  range.df$strand <- as.character(x = range.df$strand)
  range.df$strand <- ifelse(
    test = range.df$strand == "*",
    yes = "+",
    no = range.df$strand
  )
  range.df = range.df[!range.df$gene_name == "",]
  collapsed <- dplyr::group_by(range.df, gene_name) %>%
    dplyr::reframe(
      gene_id = find_mode(gene_id, na.rm = T),
      seqnames = find_mode(seqnames, na.rm = T),
      start = min(start, na.rm = T),
      end = max(end, na.rm = T),
      strand = find_mode(strand, na.rm = T),
      gene_biotype = find_mode(gene_biotype, na.rm = T)
    )
  
  collapsed = collapsed[!collapsed$gene_name == "",]
  collapsed$gene_name <- make.unique(names = collapsed$gene_name)
  gene.ranges <- makeGRangesFromDataFrame(
    df = collapsed,
    keep.extra.columns = TRUE
  )
  return(gene.ranges)
}

#Finds Mode of a Vector
#
# @param x      vector
# @param na.rm  whether to remove NAs before calculating mode
#' @export
# @return Returns a sparse matrix
find_mode <- function(x, na.rm = T) {
  u = unique(x)
  if(na.rm){
    u = u[!is.na(u)]
  }
  tab = tabulate(match(x, u))
  u[tab == max(tab)]
}

# From SIgnac 
#Find peaks near genes
#
# Find peaks that are within a given distance threshold to each gene
#
# @param peaks A GRanges object containing peak coordinates
# @param genes A GRanges object containing gene coordinates
# @param distance Distance threshold. Peaks within this distance from the gene
# will be recorded.
# @param sep Separator for peak names when creating results matrix
#
#' @importFrom GenomicRanges findOverlaps
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom Matrix sparseMatrix
#' @importFrom GenomicRanges resize
#' @export
# @return Returns a sparse matrix
DistanceToTSS <- function(
    peaks,
    genes,
    distance = 200000,
    sep = c("-", "-")
) {
  tss <- resize(x = genes, width = 1, fix = 'start')
  genes.extended <- suppressWarnings(
    expr = Extend(
      x = tss, upstream = distance, downstream = distance
    )
  )
  overlaps <- findOverlaps(
    query = peaks,
    subject = genes.extended,
    type = 'any',
    select = 'all'
  )
  hit_matrix <- Matrix::sparseMatrix(
    i = queryHits(x = overlaps),
    j = subjectHits(x = overlaps),
    x = 1,
    dims = c(length(x = peaks), length(x = genes.extended))
  )
  rownames(x = hit_matrix) <- GRangesToString(grange = peaks, sep = sep)
  colnames(x = hit_matrix) <- genes.extended$gene_name
  return(hit_matrix)
}

