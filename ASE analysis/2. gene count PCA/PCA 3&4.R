# creating a function to find PCAs 3 and 4

PCA_more = function(object, intgroup="condition", ntop=500, returnData=FALSE)
{
  # calculate the variance for each gene
  rv <- rowVars(assay(object))
  
  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  
  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select,]))
  
  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop=FALSE])
  
  # add the intgroup factors together to create a new grouping factor
  group <- if (length(intgroup) > 1) {
    factor(apply( intgroup.df, 1, paste, collapse=":"))
  } else {
    colData(object)[[intgroup]]
  }
  
  # assembly the data for the plot
  d <- data.frame(PC3=pca$x[,3], PC4=pca$x[,4], group=group, intgroup.df, name=colnames(object))
  
  if (returnData) {
    attr(d, "percentVar") <- percentVar[3:4]
    return(d)
  }
  
  ggplot(data=d, aes_string(x="PC3", y="PC4", color="group")) + geom_point(size=3) + 
    xlab(paste0("PC3: ",round(percentVar[3] * 100),"% variance")) +
    ylab(paste0("PC4: ",round(percentVar[4] * 100),"% variance")) +
    coord_fixed()
}
