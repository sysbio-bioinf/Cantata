evalNetworks <- function(folder, original)
{
  del <- c("map3k7", "map3k1", "mkk4", "ras")
  #read networks
  networks <- extractBestNetworks(folder)
  #build adjacency matrices for each ne
  adjMatrices <- createAdjacencyMatrices(networks)
  #for each possible edge in network list networks, that contain this edge
  edges <- lookUpEdges(adjMatrices,original$genes)
  #floyd warshall 
  reducedEdges <- getAllPaths(edges,del)
  #create probability matrix 
  countMatrix <- createCountMatrix(reducedEdges)
  
  rownames(countMatrix) <- original$genes
  colnames(countMatrix) <- original$genes
  
  #genes to delete
  
  
  shrinked <- countMatrix[is.na(match(rownames(countMatrix),del)),is.na(match(colnames(countMatrix),del))]
  shrinked[shrinked < 20] <- 0
  #shrinked <- t(shrinked)
  #create true matrix from orig network
  trueMatrix <- createAdjacencyMatrices(list(original))
  trueEdges <- lookUpEdges(trueMatrix,original$genes)
  reducedTrueEdges <- getAllPaths(trueEdges,del) #trueEdges #
  print(reducedTrueEdges)
  trueMatrix <- createCountMatrix(reducedTrueEdges)
  rownames(trueMatrix) <- original$genes
  colnames(trueMatrix) <- original$genes
  
  shrinkedTrueMatrix <- trueMatrix[is.na(match(rownames(trueMatrix),del)),is.na(match(colnames(trueMatrix),del))]
  #shrinkedTrueMatrix <- t(shrinkedTrueMatrix)
  print(shrinkedTrueMatrix)
  res <- compareNets(countMatrix=shrinked,trueNetMatrix=shrinkedTrueMatrix, original)
}  


extractBestNetworks <- function(folder)
{
  subdirs <- list.dirs(folder,recursive=F)
  files <- list()
  for (dir in subdirs)
  {
    files[[length(files) + 1]] <- paste(dir,"resultnet_1.txt",sep="/")
  }
  
  networks <- lapply(files, loadNetwork)
  
  return(networks)
  
}


createAdjacencyMatrices <- function(networks)
{
  dependencies <- lapply(networks, function(net)
  {
    sapply(net$interactions, function(int)
    {
      r <- rep(0, length(net$genes))
      r[int$input] <- 1
      return(r)                                   
    })
  })
  
  return(dependencies)
}


lookUpEdges <- function(adjMatrices, genes)
{
  #init matrix
  matrix <- vector('list', nrow(adjMatrices[[1]]))
  names(matrix) <- genes
  matrix <- lapply(matrix, function(x) {x <- rep(list(NA), ncol(adjMatrices[[1]]))
                                        names(x) <- genes
                                        return(x)})
  for(i in 1:length(adjMatrices))
  {
    edgeIndices <- which(adjMatrices[[i]] == 1, arr.ind=T)
    
    for(j in 1:nrow(edgeIndices))
    {
      if(all(is.na(matrix[[edgeIndices[j,1]]][[edgeIndices[j,2]]])))
        matrix[[edgeIndices[j,1]]][[edgeIndices[j,2]]] <- i
      else
        matrix[[edgeIndices[j,1]]][[edgeIndices[j,2]]] <- c(matrix[[edgeIndices[j,1]]][[edgeIndices[j,2]]], i)
    }
  }
  return(matrix)  
}

getAllPaths <- function(edges, nodes)
{
  allNodes <- names(edges)
  for(k in nodes)
  {
    for(i in allNodes)
    {
      for(j in allNodes)
      {
        if(!any(is.na(c(edges[[i]][[k]],edges[[k]][[j]]))))
          edges[[i]][[j]] <- union(edges[[i]][[j]][!is.na(edges[[i]][[j]])], intersect(edges[[i]][[k]], edges[[k]][[j]]))
      }
    }
  }
  
  return(edges)
}

createCountMatrix <- function(reducedEdges)
{
  mat <- sapply(reducedEdges, function(x) sapply(x, function(y) sum(!is.na(y))))
  for (input in c("il1a","tgfa","igf1","tnfa"))
    mat[input, input] <- 0
  mat
}

compareNets <- function(countMatrix, trueNetMatrix,trueNet,
                        plotIt=TRUE,...)
{
  #require(gmp)
  library(BoolNet)

  #g <- plotNetworkWiring(trueNet, plotIt=FALSE)
  set.seed(123124)
  
  
  
  
  
  reconstDepMat <- countMatrix
  trueDepMat <- trueNetMatrix
  
  trueNet <- list()
  trueNet$genes <- colnames(trueNetMatrix)
  
  totalFuncs <- 0
  
  
  

  
  
  #sumFuncs <- apply(reconstDepMat,1,sum)
  print(reconstDepMat)
  #print(sumFuncs)
  
  if (plotIt)
  {
    library(igraph)
    #source("plot.igraph.R")
    vertices <- trueNet$genes
    print("---------vertices--------")
    print(vertices)
    vertexCols <- c()
    edgeCols <- c()
    edgeLabels <- c()
    edgeTypes <- c()
    edgeWidths <- c()
    edges <- c()
    for (gene in trueNet$genes)
    {
        vertexCols <- c(vertexCols, "lightblue")
        factors <- unique(c(trueNet$genes[reconstDepMat[gene,] > 0],
                            trueNet$genes[trueDepMat[gene,] > 0]))
        print(factors)
      
        for (factor in factors)
        {
          edges <- rbind(edges,c(factor, gene))
          if (reconstDepMat[gene,factor] > 0 && trueDepMat[gene,factor] > 0)
            {
              edgeCols <- c(edgeCols,"green")
              edgeTypes <- c(edgeTypes, "solid")
              edgeLabels <- c(edgeLabels, paste(reconstDepMat[gene,factor],"%",sep=""))
              edgeWidths <- c(edgeWidths, max(1,
                                              round(reconstDepMat[gene,factor]/100 * 5)))
            }
          else
            if (reconstDepMat[gene,factor] > 0 &&  trueDepMat[gene,factor] == 0)
            {
              edgeCols <- c(edgeCols, "red")
              edgeTypes <- c(edgeTypes, "solid")
              edgeLabels <- c(edgeLabels, paste(reconstDepMat[gene,factor],"%",sep=""))
              edgeWidths <- c(edgeWidths, max(1,
                                              round(reconstDepMat[gene,factor]/100 * 5)))
            }
          else
            if (reconstDepMat[gene,factor] == 0 && trueDepMat[gene,factor] > 0)  
          {
            edgeCols <- c(edgeCols, "red")
            edgeTypes <- c(edgeTypes, "dashed")
            edgeLabels <- c(edgeLabels, "")
            edgeWidths <- c(edgeWidths, 1)
            print(paste("gene :", gene, " factor : ", factor))
          }
      }
    }
  }
  
  g <- graph.data.frame(as.data.frame(edges), vertices=data.frame(vertices), directed=TRUE)
  
  g <- set.edge.attribute(g, "color", value=edgeCols, index=seq_len(nrow(edges)))
  g <- set.edge.attribute(g, "label.color", value=edgeCols, index=seq_len(nrow(edges)))
  g <- set.edge.attribute(g, "label", value=edgeLabels, index=seq_len(nrow(edges)))
  g <- set.edge.attribute(g, "lty", value=edgeTypes, index=seq_len(nrow(edges)))
  g <- set.edge.attribute(g, "width", value=edgeWidths, index=seq_len(nrow(edges)))
  g <- set.vertex.attribute(g, "label", value=vertices, index=seq_along(vertices))
  g <- set.vertex.attribute(g, "color", value=vertexCols, index=seq_along(vertices))
  
  #autocurve.edges(g)
  #pdf(file="../../CantataPaper/plot.pdf",width=20,height=20)
  lay <- layout.sugiyama(g, hgap=1.5)$layout
  #par(mai=c(0,0,1,0), mar=c(2,2,2,2))
  tkplot(g,
       edge.arrow.size = 0.5,
       vertex.size=10,
       vertex.size2=7,
       vertex.label.cex=1.2,
       vertex.shape="crectangle",
       edge.label.cex=1.0,
       edge.arrow.mode=2,
       edge.arrow.width=2,
       edge.curved=F,
       layout=lay,
       ...)
  #dev.off()
  #reconstDepMat <- (reconstDepMat > 0) * 1
  #TP <- sum(trueDepMat * reconstDepMat)
  #FP <- sum((1-trueDepMat) * reconstDepMat)
  #TN <- sum((1-trueDepMat) * (reconstDepMat == 0))
  #FN <- sum(trueDepMat * (reconstDepMat == 0))
  #TN <- sum((1-trueDepMat) * (sumFuncs - reconstDepMat))
  #FN <- sum(trueDepMat * (sumFuncs - reconstDepMat))
  #prec <- TP/(TP+FP) * 100
  #recall <- TP/(TP+FN) * 100
  #spec <- TN/(TN+FP) * 100
  
  #res <- round(c(TP,FP,TN,FN, prec, recall, spec),2)
  #names(res) <- c("TP","FP","TN","FN", "Precision", "Recall/Sensitivity", "Specificity")
  
  #return(res)
  return(g)
}


