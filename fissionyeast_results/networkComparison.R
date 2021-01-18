library(BoolNet)

pairwiseDistances <- function(files)
{
  tts <- lapply(files,function(file)
          {
            n <- loadNetwork(file)
            tt <- getTransitionTable(getAttractors(n))
            return(tt[(length(n$genes)+1):(2*length(n$genes))])
          })
  dist <- matrix(ncol=length(tts),nrow=length(tts),NA)        
  for (i in 1:(length(tts)-1))
    for (j in ((i+1):length(tts)))
    {
      print(i)
      print(j)
      dist[i,j] <- sum(sapply(1:nrow(tts[[i]]),function(k)
      {
         any(tts[[i]][k,] != tts[[j]][k,])
      }))
    }
  return(dist)
}

getRelevantDependencies <- function(networkFile,exPath="Release/NetworkAnalyzer")
{
  network <- loadNetwork(networkFile)
  att1 <- getAttractors(network)
  mat <- sapply(network$genes,function(target)
        {
          sapply(network$genes,function(dependency)
          {
            system(paste(exPath,"--truncate -n", networkFile,
                         "-o tmpnet.txt -tg", target, "-dg", dependency),intern=TRUE)
            n <- loadNetwork("tmpnet.txt")
            att2 <- getAttractors(n)
            return(!compareAttractors(att1, att2))
          })
        })
  system("rm tmpnet.txt")
  return(mat)        
}

compareNetworks <- function(origFile, modifiedFiles, bodySeparator = ",", compareAttractors=TRUE, compareTransitionTables=FALSE)
{
  #source("../functionClassAnalysis.R")
  origNet <- loadNetwork(origFile)
  #origFunctionClasses <- summarizeFunctionClasses(origNet)
  
  if (length(origNet$genes) > 15)
    startStates = 10000
  else
    startStates = NULL
  
  if (compareAttractors || compareTransitionTables)
  {
    sink("/dev/null")
    origAtt <- getAttractors(origNet, startStates=startStates)
  }
  origNumDependencies <- sum(sapply(origNet$interactions,function(int)length(int$input)))
 
  modifiedNets <- lapply(modifiedFiles, function(file)loadNetwork(file, bodySeparator=bodySeparator)) 
  
  numDependencies <- sapply(modifiedNets,function(net)sum(sapply(net$interactions,function(int)length(int$input))))
  
  if (compareAttractors || compareTransitionTables)
  {
    modAtt <- lapply(modifiedNets, function(net)
                     getAttractors(net, startStates=startStates))
    sink()
    attractorsEqual <- sapply(modAtt, function(att)
                          {
                            compareAttractors(origAtt, att)
                          })
  }
  else
    attractorsEqual <- NULL
                            
  if (compareTransitionTables)
  {
    ttEqual <- sapply(modAtt, function(att)
                      {
                        identical(att$stateInfo$table, origAtt$stateInfo$table)
                      })
  }
  else
    ttEqual <- NULL
  
  if (compareAttractors)
  {
    usedAttractors <- attractorsEqual
  }
  else
    usedAttractors <- rep(TRUE, length(modifiedNets))
  
  dependencies <- lapply(modifiedNets[usedAttractors], function(net)
                      {
                          sapply(net$interactions, function(int)
                          {
                            r <- rep(0, length(origNet$genes))
                            r[int$input] <- 1
                            return(r)                                   
                          })
                      }) 
    
  insertions <- lapply(modifiedNets[usedAttractors], function(net)
                      {
                          if (length(origNet$interactions) != length(net$interactions))
                          {
                            print(origNet$interactions)
                            print(net$interactions)
                          }
                          mapply(function(int1, int2)
                          {
                            s <- setdiff(int2$input, int1$input)
                            r <- rep(0, length(origNet$genes))
                            r[s] <- 1
                            return(r)                                   
                          }, origNet$interactions, net$interactions)
                      })
                      
  deletions <- lapply(modifiedNets[usedAttractors], function(net)
                      {
                          mapply(function(int1, int2)
                          {
                            s <- setdiff(int1$input, int2$input)
                            r <- rep(0, length(origNet$genes))
                            r[s] <- 1
                            return(r)                                   
                          }, origNet$interactions, net$interactions)
                      })
                      
  #functionClasses <- lapply(modifiedNets[usedAttractors], function(net)
  #                    {
  #                        r <- summarizeFunctionClasses(net)
  #                        r[is.na(r)] <- 0
  #                        r
  #                    })                 

  numInsertions <- sapply(insertions, function(ins)
                           {
                            sum(ins)
                           })
                           
  numDeletions <- sapply(deletions, function(del)
                           {
                            sum(del)
                           })
                      
  numModifications <- mapply(function(ins,del)
                           {
                            sum(ins) + sum(del)
                           }, insertions, deletions)

  dependencies <- Reduce(function(a,b)a+b, dependencies)
  
  insertions <- Reduce(function(a,b)a+b, insertions)
  
  deletions <- Reduce(function(a,b)a+b, deletions)
  
  #functionClasses <- Reduce(function(a,b)a+b, functionClasses)
                      
  rownames(dependencies) <- origNet$genes
  
  rownames(insertions) <- origNet$genes                                            
                      
  rownames(deletions) <- origNet$genes        

  return(list(origNumDependencies = origNumDependencies,
              #origFunctionClasses = origFunctionClasses,
              numDependencies = numDependencies[usedAttractors],
              dependencies = dependencies,
              insertions = insertions, deletions = deletions, 
              #functionClasses = functionClasses,
              numInsertions = numInsertions,
              numDeletions = numDeletions,
              numModifications = numModifications,
              usedAttractors = usedAttractors,
              attractorsEqual = attractorsEqual,
              ttEqual = ttEqual))
}

summarizeComparison <- function(res, fractions=TRUE)
{
  if (!is.null(res$attractorsEqual))
    cat("Number of networks with equal attractors: ",sum(res$attractorsEqual), 
        "/", length(res$attractorsEqual), "\n", sep="")
        
  if (!is.null(res$ttEqual))
    cat("Number of networks with equal transition tables: ",sum(res$ttEqual), 
        "/", length(res$ttEqual), "\n", sep="")        

  tm <- table(res$numModifications)
  ti <- table(res$numInsertions)
  td <- table(res$numDeletions)
 
  if (fractions)
  {
    res$dependencies <- round(res$dependencies/sum(res$usedAttractors) * 100, 2)
    res$insertions <- round(res$insertions/sum(res$usedAttractors) * 100, 2)
    res$deletions <- round(res$deletions/sum(res$usedAttractors) * 100, 2)
    res$functionClasses <- round(res$functionClasses/sum(res$usedAttractors) * 100, 2)
    
    tm <- round(tm / sum(tm) * 100, 2)
    ti <- round(ti / sum(ti) * 100, 2)    
    td <- round(td / sum(td) * 100, 2)    
  
  }
  
  cat("\nNumber of dependencies in the original network:\n")
  print(mean(res$origNumDependencies))
  
  cat("\nMean number of dependencies in networks:\n")
  print(mean(res$numDependencies))
  
  cat("\nNumber of dependency modifications in networks:\n")
  print(tm)
  
  cat("\nNumber of inserted dependencies in networks:\n")
  print(ti)

  cat("\nNumber of removed dependencies in networks:\n")
  print(td)
  
  cat("\nInsertions (columns: genes, rows: insertions)\n")
  print(res$insertions)
  cat("\nDeletions (columns: genes, rows: deletions)\n")  
  print(res$deletions)    
  
  cat("\nFunction classes in the original model:\n")
  print(res$origFunctionClasses)
  
  cat("\nFunction classes in the networks:\n")
  print(res$functionClasses)
}

compareAttractors <- function(att1, att2)
{
    if (length(att1$attractors) != length(att2$attractors))
      return(FALSE)
    return(all(mapply(function(att1, att2)
                  {
                    if (ncol(att1$involvedStates) != ncol(att2$involvedStates))
                      return(FALSE)
                    else
                      return(all(att1$involvedStates == att2$involvedStates))
                  }, att1$attractors, att2$attractors)))
}

readDependencies <- function(net, removeNegations=TRUE)
{
  res <- lapply(net$interactions,function(int)
          {
            s <- gsub("&","",int$expression,fixed=TRUE)
            s <- gsub("|","",s,fixed=TRUE)
            if (removeNegations)
            {
              s <- gsub("!","",s,fixed=TRUE)
            }
            else
              s <- gsub("! ","!",s,fixed=TRUE)
            s <- gsub("(","",s,fixed=TRUE)
            s <- gsub(")","",s,fixed=TRUE)
            s <- strsplit(s," ",fixed=TRUE)
            s <- s[[1]][nchar(s[[1]]) != 0]
          })
 return(res)                
}

compareDependencies <- function(origNet, reconstructedNets, removeNegations=FALSE, summarize=TRUE)
{
  origNet <- loadNetwork(origNet)
  reconstructedNets <- lapply(reconstructedNets, loadNetwork)
  origDep <- readDependencies(origNet,removeNegations=removeNegations)
  depChanges <- lapply(reconstructedNets, function(net)
  {
    reconstDep <- readDependencies(net,removeNegations=removeNegations)
    mapply(function(gene1, gene2)
    {
      cnt <- rep(0,length(origNet$genes)*2)
      names(cnt) <- c(origNet$genes, paste0("!",origNet$genes))
      diff <- setdiff(gene2,gene1)

      if (length(diff)>0)
        cnt[diff] <- 1
      return(cnt)
    }, origDep, reconstDep, SIMPLIFY=TRUE)
  })

  if (summarize)
  {
    depChanges <- Reduce(function(a,b)a+b,depChanges)

    indices <- which(depChanges>0,arr.ind=TRUE)
    tab <- t(apply(indices, 1, function(x)
    {
      c("target"=colnames(depChanges)[x[2]],
        "factor"=rownames(depChanges)[x[1]],
        "count"=depChanges[x[1],x[2]])
    }))
    print(tab)
    rownames(tab) <- c()
    tab <- data.frame("target"=tab[,1],"factor"=tab[,2],"count"=as.integer(tab[,3]))
    tab <- tab[order(tab[,"count"],decreasing=TRUE),]    
    return(tab)
  }
  else
    return(depChanges)
}

difftable <- function(x,y,levels=unique(c(x,y)))
{
  tab1 <- table(factor(x,levels=levels))
  tab2 <- table(factor(y,levels=levels))
  return(tab1-tab2)
}

countReconstructions <- function(origNet, truncatedNet, reconstructedNets)
{
  origDep <- readDependencies(origNet)
  truncatedDep <- readDependencies(truncatedNet)

  deletions <- mapply(function(o,t)difftable(o,t,origNet$genes), origDep, truncatedDep)
  
  reconst <- mapply(function(o,t)difftable(o,t,origNet$genes), origDep, truncatedDep)
  print(deletions)
}

networkSimilarity <- function(net1, net2)
{
  genes <- c(net1$genes,paste("!",net1$genes,sep=""))
  dep1 <- lapply(readDependencies(net1, removeNegations=FALSE),function(dep)table(factor(dep,levels=genes)))
  dep2 <- lapply(readDependencies(net2, removeNegations=FALSE),function(dep)table(factor(dep,levels=genes)))
  
  intersects <- mapply(function(gene1,gene2)mapply(min,gene1,gene2),dep1,dep2)
  
  return(sum(intersects)/sqrt(sum(unlist(dep1))*sum(unlist(dep2))))
}


asLatexTable <- function(table1, table2, relevantDependencies, simpleComparison=FALSE)
{
  table1 <- t(table1)
  table2 <- t(table2)
  relevantDependencies <- t(relevantDependencies)
  
  res <- 
  paste(sep="",
        "\\begin{footnotesize}\n\\begin{tabularx}{\\linewidth}{|l|", paste(rep("AB|",ncol(table1)), collapse=""), "}\\hline\n",
        paste(paste("  & \\multicolumn{2}{c|}{",gsub("_","\\_",colnames(table1),fixed=TRUE),"}",
                    sep=""),collapse=""),"\\\\\\hline\n",  
        paste(sapply(1:nrow(table1), function(i)
        {
          paste(gsub("_","\\_",rownames(table1)[i],fixed=TRUE),
          paste(mapply(function(e1, e2, rel)
                 {
                  if (simpleComparison)
                  {
                    if (e2 == 0 && e1 != 0)
                      paste("\\cellcolor{green} ",sprintf("%.1f",e1), " & \\cellcolor{green} ", sprintf("%.1f",e2),sep="")
                    else
                    if (abs(e1) - abs(e2) > 0)
                      paste("\\cellcolor{YellowGreen} ",sprintf("%.1f",e1), " & \\cellcolor{YellowGreen} ", sprintf("%.1f",e2),sep="")
                    else
                    if (abs(e1) - abs(e2) < 0)
                      paste("\\cellcolor{RedOrange} ",sprintf("%.1f",e1), " & \\cellcolor{RedOrange} ", sprintf("%.1f",e2),sep="")
                    else
                      paste(sprintf("%.1f",e1), " & ", sprintf("%.1f",e2),sep="")
                  }
                  else
                  if (rel)
                  {
                    if (e1 == e2)#(abs(e1 - e2) < 0.2)
                    {
                      if (e1 == 0)
                        paste("\\cellcolor{cyan} ",sprintf("%.1f",e1), " & \\cellcolor{cyan} ", sprintf("%.1f",e2),sep="")
                      else
                        paste("\\cellcolor{blue} ",sprintf("%.1f",e1), " & \\cellcolor{blue} ", sprintf("%.1f",e2),sep="")
                    }
                    else
                    if (e1 < e2)
                      paste("\\cellcolor{red} ",sprintf("%.1f",e1), " & \\cellcolor{red} ", sprintf("%.1f",e2),sep="")
                    else
                    {
                      if (e2 < 0.005)
                        paste("\\cellcolor{darkgreen} ",sprintf("%.1f",e1), " & \\cellcolor{darkgreen} ", sprintf("%.1f",e2),sep="")
                      else
                        paste("\\cellcolor{green} ",sprintf("%.1f",e1), " & \\cellcolor{green} ", sprintf("%.1f",e2),sep="")
                    }
                  }
                  else
                    paste(sprintf("%.1f",e1), " & ", sprintf("%.1f",e2),sep="")
                 },table1[i,],table2[i,],relevantDependencies[i,]),collapse=" & "), sep=" & ")
        }),collapse="\\\\\\hline\n"),
        "\\\\\\hline\n\\end{tabularx}\n\\end{footnotesize}\n")


  #paste(apply(table,1,function(row)paste(,collapse=" & ")),collapse="\\\\\\hline\n")
}

plotReconstructionHeatmap <- function(truncated, reconstructed, numTruncated, numReconstructed, crucial=NULL, ncolor=100, printPercentage=FALSE, limits=30)
{
  truncated <- round(truncated/numTruncated*100,2)
  reconstructed <- round(reconstructed/numReconstructed*100,2)
  library(RColorBrewer)
  library(gplots)
  pal <- colorRampPalette(c("#D7191C","#D73027","#F46D43","#FFFFFF","#66BD63","#1A9850","#1A9641"))(ncolor-1)
    
 
  mat <- truncated - reconstructed
  cellnote <- t(matrix(nrow=nrow(mat),paste(sprintf("%.1f",truncated), "=>",sprintf("%.1f",reconstructed))))
  
  if (printPercentage)
  { 
    percentage <- t(-((truncated - reconstructed)/truncated) * 100)
    idx <- is.finite(percentage)
    cellnote[idx] <- paste(cellnote[idx],"\n(",sprintf("%+.0f",percentage[idx]),"%)")
  }
  
  cellnote[t(truncated) == 0 & t(reconstructed) == 0] <- ""
  breaks <- seq(-limits,limits,length.out=ncolor)
  print(breaks)
  print(mat)
  
  addExpr <<- function()
  {
    box(col="darkgrey")
    #rect(0.5,0.5,nrow(mat)+0.5,ncol(mat)+0.5,border="darkgrey")  
    grid(nx=nrow(mat),ny=ncol(mat), col="darkgrey", lty="solid")
    if (!is.null(crucial))
    {
      ind <- which(crucial,arr.ind=TRUE)
      ind[,2] <- nrow(mat) - ind[,2] + 1
      print(cbind(ind[,1] - 0.5, ind[,2] - 0.5, ind[,1] + 0.5, ind[,2] + 0.5))
      rect(ind[,1] - 0.5, ind[,2] - 0.5, ind[,1] + 0.5, ind[,2] + 0.5, border="yellow", lwd=2)
    }
  }
  
  try(heatmap.2(t(mat), cellnote=cellnote, notecol="black", notecex=0.75, 
      scale="none", trace="none", breaks=breaks, dendrogram="none", 
      Rowv=FALSE, Colv=FALSE, col=pal, key=FALSE, density.info="n",keysize=0, add.expr = {addExpr()}), silent=F)
      
}

preprocessGeneNames <- function(geneNames)
{
  replacements <- c("Wee1/Mik1","Cdc2/Cdc13","Cdc2/Cdc13*")
  names(replacements) <- c("Wee1_Mik1","Cdc2_Cdc13","Cdc2_Cdc13A")
  
  return(sapply(geneNames,function(gene)
  {
    if (!is.na(replacements[gene]))
      return(replacements[gene])
    else
      return(gene)
  }))  
}

plotReconstructionArcGraph <- function(res, originalNetwork="cellcycleorigin.txt", 
                                       cutAt=0.05, cutMode=c("test","percentage"))
{
   recPerc <- t(res$reconstructed$dependencies/sum(res$reconstructed$usedAttractors))
   truncPerc <- t(res$truncated$dependencies/length(res$truncated$usedAttractors))
   
   changes <- recPerc - truncPerc
   print(changes)
   
   network <- loadNetwork(originalNetwork)
   
   origDependencies <- t(sapply(network$interactions, function(int)
                          {
                            r <- rep(0, length(network$genes))
                            r[int$input] <- 1
                            return(r)                                   
                          }))
   
   edges <- switch(match.arg(cutMode),
    percentage = which(abs(changes) > cutAt, arr.ind=TRUE),
    test = 
    {
        pvals <- sapply(1:ncol(changes),function(i)
             sapply(1:nrow(changes),function(j)
             {
               fisher.test(matrix(c(res$reconstructed$dependencies[i,j],
                                  sum(res$reconstructed$usedAttractors)-
                                      res$reconstructed$dependencies[i,j],
                                  res$truncated$dependencies[i,j],
                                  length(res$truncated$usedAttractors)-
                                      res$truncated$dependencies[i,j]),ncol=2))$p.value
             }))
        pvals <- matrix(p.adjust(pvals,method="bonferroni"),nrow=nrow(pvals))
        print(pvals)
        which(pvals < cutAt, arr.ind=TRUE)
    })
    
   print(edges)
   
   vals <- apply(edges,1,function(x)changes[x[1],x[2]])
   before <- apply(edges,1,function(x)truncPerc[x[1],x[2]])
   after <- apply(edges,1,function(x)recPerc[x[1],x[2]])   
   

   cols <- apply(edges,1,function(x)
                 {
                    if (origDependencies[x[1],x[2]] == 1) 
                        "darkgreen"
                    else
                        "blue"
                 })

   network$genes <- preprocessGeneNames(network$genes)
   edgeNames <- apply(edges,1,function(x)paste(network$genes[x[2]],network$genes[x[1]],sep="->"))
   
   edgeList <- cbind(1:length(vals),1:length(vals) + length(vals))
   #g <- graph.data.frame(edgeList,vertices=as.data.frame(1:(2*length(vals))))
   #print(get.edgelist(g))
   
   #nodeNames <- unlist(lapply(edgeNames,rep,2))
   #nodeCex <- t(apply(edges,1,function(x)
   #           {
   #             round(c(max(truncPerc[x[1],x[2]] * 5,0.5), max(recPerc[x[1],x[2]] * 5,0.5)))
   #           }))
   
   nodeRadius <- t(apply(edges,1,function(x)
              {
                sqrt(c(truncPerc[x[1],x[2]], recPerc[x[1],x[2]])/(4*pi)) * 0.5
              }))
   print(nodeRadius)
   
   nodeCols <- t(sapply(cols,rep,2))
   
   #library(arcdiagram)
   
   par(mar=c(5, 15, 4, 15)+0.1)
   plot(NA,type="n",xlim=c(-0.5,11.5),ylim=c(length(vals),0),xaxt="n",yaxt="n",bty="n",xlab="",ylab="",
        main="Insertion/deletion of regulatory dependencies in the reconstruction")
   axis(side=2,at=1:length(vals),labels=edgeNames, las=1, cex=0.5)
   axis(side=4,at=1:length(vals),labels=edgeNames, las=1, cex=0.5)   
   
   #points(x=c(rep(1,length(vals)),rep(10,length(vals))),
   #       y=rep(1:length(vals),2),pch=19,
   #       cex=as.numeric(nodeCex),
   #       col=as.character(nodeCols))
   symbols(x=c(rep(1,length(vals)),rep(10,length(vals))),
           y=rep(1:length(vals),2),
           circles=as.numeric(nodeRadius),
           fg=as.character(nodeCols),
           bg=as.character(nodeCols),
           inches=FALSE,
           add=TRUE)

   for (i in 1:length(vals))
   {
      #lines(x=c(2,8),y=c(i,i),col=cols[i],lwd=round(vals[i] * 50))
      if (vals[i] > 0)
        polygon(x=c(1.5,9.5,9.5), 
                y=c(i, i + (0.5 * abs(vals[i])), i - (0.5 * abs(vals[i]))),
                col=cols[i],
                border=cols[i])
      else                
        polygon(x=c(1.5,1.5,9.5), 
                y=c(i + (0.5 * abs(vals[i])), i - (0.5 * abs(vals[i])), i),
                col=cols[i],
                border=cols[i])
                
    text(4.5,i-0.5,sprintf("%+.1f%%",vals[i]*100),col=cols[i],cex=0.75)
    text(-0.5,i,sprintf("%.1f%%",before[i]*100),col=cols[i],cex=0.75,pos=4)
    text(11.5,i,sprintf("%.1f%%",after[i]*100),col=cols[i],cex=0.75,pos=2)        
                
   }
   mtext(at=-1,adj=1,"Input networks")
   mtext(at=12,adj=0,"Reconstructed networks")
      
   return(pvals)
   #arcplot(edgeList, 
   #        sorted=TRUE,
   #        labels=nodeNames,
   #        lwd.arcs=round(vals * 50),
   #        col.arcs=cols,
   #        col.nodes=nodeCols,
   #        cex.nodes=nodeCex,
   #        show.nodes=TRUE, pch.nodes=21, bg.nodes="gray90", lwd.nodes=2)
   
   
}

plotDependencies <- function(dependencies, original, comparison, numReconst, numComparison, 
                             doTest=(!missing(comparison)), crucial, colors, cex=0.75, cex.axis=0.75, ...)#(!missing(comparison)))
{
  dependencies <- t(dependencies)
  colnames(dependencies) <- preprocessGeneNames(colnames(dependencies))
  rownames(dependencies) <- preprocessGeneNames(rownames(dependencies))
  
  if (!missing(comparison))
    comparison <- t(comparison)
    
  if (missing(colors))
  {
    if (length(find("pal")) == 0)
    {
      library(RColorBrewer)
      #ncolor <- 201
      #pal <- colorRampPalette(c("#D7191C","#FFFFFF","#1A9641"))(ncolor)
      pal_neg <- unique(colorRampPalette(c("white","#1a4396"))(1000))
      pal_pos <- unique(colorRampPalette(c("white","#18943f"))(1000))
      source("../colorspace.R")
      pal <<- opt_color_space(pal_pos,pal_neg,n=100,distfun=d_lab)
    }
  }
  else
    pal <- colors

  plot(c(),c(),xlim=c(0,ncol(dependencies)),
       ylim=c(0,nrow(dependencies)),xlab="",ylab="",
       axes=FALSE,xaxs="i",yaxs="i", ...)


  axis(1,(1:nrow(dependencies))-0.5,colnames(dependencies), yaxt='s', las=2, tick=FALSE, cex.axis=cex.axis)
  axis(2,(1:nrow(dependencies))-0.5,rownames(dependencies), yaxt='s', las=2, tick=FALSE, cex.axis=cex.axis)

  if (doTest)
  {
    pvals <- sapply(1:ncol(dependencies),function(i)
             sapply(1:nrow(dependencies),function(j)
             {
                #wilcox.test(c(rep(1,dependencies[j,i]),rep(0,numReconst-dependencies[j,i])),
                #            c(rep(1,comparison[j,i]),rep(0,numComparison-comparison[j,i])))$p.value
                fisher.test(matrix(c(dependencies[j,i],numReconst-dependencies[j,i],
                                     comparison[j,i],numComparison-comparison[j,i]),ncol=2))$p.value
                #fisher.test(matrix(c(round(dependencies[j,i]/numReconst*10000),
                #                     round((numReconst-dependencies[j,i])/numReconst*10000),
                #                     round(comparison[j,i]/numComparison*10000),
                #                     round((numComparison-comparison[j,i])/numComparison*10000)),
                #                     ncol=2,byrow=TRUE))$p.value
                #prop.test(matrix(c(dependencies[j,i],numReconst-dependencies[j,i],
                #                   comparison[j,i],numComparison-comparison[j,i]),ncol=2,byrow=TRUE))$p.value
                #p <- prop.test(x=c(comparison[j,i],dependencies[j,i]),
                #               n=c(numComparison, numReconst))$p.value
                #if (is.nan(p))
                #  p <- 1
                  
                #p  
                #fisher.test(x=c(dependencies[j,i],numReconst-dependencies[j,i]),
                #            y=c(comparison[j,i],numComparison-comparison[j,i]))$p.value
             }))
    cat("pvals:\n")
    print(pvals)             
    pvals <- matrix(p.adjust(pvals,method="bonferroni"),nrow=nrow(pvals))
  }
  

  # plot active and inactive states
  for(i in 1:ncol(dependencies))
    for(j in 1:nrow(dependencies))
    {
      val <- dependencies[j,i]/numReconst * 100

      colIdx <- round(val*((length(pal)%/%2) / 100))
            
      if  (i %in% original$interactions[[j]]$input)
        rect(i-1,j-1,i,j,col=pal[colIdx+length(pal)/2],border=NA,lwd=2)
      else
      {
        rect(i-1,j-1,i,j,col=pal[-colIdx+length(pal)/2],border=NA)
      }
      if (missing(comparison))
        text(i-0.5,j-0.5,sprintf("%.1f",val), cex=cex)
      else
      {
        t <- sprintf("%.1f/(%+.1f)",val,val-(comparison[j,i]/numComparison*100))
        if (doTest)
        {
          p <- pvals[j,i]
          
          #if (is.nan(p))
          #  p <- 1
          if(p == 1)
            t2 <- "1"
          else 
            t2 <- formatC(p, format = "e", digits = 2)
          
          t2 <- paste("p=", t2, sep="")
          if (p < 0.01)
            t2 <- paste(t2,"**",sep="")
          else
          if (p < 0.05)
            t2 <- paste(t2,"*",sep="")
        }
        text(i-0.5, j, t, pos=1, cex=cex, font=2)
        text(i-0.5, j-0.5, t2, pos=1, cex=0.8*cex)
      }
    }
    
  abline(h=0:nrow(dependencies), col="darkgrey", lty="solid")
  abline(v=0:ncol(dependencies), col="darkgrey", lty="solid")
  
  #for(i in 1:ncol(dependencies))
  #  for(j in 1:nrow(dependencies))
  #  {
  #    if ((missing(crucial) &&  i %in% original$interactions[[j]]$input) ||
  #        (!missing(crucial) && crucial[i,j]))
  #        rect(i-1,j-1,i,j,col=NA,border="blue",lwd=2)
  #  }
  
}

