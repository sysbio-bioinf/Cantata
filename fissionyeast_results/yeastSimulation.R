library(BoolNet)
library(snowfall)

cantataPath <- "../bin/cantata"

reseed <- function()
{
  seeds <<- round(runif(min=0,max=1000000,n=100))
}

set.seed(123456)
reseed()

source("networkComparison.R")

compareProperties <- function(net)
{
  a <- getAttractors(net)
  load("relevantAttractors.dat")
  load("path.dat")
  attEqual <- sapply(atts, function(att1)
                  {
                    any(sapply(a$attractors, function(att2)
                    {
                      if (ncol(att1$involvedStates) != ncol(att2$involvedStates))
                        return(FALSE)
                      else
                        return(all(att1$involvedStates == att2$involvedStates))
                    }))
                  })
  chainEqual <- identical(path,getPathToAttractor(network=net,state=c(1,0,0,1,1,0,0,1,0,0),includeAttractorStates="first"))
  return(list(attEqual=attEqual,chainEqual=chainEqual))
}

batchStart <- function(resultFolder="result",baseFileFolder="truncatedNetworks",fileName="fissionyeast_trunc",subfolders="3",additionalOptions="",...)
{
  for (folder in subfolders)
  {
    i <- 1
    while(TRUE)
    {
      if (file.exists(paste(resultFolder,"/",folder,"_",i,sep="")))
        i <- i+1
      else
        break
    }
    reseed()
    workFolder <- paste(resultFolder,"/",folder,"_",i,sep="")
    system(paste("mkdir",workFolder))
    system(paste("cp ",baseFileFolder,"/",folder,"/",fileName,"* ",workFolder,sep=""))
    save(seeds,file=paste(workFolder,"/seeds.dat",sep=""))
    runYeastTests(folder=workFolder, additionalOptions=additionalOptions, ...)
  }
}

getRelevantDependencies <- function(networkFile,exPath=cantataPath)
{
  network <- loadNetwork(networkFile,lowercaseGenes=TRUE)
  att1 <- getAttractors(network)
  mat <- sapply(network$genes,function(target)
        {
          sapply(network$genes,function(dependency)
          {
            system(paste(exPath,"--truncate -n", networkFile,
                         "-o tmpnet.txt -tg", target, "-dg", dependency),intern=TRUE)
            n <- loadNetwork("tmpnet.txt",lowercaseGenes=TRUE) 
            n <- fixGenes(n,n$genes,-1)           
            return(!all(unlist(compareProperties(n))))
          })
        })
  system("rm tmpnet.txt")
  return(mat)        
}

delayedBatch <- function(waitFor=cantataPath,
                         waitInterval=60,
                         resultFolder="result",
                         baseFileFolder="perturbedNetworks",fileName="fissionyeast_trunc",
                         subfolders,additionalOptions="")
{
  cat("Waiting for process ",waitFor,"...\n",sep="")
  while(TRUE)
  {
    proc <- system("ps aux",intern=TRUE)
    if (length(proc[grep(waitFor,proc,fixed=TRUE)]) == 0)
      break;
    Sys.sleep(waitInterval)
  }
  cat("Starting batch process...\n")
  batchStart(resultFolder=resultFolder,
             baseFileFolder=baseFileFolder,
             fileName=fileName,
             subfolders=subfolders,
             additionalOptions=additionalOptions)
}

generateYeastNets <- function(nets=1:100,  numModifications=3, insertionProb=0.0, allowEmpty=FALSE, fullDeletion=FALSE, folder="result")
{

    for (i in nets)
    {
      print(i)
      while(TRUE)
      {
        if (allowEmpty)
          system(paste(cantataPath, " --truncate -n fissionyeast_dnf.txt -o ",folder,"/fissionyeast_trunc_",i,".txt -nd ",
                     numModifications," -rs ",round(runif(n=1,min=0,max=10000)),
                     " -ip ",insertionProb,
                     " -ae", 
                     if (fullDeletion){" -fd"}else{""},                     
                     sep=""), intern=T)
        else
          system(paste(cantataPath, " --truncate -n fissionyeast_dnf.txt -o ",folder,"/fissionyeast_trunc_",i,".txt -nd ",
                     numModifications," -rs ",round(runif(n=1,min=0,max=10000)), 
                     " -ip ",insertionProb,
                     if (fullDeletion){" -fd"}else{""},
                     sep=""), intern=T)
        
        net <- loadNetwork(paste(folder,"/fissionyeast_trunc_",i,".txt",sep=""),lowercaseGenes=TRUE)
        net <- fixGenes(net, net$genes, -1)
        eq <- compareProperties(net)
        if (all(unlist(eq) == FALSE))
          break
      }
    }
}

runYeastTests <- function(folder="result", nets=1:100, runs=5, parallel=TRUE, skipExisting=TRUE, additionalOptions="")
{
  load(paste0(folder,"/seeds.dat"))
  if (skipExisting)
  {
    existing <- sapply(nets,function(i)
                       file.exists(paste(folder,"/result_fissionyeast_trunc_",i,".txt",sep="")))
    if (sum(existing) > 0)
    {                       
      cat("Skipping networks with result files:\n")
      print(nets[existing])
      nets <- nets[!existing]
    }             
  }
  sfInit(parallel=parallel, cpus=4)
  sfExport("seeds")
  sfLapply(1:length(nets),function(i)
  {
    print(nets[i])
    system(paste(cantataPath, " --optimize -n ",folder,"/fissionyeast_trunc_",nets[i],
                 ".txt -o ",folder,"/result_fissionyeast_trunc_",nets[i],".txt -r fissionyeast-rules.txt -on ",
                 folder,"/resultnet_fissionyeast_trunc_",nets[i],"_%d.txt -ni 1000 -ns ",runs," -rs ",seeds[i], additionalOptions, sep=""), intern=T)
  })
}

runSingleTest <- function(number)
{
  system(paste(cantataPath, " --optimize -n ",folder,"/fissionyeast_trunc_",number,
                 ".txt -o ",folder,"/result_fissionyeast_trunc_",number,".txt -r fissionyeast-rules.txt -on ",
                 folder,"/resultnet_fissionyeast_trunc_",number,"_%d.txt -ni 1000 -ns 5 -rs ",seeds[number], sep=""), intern=T)
}

countReconstructionsForNets <- function(indices, folder="result")
{
  origFile <- "fissionyeast_dnf.txt"
  origNet <- loadNetwork(origFile,lowercaseGenes=TRUE)
  
  res <- Reduce("+",lapply(indices,function(index)
  {
    truncatedFile <- paste(folder,"/fissionyeast_trunc_",index,".txt",sep="")
    truncatedNet <- loadNetwork(truncatedFile,lowercaseGenes=TRUE)
  
    reconstructedFiles <- list.files(folder,
                                     paste("resultnet_fissionyeast_trunc_",index,"_[0-9]*.txt",sep=""), 
                                     full.names=TRUE)
    
    reconstructedNets <- lapply(reconstructedFiles,function(x)loadNetwork(x,lowercaseGenes=TRUE))

    countReconstructions(origNet, truncatedNet, reconstructedNets)
  }))
  return(res)
}

analyzeYeastResults <- function(numNets=100, folder="result", partial=TRUE)
{
  if (!partial)
    indices <- 1:numNets
  else
  {
    cat("Examined networks:\n")
    indices <- (1:numNets)[sapply(1:numNets,function(i)
                                  file.exists(paste(folder,"/result_fissionyeast_trunc_",i,".txt",sep="")))]
    print(indices)
  }
  cat("Changes in mutated networks:\n")
  files <- paste(folder,"/fissionyeast_trunc_",indices,".txt",sep="")  
  res1 <- compareNetworks("fissionyeast_dnf.txt", files, compareAttractors=FALSE)
  summarizeComparison(res1)
  roundedOrigDependencies <- round(res1$dependencies/numNets * 100, 2)   
  roundedOrigDeletions <- round(res1$deletions/numNets * 100, 2) 
  roundedOrigInsertions <- round(res1$insertions/numNets * 100, 2)  
  
  files <- lapply(indices,function(i)list.files(folder,
                                                  paste("resultnet_fissionyeast_trunc_",i,"_[0-9]*.txt",sep=""), 
                                                  full.names=TRUE))
  
  cat("Networks that could not be reconstructed:\n")
  numReconst <- sapply(files,length)
  print(which(numReconst == 0))
  
  cat("Number of networks that could not be reconstructed:\n")
  print(sum(numReconst == 0))
  
  cat("Networks matching the conditions:\n")
  files <- unlist(files)
  
  matchingParts <- t(sapply(files,function(file)
        {
          net <- loadNetwork(file,lowercaseGenes=TRUE)
          net <- fixGenes(net,1:length(net$genes),-1)
          unlist(compareProperties(net))
        }))
  
  #print(matchingParts)
  matching <- apply(matchingParts,1,function(m)all(unlist(m)))
  
  cat(sum(matching),"/",length(files),"\n",sep="")
  
  cat("Average number of candidates (matching conditions exactly):\n")
  print(sum(matching)/sum(numReconst != 0))
  
  cat("Average number of candidates (all):\n")
  print(length(files)/sum(numReconst != 0))
  
  #cat("Networks not matching the conditions:\n")
  #print(files[!matching])

  cat("Changes in repaired networks:\n")                                                          

  res2 <- compareNetworks("fissionyeast_dnf.txt", files, compareTransitionTables=TRUE, compareAttractors=FALSE)
  
  summarizeComparison(res2)
 
  roundedReconstDependencies <- round(res2$dependencies/sum(res2$usedAttractors) * 100, 2)  
  roundedReconstDeletions <- round(res2$deletions/sum(res2$usedAttractors) * 100, 2) 
  roundedReconstInsertions <- round(res2$insertions/sum(res2$usedAttractors) * 100, 2) 
  relevantDependencies <- getRelevantDependencies("fissionyeast_dnf.txt",cantataPath)
 
  delDiff <- (roundedOrigDeletions - roundedReconstDeletions)

  recovered <- c(sum(roundedReconstDeletions[as.logical(relevantDependencies)] < 0.0005 & 
                   roundedOrigDeletions[as.logical(relevantDependencies)] > 0.0005),
               sum(roundedReconstDeletions[as.logical(relevantDependencies)] >= 0.0005 & 
                   delDiff[as.logical(relevantDependencies)] > 0.0),
               sum(roundedOrigDeletions[as.logical(relevantDependencies)] > 0 & 
                   delDiff[as.logical(relevantDependencies)] == 0),
               sum(delDiff[as.logical(relevantDependencies)] < 0.0))                        

  names(recovered) <- c("Always","Partially","No change","Deletion")
  print(recovered)                      
  
  cat("Deletions:\n")
  cat(asLatexTable(roundedOrigDeletions,roundedReconstDeletions,relevantDependencies))
  cat("Insertions:\n")
  disabledRel <- matrix(FALSE, nrow=nrow(relevantDependencies), ncol=ncol(relevantDependencies))
  cat(asLatexTable(roundedOrigInsertions,roundedReconstInsertions,disabledRel, simpleComparison=T))
  cat("Overall:\n")
  cat(asLatexTable(roundedOrigInsertions - roundedOrigDeletions,roundedReconstInsertions - roundedReconstDeletions,disabledRel, simpleComparison=T))
  return(invisible(list(truncated=res1,reconstructed=res2, crucial=relevantDependencies)))
}

writeSummaries <- function(folders, ...)
{
  #source("../cellcycle/calculateNumberOfCandidates.R")
  return(lapply(folders, function(folder)
  {
    sink(paste(folder,"/summary.txt",sep=""))
    r <- analyzeYeastResults(folder=folder)
    sink()
    
    pdf(paste(folder,"/deletions.pdf",sep=""),width=7,height=5)
    plotReconstructionHeatmap(truncated=r$truncated$deletions, 
                              reconstructed=r$reconstructed$deletions, 
                              crucial=r$crucial,
                              numTruncated=length(r$truncated$usedAttractors), 
                              numReconstructed=sum(r$reconstructed$usedAttractors))
    dev.off()
    
    pdf(paste(folder,"/insertions.pdf",sep=""),width=7,height=5)
    plotReconstructionHeatmap(truncated=r$truncated$insertions,
                              reconstructed=r$reconstructed$insertions, 
                              numTruncated=length(r$truncated$usedAttractors), 
                              numReconstructed=sum(r$reconstructed$usedAttractors), 
                              limit=10)
    dev.off()

    pdf(paste(folder,"/depchanges.pdf",sep=""),width=11,height=12)
    pvals <- plotReconstructionArcGraph(r, originalNetwork="fissionyeast_dnf.txt", 
                                       cutAt=0.05, cutMode="test")
    dev.off()
    
    origFile <- "fissionyeast_dnf.txt"
    origNet <- loadNetwork(origFile,lowercaseGenes=TRUE)
    pdf(paste(folder,"/overview.pdf",sep=""), width=9, height=6)
    par(mar = c(6,5,1,1) + 0.1)
    #par(mfrow=c(1,2))
    #load("res_bestfit_sameIndegree_1000000.dat")
    #plotDependencies(r$truncated$dependencies, original=origNet, 
    #                 numReconst=length(r$truncated$attractorsEqual),
    #                 crucial=r$crucial, 
    #                 main="Perturbed networks", ...)
    plotDependencies(r$reconstructed$dependencies, original=origNet,
                     numReconst=sum(r$reconstructed$usedAttractors),
                     comparison=r$truncated$dependencies,
                     numComparison=length(r$truncated$usedAttractors),
                     crucial=r$crucial,                    
                     ...)
    dev.off()                     
    #pdf("overview_bestfit_1000000.pdf", width=9, height=6)
    #par(mar = c(6,5,1,1) + 0.1)
    #buildAndPlotDependencies(res_bestfit, 
    #                         crucial=r$crucial,...)
    #dev.off()
    
    return(list(res= r, pvals = pvals))
  }))
}

findNetworksWithModification <- function(numNets=100,geneIdx,dependencyIdx,
                                         type=c("deletion","insertion","existence"),
                                         resultFolder="result",silent=FALSE)
{
  origFile <- "fissionyeast_dnf.txt"
  origNet <- loadNetwork(origFile)
  
  if (is.character(geneIdx))
    geneIdx <- which(origNet$genes == geneIdx)

  if (is.character(dependencyIdx))
    dependencyIdx <- which(origNet$genes == dependencyIdx)
  
  modifiedFiles <- unlist(lapply(1:numNets,function(i)list.files(resultFolder,
                                                  paste("resultnet_fissionyeast_trunc_",i,"_[0-9]*.txt",sep=""), 
                                                  full.names=TRUE)))

  modifiedNets <- lapply(modifiedFiles, function(file)loadNetwork(file))
  deletions <- sapply(modifiedNets, function(net)
                      {
                        int1 <- origNet$interactions[[geneIdx]]
                        int2 <- net$interactions[[geneIdx]]
                        s <- switch(type,
                          deletion = setdiff(int1$input, int2$input),
                          insertion = setdiff(int2$input, int1$input),
                          existence = int2$input)
                        return(dependencyIdx %in% s)                          
                      })
  if (!silent)
  {
    res <- compareNetworks("fissionyeast_dnf.txt", modifiedFiles[deletions], compareAttractors=FALSE)
    summarizeComparison(res)

    print(sapply(modifiedNets[deletions],function(net)net$interactions[[geneIdx]]$expression))
  }                    
  return(modifiedFiles[deletions])
}
