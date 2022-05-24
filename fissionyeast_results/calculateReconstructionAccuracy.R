evaluateReconstructionAccuracy <- function(index, folder = "result", originalNet = "fissionyeast_dnf.txt")
{
  original <- loadNetwork(originalNet)
  
  #read reconstructed results
  reconstructedFils <- list.files(folder,
                                  paste("resultnet_fissionyeast_trunc_",index,"_[0-9]*.txt",sep=""), 
                                  full.names=TRUE)
  reconstructedNetworks <- lapply(reconstructedFils, 
                                  function(f) loadNetwork(f))

  #compute accuracy of differences
  poss <- 1:length(original$genes)
  
  #unify results
  unionReconstructed <- Reduce(f = function(new, accum){
    res <- new
    res$interactions <- mapply(function(a,b) {a$input <- union(a$input,b$input); a}, new$interactions, accum$interactions)}
    res)
  
  accuracies <- mapply(function(pert, rec){
    
      TP <- length(intersect(rec$input,pert$input))
      FP <- length(setdiff(rec$input,pert$input))
      TN <- length(intersect(poss, intersect(setdiff(poss,rec$input), setdiff(poss,pert$input))))
      FN <- length(intersect(setdiff(poss,rec$input),orig$input))
    
      return((TP+TN)/(TP+TN+FP+FN))
    #return((FP/(FP+TN)))
    }, original$interactions, unionReconstructed$interactions)
  return(accuracy)                            
}

analyzeYeastReconstructionAccuracy <- function(numNets=100, folder="result", partial=TRUE)
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
  
  return(lapply(indices, function(i) evaluateReconstructionAccuracy(i,folder = folder)))
}


recAccuracy05 <- analyzeYeastReconstructionAccuracy(folder="result/5_depchange_0.25_1/")
print(paste("Mean Accuracy, SD, 05 perturbation experiment:", mean(unlist(recAccuracy05)), ", ", sd(unlist(recAccuracy05)), sep = ""))
recAccuracy07 <- analyzeYeastReconstructionAccuracy(folder="result/7_depchange_0.25_1/")
print(paste("Mean Accuracy, SD, 07 perturbation experiment:", mean(unlist(recAccuracy07)), ", ", sd(unlist(recAccuracy07)), sep = ""))


subtractNetworkDependencies <- function(origNet, newNet)
{
  return(mapply(function(origN, newN) {setdiff(union(origN$input,newN$input), intersect(origN$input,newN$input))}, origNet$interactions, newNet$interactions))
}

evaluateNetworkReconstruction <- function(index, folder = "result", originalNet = "fissionyeast_dnf.txt")
{
  original <- loadNetwork(originalNet)
  
  #load perturbed network
  perturbedFils <- list.files(folder,
                              pattern = paste("^fissionyeast_trunc_",index,".txt",sep=""),
                              full.names=TRUE)
  perturbedNetwork <- lapply(perturbedFils, 
                             function(f) loadNetwork(f))[[1]]
  
  #compute differing edges between perturbed and original
  differenceSetPerturbed <- subtractNetworkDependencies(original, perturbedNetwork)
  
  #read reconstructed results
  reconstructedFils <- list.files(folder,
                                  paste("resultnet_fissionyeast_trunc_",index,"_[0-9]*.txt",sep=""), 
                                  full.names=TRUE)
  if(is.null(reconstructedFils))
    return(NA)
  reconstructedNetworks <- lapply(reconstructedFils, 
                                  function(f) loadNetwork(f))
  
  #measure differences in recon networks to original
  differenceSetsReconstructed <- lapply(reconstructedNetworks, function(n) subtractNetworkDependencies(perturbed, n))
  #extract common differences accros results
  differenceSetReconstructed <- Reduce(f= function(orig, accum) {mapply(intersect, accum, orig)}, differenceSetsReconstructed)
  
  #compute accuracy of differences
  poss <- 1:length(original$genes)
  accuracies <- mapply(function(pert, rec){
    
    TP <- length(intersect(rec,pert))
    FP <- length(setdiff(rec,pert))
    TN <- length(intersect(poss, intersect(setdiff(poss,rec), setdiff(poss,pert))))
    FN <- length(intersect(setdiff(poss,rec),orig))
    
    return((TP+TN)/(TP+TN+FP+FN))
    #return((FP/(FP+TN)))
  }, differenceSetPerturbed, differenceSetReconstructed)
  
  return(accuracies)                            
}

analyzeYeastReconstructionDifferenceAccuracy <- function(numNets=100, folder="result", partial=TRUE)
{
  if (!partial)
    indices <- 1:numNets
  else
  {
    cat("Examined networks:\n")
    indices <- (1:numNets)[sapply(1:numNets,function(i)
      file.exists(paste(folder,"resultnet_fissionyeast_trunc_",i,"_1.txt",sep="")))]
    print(indices)
  }
  
  files <- paste(folder,"/fissionyeast_trunc_",indices,".txt",sep="")  
  
  return(lapply(indices, function(i) evaluateNetworkReconstruction(i,folder = folder)))
}


diffAccuracy05 <- analyzeYeastReconstructionDifferenceAccuracy(folder="result/5_depchange_0.25_1/")
print(paste("Mean Accuracy, SD, 05 perturbation experiment, based on reverting differences between perturbed and original:", mean(unlist(diffAccuracy05)), ", ", sd(unlist(diffAccuracy05)), sep = ""))
diffAccuracy07 <- analyzeYeastReconstructionDifferenceAccuracy(folder="result/7_depchange_0.25_1/")
print(paste("Mean Accuracy, SD, 07 perturbation experiment, based on reverting differences between perturbed and original:", mean(unlist(diffAccuracy07)), ", ", sd(unlist(diffAccuracy07)), sep = ""))


