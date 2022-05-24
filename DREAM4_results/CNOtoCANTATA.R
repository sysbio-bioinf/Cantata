library(CellNOptR)
#source("CellNOptR/R/makeCNOlist.R")

buildNet_extended <- function(model, opt, genes, stimuli=c(), inhibitors=c())
{
  disjunctions <- model$reacID[opt$bString == 1]
  disjunctions <- lapply(disjunctions,function(x)strsplit(x,"=",fixed=TRUE)[[1]])
  targets <- sapply(disjunctions,function(x)x[2])
  factors <- lapply(disjunctions,function(x)
             {
                lit <- strsplit(x[1],"+",fixed=TRUE)[[1]]
                #return(lit)
                gen <- gsub("!","",lit,fixed=TRUE)
                idx <- grep("!",lit,fixed=TRUE)
                neg <- rep(FALSE,length(gen))
                if (length(idx) > 0)
                  neg[idx] <- TRUE
                
                mapply(function(gene,negated)
                {
                  if (gene %in% inhibitors)
                    negated <- !negated
                    
                  if (negated)
                    paste("!",gene,sep="")
                  else
                    gene  
                }, gen, neg, SIMPLIFY=TRUE)
             })

  if (missing(genes))
    genes <- c()
  genes <- unique(c(genes,targets,gsub("!","",unlist(factors),fixed=TRUE)))
  rules <- sapply(genes, function(gene)
  {
    matchingRules <- sapply(factors[targets == gene],paste,collapse=" & ")
    if (length(matchingRules) == 0)
    #  if (gene %in% inhibitors)
    #   paste(gene,", 0",sep="")
    #  else
        paste(gene,", ",gene,sep="")#", 0",sep="")
    else
    #if (gene %in% stimuli)
    #  paste(gene,", ", gene, " | (",paste(matchingRules,collapse=") | ("),")",sep="")
    #else
    if (gene %in% inhibitors)
      paste(gene,", !((",paste(matchingRules,collapse=") | ("),"))",sep="")
    #  paste(gene, ", (t1  & !",gene,") | (",paste(matchingRules,collapse=") | ("),")",sep="")
    else
      paste(gene,", (",paste(matchingRules,collapse=") | ("),")",sep="")
  })
  #rules <- c(rules, paste(stimuli,"_in, 1",sep=""), paste(inhibitors,"_in, 1",sep=""))
  return(rules)
}

MIDASToRules <- function(midasFile, ruleFile, normalize=FALSE, allInputs=FALSE)
{
  #library(Binarization)
  #source("CellNOptR/R/makeCNOlist.R")
  if (is.character(midasFile))
  {
    midas <- readMIDAS(midasFile, verbose=T)
    list <- makeCNOlist(midas, subfield=FALSE)
  }  
  else
  {
    list <- midasFile
  }

  #binThresholds <- apply(midas$dataMatrix[,midas$DVcol,drop=FALSE],2,
  #                      function(x)binarize.BASC(x)@threshold)
  #binarizedMeasurements <- lapply(list$valueSignals,function(m)t(t(m) > binThresholds))
  
  if (normalize)
    list <- normaliseCNOlist(list)
  else
  if ((class(list)=="CNOlist")==FALSE){
        list = CellNOptR::CNOlist(list)
  } 
  binarizedMeasurements <- lapply(list@signals,function(m)m > 0.5)
  #binThresholds <- apply(midas$dataMatrix[,midas$DVcol,drop=FALSE],2,
  #                      function(x)binarize.BASC(x)@threshold)
  #binarizedMeasurements <- lapply(list@signals,function(m)t(t(m) > binThresholds))

  sink(ruleFile)
  for (i in 1:nrow(list@cues))
  {
    cat("Attractor:\nInitial condition:\n")
    stimuli <- list@stimuli[i,]

    stimuli <- stimuli[!is.na(stimuli) & !is.nan(stimuli)]
    
    inhibitors <- list@inhibitors[i,]

    inhibitors <- 1 - inhibitors[!is.na(inhibitors) & !is.nan(inhibitors)]
    
    cues <- c(stimuli, inhibitors)
    
    if (allInputs)
    {
        row <- binarizedMeasurements[[1]][i,]
        row <- row[!is.na(row) & !is.nan(row)]
        
        drop <- c(sapply(names(row),function(n)
                          ((n %in% names(stimuli)) & (stimuli[n] != 0) & (row[n] == 0)) ||
                          ((n %in% names(inhibitors)) & (inhibitors[n] == 0) & (row[n] != 0))))
        if (sum(drop) > 0)
        {
          warning(paste("Dropping inconsistent measurements for genes:",
                        paste(names(drop)[drop],collapse=", ")))
          row <- row[!drop]
        }
        inputs <- c(cues, row)
    }
    else
      inputs <- cues
    
    cat(paste(mapply(function(gene,active)
                {
                  if (active)
                    gene
                  else
                    paste("!",gene,sep="")
                },gsub("-","_",names(inputs),fixed=TRUE), inputs != 0),
                collapse=" "),"\n",sep="")
    cat("Fixed genes:\n")
    cat(paste(c(gsub("-","_",names(stimuli),fixed=TRUE)[stimuli != 0],
                gsub("-","_",paste("!",names(inhibitors),sep=""),fixed=TRUE)[inhibitors == 0]),
                collapse=" "),"\n",sep="")
    cat("State specifications:\n")
    for (meas in binarizedMeasurements[-1])
    {
        row <- meas[i,]
        row <- row[!is.na(row) & !is.nan(row)]
        
        drop <- c(sapply(names(row),function(n)
                          ((n %in% names(stimuli)) & (stimuli[n] != 0) & (row[n] == 0)) ||
                          ((n %in% names(inhibitors)) & (inhibitors[n] == 0) & (row[n] != 0))))
        names(drop) <- names(row)
        if (sum(drop) > 0)
        {
          warning(paste("Dropping inconsistent measurements for genes:",
                        paste(names(drop)[drop],collapse=", ")))
          row <- row[!drop]
        }
        cat(paste(mapply(function(gene,active)
        {
          if (active)
            gene
          else
            paste("!",gene,sep="")
        }, gsub("-","_",names(row),fixed=TRUE), row),
        collapse=" "),"\n",sep="")
    }
    cat("\n")
  }
  sink()
  return(invisible(list))
}

SIFToBoolNet <- function(sifFile, boolnetFile, CNOlist, fixInputs=TRUE, preprocess=TRUE, ignoreAnds=TRUE)
{
  if (preprocess)
  {
    sif <- readSIF(sifFile)
    sif <- preprocessing(data=CNOlist, model=sif)
    system("rm tmp.sif")
    writeSIF(sif,file="tmp.sif")
    sifFile <- "tmp.sif"
  }
  
  
  if ((class(CNOlist)=="CNOlist")==FALSE){
      CNOlist = CellNOptR::CNOlist(CNOlist)
  } 

  sif <- read.table(sifFile,sep="\t",header=FALSE,stringsAsFactors=FALSE)
  if (ncol(sif) != 3)
    sif <- read.table(sifFile,sep=" ",header=FALSE,stringsAsFactors=FALSE)
  
  genes <- unique(c(sif[,1],sif[,3],colnames(CNOlist@cues)))
  genes <- gsub("-","_",genes,fixed=TRUE)
  sif[,1] <- gsub("-","_",sif[,1],fixed=TRUE)
  sif[,3] <- gsub("-","_",sif[,3],fixed=TRUE)  
  #missingGenes <- sort(setdiff(genes,unique(c(cues,signals))))
  #if (length(missingGenes) > 0)
  #  warning(paste("The following network genes were not present in the MIDAS file:",paste(missingGenes,collapse=",")))
  #missingGenes <- sort(setdiff(unique(c(cues,signals)),genes))
  #if (length(missingGenes) > 0)
  #  warning(paste("The following MIDAS genes were not present in the network:",paste(missingGenes,collapse=",")))
  
  #cat("Both cues and signals:\n")
  #print(intersect(cues,signals))
  
  andIdx <- grep("and",genes,fixed=TRUE)
  andIdx <- andIdx[order(genes[andIdx],decreasing=TRUE)]
  ands <- genes[andIdx]

  andStrings <- sapply(ands, function(gene)
  {
    facIdx <- which(sif[,3] == gene)
    if (ignoreAnds && length(facIdx) > 1) 
      return(NA)
      
    if (length(facIdx) > 0)
      paste("(",paste(mapply(function(gene,active)
              {
                if (active == 1)
                  gene
                else
                  paste("!",gene,sep="")
              }, sif[facIdx,1], sif[facIdx,2]), collapse=" & "),")",sep="")
  })
  
  if (ignoreAnds)
    ignoredAnds <- ands[is.na(andStrings)]
  else
    ignoredAnds <- c()
    
  geneStrings <- sapply(genes[-andIdx], function(gene)
  {
    facIdx <- which(sif[,3] == gene & !(sif[,1] %in% ignoredAnds))    
      
    if (length(facIdx) > 0)
      paste(mapply(function(gene,active)
              {
                if (active == 1)
                  gene
                else
                  paste("!",gene,sep="")
              }, sif[facIdx,1], sif[facIdx,2]), collapse=" | ")
    else
    if (fixInputs && (gene %in% colnames(CNOlist@cues)) && 
        !(gene %in% colnames(CNOlist@signals[[1]])))
      paste(gene, "[1.0]", sep="")
    else    
      gene
  })
  
  geneStrings <- sapply(geneStrings, function(s)
  {
    
    for (a in 1:length(andIdx))
    {
      s <- gsub(ands[a], andStrings[a], s, fixed=TRUE)
    }
    s
  })
  genes <- genes[-andIdx]
  print(genes)
  #print(outputStrings)
  sink(boolnetFile)
  cat("targets, factors\n")  
  for (i in 1:length(genes))
  {
    cat(genes[i],", ",geneStrings[i],"\n",sep="")
  }    
  sink()
}

SIFToBoolNet.old <- function(sifFile, boolnetFile, CNOlist, fixInputs=TRUE, preprocess=TRUE)
{
  sif <- read.table(sifFile,sep="\t",header=FALSE,stringsAsFactors=FALSE)
  if (ncol(sif) != 3)
    sif <- read.table(sifFile,sep=" ",header=FALSE,stringsAsFactors=FALSE)
  
  genes <- unique(c(sif[,1],sif[,3]))
  
  if (preprocess)
  {
    genes <- unique(c(preprocessing(data=CNOlist, model=readSIF(sifFile))$namesSpecies,
                      colnames(CNOlist@cues)))
    dropIdx <- which(apply(sif,1,function(row)!(row[1] %in% genes) && (row[3] %in% genes)))
    if (length(dropIdx) > 0)
      warning(paste("Eliminated dependencies:\n",
              paste(apply(sif[dropIdx,,drop=FALSE],1,function(row)paste(row[1],row[3],sep="->")),
                    collapse="\n"),sep=""))
    sif <- sif[-dropIdx,,drop=FALSE]                    
  }
  genes <- gsub("-","_",genes,fixed=TRUE)
  sif[,1] <- gsub("-","_",sif[,1],fixed=TRUE)
  sif[,3] <- gsub("-","_",sif[,3],fixed=TRUE)  
  #missingGenes <- sort(setdiff(genes,unique(c(cues,signals))))
  #if (length(missingGenes) > 0)
  #  warning(paste("The following network genes were not present in the MIDAS file:",paste(missingGenes,collapse=",")))
  #missingGenes <- sort(setdiff(unique(c(cues,signals)),genes))
  #if (length(missingGenes) > 0)
  #  warning(paste("The following MIDAS genes were not present in the network:",paste(missingGenes,collapse=",")))
  
  #cat("Both cues and signals:\n")
  #print(intersect(cues,signals))
    
  sink(boolnetFile)
  cat("targets, factors\n")
  for (gene in genes)
  {
    facIdx <- which(sif[,3] == gene)
    if (length(facIdx) > 0)
      cat(gene,", ",paste(mapply(function(gene,active)
              {
                if (active == 1)
                  gene
                else
                  paste("!",gene,sep="")
              }, sif[facIdx,1], sif[facIdx,2]), collapse=" | "), "\n", sep="")
    else
    if (fixInputs && (gene %in% colnames(CNOlist@cues)) && 
        !(gene %in% colnames(CNOlist@signals[[1]])))
      cat(gene, ", ", gene, "[1.0]\n", sep="")
    else    
      cat(gene, ", ", gene, "\n", sep="")
  }
  sink()
}

simulate_CNO <- function(model, bString, simList=NULL, CNOlist, indexList=NULL)
{

    if ((class(CNOlist)=="CNOlist")==FALSE){
        CNOlist = CellNOptR::CNOlist(CNOlist)
    } 
    if (is.null(simList)==TRUE){
        simList <- prep4sim(model)
    }
    if (is.null(indexList)==TRUE){
        indexList = indexFinder(CNOlist, model)
    }

    # keep simList and indxList for back compatibility ?
    modelCut <- cutModel(model, bString)
    simListCut <- cutSimList(simList, bString)

    # t0
    Sim0 <- simulatorT0(CNOlist=CNOlist, model=modelCut, simList=simListCut, indexList=indexList)
    #Sim0[,indexList$inhibited] <- 1 - Sim0[,indexList$inhibited]
    simRes0 <- as.matrix(Sim0[,indexList$signals])

    # t1
    Sim <- simulatorT1(CNOlist=CNOlist, model=modelCut, simList=simListCut, indexList=indexList)
    #Sim[,indexList$inhibited] <- 1 - Sim[,indexList$inhibited]
    
    colnames(Sim) <- model$namesSpecies
    simRes <- as.matrix(Sim[,indexList$signals])

    sig <- #as.matrix(Sim[,c(indexList$stimulated,indexList$inhibited)])
              as.matrix(Sim[,c(indexList$stimulated,indexList$inhibited)])
    
    simResults <- list(input=CNOlist@cues,
                       t0=simRes0, t1=simRes, trueSig=CNOlist@signals[[2]])

    return(simResults)
}

simulate_CANTATA <- function(CNOlist, network)
{
 if ((class(CNOlist) == "CNOlist") == FALSE) {
        CNOlist = CellNOptR::CNOlist(CNOlist)
    }
  input <- cbind(CNOlist@stimuli,1 - CNOlist@inhibitors)
  colnames(input) <- c(colnames(CNOlist@stimuli), colnames(CNOlist@inhibitors))
    
  t1 <- t(apply(input,1,function(inp)
  {
        
    startState <- rep(0,length(network$genes))
    names(startState) <- network$genes
    startState[names(inp[inp == 1])] <- 1
    
    stim <- inp[colnames(CNOlist@stimuli)]
    inh <- inp[colnames(CNOlist@inhibitors)]

    n <- fixGenes(network,names(stim[stim == 1]),1)
    n <- fixGenes(n,names(inh[inh == 0]),0)
    sequence <- getPathToAttractor(n, startState)
        
    endSeq <- min(nrow(sequence), floor(1.2*length(network$genes)))
    return(as.integer(sequence[endSeq,c(colnames(CNOlist@signals[[1]]),colnames(CNOlist@cues))]))

  }))
  colnames(t1) <- c(colnames(CNOlist@signals[[1]]),colnames(CNOlist@cues))
  return(list(input=cbind(CNOlist@stimuli,CNOlist@inhibitors), 
              t0 = CNOlist@signals[[1]], 
              t1 = t1[,colnames(CNOlist@signals[[1]])],
              trueSig=CNOlist@signals[[2]]))
}

optimizeAndSimulate <- function(midasFile,sifFile,maxGens=5000,stallGenMax=Inf,
                                popSize=200,elitism=popSize/10,sizeFac=1e-04, maxTime=Inf)
{
  model_orig <- readSIF(sifFile)
  list <- makeCNOlist(readMIDAS(midasFile, verbose=TRUE), subfield=FALSE)

  list$namesCues <- gsub("-","_",list$namesCues,fixed=TRUE)
  list$namesSignals <- gsub("-","_",list$namesSignals,fixed=TRUE)  
  
  model_orig$reacID <- gsub("-","_",model_orig$reacID,fixed=TRUE)
  model_orig$namesSpecies <- gsub("-","_",model_orig$namesSpecies,fixed=TRUE)
  rownames(model_orig$interMat) <- gsub("-","_",rownames(model_orig$interMat),fixed=TRUE)
  colnames(model_orig$interMat) <- gsub("-","_",colnames(model_orig$interMat),fixed=TRUE)
  rownames(model_orig$notMat) <- gsub("-","_",rownames(model_orig$notMat),fixed=TRUE)
  colnames(model_orig$notMat) <- gsub("-","_",colnames(model_orig$notMat),fixed=TRUE)
  
  model <- preprocessing(list, model_orig)
  initBstring<-rep(1,length(model$reacID))
  t <- system.time(opt<-gaBinaryT1(
           CNOlist=list,
           model=model,
           sizeFac=sizeFac,
           initBstring=initBstring,
           maxGens=maxGens, popSize=popSize,  elitism=elitism, 
           stallGenMax=stallGenMax,
           maxTime=maxTime, verbose=TRUE))

  simResults <- simulate_CNO(model=model,#_orig,
                             CNOlist=list,
                             bString=opt$bString)
  cat("Elapsed time:\n")
  print(t)                             
  return(list(model_orig=model_orig, model=model, opt=opt, CNOlist=list, simulation=simResults))
}

evalLiverDREAM <- function()
{
  for (i in 1:5)
  {
    res <- optimizeAndSimulate(midasFile="MD-LiverDREAM.csv",
                               sifFile="PKN-LiverDREAM.sif",
                               maxGens=1000)
    save(res,file=paste("cno_result/LiverDREAM/cnoresult_",i,".dat",sep=""))
    res <- optimizeAndSimulate(midasFile="MD-LiverDREAM.csv",
                               sifFile="PKN-LiverDREAM.sif",
                               maxGens=1000,
                               sizeFac=0)
    save(res,file=paste("cno_result/LiverDREAM/cnoresult_sizefac0_",i,".dat",sep=""))
  }
}

findBestRatedNetwork <- function(folder="cno_result/LiverDREAM",mask="cnoresult_[0-9]*.dat",
                                 sizeFac = 1e-04,NAFac = 1)
{
  f <- list.files(folder,
                  mask,
                  full.names=TRUE)
  
  results <- lapply(f,function(n)
                    {
                      load(n)
                      return(res)
                    })
  scores <- sapply(results, int_getfit, sizeFac=sizeFac, NAFac=NAFac)
  bestIdx <- which.min(scores)
  cat(sprintf("Best score: %.4f for file %s\n", scores[bestIdx],f[bestIdx]))
  return(results[[bestIdx]])
}

int_getfit <- function(CNOresult,sizeFac = 1e-04,NAFac = 1)
{
  model <- cutModel(CNOresult$model, CNOresult$opt$bString)
  simList<-prep4sim(CNOresult$model)
  simList <- cutSimList(simList, CNOresult$opt$bString)
  indices <-indexFinder(CNOresult$CNOlist,CNOresult$model,verbose=FALSE)
  simResults<-simulatorT1(
         CNOlist=CNOresult$CNOlist,
         model=model,
         simList=simList,
         indexList=indices)
    Score<-getFit(
         simResults=simResults,
         CNOlist=CNOresult$CNOlist,
         model=model,
         indexList=indices,
         timePoint="t1",
         sizeFac=sizeFac,
         NAFac=NAFac, 
     nInTot=length(which(model$interMat == -1)))
  return(Score)
}

getFit1 <- function (simResults, CNOlist, model, indexList, timePoint = c("t1", 
    "t2"), sizeFac = 1e-04, NAFac = 1, nInTot, simResultsT0 = NA) 
{
    if ((class(CNOlist) == "CNOlist") == FALSE) {
        CNOlist = CellNOptR::CNOlist(CNOlist)
    }
    simResults <- simResults[, indexList$signals]
    if (timePoint == "t1") {
        tPt <- 2
    }
    else {
        if (timePoint == "t2") {
            tPt <- 3
        }
        else {
            tPt <- timePoint
        }
    }
    if (tPt == 2 && is.na(simResultsT0) == FALSE) {
        Diff0 <- simResultsT0[, indexList$signals] - CNOlist@signals[[1]]
        Diff <- simResults - CNOlist@signals[[tPt]]
        r0 <- Diff0^2
        r <- Diff^2
        r <- rbind(r0, r)
        deviationPen <- sum(r[!is.na(r)])/2
    }
    else {

        Diff <- simResults - CNOlist@signals[[tPt]]
        r <- Diff^2
        deviationPen <- sum(r[!is.na(r)])
    }
    print(mean(Diff ^ 2,na.rm=T))
    print(deviationPen)
    NAPen <- NAFac * length(which(is.na(simResults)))
    nDataPts <- dim(CNOlist@signals[[tPt]])[1] * dim(CNOlist@signals[[tPt]])[2]
    nInputs <- length(which(model$interMat == -1))
    sizePen <- (nDataPts * sizeFac * nInputs)/nInTot
    score <- deviationPen + NAPen + sizePen
    return(score)
}

prepareValidationSet <- function(trainSet="MSB2009-MIDAS-RawDataTraining.csv", validationSet="MSB2009-MIDAS-RawDataValidation.csv")
{
  trainCNOlist <- normaliseCNOlist(makeCNOlist(readMIDAS(trainSet),subfield=F))
  
  valCNOlist <- normaliseCNOlist(makeCNOlist(readMIDAS(validationSet),subfield=F))
  
  trainCueNames <- colnames(trainCNOlist@cues)
  valCueNames <- colnames(valCNOlist@cues)
  
  trainSignalNames <- colnames(trainCNOlist@signals[[1]])
  valSignalNames <- colnames(valCNOlist@signals[[1]])
    
  removeCueCols <- setdiff(valCueNames,trainCueNames)
  #removeSignalCols <- setdiff(valSignalNames,trainSignalNames)
  commonCues <- intersect(trainCueNames,valCueNames)
  
  commonSignals <- intersect(trainSignalNames,valSignalNames)
  
  removeRows <- which(apply(valCNOlist@cues[,removeCueCols,drop=FALSE],1,function(x)any(!is.na(x) & x>0)))
  
  newCues <- valCNOlist@cues[-removeRows,commonCues,drop=FALSE]
  valCNOlist@cues <- matrix(nrow=nrow(valCNOlist@cues) - length(removeRows),ncol=length(trainCueNames),0)
  colnames(valCNOlist@cues) <- trainCueNames
  valCNOlist@cues[,commonCues] <- newCues
  
  valCNOlist@inhibitors <- valCNOlist@cues[,colnames(trainCNOlist@inhibitors),drop=FALSE]
  valCNOlist@stimuli <- valCNOlist@cues[,colnames(trainCNOlist@stimuli),drop=FALSE]  
  
  for (i in 1:length(valCNOlist@signals))
  {
    sig <- matrix(nrow=nrow(valCNOlist@signals[[i]]) - length(removeRows),ncol=length(trainSignalNames),NA)
    colnames(sig) <- trainSignalNames
    sig[,commonSignals] <- valCNOlist@signals[[i]][removeRows,commonSignals,drop=FALSE]
    valCNOlist@signals[[i]] <- sig
  }
  return(valCNOlist)
}

crossValidate <- function(sifFile, midasFile, func, ntimes = 10, nfold = 10, cpus=4, outputDir, preprocess=TRUE, fixInputs=TRUE, seed=1234567, ...)
{
  system(paste("mkdir -p ",outputDir,sep=""))
  sink(paste(outputDir,"/seed.txt",sep=""))
  cat(seed)
  sink()
  library(CellNOptR)
  CNOl <- makeCNOlist(readMIDAS(midasFile),subfield=FALSE)
  if ((class(CNOl)=="CNOlist")==FALSE){
    CNOl = CellNOptR::CNOlist(CNOl)
  }
  model_orig <- readSIF(sifFile)
  SIFToBoolNet(sifFile, "model_draft.txt", CNOl, fixInputs=fixInputs, preprocess=preprocess)
  SIFToBoolNet(sifFile, paste(outputDir,"/model_draft.txt",sep=""), 
               CNOl, fixInputs=FALSE, preprocess=preprocess)
  
  model_orig$reacID <- gsub("-","_",model_orig$reacID,fixed=TRUE)
  model_orig$namesSpecies <- gsub("-","_",model_orig$namesSpecies,fixed=TRUE)
  rownames(model_orig$interMat) <- gsub("-","_",rownames(model_orig$interMat),fixed=TRUE)
  colnames(model_orig$interMat) <- gsub("-","_",colnames(model_orig$interMat),fixed=TRUE)
  rownames(model_orig$notMat) <- gsub("-","_",rownames(model_orig$notMat),fixed=TRUE)
  colnames(model_orig$notMat) <- gsub("-","_",colnames(model_orig$notMat),fixed=TRUE)  
  
  set.seed(seed)
  numSamples <- nrow(CNOl@cues)
  splits <- lapply(1:ntimes,function(run)
  # for each run
  {
    # calculate folds
    permut <- sample(1:numSamples, numSamples,replace=FALSE)
    indices <- lapply(1:nfold, function(i)
    {
	    # split the samples in nfold groups
	    permut[seq(i, numSamples, nfold)]
    })
  })
  library(snowfall)
  sfInit(parallel=(cpus>1),cpus=cpus)
  sfExportAll()
  seeds <- round(runif(0,.Machine$integer.max,n=ntimes*nfold))
  
  if (cpus > 1)
    sfExport("splits","CNOl","model_orig")
  
  sfLibrary(CellNOptR)
  sfLibrary(BoolNet)
  for (i in 1:ntimes)
  {
    for (j in 1:nfold)
    {
      system(paste("mkdir -p ",outputDir,"/",i,"_",j,sep=""))
    }
  }
  partitions <- expand.grid(1:ntimes,1:nfold)
  partitions <- cbind(partitions,seeds)
  partitions <- as.list(as.data.frame(t(partitions)))

  res <- sfClusterApplyLB(partitions,function(part,...)
  {
    testIndices <- splits[[part[1]]][[part[2]]]
    trainIndices <- (1:numSamples)[-testIndices]
    CNOlistTrain <- CNOl
    CNOlistVal <- CNOl
    CNOlistTrain@cues <- CNOlistTrain@cues[trainIndices,,drop=FALSE]
    CNOlistTrain@stimuli <- CNOlistTrain@stimuli[trainIndices,,drop=FALSE]
    CNOlistTrain@inhibitors <- CNOlistTrain@inhibitors[trainIndices,,drop=FALSE]
    
    CNOlistVal@cues <- CNOlistVal@cues[testIndices,,drop=FALSE]
    CNOlistVal@stimuli <- CNOlistVal@stimuli[testIndices,,drop=FALSE]
    CNOlistVal@inhibitors <- CNOlistVal@inhibitors[testIndices,,drop=FALSE]

    for (i in 1:length(CNOlistTrain@signals))
    {
      CNOlistTrain@signals[[i]] <- CNOlistTrain@signals[[i]][trainIndices,,drop=FALSE]
      CNOlistVal@signals[[i]] <- CNOlistVal@signals[[i]][testIndices,,drop=FALSE]
    }  
    func(resultFolder=paste(outputDir,"/",part[1],"_",part[2],sep=""),
         CNOlistTrain=CNOlistTrain,
         CNOlistVal=CNOlistVal,
         model=model_orig,
         seed=part[3],
         ...)
  },...)
  save(splits,file=paste(outputDir,"/splits.dat",sep=""))
  save(res,file=paste(outputDir,"/result.dat",sep=""))
  return(res)
}

cv_CNO <- function(resultFolder,CNOlistTrain,CNOlistVal,model,seed,
                   maxGens=1000,stallGenMax=Inf,
                   popSize=200,elitism=popSize/10,sizeFac=1e-04, NAFac=1, maxTime=Inf, 
                   numStarts = 5, ...)
{
  if (!file.exists(paste(resultFolder,"/cnoresult.dat",sep="")))
  {
    set.seed(seed)
    model_cut <- preprocessing(CNOlistTrain, model)
    initBstring<-rep(1,length(model_cut$reacID))
    res <- c()
    fit <- Inf
    for (i in 1:numStarts)
    {
      opt<-gaBinaryT1(
               CNOlist=CNOlistTrain,
               model=model_cut,
               initBstring=initBstring,
               maxGens=maxGens, 
               popSize=popSize, elitism=elitism, 
               stallGenMax=stallGenMax,
               sizeFac=sizeFac,
               NAFac=NAFac,
               maxTime=maxTime, 
               verbose=TRUE,
               ...)

      simResults_train <- simulate_CNO(model=model_cut,
                                       CNOlist=CNOlistTrain,
                                       bString=opt$bString)
                                                 
      r <- list(model_orig=model, model=model_cut, opt=opt, 
                CNOlist=CNOlistTrain, 
                simulation=simResults_train)
      f <- int_getfit(r, sizeFac=sizeFac, NAFac=NAFac)
      if (f < fit)
      {
        simResults_val <- simulate_CNO(model=model_cut,
                                       CNOlist=CNOlistVal,
                                       bString=opt$bString)
        res <- list(model_orig=model, model=model_cut, opt=opt, 
                    CNOlistTrain=CNOlistTrain, 
                    CNOlistVal=CNOlistVal,  
                    simulation=simResults_val)
        fit <- f
      }
    }  
    save(res,file=paste(resultFolder,"/cnoresult.dat",sep=""))
  }
  else
    load(paste(resultFolder,"/cnoresult.dat",sep=""))
  return(res$simulation)            
}

cv_CANTATA <- function(resultFolder,CNOlistTrain,CNOlistVal,model,seed,
                       popSize=100,numIter=1000,maxStates=20,maxTransitions=1000,numStarts=5,
                       normalize=FALSE,allInputs=FALSE,additionalParams="",...)
{
  if (!file.exists(paste(resultFolder,"/result.txt",sep="")))
  {
    if(!file.exists( paste(resultFolder,"/rules.txt",sep="")))
    {
        MIDASToRules(CNOlistTrain, paste(resultFolder,"/rules.txt",sep=""), 
                                  normalize=normalize, allInputs=allInputs)
    }
    system(paste("../Release/NetworkAnalyzer --optimize -n model_draft.txt -r ",resultFolder,"/rules.txt ",
                 "-ps ",popSize," -ni ",numIter," -mt ",maxTransitions," -ms ",maxStates," -ns ",numStarts,
                 " -o ",resultFolder,"/result.txt -on ",
                 resultFolder,"/resultnet_%d.txt -rs ",seed," -me 1",additionalParams,sep=""))
  }
  network <- loadNetwork(paste(resultFolder,"/resultnet_1.txt",sep=""))
  network <- fixGenes(network, network$genes, -1)
  simResults <- simulate_CANTATA(CNOlistVal, network)
  res <- list(network=network, 
              CNOlistTrain=CNOlistTrain, 
              CNOlistVal=CNOlistVal,  
              simulation=simResults)
  save(res,file=paste(resultFolder,"/cantataresult.dat",sep=""))
  return(simResults)
}

evaluateCV <- function(folder, mode=c("mse","roundmse"), fillNAs=TRUE, training=FALSE, reparse=FALSE)
{
  if (file.exists(paste(folder,"/result.dat",sep="")) && !reparse && !training)
  {
    load(paste(folder,"/result.dat",sep=""))
  }
  else
  if (file.exists(paste(folder,"/result_train.dat",sep="")) && !reparse && training)
  {
    load(paste(folder,"/result_train.dat",sep=""))
  }
  else
  {
    folders <- list.files(path=folder, include.dirs=TRUE, full.names=TRUE)
    folders <- folders[file.info(folders)$isdir]
    res <- lapply(folders, function(folder)
    {
      if (file.exists(paste(folder,"/cantataresult.dat",sep="")))
      {
        load(paste(folder,"/cantataresult.dat",sep=""))
        
        if (!training)
          sim <- res$simulation
        else
          sim <- simulate_CANTATA(res$CNOlistTrain, res$network)
      }
      else
      if (file.exists(paste(folder,"/cnoresult.dat",sep="")))
      {
        load(paste(folder,"/cnoresult.dat",sep=""))
        
        if (!training)
          sim <- res$simulation
        else
          sim <- simulate_CNO(CNOlist=res$CNOlistTrain, 
                              model=res$model, 
                              bString=res$opt$bString)
      }
      else
        return(NULL)
      
      return(sim)

    })
    idx <- sapply(res,is.null)
    cat("Missing results: \n")
    print(folders[idx])
    res <- res[!idx]
  }
  if (!file.exists(paste(folder,"/result_train.dat",sep="")) && training)
    save(res, file=paste(folder,"/result_train.dat",sep=""))
  return(evaluateCVresults(res, mode=mode, fillNAs=fillNAs))
}

evaluateCVresults <- function(res, mode=c("mse","roundmse"), fillNAs=TRUE)
{
  errors <- unlist(sapply(res,function(x)
            {
              if (fillNAs)
              {
                naPos <- which(is.na(x$t1))
                x$t1[naPos] <- sample(c(0,1),size=length(naPos),replace=TRUE)
              }
              switch(match.arg(mode, c("mse","roundmse")),
                mse = as.numeric((x$t1-x$trueSig)^2),
                roundmse = as.numeric((x$t1-round(x$trueSig))^2)
              )
            }))
  return(list(error=errors,
              mean=mean(errors,na.rm=T),
              median=median(errors,na.rm=T)))
}

evaluateDependencies <- function(folder,excludeSelf=FALSE)
{
  #source("../networkComparison.R")
  folders <- list.files(path=folder, include.dirs=TRUE, full.names=TRUE)
  folders <- folders[file.info(folders)$isdir]
  res <- lapply(folders, function(folder)
  {
    if (file.exists(paste(folder,"/cantataresult.dat",sep="")))
    {
       load(paste(folder,"/cantataresult.dat",sep=""))
      net <- loadNetwork(paste(folder,"/resultnet_1.txt",sep=""))
    }
    else
    if (file.exists(paste(folder,"/cnoresult.dat",sep="")))
    {
      load(paste(folder,"/cnoresult.dat",sep=""))
      sink("tmpnet.txt")
      cat("targets,factors\n")
      cat(paste(buildNet_extended(res$model,res$opt,
                                stimuli=colnames(res$CNOlistTrain@stimuli),
                                inhibitors=colnames(res$CNOlistTrain@inhibitors)),collapse="\n"),
                                "\n",sep="")
      sink()                                
      net <- loadNetwork("tmpnet.txt")                    
    }
    else
      return(NULL)
    
    if (excludeSelf)
      keep <- !sapply(1:length(net$interactions),function(i)
                     {
                       int <- net$interactions[[i]]
                       return(length(int$input) == 1 && int$input[1] == i)
                     })
    else
      keep <- rep(TRUE,length(net$genes))
    
    return(c("genes"=length(net$interactions),
             "dependencies"=sum(sapply(net$interactions[keep],function(x)length(x$input)))))
  })
  res <- res[!sapply(res,is.null)]
  return(c("genes"=mean(sapply(res,function(x)x["genes"])),
           "dependencies"=mean(sapply(res,function(x)x["dependencies"]))))
}

calculateError <- function(simResults, CNOlist, rounded=FALSE, fillNAs=TRUE)
{
  if (missing(CNOlist))
    signals <- simResults$trueSig
  else
  {
    if ((class(CNOlist)=="CNOlist")==FALSE){
      CNOlist = CellNOptR::CNOlist(CNOlist)
    } 
    signals <- CNOlist@signals[[2]]
  }
  
  if (fillNAs)
  {
    naPos <- which(is.na(simResults$t1))
    simResults$t1[naPos] <- sample(c(0,1),size=length(naPos),replace=TRUE)
  }
  
  if (rounded)
    error <- (simResults$t1 -  round(signals[,colnames(simResults$t1)]))^2
  else    
    error <- (simResults$t1 -  signals[,colnames(simResults$t1)])^2
  rowwiseError <- apply(error,1,mean,na.rm=T)
  meanError <- mean(error,na.rm=T)
  return(list(error=error,rowwiseError=rowwiseError,
              meanError=meanError))
}

calculateError_CANTATA <- function(network, CNOlist, rounded=FALSE)
{
  sim <- simulate_CANTATA(CNOlist=CNOlist, network=network)
  return(calculateError(sim,CNOlist,rounded))
}

calculateError_CNO <- function(CNOresult, CNOlist, rounded=FALSE)
{
   sim <- simulate_CNO(res$model, res$opt$bString, CNOlist=CNOlist)
   return(calculateError(sim,CNOlist,rounded))
}

summarizePredictions <- function(cvRes, CNOlist)
{
  if ((class(CNOlist)=="CNOlist")==FALSE)
      CNOlist = CellNOptR::CNOlist(CNOlist)
  keys <- apply(CNOlist@cues[,colnames(cvRes[[1]]$input),drop=FALSE], 1, paste, collapse="")

  values <- lapply(keys, function(key)
             {
                res <- c()
                for (i in 1:length(cvRes))
                {
                  for (j in 1:nrow(cvRes[[i]]$input))
                  {
                    if (paste(cvRes[[i]]$input[j,],collapse="") == key)
                      res <- rbind(res, cvRes[[i]]$t1[j,,drop=FALSE])
                  }
                }
                return(res)                
             })
       
  values <- t(sapply(values,function(mat)apply(mat,2,mean)))
  return(list(input=CNOlist@cues,t1=values,trueSig=CNOlist@signals[[2]]))
}

plotPredictions <- function(simResults, CNOlist)
{
  error <- abs(simResults$t1 - simResults$trueSig)
  error[is.na(error) | is.nan(error)] <- 0

  if (!missing(CNOlist))
  {
    if ((class(CNOlist)=="CNOlist")==FALSE)
      CNOlist = CellNOptR::CNOlist(CNOlist)
    
    simResults$input <- simResults$input[,c(colnames(CNOlist@stimuli),
                                            colnames(CNOlist@inhibitors)),drop=FALSE]
    
  }

  library(RColorBrewer)
  library(gplots)

  pal <- colorRampPalette(c("#FFFFFF","#FFB5B5","#FF0000"))(101)
         #colorRampPalette(c("#FFFFFF","#FFB5B5","#FF0000"))(101)
         #colorRampPalette(c("#FFFFFF","#FF0000"))(101)
  pal <- sapply(pal, function(col)
         {
          rgb <- col2rgb(col)
          return(rgb(rgb["red",1],rgb["green",1],rgb["blue",1],50,maxColorValue=255))
         })

  plot(NA,NA,type="n",xlim=c(0,ncol(simResults$input) + ncol(simResults$t1) + 0.5), 
                      ylim=c(0,nrow(simResults$input)),xaxt="n",yaxt="n",xaxs="i",yaxs="i",
                      xlab="",ylab="",bty="n",frame.plot=FALSE,ann=FALSE)

  for (i in 1:ncol(simResults$input))
  {
    for (j in 1:nrow(simResults$input))
    {
      rect(i-1,j-1,i,j,col=c("black","white")[as.integer(simResults$input[j,i])],border="gray")
    }
    text(y=-0.2,x=(1:ncol(simResults$input)) - 0.5, labels=colnames(simResults$input), 
         srt=90, adj=1, xpd=TRUE, cex=0.75)
  }
  
  offset <- ncol(simResults$input) + 0.5
  
  if (!missing(CNOlist))
  {
    text(y=nrow(simResults$input)+0.2,x=c(0,ncol(CNOlist@stimuli),offset),
         labels=c("Stimuli","Inhibitors","Simulated vs. measured steady states"),
         xpd=TRUE,adj=c(0,0),cex=0.75)
    abline(v=ncol(CNOlist@stimuli),lwd=2,col="gray")
  }
  else
    text(y=nrow(simResults$input)+0.2,x=c(0,offset),
         labels=c("Perturbations","Simulated vs. measured steady states"),
         xpd=TRUE,adj=c(0,0),cex=0.75)
  for (i in 1:ncol(simResults$t1))
  {
    for (j in 1:nrow(simResults$t1))
    {
      rect(offset + i - 1, j - 1, offset + i, j,col=pal[round(error[j,i] * 100)+1],border=NA)
      rect(offset + i - 0.9, j - 0.9, offset + i - 0.5, j - 0.9 + simResults$t1[j,i]*0.8,
           col="gray40",border=NA)
      rect(offset + i - 0.5, j - 0.9, offset + i - 0.1, j - 0.9 + simResults$trueSig[j,i]*0.8,
           col="gray60",border=NA)
      rect(offset + i - 1, j - 1, offset + i, j,col=NA,border="darkgrey")
      #rect(offset + i - 0.95, j - 0.95, offset + i - 0.05, j - 0.05,col=NA,
      #    border=pal[round(error[j,i] * 100)+1], lwd=2)
    }
    text(y=-0.2,x=(1:ncol(simResults$t1)) + offset - 0.5, labels=colnames(simResults$t1), 
         srt=90, adj=1, xpd=TRUE, cex=0.75)
  }
  rect(0,0,
       ncol(simResults$input) + ncol(simResults$t1) + 0.5,
       nrow(simResults$input),lwd=2,border="gray",col=NA,xpd=TRUE)
}


plotCV <- function(folder, CNOlist, training=FALSE)
{
  if (training)
    load(paste(folder,"/result_train.dat",sep=""))
  else
    load(paste(folder,"/result.dat",sep=""))
  sumRes <- summarizePredictions(res, CNOlist)
  
  if (training)
    pdf(paste(folder,"/predictions_train.pdf",sep=""))
  else
    pdf(paste(folder,"/predictions.pdf",sep=""))
  plotPredictions(sumRes,CNOlist)
  dev.off()
}

plotCVDependencies <- function(folder, origNet=paste(folder,"/model_draft.txt",sep=""))
{
  library(BoolNet)
  source("../networkComparison.R")
  if (is.character(origNet))
    origNet <- loadNetwork(origNet)
  folders <- list.files(path=folder, include.dirs=TRUE, full.names=TRUE)
  folders <- folders[file.info(folders)$isdir]
  dependencies <- lapply(folders, function(f)
  {
    net <- loadNetwork(paste(f,"/resultnet_1.txt",sep=""))
    sapply(net$interactions, function(int)
    {
      r <- rep(0, length(origNet$genes))
      r[int$input] <- 1
      return(r)                                   
    })
  })
  dependencies <- Reduce(function(a,b)a+b, dependencies)
  rownames(dependencies) <- origNet$genes
  colnames(dependencies) <- origNet$genes
  pdf(paste(folder,"/overview.pdf",sep=""))
  plotDependencies(dependencies, original=origNet, numReconst=length(folders), cex.axis=0.5, cex=0.5)
  dev.off()
}

findNetworksWithModification <- function(geneIdx,dependencyIdx, 
                                         type=c("deletion","insertion","existence"),
                                         resultFolder="LiverDREAM",silent=FALSE)
{
  source("../networkComparison.R")
  origFile <- paste(resultFolder,"/model_draft.txt",sep="")
  origNet <- loadNetwork(origFile)
  
  if (is.character(geneIdx))
    geneIdx <- which(origNet$genes == geneIdx)

  if (is.character(dependencyIdx))
    dependencyIdx <- which(origNet$genes == dependencyIdx)
  
  folders <- list.files(path=resultFolder, include.dirs=TRUE, full.names=TRUE)
  folders <- folders[file.info(folders)$isdir]
  modifiedFiles <- paste(folders,"/resultnet_1.txt",sep="")

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
    res <- compareNetworks(paste(resultFolder,"/model_draft.txt",sep=""), 
                           modifiedFiles[deletions], compareAttractors=FALSE)
    summarizeComparison(res)

    #print(sapply(modifiedNets[deletions],function(net)net$interactions[[geneIdx]]$expression))
    print(sapply(modifiedNets[deletions],function(net)
          {
            c("expression"=net$interactions[[geneIdx]]$expression,
              "deleted"=paste(origNet$genes[setdiff(origNet$interactions[[geneIdx]]$input,
                                       net$interactions[[geneIdx]]$input)], collapse=", "),
              "inserted"=paste(origNet$genes[setdiff(net$interactions[[geneIdx]]$input,
                                      origNet$interactions[[geneIdx]]$input)], collapse=", "))
          }))
  }                    
  return(modifiedFiles[deletions])
}



start_liverdream <- function()
{
  #res_dream_cant <<- crossValidate(sifFile="PKN-LiverDREAM.sif",
  #                            midasFile="MD-LiverDREAM.csv",
  #                            outputDir="LiverDREAM/cv_cant",
  #                            ntimes=10,nfold=10,cpus=6,preprocess=TRUE,
  #                            numStarts=1,maxStates=100,func=cv_CANTATA,fixInputs=FALSE,seed=42811)
  #res_dream_cno <<- crossValidate(sifFile="PKN-LiverDREAM.sif",
  #                            midasFile="MD-LiverDREAM.csv",
  #                            outputDir="LiverDREAM/cv_cno",
  #                            fixInputs=FALSE,
  #                            ntimes=10,nfold=10,cpus=4,preprocess=TRUE,
  #                            numStarts=1,func=cv_CNO,seed=42811)
  res_dream_cant <<- crossValidate(sifFile="PKN-LiverDREAM.sif",
                              midasFile="MD-LiverDREAM.csv",
                              outputDir="LiverDREAM/cv_cant_3",
                              ntimes=10,nfold=10,cpus=6,preprocess=TRUE,
                              numStarts=1,maxStates=200,func=cv_CANTATA,
                              fixInputs=TRUE,additionalParams=" -tw 0.25/0.5/0.25",
                              seed=113224)
  #res_dream_cno <<- crossValidate(sifFile="PKN-LiverDREAM.sif",
  #                            midasFile="MD-LiverDREAM.csv",
  #                            outputDir="LiverDREAM/cv_cno_2",
  #                            fixInputs=FALSE,
  #                            ntimes=10,nfold=10,cpus=4,preprocess=TRUE,
  #                            numStarts=1,func=cv_CNO,seed=113224)
}

preprocess_liverdream_raw <- function(reorder=FALSE,dropTP=180,...)
{
  CNOl <- makeCNOlist(readMIDAS("LiverDataDREAMMIDAS.csv"), subfield=FALSE)
  CNOl$namesStimuli <- tolower(CNOl$namesStimuli)
  CNOl$namesInhibitors <- tolower(CNOl$namesInhibitors)  
  CNOl$namesCues <- tolower(CNOl$namesCues)
  CNOl$namesSignals <- tolower(CNOl$namesSignals)
  
  if (length(dropTP) > 0)
  {
    keepIdx <- !(CNOl$timeSignals %in% dropTP)
    CNOl$timeSignals <- CNOl$timeSignals[keepIdx]
    CNOl$valueSignals <- CNOl$valueSignals[keepIdx]
  }
  CNOl <- normaliseCNOlist(CNOl, ...)
  
  if (reorder)
  {
    data(CNOlistDREAM)
    orderStrings1 <- apply(CNOl@cues[,CNOlistDREAM$namesCues],1,paste,collapse="")
    orderStrings2 <- apply(CNOlistDREAM$valueCues,1,paste,collapse="")
    ordering <- sapply(orderStrings2, function(x)which(orderStrings1 == x))
    CNOl@cues <- CNOl@cues[ordering,CNOlistDREAM$namesCues]
    CNOl@stimuli <- CNOl@stimuli[ordering,CNOlistDREAM$namesStimuli]
    CNOl@inhibitors <- CNOl@inhibitors[ordering,CNOlistDREAM$namesInhibitors]    
    for (i in 1:length(CNOl@signals))
    {
      CNOl@signals[[i]] <- CNOl@signals[[i]][ordering,CNOlistDREAM$namesSignals]
      CNOl@variances[[i]] <- CNOl@variances[[i]][ordering,CNOlistDREAM$namesSignals]
    }
  }
  return(CNOl)
}

cutLiverMSB <- function()
{
  CNOl <- CNOlist("MSB2009-MIDAS-RawDataTraining.csv")
  colnames(CNOl@cues) <- gsub("mek","mek12",fixed=T,tolower(colnames(CNOl@cues)))
  colnames(CNOl@stimuli) <- gsub("mek","mek12",tolower(colnames(CNOl@stimuli)))
  colnames(CNOl@inhibitors) <- gsub("mek","mek12",tolower(colnames(CNOl@inhibitors)))    

  for (i in 1:length(CNOl@signals))
    colnames(CNOl@signals[[i]]) <- gsub("mek","mek12",fixed=T,gsub("jnk","jnk12",fixed=T,gsub("erk","erk12",fixed=T,tolower(colnames(CNOl@signals[[i]])))))
    
  data(CNOlistDREAM)
  CNOlistDREAM <- CNOlist(CNOlistDREAM)

  conditions <- c()
  for (i in 1:nrow(CNOlistDREAM@cues))
  {
    
    for (j in 1:nrow(CNOl@cues))
    {

     if (all(CNOl@cues[j,colnames(CNOlistDREAM@cues)] == CNOlistDREAM@cues[i,]))
       conditions <- c(conditions,j)
    }
  }
  print(setdiff(colnames(CNOlistDREAM@signals[[2]]), colnames(CNOl@signals[[2]])))
  return(list(cues=CNOl@cues[conditions,colnames(CNOlistDREAM@cues)], 
              signals = CNOl@signals[[2]][,colnames(CNOlistDREAM@signals[[2]])]))
}

start_liverdream_raw <- function()
{

}

start_cno <- function()
{
  res_dream <<- crossValidate(sifFile="PKN-LiverDREAM.sif",
                              midasFile="MD-LiverDREAM.csv",
                              outputDir="LiverDREAM/cv_cno",
                              ntimes=10,nfold=10,cpus=10,preprocess=TRUE,func=cv_CNO,seed=1234321)
                              #ntimes=2,nfold=2,maxGens=10,cpus=4,preprocess=TRUE,func=cv_CNO,seed=1234321)
  res_msb <<- crossValidate(sifFile="PKN-LiverMSB2009.sif",
                              midasFile="MSB2009-MIDAS-RawDataTraining.csv",
                              outputDir="LiverMSB/cv_cno",
                              ntimes=10,nfold=10,cpus=10,preprocess=TRUE,func=cv_CNO,seed=1234567890)
}

start_cantata <- function()
{
  #res_dream <<- crossValidate(sifFile="PKN-LiverDREAM.sif",
  #                            midasFile="MD-LiverDREAM.csv",
  #                            outputDir="LiverDREAM/cv_cant",
  #                            ntimes=10,nfold=10,cpus=10,preprocess=TRUE,func=cv_CANTATA,seed=1234321)
  res_msb <<- crossValidate(sifFile="PKN-LiverMSB2009.sif",
                              midasFile="MSB2009-MIDAS-RawDataTraining.csv",
                              outputDir="LiverMSB_1/cv_cant",
                              ntimes=10,nfold=10,cpus=6,maxStates=100,numStarts=1,
                              preprocess=TRUE,
                              func=cv_CANTATA,seed=12345678)
  res_msb <<- crossValidate(sifFile="PKN-LiverMSB2009.sif",
                              midasFile="MSB2009-MIDAS-RawDataTraining.csv",
                              outputDir="LiverMSB/cv_cant",
                              ntimes=10,nfold=10,cpus=6,maxStates=100,preprocess=TRUE,
                              func=cv_CANTATA,seed=1234567890)                              
}
