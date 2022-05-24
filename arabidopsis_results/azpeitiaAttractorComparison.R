#compare attractors azpeitia

cantata_azpeitia <- loadNetwork("./Azpeitia2013RSCN_Net_noRed.txt")

attr <- getAttractors(cantata_azpeitia)
binAttr <- lapply(attr$attractors, function(a) {r <- BoolNet:::dec2bin(a$involvedStates, 11); names(r) <- attr$stateInfo$genes; return(r)})

attractorMatrix <- matrix(data = c(1,0,0,0,1,0,0,1,0,0,0,1,1,0,0,0,0,0,1,0,0,0,1,1,1,1,0,1,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,1,1,1,1,0,0,0,0,1,0,0,0,0,1,1,1,1,0,1,1,0,0,0,0,0,1,0,0,0,0,1,0,0,1,1,1,1,1,0,0,1,1,0,1,0,0), byrow = T, ncol = 11)
attractorMatrix <- apply(attractorMatrix, MARGIN = c(1,2), as.logical)
colnames(attractorMatrix) <- c("SHR", "miRNA165", "JKD", "MGP", "PHB", "SCR", "IAA5", "Auxin", "WOX5", "CLE", "ACR")

apply(attractorMatrix, MARGIN = 1, function(a) any(sapply(binAttr, function(b) all(b[colnames(attractorMatrix)], a))))

rownames(attractorMatrix) <- c("CVC", "PVC", "End", "Cor", "LCC", "VI", "CEI", "CLEI", "QC")

library("ggplot")
library("reshape2")

attractorMatrix <- data.frame(attractorMatrix)
attractorTable <- data.frame(gene = as.vector(replicate(9, colnames(attractorMatrix))), 
                             group = as.vector(sapply(rownames(attractorMatrix), function(name) replicate(11,name))), 
                             value = as.numeric(as.vector(t(attractorMatrix))))
attractorTable$gene <- factor(attractorTable$gene, levels = c("SHR", "miRNA165", "JKD", "MGP", "PHB", "SCR", "IAA5", "Auxin", "WOX5", "CLE", "ACR"))
attractorTable$group <- factor(attractorTable$group, levels = c("CVC", "PVC", "End", "Cor", "LCC", "VI", "CEI", "CLEI", "QC"))
attractorTable$value <- factor(attractorTable$value, levels = c("0","1"))


ggplot(data = attractorTable, aes(x = gene, y = group, fill = value)) + geom_tile(color= "white", linetype=1, lwd = 2) + facet_wrap(~group, scales = "free_y", nrow = 9) + scale_x_discrete(position = "top") + scale_fill_manual(breaks=c(0,1), values=c(rgb(58,58,58,1, max=255), "green")) + xlab("Gene") + ylab("Attractor (cell type)") + theme_minimal() + theme(
  strip.background = element_blank(),
  strip.text.x = element_blank()
)



