create450kbg <- function(probes, groups=NULL){
  stopifnot(is.character(probes))
  if (!is.null(groups)){
    stopifnot(all(groups %in% c("TSS1500", "TSS200", "Body", "1stExon", "5'UTR", "3'UTR")))
  }
  #Dummy
  dummy <- matrix(nrow=length(probes), ncol=2)
  rownames(dummy) <- probes
  RSanno <- RatioSet(dummy, annotation=c(array="IlluminaHumanMethylation450k", annotation="ilmn12.hg19"))
  anno <- getAnnotation(RSanno)
  anno <- as.data.frame(anno)
    
  if (!is.null(groups)){
    message("Subgrouping detected\n")
    groupgrep <- paste(groups, collapse="|")
    splitanno <- anno[, c("UCSC_RefGene_Name", "UCSC_RefGene_Group")]
    splitanno <- splitanno[grep(groupgrep, splitanno$UCSC_RefGene_Group),]
    probes <- rownames(splitanno)
    splitanno <- as.matrix(splitanno)
    message("Splitting annotation...\n")
    splitanno <- apply(splitanno, 1, function(x) cbind(unlist(strsplit(x[1], ";")), unlist(strsplit(x[2], ";"))))
    splitanno <- sapply(splitanno, function (x) apply(as.matrix(x), 1, function (y) paste(y, collapse=".")))
    message("Annealing annotation...")
    groupprobes <- lapply(splitanno, function (x) x[grep(groupgrep, x)])
    groupgenes <- lapply(groupprobes, function (x) unique(gsub("\\..*", "", (x))))
    bg0 <- table(unlist(groupgenes))
    bg <- as.integer(bg0)
    names(bg) <- names(bg0)
    res <- c(bg, probes, groups)
    res <- list(bg, probes, groups)
    names(res) <- c("background", "probes", "groups")
      
    
  } else {
    uniquepp <- lapply(strsplit(as.character(anno$UCSC_RefGene_Name), split=';'), function(x) unique(unlist(x)))
    genenames <- unique(unlist(uniquepp))
    bg0 <- table(unlist(uniquepp))
    bg <- as.integer(bg0)
    names(bg) <- names(bg0)
    res <- list(bg, probes, groups)
    names(res) <- c("background", "probes", "groups")
  }
  class(res) <- "background"
  return(res)  
}