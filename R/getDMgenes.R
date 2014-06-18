getDMgenes <- function(dmprobes, background){
  stopifnot(is(background, "background"))
  stopifnot(is.character(dmprobes))
  
  dmprobes <- dmprobes[dmprobes %in% background$probes]
  #Dummy
  dummy <- matrix(nrow=length(dmprobes), ncol=2)
  rownames(dummy) <- dmprobes
  RSanno <- RatioSet(dummy, annotation=c(array="IlluminaHumanMethylation450k", annotation="ilmn12.hg19"))
  anno <- getAnnotation(RSanno)
  anno <- as.data.frame(anno)
    
  if (!is.null(background$groups)){
    message("Subgrouping detected\n")
    groupgrep <- paste(background$groups, collapse="|")
    splitanno <- anno[, c("UCSC_RefGene_Name", "UCSC_RefGene_Group")]
    splitanno <- splitanno[grep(groupgrep, splitanno$UCSC_RefGene_Group),]
    splitanno <- as.matrix(splitanno)
    message("Splitting annotation...\n")
    splitanno <- apply(splitanno, 1, function(x) cbind(unlist(strsplit(x[1], ";")), unlist(strsplit(x[2], ";"))))
    splitanno <- sapply(splitanno, function (x) apply(as.matrix(x), 1, function (y) paste(y, collapse=".")))
    message("Annealing annotation...")
    groupprobes <- lapply(splitanno, function (x) x[grep(groupgrep, x)])
    groupgenes <- lapply(groupprobes, function (x) unique(gsub("\\..*", "", (x))))
    genenames <- unique(unlist(groupgenes))
    
        
  } else {
    uniquepp <- lapply(strsplit(as.character(anno$UCSC_RefGene_Name), split=';'), function(x) unique(unlist(x)))
    genenames <- unique(unlist(uniquepp))
    
  }
  
  return(genenames)    

}