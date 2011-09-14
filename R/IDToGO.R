IDToGO <-function(ids, universe=NULL, ontology=c("BP", "MF", "CC"), ID=c("entrez", "genbank", "alias", "ensembl", "symbol", "genename", "unigene"), speciesdb=c("org.Mm.eg.db", "org.Hs.eg.db"), n=NULL, p.adj=p.adjust.methods, p.cutoff=0.05, orderBy=c("classic", "elim", "parentchild")) {

    ontology <- match.arg(ontology)
    ID <- match.arg(ID)
    speciesdb <- match.arg(speciesdb)
    p.adj <- match.arg(p.adj)
    orderBy <- match.arg(orderBy)

    ID2GO <- annFUN.org(ontology, feasibleGenes=universe, speciesdb, ID=ID)
    IDnames <- names(inverseList(ID2GO))
    genelist <- factor(as.integer(IDnames %in% ids))
    if (any(is.na(genelist))) warning("Only ", sum(ids %in% IDnames), " out of ", length(ids), " ids mapped!")
    names(genelist) <- IDnames

    GOdata <- new("topGOdata", ontology=ontology, allGenes=genelist, annot=annFUN.GO2genes, GO2genes=ID2GO, description="Genes")
    GOtest <- list("classic"=runTest(GOdata, algorithm="classic", statistic="fisher"),
                      "elim"=runTest(GOdata, algorithm="elim", statistic="fisher"),
               "parentchild"=runTest(GOdata, algorithm="parentchild", statistic="fisher"))
    if (is.null(n)) { #use pvalue cutoff
        p <- topGO::score(GOtest[[orderBy]])
        p <- p.adjust(p, p.adj)

        #No significant terms
        if (sum(p<p.cutoff)==0) return(data.frame())

        tmp <- GenTable(GOdata, classic=GOtest[[1]], elim=GOtest[[2]],  parentchild=GOtest[[3]], topNodes=sum(p<p.cutoff))
        tmp <- tmp[order(tmp[[orderBy]]),]
        tmp[[paste(orderBy, "adj")]] <- p[tmp$GO.ID]
        return(tmp)
    } else { #return set number of terms
        if (n==0) n <- length(topGO::score(GOtest[[orderBy]])) #return ALL terms
        return(GenTable(GOdata, classic=GOtest[[1]], elim=GOtest[[2]], parentchild=GOtest[[3]], orderBy=orderBy, topNodes=n))
    }
}

