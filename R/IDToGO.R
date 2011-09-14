IDToGO <-function(ids, universe=NULL, ontology=c("BP", "MF", "CC"), ID=c("entrez", "genbank", "alias", "ensembl", "symbol", "genename", "unigene"), speciesdb=c("org.Mm.eg.db", "org.Hs.eg.db"), n=NULL, p.adj=NULL, p.cutoff=0.05, orderBy=c("classic", "elim", "parentchild"), tests=NULL) {

    ontology <- match.arg(ontology)
    ID <- match.arg(ID)
    speciesdb <- match.arg(speciesdb)
    orderBy <- match.arg(orderBy)
    tests <- unique(c(orderBy, tests))
    names(tests) <- tests
    if (any(!tests %in% c("classic", "elim", "parentchild")))
        stop("One or more of ", paste(tests, collapse=" "), "not recognised - must be classic, elim or parentchild")
    if (!is.null(p.adj)) if (!p.adj %in% p.adjust.methods)
        stop("If p.adj is specified, it must be one of p.adjust.methods")

    ID2GO <- annFUN.org(ontology, feasibleGenes=universe, speciesdb, ID=ID)
    IDnames <- names(inverseList(ID2GO))
    genelist <- factor(as.integer(IDnames %in% ids))
    if (any(is.na(genelist))) warning("Only ", sum(ids %in% IDnames), " out of ", length(ids), " ids mapped!")
    names(genelist) <- IDnames

    GOdata <- new("topGOdata", ontology=ontology, allGenes=genelist, annot=annFUN.GO2genes, GO2genes=ID2GO, description="Genes")

    GOtest <- lapply(tests, function(x) runTest(GOdata, algorithm=x, statistic="fisher"))

    #gimme some pvalue lovin
    p <- topGO::score(GOtest[[orderBy]])
    if (!is.null(p.adj)) p <- p.adjust(p, p.adj) #we want to adjust dat bad boi

    if (is.null(n)) { #determine n from the pvalue cutoff
        if (sum(p<p.cutoff)==0) return(data.frame()) else n <- sum(p<p.cutoff)
    } else if (n==0) n <- length(topGO::score(GOtest[[orderBy]])) #return ALL terms

    tmp <- do.call(GenTable, c(object=GOdata, GOtest, orderBy=orderBy, topNodes=n))
    if (!is.null(p.adj)) tmp[[paste(orderBy, "adj")]] <- p[tmp$GO.ID] #tack on the adjusted
    return(tmp)
}

