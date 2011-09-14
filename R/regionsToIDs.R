regionsToIDs <- function(x, txdb, where=c("TSS", "gene"), within=100) {
    where <- match.arg(where)
    tx <- transcripts(txdb, columns=c("gene_id"))
    if (where=="TSS") tx <- resize(resize(tx, 1), within*2, fix="center")
    as.character(values(subsetByOverlaps(tx, x))$gene_id)
}


