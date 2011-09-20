makeGOTopTable <- function(ids, outfile="topGOout.tex", main="Significant GO Terms", n.print=25, ...) {
    resTable <- IDToGO(ids, ...)
    writeLines(c("\\documentclass[a4paper,12pt]{article}",
                 "\\usepackage[landscape]{geometry}",
                 "\\usepackage{fullpage}",
                 "\\usepackage{longtable}",
                 "\\begin{document}",
                 paste("\\section*{\\centering{", main, "}}", sep="")), con=outfile)
    if (nrow(resTable)==0) warning("zomg! Results table is empty!") else print(xtable(resTable[1:min(n.print, nrow(resTable)),]), type="latex", floating=FALSE, tabular.environment="longtable", size="small", file=outfile, append=TRUE)
    f1 <- file(outfile, open="at")
    writeLines("\\end{document}", con=f1)
    close(f1)
    tools::texi2dvi(outfile, pdf=TRUE, clean=TRUE)
    invisible(resTable)
}

