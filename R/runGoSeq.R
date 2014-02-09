runGoSeq<-function(de.genes,all.genes,onto,genome,ID,rep=1000,p.adj="BH",adj.p.val=.05)
{
require(goseq)
if (length(onto) != 1) {stop("ontology must be one of GO:BP GO:MF or GO:CC")} 
	
	    cat("filter for genes annotated to root.term \n")
	
		if (onto=="GO:BP") root.term<-"GO:0008150"

		if (onto=="GO:MF") root.term<-"GO:0003674"

		if (onto=="GO:CC") root.term<-"GO:0005577"
	
		cat(paste(onto,root.term,"\n"))
	
	out<-getgo(all.genes,genome,ID)
	#filter for genes not in db
	out2 <- out[!is.na(names(out))]
	#filter for genes not annotated to root.term
	out3 <- out2[!is.na(sapply(out2,function(x) match(root.term,x)))]

	cat("create gene vector \n")
	
	#make a vector of 1 for DE 0 for not, names are genes.
	assayed.genes<-names(out3)
	#this is the background for gene ontology.

	both<-as.character(na.omit(assayed.genes[match(de.genes,assayed.genes)]))
	#these are DE genes.
	
	gene.vector.both<-as.integer(assayed.genes%in%both)
	names(gene.vector.both)<-assayed.genes

	cat( paste(onto,root.term,"background =",length(gene.vector.both),"genes","\n"))
	cat( paste(onto,root.term,"de.genes =",sum(gene.vector.both),"genes","\n"))	

	cat("make pwf \n")

	pdf(file=paste(gsub("GO:","",onto),"pwf.pdf",sep=""))

	pwf_both<-nullp(gene.vector.both,genome,ID)

	dev.off()

	cat("go seq \n")

	GO.both<-goseq(pwf_both,genome,ID,method="Sampling",test.cats=c(onto),repcnt=rep)

	# pvalue correction
    cat("\n p value adjust \n")
	enriched<-as.integer(p.adjust(GO.both$over_represented_pvalue,method=p.adj)< adj.p.val)
    enriched_adjusted<-p.adjust(GO.both$over_represented_pvalue,method=p.adj)
	depleted<-as.integer(p.adjust(GO.both$under_represented_pvalue,method=p.adj)< adj.p.val)
    depleted_adjusted<-p.adjust(GO.both$under_represented_pvalue,method=p.adj)
	
	significant.GO<-cbind(GO.both,"enriched"=enriched,"enriched_adjusted"=enriched_adjusted,"depleted"=depleted,"depleted_adjusted"=depleted_adjusted)
	
	return(significant.GO)
	cat("done \n")

}
