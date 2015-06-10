
mut2tran_list <- commandArgs(trailingOnly = TRUE)[1]
sample_list <- commandArgs(trailingOnly = TRUE)[2]
gene_list <- commandArgs(trailingOnly = TRUE)[3]
output_fig <- commandArgs(trailingOnly = TRUE)[4]


samples <- read.table(sample_list, header = FALSE, stringsAsFactors = FALSE)[,1]
genes <- read.table(gene_list, header = FALSE, stringsAsFactors = FALSE)[,1]

sampleNum <- length(samples)
geneNum <- length(genes)

xspace <- 0.05;
yspace <- 0.05;



jpeg(filename = output_fig, width = (sampleNum * (1 + xspace) - xspace) * 40 + 240, height = (geneNum * (1 + yspace) - yspace) * 40 + 240, quality = 100);

plot.new()
par(mar=c(0.1, 6.1, 6.1, 0.1))
plot.window(xlim=c(0, sampleNum * (1 + xspace) - xspace), ylim=rev(c(0, geneNum * (1 + yspace) - yspace )))


for (i in 1:sampleNum) {
	for (j in 1:geneNum) {
		rect( (i - 1) * (1 + xspace) , (j - 1) * (1 + yspace), (i - 1) * (1 + xspace) + 1, (j - 1) * (1 + yspace) + 1, col="gray96", border=F);
	}
}

for (i in 1:sampleNum) {
 	mtext(samples[i], side=3, line= -1, at=(i - 1)*(1+ xspace) + 0.5, las=2, adj = 0, cex = 1.5)
}

for (i in 1:geneNum) {
	mtext(genes[i], side=2, line= -1, at=(i - 1) * (1 + yspace) + 0.5, las=1, adj = 1, cex = 1.5)
}


result <- read.table(mut2tran_list, sep="\t", header=F, stringsAsFactors = FALSE);

for (i in 1:nrow(result)) {
	if (result[i, 1] %in% genes)  {

		gind <- which(genes == result[i,1] );
		sind <- which(samples == result[i,2]);
	
		tcol <- "black";
        pflag <- 0
		if (result[i,3] == "point mutation,indel") {
			tcol <- "yellow";
            pflag <- 1;
		} 	
		if (result[i,3] == "structural variation") {
			tcol <- "purple";
            pflag <- 1;
		} 	
		if (result[i,3] == "HBV integration") {
			tcol <- "cyan";
            pflag <- 1;
		} 	

        if (pflag == 1) {
            tx <- c((sind - 1) * (1 + xspace), (sind - 1) * (1 + xspace) + 1, (sind - 1) * (1 + xspace) );
            ty <- c((gind - 1) * (1 + yspace), (gind - 1) * (1 + yspace), (gind - 1) * (1 + yspace) + 1);
            polygon( tx, ty,  col=tcol, border=F);
        }


        pflag <- 0;
        if (result[i,4] == "splicing aberration") {
            tcol <- "green";
            pflag <- 1;
        }
        if (result[i,4] == "gene fusion") {
            tcol <- "red";
            pflag <- 1;
        }
		if (result[i,4] == "HBV fusion") {
			tcol <- "blue";
			pflag <- 1;
		} 	
		if (result[i,4] == "over expression") {
			tcol <- "chocolate";
			pflag <- 1;
		} 

		if (pflag == 1) {
			tx <- c((sind - 1) * (1 + xspace), (sind - 1) * (1 + xspace) + 1, (sind - 1) * (1 + xspace) + 1 );
			ty <- c((gind - 1) * (1 + yspace) + 1, (gind - 1) * (1 + yspace) + 1, (gind - 1) * (1 + yspace));
		    polygon( tx, ty,  col=tcol, border=F);
        }

	}
}










