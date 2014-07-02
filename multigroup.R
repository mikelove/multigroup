pheno <- read.delim("GTEx_Analysis_Annotations_Sample_DS__Pilot_2013_01_31.txt")
genecounts <- read.delim("GTEx_Analysis_RNA-seq_RNA-SeQCv1.1.8_gene_reads__Pilot_2013_01_31_patch1.gct",skip=2)

pheno$cleanid <- make.names(pheno$SAMPID)

# which samples correspond to RNA-Seq (some are exon array)
sequenced <- pheno$SMGEBTCHT == "TrueSeq.v1"
table(sequenced)

# which samples are in the gene counts table
ingenecounts <- pheno$cleanid %in% colnames(genecounts)
table(ingenecounts)

# build a list of samples for tissues
tissuecompare <- list()
n <- 5 # samples per tissue

tissues <- c("Blood","Brain","Heart","Lung","Muscle","Skin")

for (i in seq_along(phenolist)) {
  tissuecompare[[i]] <- as.character(head(pheno[pheno$SMTS == tissues[i]
                                                  & sequenced & ingenecounts,"cleanid"],n))
}

# the first two rows of the gene count is the name and description
idx <- match(unlist(tissuecompare), colnames(genecounts))

# extract only these columns specified above
genecounts.sub <- genecounts[,idx]
rownames(genecounts.sub) <- genecounts$Name
tissue <- factor(pheno$SMTS[match(colnames(genecounts.sub), pheno$cleanid)])

rs <- rowSums(genecounts.sub)
hist(log10(rs+1),col="grey",breaks=40)

rowidx <- rs > 10
table(rowidx)

##

library(DESeq2)
dds <- DESeqDataSetFromMatrix(genecounts.sub[rowidx,],
                              DataFrame(sample=colnames(genecounts.sub),
                                        t=tissue),
                              ~ t)
mcols(dds)$description <- genecounts$Description[rowidx]

dds <- dds[1:1000,]
dds <- DESeq(dds)
res <- results(dds,name="tBlood")
res <- results(dds,contrast=c("t","Brain","Blood"))
res <- results(dds,contrast=list(c("tBlood","tBrain"),
                     c("tHeart","tSkin","tLung","tMuscle")),
               listValues=c(1/2,-1/4))
resSort <- res[order(-res$stat),]

##

gene <- rownames(resSort)[1]
intercept <- coef(dds)[gene,1]
cond <- coef(dds)[gene,-1]
condSE <- coef(dds,SE=TRUE)[gene,-1]
plotCounts(dds, gene, "t", transform=TRUE,
           main=mcols(dds,use.names=TRUE)[gene,"description"])

# draw the intercept and tissue specific fold changes
abline(h=intercept, lwd=2, col="green", lty=2)
x <- 1:6 + .1
qn <- qnorm(.975)
segments(x, intercept, x, intercept+cond, col="red", lwd=2)
points(x, intercept+cond, col="red", lwd=2)
arrows(x, intercept+cond, x, intercept+cond+qn*condSE, col="red", lwd=2, angle=90, length=.075)
arrows(x, intercept+cond, x, intercept+cond-qn*condSE, col="red", lwd=2, angle=90, length=.075)

