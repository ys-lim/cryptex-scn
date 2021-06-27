############
## DEXSEQ ##
############

## This part analyses the pycount files generated during identification step
## We want to test for differential exon usage

library("DEXSeq")
library(BiocParallel)

# Count files from identification step
STAR_scn_countFiles = list.files(path=".",pattern=".txt$", full.names=TRUE)
basename(STAR_scn_countFiles)

# all_exons_and_cryptics.gff
STAR_scn_flattenedFile = list.files(path=".", pattern="gff$", full.names=TRUE)
basename(STAR_scn_flattenedFile)

STAR_scn_sampleTable = data.frame(
  row.names = c( "cko1", "cko2", "cko3","ctr1", "ctr2", "ctr3"),
  condition = c("knockout", "knockout", "knockout","control", "control", "control"),
  libType = c( "pair-end", "pair-end", "pair-end", 
               "pair-end", "pair-end", "pair-end") )
STAR_scn_sampleTable #needs to be in same order as basename(countFiles) above

STAR_BPPARAM = SnowParam(workers=4)

# This does the testing of DE usage
STAR_scn_DexSeqExons.loc <- DEXSeqDataSetFromHTSeq(STAR_scn_countFiles,
                                              sampleData = STAR_scn_sampleTable,
                                              design = ~ sample + exon + condition * exon,
                                              flattenedfile = STAR_scn_flattenedFile)
                                              
# From actual script: formula1 <-  ~ sample + type * exon + condition * exon
## TEST: using ~ sample + type * exon + condition * exon formula from cryptex script
test_scn_DexSeqExons.loc <- DEXSeqDataSetFromHTSeq(STAR_scn_countFiles,
                                                   sampleData = STAR_scn_sampleTable,
                                                   design = ~ sample + type * exon + condition * exon,
                                                   flattenedfile = STAR_scn_flattenedFile)
# Error in DEXSeqDataSet(dcounts, sampleData, design, exons, genesrle, exoninfo[matching],  : 
# the variables 'type' of the parameter 'design' are not specified in the columns of the sampleData
# ^^ if have error downstream, go back and try this again

colData(STAR_scn_DexSeqExons.loc)
resultsNames(STAR_scn_DexSeqExons.loc)

# Other statistical stuff
STAR_scn_DexSeqExons.loc <- estimateSizeFactors(STAR_scn_DexSeqExons.loc) #10.34am-10.35am
STAR_scn_DexSeqExons.loc <- DEXSeq::estimateDispersions(STAR_scn_DexSeqExons.loc) # takes very long! about 1 hour++
STAR_scn_DexSeqExons.loc <- DEXSeq::testForDEU(STAR_scn_DexSeqExons.loc, BPPARAM=STAR_BPPARAM) # 3.19pm-4.20pm
STAR_scn_DexSeqExons.loc <- DEXSeq::estimateExonFoldChanges(STAR_scn_DexSeqExons.loc, BPPARAM=STAR_BPPARAM) # about 1 hour?
rowData(STAR_scn_DexSeqExons.loc)

# Cleaning up dexseq output
STAR_scn_res_new <- DEXSeq::DEXSeqResults(STAR_scn_DexSeqExons.loc) #few min 
head(STAR_scn_res_new)
STAR_logname <- grep(names(STAR_scn_res_new), pattern = 'log2fold', value = TRUE)
STAR_scn_res.clean <- as(STAR_scn_res_new[, c('groupID', 'featureID', 'exonBaseMean', logname, 'dispersion', 'stat', 'pvalue')], 'data.frame')
head(STAR_scn_res.clean)
names(STAR_scn_res.clean)<- c("EnsemblID", "exonID", "meanBase", "log2FoldChange", "dispersion", "stat", "pvalue")
STAR_scn_res.clean$FDR <- p.adjust(STAR_scn_res.clean$pvalue, method = 'fdr')
STAR_scn_res.clean$chromosome <- as.character(seqnames(STAR_scn_res_new$genomicData))
STAR_scn_res.clean$exon.start <- start(STAR_scn_res_new$genomicData)
STAR_scn_res.clean$exon.end <- end(STAR_scn_res_new$genomicData)

library("biomaRt")
listEnsembl(version=94)
STAR_scn_ensembl <- useEnsembl(biomart = "ensembl")
STAR_scn_ensembl
STAR_scn_datasets <- listDatasets(STAR_scn_ensembl)
head(STAR_scn_datasets)
searchDatasets(mart = STAR_scn_ensembl, pattern = "mmusculus_gene_ensembl")
STAR_scn_mouse_annotation <- useDataset(dataset = "mmusculus_gene_ensembl", mart = STAR_scn_ensembl)
listAttributes(STAR_scn_mouse_annotation)
searchAttributes(STAR_scn_mouse_annotation, pattern="ensembl_gene_id")
library(dplyr)
# added start & end position & chr name (needed for downstream functions)
STAR_scn_mouse_annotation_infos <- getBM(attributes=c('ensembl_gene_id','external_gene_name', 'start_position', 'end_position','strand'),
                                    mart = STAR_scn_mouse_annotation,useCache = FALSE)
head(STAR_scn_mouse_annotation_infos)

STAR_scn_res.clean$external_gene_id <- STAR_scn_mouse_annotation_infos$ensembl_gene_id[ match(STAR_scn_res.clean$EnsemblID, table = STAR_scn_mouse_annotation_infos$ensembl_gene_id) ]
head(STAR_scn_res.clean)
tail(STAR_scn_res.clean)

STAR_scn_res.clean <- STAR_scn_res.clean[, c('external_gene_id', "EnsemblID", "exonID", "meanBase", "log2FoldChange", "dispersion", "stat", "pvalue", "FDR", "chromosome", "exon.start", "exon.end")]  ### reorder the names nicely

## add strand if available
if ('strand' %in% names(STAR_scn_mouse_annotation_infos)) STAR_scn_res.clean$strand <- STAR_scn_mouse_annotation_infos$strand[ match(STAR_scn_res.clean$EnsemblID, table = STAR_scn_mouse_annotation_infos$ensembl_gene_id) ]

STAR_scn_res.clean <- STAR_scn_res.clean[ order(STAR_scn_res.clean$pvalue),]  ##reorder the rows
write.csv(x = STAR_scn_res.clean,
          file="24may_scn_SignificantExons_STAR.csv",
          row.names = FALSE)

# getting cryptic exons
STAR_scn_res.clean.cryptics <- subset(STAR_scn_res.clean,STAR_scn_res.clean$FDR < 0.01 & grepl("i",STAR_scn_res.clean$exonID))
head(STAR_scn_res.clean.cryptics)
STAR_scn_res.clean.cryptics.up <- subset(STAR_scn_res.clean.cryptics, STAR_scn_res.clean.cryptics$log2FoldChange > 0)
head(STAR_scn_res.clean.cryptics.up)
STAR_scn_res.clean.cryptics.down <- subset(STAR_scn_res.clean.cryptics, STAR_scn_res.clean.cryptics$log2FoldChange < 0)
head(STAR_scn_res.clean.cryptics.down)
STAR_scn_res.clean.cryptics.NA <- subset(STAR_scn_res.clean.cryptics, is.na(STAR_scn_res.clean.cryptics$log2FoldChange))
head(STAR_scn_res.clean.cryptics.NA)
nrow(STAR_scn_res.clean.cryptics.NA) #7 rows

write.table(x = STAR_scn_res.clean.cryptics,
            file="24may_scn_CrypticExons_STAR.tab",
            row.names = FALSE)

STAR_scn_codes <- c("Significant Cryptic events (FDR < 0.01):", "Up-going:", "Down-going:", "NA (possible error):")
STAR_scn_counts <- c(dim(STAR_scn_res.clean.cryptics)[1], dim(STAR_scn_res.clean.cryptics.up)[1], dim(STAR_scn_res.clean.cryptics.down)[1], dim(STAR_scn_res.clean.cryptics.NA)[1] )
STAR_scn_report <- data.frame(STAR_scn_codes,STAR_scn_counts)
STAR_scn_report

write.table(x = STAR_scn_report,
            file="24may_scn_Cryptic_Report_STAR.tab",
            row.names = F, quote = F)
write.table(x = STAR_scn_res.clean.cryptics.up[,c(10:12,1,3)],
            file="24may_scn_Cryptics_UP_STAR.bed",
            quote=F, row.names=F, col.names=F, sep="\t")
write.table(x = STAR_scn_res.clean.cryptics.down[,c(10:12,1,3)],
            file="24may_scn_Cryptics_DOWN_STAR.bed",
            quote=F, row.names=F, col.names=F, sep="\t")

STAR_scn_n.sig <- sum(STAR_scn_res.clean$FDR < 0.01, na.rm = TRUE)
if (STAR_scn_n.sig <= 50) {
  STAR_scn_res.cleanSigs <- subset(STAR_scn_res.clean, FDR<0.01 & grepl("i",STAR_scn_res.clean$exonID))
} else STAR_scn_res.cleanSigs <- subset(STAR_scn_res.clean, grepl("i",STAR_scn_res.clean$exonID))[1:50,]

STAR_scn_genes.to.plot <- unique(STAR_scn_res.cleanSigs$EnsemblID)
STAR_scn_pretty.gene.names <- as.character(STAR_scn_mouse_annotation_infos$external_gene_name[ match(STAR_scn_genes.to.plot, table = STAR_scn_mouse_annotation_infos$ensembl_gene_id) ])

for (i in 1:length(STAR_scn_genes.to.plot)) {
  STAR_scn_gene <- as.character(STAR_scn_genes.to.plot[i])
  
  if (!is.na(STAR_scn_pretty.gene.names[ i ])) {
    STAR_scn_gene.pretty <- as.character(STAR_scn_pretty.gene.names[ i ])
    
    message(i, ' ', STAR_scn_gene, ' ', STAR_scn_gene.pretty)
    
    STAR_scn_output.pdf <- paste(STAR_scn_gene.pretty, '.pdf', sep = '')
    pdf(STAR_scn_output.pdf, width = 16, height = 9.8)
    plotDEXSeq(STAR_scn_res_new,
               geneID = STAR_scn_gene,  ##I suspect it has to be gene, otherwise it crashes
               cex.axis = 1.2,
               cex=1.3,
               lwd=2,
               legend=TRUE,
               displayTranscripts = TRUE,
               names = TRUE,
               main = STAR_scn_gene.pretty)
    dev.off()
    print(STAR_scn_output.pdf)
  }
}


STAR_nfasc <- "ENSMUSG00000026442"
STAR_nfasc.gene.name <- as.character("Nfasc")
STAR_nfasc.pdf <- paste('Nfasc.pdf')
pdf(STAR_nfasc.pdf, width = 16, height = 9.8)
plotDEXSeq(STAR_scn_res_new,
           geneID = STAR_nfasc,  ##I suspect it has to be gene, otherwise it crashes
           cex.axis = 1.2,
           cex=1.3,
           lwd=2,
           legend=TRUE,
           displayTranscripts = TRUE,
           names = TRUE,
           main = STAR_nfasc.gene.name)
dev.off()
print(STAR_nfasc.pdf)

# Plot these graphs which I have no idea what they mean
pdf(file=paste('24MayDEXSeq-MeanVsDispPoints_STAR.pdf', sep = ''))
plotDispEsts (STAR_scn_DexSeqExons.loc)
dev.off()


pdf(file=paste('24MayDEXSeq-MeanVsDispCircles_STAR.pdf', sep = ''))
plotMA(data.frame(baseMean = STAR_scn_res.clean[,6],log2FoldChange = STAR_scn_res.clean[,7], padj = STAR_scn_res.clean[,5] < 0.1),
       ylim=c(-4,4), cex=0.8)
dev.off()
