# http://seqanswers.com/forums/archive/index.php/t-42420.html 
# Using dexseq to detect intron retention events in diseased samples
# Count how many reads map to each intron then use dexseq to do comparison
# Identification of differential intron retention using dexseq --> need to prepare intron count & annotation files

# INPUT: an annotation GFF file prepared by dexseq
# PROCESS: adds in "intronic_part" records with associated gene_ids
# OUTPUT: 1 file containing = original exonic bins + newer intronic bins

library(GenomicRanges)
library(rtracklayer)
library(parallel)

# import gtf file & obtain metadata
scn_gtf <- import.gff2("Mus_musculus.GRCm38.94.dexseq.gtf")
scn_meta <- elementMetadata(scn_gtf)

# add "intronic_part" as a factor in metadata
elementMetadata(scn_gtf)$type <- factor(elementMetadata(scn_gtf)$type, levels=c(levels(elementMetadata(scn_gtf)$type), "intronic_part"))
# filter out NA values in exonic_part_number column
scn_USE <- which(!is.na(elementMetadata(scn_gtf)$exonic_part_number))

# view current exonic_part_number column
elementMetadata(scn_gtf)$exonic_part_number
# format exonic_part_number and filter out NA values
scn_exonic_parts <- sprintf("%03s", elementMetadata(scn_gtf)$exonic_part_number[scn_USE])
# view filtered exonic_part_number column
scn_exonic_parts
elementMetadata(scn_gtf)$exonic_part_number <- as.character(elementMetadata(scn_gtf)$exonic_part_number)
# replace exonic_part_number column with filtered version
elementMetadata(scn_gtf)$exonic_part_number[scn_USE] <- scn_exonic_parts

# GRangesList
scn_grl <- split(scn_gtf, elementMetadata(scn_gtf)$gene_id)
scn_grl

# function that takes in a GRangesList and outputs GRangesList with intronic parts incorporated between the exonic regions
add_introns <- function(gr) {
  exons <- gr[which(elementMetadata(gr)$type=="exonic_part"),] # load in exonic parts types
  if(length(exons) > 1) {# exons=a list of exons 
    seqname <- seqnames(exons)[-1]# all but the 1st element
    starts <- end(exons)+1# end(x) grabs the last value of an interval (we are starting at the end of exons) Eg. 124321 + 1
    starts <- starts[-length(starts)] # gets the length of starts (1)
    ends <- start(exons)-1# get the start value of exon interval +1
    ends <- ends[-1]# all but the first element
    bounds <- IRanges(start=starts, end=ends) # defining integer ranges
    strand <- strand(exons)[-1] # all but the first element
    # GRanges object is a collection of genomic features 
    introns <- GRanges(seqnames=seqname, ranges=bounds, strand=strand) # constructing a new GRanges object
    intron_ids <- sprintf("%03i", c(1:length(introns))) # naming the introns (intron_ids) 
    
    # Remove 0-width introns based on their ranges
    DISCARD <- which(width(introns) <= 0) 
    # which() returns a vector of introns that contains width <= 0
    # the vector of introns to discard
    # QNS: WHY will there be introns that are negative width?
    if(length(DISCARD) > 0) { # if there are things to discard
      introns <- introns[-DISCARD] # grabs all but the DISCARD elements
      intron_ids <- intron_ids[-DISCARD] # grabs all but the DISCARD id
    }
    if(length(introns) > 0) {# if there are introns to create
      # create the meta-data
      df <- as.data.frame(elementMetadata(exons))# load in the current exons as a meta-dataframe where exons <- gr[which(elementMetadata(gr)$type=="exonic_part"),]
      nrows <- length(introns) # number of new rows to add for introns
      
      # i think concatenating df (exons) with empty rows for introns
      metadf <- df[1:nrows,] # does this need to deal with gene_id and transcripts differently?
      
      # basically just transforming the gene_id and transcripts cols to character
      metadf <- transform(metadf, gene_id=as.character(gene_id), transcripts=as.character(transcripts))
      
      # replicate "NA" nrows times (i think this is added as rows)
      metadf$transcripts <- as.character(c(rep(NA, nrows)))
      
      # replicate "intronic_part" nrows times under the $type column --> then factor the $type column in metadf
      metadf$type <- factor(c(rep("intronic_part", nrows)), levels=levels(metadf$type))
      
      # appending intron_ids to $exonic_part_number column in metadf
      metadf$exonic_part_number <- intron_ids # we created the intron_id earlier
      
      # similar to mcols function (get or set metadata columns)
      #i think i load in metadf into introns metadata
      elementMetadata(introns) <- metadf
      
      # Merge the GRanges containing exons and introns
      # introns is its own metadata
      gr <- append(gr, introns)
      gr <- gr[order(start(gr), elementMetadata(gr)$type),] # re-sort from start of gr, according to type (exonic_part or intronic_part)
    }
  }
  return(gr)
}

# endoapply takes quite long (around 1 hour) # 10.23am-11.17am
scn_with_introns <- endoapply(scn_grl, add_introns)
# should return a GRangesList with intronic regions incorporated
scn_with_introns

# chrom number of each sequence
scn_chroms <- sapply(scn_with_introns, function(x) as.factor(seqnames(x))[1])
# starting number of interval of each sequence
scn_starts <- sapply(scn_with_introns, function(x) start(x)[1])
# re-sort the sequence according to chromosome number and start position
scn_o <- order(scn_chroms, scn_starts)
# new spc_with_introns ordered according to chr number and start position
scn_with_introns_ordered <- scn_with_introns[scn_o]
scn_with_introns_ordered

# function that takes in a GRangesList and outputs a GTF file
asGFF2 <- function(x) {
  df <- as.data.frame(x) # take in a GRangesList, x
  aggregates <- which(df$type == "aggregate_gene") # returns TRUE for the positions that are aggregate genes
  
  # character creates a character vector of the specified length
  meta <- character(nrow(df)) # meta is a blank character vector with length = no of rows in df
  # meta has a column called aggregates containing all the gene_id of the aggregates in df
  meta[aggregates] <- sprintf("gene_id \"%s\"", df$gene_id[aggregates])
  # This gives introns a transcript "NA" field, which may not be ideal
  # keeps everything but aggregates
  meta[-aggregates] <- sprintf("transcripts \"%s\"; exonic_part_number \"%s\"; gene_id \"%s\"", df$transcripts[-aggregates], df$exonic_part_number[-aggregates], df$gene_id[-aggregates])
  # paste():Concatenate vectors after converting to character. converts its arguments (via as.character) to character strings
  paste(df$seqnames, "dexseq_prepare_annotation.py", df$type, df$start, df$end, ".", df$strand, ".", meta, sep="\t")
}

# asGFF2 function takes a while too 11.48am
scn_final_output <- unlist(lapply(scn_with_introns_ordered, asGFF2))
scn_final_output
# write out into a gtf file that can be used
write.table(scn_final_output, file="Mus_musculus.GRCm38.94.dexseq.introns.gtf", row.names=F, col.names=F, quote=F)
