# LOAD IN DEXSEQ RESULTS

###=========
## STRATEGY
###=========
# Each cryptic exon sits within an intron but the start and end of that intron cannot be inferred from the GFF file alone. Instead the canonical junction should be the junction that falls outside of the cryptic exons's coordinates with the most junction reads.

########################
# DECLARING FUNCTIONS #
########################
library(data.table,quietly=T)
library(dplyr)
library(stringr)
library(GenomicRanges,quietly=T)
library(ggplot2)

# This function reads in all the SJ.tab files given in the list and merges them all together to give a list of unique SJs with total numbers of occurences across all the datasets in the list.



########
## 1. ## 
########
# This function merges all SJ.tab.out files within the given SJ.tab.list to produce a list of unique SJs with total numbers of occurences across all the datasets in the list
merge.SJ.files <- function(SJ.tab.list){
	# First read in all the files in the list
	SJ.file.list <- list()
	for(i in 1:length(SJ.tab.list)){
        	SJ <- as.data.frame(fread(SJ.tab.list[i]))
        	SJ.file.list[[i]] <- SJ
        }
	# rbind all together
	SJ.merge <- rbindlist(SJ.file.list)
	SJ.merge <- as.data.frame(SJ.merge)
	#Use dplyr to group the total list of splice junctions by unique SJ and then sum the total of the cases.
	by_coord <- group_by(SJ.merge, paste(V1,V2,V3))
	SJ.summary <- summarise(by_coord,
        	count.unique = sum(V7),
        	count.multi = sum(V8), # count of multi-mapping reads
        	strand = mean(V4),
			intron.motif = mean(V5)
        	)
	names(SJ.summary)[1] <- "coord"
	#split the coordinate reference into columns
	SJ.summary$chr <- str_split_fixed(SJ.summary$coord, " ", 3)[,1]
	SJ.summary$start <- str_split_fixed(SJ.summary$coord, " ", 3)[,2]
	SJ.summary$end <- str_split_fixed(SJ.summary$coord, " ", 3)[,3]
	#this is where columns are dropped?
	SJ.summary <- SJ.summary[,c(6,7,8,2,3,4,5)]
	SJ.summary <- as.data.frame(SJ.summary)
	SJ.summary$start <- as.numeric(SJ.summary$start)
	SJ.summary$end <- as.numeric(SJ.summary$end)
	# make strand readable
	SJ.summary[SJ.summary$strand == 1,]$strand <- "+"
	SJ.summary[SJ.summary$strand == 2,]$strand <- "-"
    SJ.summary[SJ.summary$strand == 0,]$strand <- "*"
    #make intron motif readable
	motif.list <- c("non-canonical","GT/AG","CT/AC","GC/AG","CT/GC","AT/AC","GT/AT")
	for(i in 1:7){
		if(length(SJ.summary[SJ.summary$intron.motif == (i - 1),]$intron.motif) > 0){
			SJ.summary[SJ.summary$intron.motif == (i - 1),]$intron.motif <- motif.list[i]
		} 
	}
	#unique SJs will become "score" field in GRanges object
	names(SJ.summary)[4] <- "score"
	return(SJ.summary)
}

########
## 2. ## 
########

# this function below allows you to specify which GRange object to query
# this function outputs the canonical junction flanking the cryptic tag (since there could be multiple)
canonical_junction_query <- function(CE.chr,CE.start,CE.end, SJ.GRange){ # iterate each SJ through the GRanges of SJs
        junction <- SJ.GRange[seqnames(SJ.GRange)==CE.chr & start(SJ.GRange) <= as.numeric(CE.start)+1 & end(SJ.GRange) >= as.numeric(CE.end)-1] # obtain that particular SJ
        junction <- head(junction[order(score(junction),decreasing=T)],1) # we only want the most abundantly supported SJ
        return(junction)
}

########
## 3. ## 
########
# this function looks for junctions in the case dataset that exactly match those canonical junctions discovered in the control dataset
canonical_junction_replication <- function(CE.chr,canonical.start,canonical.end,SJ.GRange){
	junction <- SJ.GRange[seqnames(SJ.GRange) == CE.chr & start(SJ.GRange) == as.numeric(canonical.start) & end(SJ.GRange) == as.numeric(canonical.end)]
	return(junction)
}

########
## 4. ## 
########
#this function looks for SJs that span from the upstream end of the canonical intron to the 5' end of the cryptic exon. Only the most abundant SJ will be reported.
#does this have to be canonical.start or canonical.start +1?
# This is where the EXACT cryptic exon 5' site is discovered
upstream_junction_query <- function(CE.chr,CE.start,CE.end,canonical.start,SJ.GRange){
	junction <- SJ.GRange[seqnames(SJ.GRange) == CE.chr & start(SJ.GRange) == as.numeric(canonical.start) & end(SJ.GRange) >= as.numeric(CE.start) - 1 & end(SJ.GRange) < as.numeric(CE.end)]
	junction <- head(junction[order(score(junction),decreasing=T)],1)
    return(junction)
}

########
## 5. ## 
########
# Similar to upstream_junction_query, but for 3' site
downstream_junction_query <- function(CE.chr,CE.start,CE.end,canonical.end,SJ.GRange){
	junction <- SJ.GRange[seqnames(SJ.GRange) == CE.chr & start(SJ.GRange) > as.numeric(CE.start) & start(SJ.GRange) <= as.numeric(CE.end) + 1 & end(SJ.GRange) == as.numeric(canonical.end)]
	junction <- head(junction[order(score(junction),decreasing=T)],1)
    return(junction)
}

########
## 6. ## 
########
#applying the above functions to each cryptic tag and cleaning up the result
canonical_junction_detector <- function(SJ.summary,results.df,mode="discovery"){
	GRanges_object <-  makeGRangesFromDataFrame(SJ.summary,keep.extra.columns=T) # convert the combined SJ file into a GRanges obj
	if(mode == "discovery"){
		junctions.list <- apply(results.df, MAR=1,FUN=function(x) canonical_junction_query(x[10],x[11],x[12], GRanges_object)) # this is for control
		}
	if(mode == "replication"){
		junctions.list <- apply(results.df, MAR=1,FUN=function(x) canonical_junction_replication(x[10],x[15],x[16], GRanges_object)) # this is for case
	}
	#output is a list of GRange objects - unuseable.
	junctions.list <- unlist(GRangesList(junctions.list)) # list of canonical junctions in results.df, referenced from SJ.summary GRanges obj
	#convert into a dataframe, extracting the relevent information from the GRanges object.
	#names(GRanges) is a vector of rownames, confusingly.
	canonical.df <- data.frame(row.names=names(junctions.list),
			canonical.chr=seqnames(junctions.list),
			canonical.start=start(junctions.list),
			canonical.end=end(junctions.list),
			canonical.unique.count = score(junctions.list),
			canonical.strand = strand(junctions.list),
			intron.motif = mcols(junctions.list)[3])
	return(canonical.df)
}

# need to find the counts of canonical junctions that were assigned by canonical_junction_detector in the case dataset.

########
## 7. ## 
########
#This function assumes that the results file has been appended with the canonical start and end coordinates at positions XXX and YYY respectively.
# This function combines downstream_junction_query and upstream_junction_query to find the EXACT splice sites of cryptic exons
bridging_junction_finder <- function(SJ.summary, results.df, query.type){
	GRanges_object <-  makeGRangesFromDataFrame(SJ.summary,keep.extra.columns=T)
	if(query.type == "downstream"){ # 3'
		junctions.list <- apply(results.df, MAR=1,FUN=function(x) downstream_junction_query(x[10],x[11],x[12], x[16], GRanges_object))
	}
	if(query.type == "upstream"){ # 5'
		junctions.list <- apply(results.df, MAR=1,FUN=function(x) upstream_junction_query(x[10],x[11],x[12], x[15], GRanges_object))
	}
	#output is a list of GRange objects - unuseable.
	junctions.list <- unlist(GRangesList(junctions.list))
	#convert into a dataframe, extracting the relevent information from the GRanges object.
	#names(GRanges) is a vector of rownames, confusingly.
	if(query.type == "downstream"){
		bridging_junctions.df <- data.frame(row.names=names(junctions.list),
			chr=seqnames(junctions.list),
			cryptic.3prime=start(junctions.list),
			canonical.end=end(junctions.list),
			downstream.unique.count = score(junctions.list),
			downstream.strand = strand(junctions.list),
			intron.motif = mcols(junctions.list)[3])
	}
	if(query.type == "upstream"){
		bridging_junctions.df <- data.frame(row.names=names(junctions.list),
			chr=seqnames(junctions.list),
			canonical.start=start(junctions.list),
			cryptic.5prime=end(junctions.list),
			upstream.unique.count = score(junctions.list),
			upstream.strand = strand(junctions.list),
			intron.motif = mcols(junctions.list)[3])
	}
	return(bridging_junctions.df)
}

########
## 8. ## 
########
fix.gene.names <- function(results.df,annotation){
	# This is a function to find the most likely gene name for a given genomic range.
	# If multiple genes overlap the query range then output the names appended together with "+".
	# Get the strand information as well. If multiple genes are returned and their strands agree then output that strand.
	# If multiple genes are returned and the strands do not agree then output NA.
	gene.names.query <- function(chromosome, canonical.start, canonical.end, anno.GRange){
		gene.name <- anno.GRange[seqnames(anno.GRange) == chromosome & start(anno.GRange) <= as.numeric(canonical.start) & end(anno.GRange) >= as.numeric(canonical.end)]
		#if multiple gene names come up then collapse all the gene names into one string separated by "+".
		if(length(gene.name) > 1){
			mcols(gene.name)[,2] <- paste(mcols(gene.name)[,2],collapse="+")	
		# and strand as well!
		# if the different genes have differing strands then assign strand as NA.
			if(length(unique(as.list(strand(gene.name)))) != 1){
				strand(gene.name) <- NA
			}
			gene.name <- gene.name[1]
		}
		return(gene.name)
	}
	annotation <- as.data.frame(fread(annotation,header=T,stringsAsFactors=F))
	#sort out annotation and turn into a GRanges object. remove gm genes
	names(annotation)[4:5] <- c("start","end")
	anno.GRange<-  makeGRangesFromDataFrame(annotation,keep.extra.columns=T)
	fixed.gene.names <- apply(results.df, MAR=1,FUN=function(x) gene.names.query(x[10],x[15],x[16],anno.GRange)) 
	fixed.gene.names <- unlist(GRangesList(fixed.gene.names))
	fixed.gene.names <- data.frame(row.names=names(fixed.gene.names),gene <- mcols(fixed.gene.names)[,2], fixed.strand=strand(fixed.gene.names))
	names(fixed.gene.names)[1] <- "fixed.gene.id"
	#fixed.gene.names$fixed.strand <- gsub("+","1",fixed.gene.names$fixed.strand,fixed=T)
	#fixed.gene.names$fixed.strand <- gsub("-","-1",fixed.gene.names$fixed.strand,fixed=T)
	return(fixed.gene.names)
}

########
## 9. ## 
########
# Classify/annotate the tags with small fold changes or low read depths (of the number of canonical junction reads)

cryptic.classifier <- function(crypt.counts){
	classer <- function(subset.df,classification){ # helper function to annotate with classification of read depth/fold change
		if(dim(subset.df)[1] > 0){
			subset.df$class <- classification
			return(subset.df$class)
		}
	}
	FEW.READS.UP <- subset(crypt.counts, (as.numeric(log2FoldChange)) > 0  & canonical.control.mean.SJ < min.canonical.control.SJs ) # actual filter
	FEW.READS.UP$class <- classer(FEW.READS.UP,"FEW.READS.UP")

	FEW.READS.DOWN <- subset(crypt.counts, (as.numeric(log2FoldChange)) < 0  & canonical.control.mean.SJ < min.canonical.control.SJs ) # actual filter
	FEW.READS.DOWN$class <- classer(FEW.READS.DOWN,"FEW.READS.DOWN")

	SMALL.FOLDCHANGE.UP <- subset(crypt.counts, as.numeric(log2FoldChange) < 0.6 & as.numeric(log2FoldChange) > 0 ) # actual filter
	SMALL.FOLDCHANGE.UP$class <- classer(SMALL.FOLDCHANGE.UP, "SMALL.FOLDCHANGE.UP")

	SMALL.FOLDCHANGE.DOWN <- subset(crypt.counts, as.numeric(log2FoldChange) > -0.6 & as.numeric(log2FoldChange) < 0 ) # actual filter
	SMALL.FOLDCHANGE.DOWN$class <- classer(SMALL.FOLDCHANGE.DOWN, "SMALL.FOLDCHANGE.DOWN")
	
	CLEAN.UP <- subset(crypt.counts, as.numeric(log2FoldChange) >= 0.6  & canonical.control.mean.SJ >= min.canonical.control.SJs )
	
	SJ.UNSUPPORTED.UP <- subset(CLEAN.UP,
		(CLEAN.UP$upstream.case.mean.SJ < 1 &
		CLEAN.UP$downstream.case.mean.SJ < 1 ) | 
		CLEAN.UP$PSI.class == "TOO.SMALL.PSI" )
	SJ.UNSUPPORTED.UP$class <- classer(SJ.UNSUPPORTED.UP,"SJ.UNSUPPORTED.UP")

	SJ.SUPPORTED.UP <- subset(CLEAN.UP,
		(CLEAN.UP$upstream.case.mean.SJ >= 1 |
		CLEAN.UP$downstream.case.mean.SJ >= 1) &
		CLEAN.UP$PSI.class != "TOO.SMALL.PSI" )
	SJ.SUPPORTED.UP$class <- classer(SJ.SUPPORTED.UP,"SJ.SUPPORTED.UP")
	
	CLEAN.DOWN <- subset(crypt.counts, as.numeric(log2FoldChange) <= -0.6 & canonical.control.mean.SJ >= min.canonical.control.SJs )
	SJ.UNSUPPORTED.DOWN <- subset(CLEAN.DOWN,
	  ( CLEAN.DOWN$upstream.case.mean.SJ < 1 &
		CLEAN.DOWN$downstream.case.mean.SJ < 1 ) |
		CLEAN.DOWN$PSI.class == "TOO.SMALL.PSI" )
	SJ.UNSUPPORTED.DOWN$class <- classer(SJ.UNSUPPORTED.DOWN,"SJ.UNSUPPORTED.DOWN")

	SJ.SUPPORTED.DOWN <- subset(CLEAN.DOWN,
		(CLEAN.DOWN$upstream.case.mean.SJ >= 1 |
		CLEAN.DOWN$downstream.case.mean.SJ >= 1) &
		CLEAN.DOWN$PSI.class != "TOO.SMALL.PSI" )
	SJ.SUPPORTED.DOWN$class <- classer(SJ.SUPPORTED.DOWN,"SJ.SUPPORTED.DOWN")

	crypt.counts.classified <- rbind(SJ.SUPPORTED.UP,SJ.SUPPORTED.DOWN,SJ.UNSUPPORTED.UP,SJ.UNSUPPORTED.DOWN,SMALL.FOLDCHANGE.UP,SMALL.FOLDCHANGE.DOWN,FEW.READS.UP,FEW.READS.DOWN)

	return(crypt.counts.classified)
}

#########
## 10. ## 
#########
# This function, based on difference in delta PSI, classifies whether it is a 5'/3' bias or cassette exon
cryptic.PSI.classifier <- function(cryptic.exons){
  classer <- function(subset.df,classification){ # another helper function to print the classification
    if(dim(subset.df)[1] > 0){
      subset.df$class <- classification
      return(subset.df$class)
    }
  }
# Classify the tags with small fold changes or low read depths

  FIVEPRIME.BIAS <- subset(cryptic.exons, (as.numeric(upstream_delta_psi) >= PSI.threshold) & (as.numeric(downstream_delta_psi) < PSI.threshold) )
  FIVEPRIME.BIAS$PSI.class <- classer(FIVEPRIME.BIAS,"FIVEPRIME.BIAS")
  
  THREEPRIME.BIAS <- subset(cryptic.exons, (as.numeric(upstream_delta_psi) < PSI.threshold) & (as.numeric(downstream_delta_psi) >= PSI.threshold) )
  THREEPRIME.BIAS$PSI.class <- classer(THREEPRIME.BIAS,"THREEPRIME.BIAS")
  
  CASSETTE.LIKE <- subset(cryptic.exons, (as.numeric(upstream_delta_psi) >= PSI.threshold ) & (as.numeric(downstream_delta_psi) >= PSI.threshold) )
  CASSETTE.LIKE$PSI.class <- classer(CASSETTE.LIKE,"CASSETTE.LIKE")
  
  TOO.SMALL.PSI <- subset(cryptic.exons, (as.numeric(upstream_delta_psi) <= PSI.threshold )& (as.numeric(downstream_delta_psi) <= PSI.threshold) )
  TOO.SMALL.PSI$PSI.class <- classer(TOO.SMALL.PSI,"TOO.SMALL.PSI")
  
  cryptic.exons.PSI.classified <- rbind(FIVEPRIME.BIAS,THREEPRIME.BIAS,CASSETTE.LIKE,TOO.SMALL.PSI)
  return(cryptic.exons.PSI.classified)
}
         
				  
## Begin analysis
#############
## CONTROL ##
#############
# we start with controls first, to discover the canonical SJs flanking each cryptic tag.

library("dplyr")
library("data.table")
library("stringr")

# merge SJ tab files (from STAR)
STAR_control_SJ_list <- c("ctr1.SJ.out.tab","ctr2.SJ.out.tab","ctr3.SJ.out.tab") 
STAR_control_total_SJ_counts <- merge.SJ.files(STAR_control_SJ_list)
head(STAR_control_total_SJ_counts)
write.table(STAR_control_total_SJ_counts, "STAR_tdp43_cnp_scn_SJs_control.tab", quote=F,sep="\t",row.names=F)
				  
# Load in dexseq results containing cryptic tag information
STAR_dexseq.res <- "24may_scn_SignificantExons_STAR.csv"
STAR_crypt.res <- as.data.frame(fread(STAR_dexseq.res))
nrow(STAR_crypt.res)

# obtaining cryptic exons with FDR < 0.05
STAR_crypt.res <- subset(STAR_crypt.res,grepl("i",exonID) & FDR < 0.05)
nrow(STAR_crypt.res) #492 rows

# Finding canonical junctions flanking each cryptic tag
STAR_canonical_results_control <- canonical_junction_detector(SJ.summary = STAR_control_total_SJ_counts, 
                                                         results.df = STAR_crypt.res, 
                                                         mode = "discovery") 

head(STAR_canonical_results_control)
head(STAR_crypt.res)

# Add the CANONICAL splice site coordinates to the cryptic exon results
STAR_crypt.res$canonical.chr <- STAR_canonical_results_control$canonical.chr[match(rownames(STAR_crypt.res),rownames(STAR_canonical_results_control))]
STAR_crypt.res$canonical.start <- STAR_canonical_results_control$canonical.start[match(rownames(STAR_crypt.res),rownames(STAR_canonical_results_control))]
STAR_crypt.res$canonical.end <- STAR_canonical_results_control$canonical.end[match(rownames(STAR_crypt.res),rownames(STAR_canonical_results_control))]
STAR_crypt.res$canonical.strand <- STAR_canonical_results_control$canonical.strand[match(rownames(STAR_crypt.res),rownames(STAR_canonical_results_control))]
STAR_crypt.res$canonical.intron.motif <- STAR_canonical_results_control$intron.motif[match(rownames(STAR_crypt.res),rownames(STAR_canonical_results_control))]
STAR_crypt.res$canonical.control.mean.SJ <- STAR_canonical_results_control$canonical.unique.count[match(rownames(STAR_crypt.res),rownames(STAR_canonical_results_control))]
STAR_crypt.res$canonical.control.mean.SJ <- STAR_crypt.res$canonical.control.mean.SJ / length(STAR_control_SJ_list)
nrow(STAR_crypt.res) #492
STAR_crypt.res <- subset(STAR_crypt.res,!is.na(STAR_crypt.res$canonical.start)) # 486 left, removed 6 rows
length(STAR_control_SJ_list)
STAR_NA_crypt.res <- subset(STAR_crypt.res,is.na(STAR_crypt.res$canonical.chr))


## COUNTING UPSTREAM JUNCTIONS IN CONTROLS
# Why is this necessary? We are assuming there is a cryptic exon within control too?
# Anyway this finds the exact 5' cryptic splice sites for control cryptic tags, based on SJ splicing into the cryptic tag
STAR_upstream_results_control <- bridging_junction_finder(SJ.summary = STAR_control_total_SJ_counts,
                                                     results.df = STAR_crypt.res,
                                                     query.type = "upstream")
head(STAR_upstream_results_control)
nrow(STAR_upstream_results_control) #158 rows

head(STAR_crypt.res)
nrow(STAR_crypt.res) # 486 rows

# Add the CRYPTIC 5' splice site coordinates to cryptic exon df
STAR_crypt.res$upstream.control.canonical.start <- STAR_upstream_results_control$canonical.start[match(rownames(STAR_crypt.res),rownames(STAR_upstream_results_control))]
STAR_crypt.res$upstream.control.cryptic.5prime <- STAR_upstream_results_control$cryptic.5prime[match(rownames(STAR_crypt.res),rownames(STAR_upstream_results_control))]
STAR_crypt.res$upstream.control.strand <- STAR_upstream_results_control$upstream.strand[match(rownames(STAR_crypt.res),rownames(STAR_upstream_results_control))]
STAR_crypt.res$upstream.control.intron.motif <- STAR_upstream_results_control$intron.motif[match(rownames(STAR_crypt.res),rownames(STAR_upstream_results_control))]
STAR_crypt.res$upstream.control.mean.SJ <- STAR_upstream_results_control$upstream.unique.count[match(rownames(STAR_crypt.res),rownames(STAR_upstream_results_control))]
STAR_crypt.res$upstream.control.mean.SJ <- STAR_crypt.res$upstream.control.mean.SJ / length(STAR_control_SJ_list)


## COUNTING DOWNSTREAM JUNCTIONS IN CONTROLS
# Repeat the same for 3' cryptic site
STAR_downstream_results_control <- bridging_junction_finder(SJ.summary = STAR_control_total_SJ_counts, 
                                                       results.df = STAR_crypt.res, 
                                                       query.type = "downstream")
head(STAR_downstream_results_control)
nrow(STAR_downstream_results_control) #159 rows

head(STAR_crypt.res)
nrow(STAR_crypt.res) # 486 rows

# Add the CRYPTIC 3' splice site coordinates to cryptic exon df
STAR_crypt.res$downstream.control.cryptic.3prime <- STAR_downstream_results_control$cryptic.3prime[match(rownames(STAR_crypt.res),rownames(STAR_downstream_results_control))]
STAR_crypt.res$downstream.control.canonical.end <- STAR_downstream_results_control$canonical.end[match(rownames(STAR_crypt.res),rownames(STAR_downstream_results_control))]
STAR_crypt.res$downstream.control.strand <- STAR_downstream_results_control$downstream.strand[match(rownames(STAR_crypt.res),rownames(STAR_downstream_results_control))]
STAR_crypt.res$downstream.control.intron.motif <- STAR_downstream_results_control$intron.motif[match(rownames(STAR_crypt.res),rownames(STAR_downstream_results_control))]
STAR_crypt.res$downstream.control.mean.SJ <- STAR_downstream_results_control$downstream.unique.count[match(rownames(STAR_crypt.res),rownames(STAR_downstream_results_control))]
STAR_crypt.res$downstream.control.mean.SJ <- STAR_crypt.res$downstream.control.mean.SJ / length(STAR_control_SJ_list)

# Once Controls are done, we move on to Case

##########
## CASE ##
##########

# merge SJ tab files (from STAR)
STAR_case_SJ_list <- c("cKO1.SJ.out.tab","cKO2.SJ.out.tab","cKO3.SJ.out.tab") 
STAR_case_total_SJ_counts <- merge.SJ.files(STAR_case_SJ_list)

head(STAR_case_total_SJ_counts)
write.table(STAR_case_total_SJ_counts, "STAR_tdp43_cnp_scn_SJs_case.tab", quote=F,sep="\t",row.names=F)

# STAR_crypt.res is now augmented with the coordinates discovered previously for the canonical sites
# This function looks for canonical junctions in the case dataset that exactly match those discovered in the control dataset.
# Why do we need to do this? --> Because we want to find the SJ splicing from the canonical splice sites into the cryptic tag region

head(STAR_case_total_SJ_counts)
head(STAR_crypt.res[10])
head(STAR_crypt.res[11])
head(STAR_crypt.res[12])

# We want to find junctions in the merged SJ summary file that matches exactly the canonical splice junctions discovered earlier in Controls analysis
STAR_canonical_results_case <- canonical_junction_detector(SJ.summary = STAR_case_total_SJ_counts, 
                                                      results.df = STAR_crypt.res, 
                                                      mode = "replication")
head(STAR_canonical_results_case)
# Add numbers of canonical splicing in cases to results
STAR_crypt.res$canonical.case.mean.SJ <- STAR_canonical_results_case$canonical.unique.count[match(rownames(STAR_crypt.res),rownames(STAR_canonical_results_case))]
STAR_crypt.res$canonical.case.mean.SJ <- STAR_crypt.res$canonical.case.mean.SJ / length(STAR_control_SJ_list)

# Discovering 5' EXACT cryptic exon sites in case
STAR_upstream_results_case <- bridging_junction_finder(SJ.summary = STAR_case_total_SJ_counts, 
                                                  results.df = STAR_crypt.res, 
                                                  query.type = "upstream")
head(STAR_upstream_results_case)
head(upstream_results_case)

# Annotate cryptic exon dataset with this new information about cryptic sites
STAR_crypt.res$upstream.case.canonical.start <- STAR_upstream_results_case$canonical.start[match(rownames(STAR_crypt.res),rownames(STAR_upstream_results_case))]
STAR_crypt.res$upstream.case.cryptic.5prime <- STAR_upstream_results_case$cryptic.5prime[match(rownames(STAR_crypt.res),rownames(STAR_upstream_results_case))]
STAR_crypt.res$upstream.case.strand <- STAR_upstream_results_case$upstream.strand[match(rownames(STAR_crypt.res),rownames(STAR_upstream_results_case))]
STAR_crypt.res$upstream.case.intron.motif <- STAR_upstream_results_case$intron.motif[match(rownames(STAR_crypt.res),rownames(STAR_upstream_results_case))]
STAR_crypt.res <- STAR_crypt.res[,1:35] # had to subset this and re-run since I forgot to run certain columns previously
head(STAR_crypt.res)
STAR_crypt.res$upstream.case.mean.SJ <- STAR_upstream_results_case$upstream.unique.count[match(rownames(STAR_crypt.res),rownames(STAR_upstream_results_case))]
STAR_crypt.res$upstream.case.mean.SJ <- STAR_crypt.res$upstream.case.mean.SJ / length(STAR_case_SJ_list)

# Discovering 3' EXACT cryptic exon sites in case
# need to find the counts of canonical junctions that were assigned by canonical_junction_detector in the case dataset.
#This function assumes that the results file has been appended with the canonical start and end coordinates at positions XXX and YYY respectively.
STAR_downstream_results_case <- bridging_junction_finder(SJ.summary = STAR_case_total_SJ_counts, 
                                                    results.df = STAR_crypt.res, 
                                                    query.type = "downstream")
head(STAR_downstream_results_case)
head(downstream_results_case)
head(STAR_crypt.res)


# Annotate cryptic exon dataset with this new information about cryptic sites			  
STAR_crypt.res$downstream.case.cryptic.3prime <- STAR_downstream_results_case$cryptic.3prime[match(rownames(STAR_crypt.res),rownames(STAR_downstream_results_case))]
STAR_crypt.res$downstream.case.canonical.end <- STAR_downstream_results_case$canonical.end[match(rownames(STAR_crypt.res),rownames(STAR_downstream_results_case))]
STAR_crypt.res$downstream.case.strand <- STAR_downstream_results_case$downstream.strand[match(rownames(STAR_crypt.res),rownames(STAR_downstream_results_case))]
STAR_crypt.res$downstream.case.intron.motif <- STAR_downstream_results_case$intron.motif[match(rownames(STAR_crypt.res),rownames(STAR_downstream_results_case))]
STAR_crypt.res$downstream.case.mean.SJ <- STAR_downstream_results_case$downstream.unique.count[match(rownames(STAR_crypt.res),rownames(STAR_downstream_results_case))]
STAR_crypt.res$downstream.case.mean.SJ <- STAR_crypt.res$downstream.case.mean.SJ / length(STAR_case_SJ_list)

######################################
## CALCULATION BASED ON JUNCTION COUNTS
######################################
# create ratios of cryptic to canonical splicing
# convert all NAs to zeros
ncol(STAR_crypt.res)
head(STAR_crypt.res)
STAR_crypt.res[is.na(STAR_crypt.res)] <- "0"
head(STAR_crypt.res[,c(20,25,30,31,36,41)])
STAR_crypt.res[,c(20,25,30,31,36,41)] <- as.numeric(unlist(STAR_crypt.res[,c(20,25,30,31,36,41)]))

# calculate ratios
STAR_crypt.res$control.upstream.ratio <- STAR_crypt.res$upstream.control.mean.SJ / STAR_crypt.res$canonical.control.mean.SJ
STAR_crypt.res$control.downstream.ratio <- STAR_crypt.res$downstream.control.mean.SJ / STAR_crypt.res$canonical.control.mean.SJ

STAR_crypt.res$case.upstream.ratio <- STAR_crypt.res$upstream.case.mean.SJ / STAR_crypt.res$canonical.case.mean.SJ
STAR_crypt.res$case.downstream.ratio <- STAR_crypt.res$downstream.case.mean.SJ / STAR_crypt.res$canonical.case.mean.SJ

head(STAR_crypt.res)
nrow(STAR_crypt.res)
ncol(STAR_crypt.res)
head(crypt.res)
head(STAR_fixed.gene.names)
				  
# Annotate cryptic exons with gene names
STAR_fixed.gene.names <- fix.gene.names(STAR_crypt.res,"mouse_annotation_infos_biomart.txt")
STAR_crypt.res$fix.gene.names <- as.character(STAR_fixed.gene.names$fixed.gene.id[match(row.names(STAR_crypt.res),row.names(STAR_fixed.gene.names))])
STAR_crypt.res$fix.gene.names <- ifelse(test = is.na(STAR_crypt.res$fix.gene.names),yes = STAR_crypt.res$external_gene_id,no = STAR_crypt.res$fix.gene.names)

STAR_crypt.res$fix.strand <- as.character(STAR_fixed.gene.names$fixed.strand[match(row.names(STAR_crypt.res),row.names(STAR_fixed.gene.names))])
STAR_crypt.res$fix.strand <- ifelse(test=is.na(STAR_crypt.res$fix.strand),yes=STAR_crypt.res$strand,no=as.character(STAR_crypt.res$fix.strand))
STAR_crypt.res <- subset(STAR_crypt.res, !is.na(STAR_crypt.res$fix.gene.names))

# Calculating upstream delta psi
STAR_crypt.res$upstream_delta_psi <- ifelse(STAR_crypt.res$fix.strand == "+",
                                       yes = (STAR_crypt.res$upstream.case.mean.SJ / (STAR_crypt.res$upstream.case.mean.SJ + STAR_crypt.res$canonical.case.mean.SJ) )  - (STAR_crypt.res$upstream.control.mean.SJ / (STAR_crypt.res$upstream.control.mean.SJ + STAR_crypt.res$canonical.control.mean.SJ) ),
                                       no = (STAR_crypt.res$downstream.case.mean.SJ / (STAR_crypt.res$downstream.case.mean.SJ + STAR_crypt.res$canonical.case.mean.SJ) ) - (STAR_crypt.res$downstream.control.mean.SJ / (STAR_crypt.res$downstream.control.mean.SJ + STAR_crypt.res$canonical.control.mean.SJ) ) )

head(STAR_crypt.res)

# Calculating downstream delta psi
STAR_crypt.res$downstream_delta_psi <- ifelse(STAR_crypt.res$fix.strand == "+",
                                         yes = (STAR_crypt.res$downstream.case.mean.SJ / (STAR_crypt.res$downstream.case.mean.SJ + STAR_crypt.res$canonical.case.mean.SJ) ) - (STAR_crypt.res$downstream.control.mean.SJ / (STAR_crypt.res$downstream.control.mean.SJ + STAR_crypt.res$canonical.control.mean.SJ) ), 
                                         no =  (STAR_crypt.res$upstream.case.mean.SJ / (STAR_crypt.res$upstream.case.mean.SJ + STAR_crypt.res$canonical.case.mean.SJ) )  - (STAR_crypt.res$upstream.control.mean.SJ / (STAR_crypt.res$upstream.control.mean.SJ + STAR_crypt.res$canonical.control.mean.SJ) ) )


PSI.threshold <- 0.05
min.canonical.control.SJs <- 5
FDR_threshold <- 0.05
easy_threshold <- FALSE

# apply delta PSI classification
STAR_crypt.res <- cryptic.PSI.classifier(STAR_crypt.res)

# apply classification function
STAR_crypt.res.classified <- cryptic.classifier(STAR_crypt.res)
nrow(STAR_crypt.res.classified)
STAR_crypt.res.classified <- subset(STAR_crypt.res.classified, !is.na(class))
nrow(STAR_crypt.res)
head(STAR_crypt.res.classified)

# for graphing, slice out whether it is cryptic exon, intron retention, possible or noise.
library("stringr")
STAR_crypt.res.classified$family <- str_split_fixed(STAR_crypt.res.classified$class, pattern = "[.]UP|[.]DOWN",3)[,1]
head(STAR_crypt.res.classified)

## graphing
STAR_graph1 <- paste0("28may_sciatic_nerve_classification.pdf")
STAR_graph1.5 <- paste0("28may_sciatic_nerve_classification_counts.pdf")
STAR_graph2 <- paste0("28may_sciatic_nerve_splice_junctions.pdf")
STAR_title.code <- "Chang_et_al_2021_TDP43_cKO_sciatic_nerve"
min.canonical.control.SJs <- 5
STAR_graph_title <- paste0( gsub("_"," ", title.code),"\n(","Mouse sciatic nerve",")\nPSI.threshold = ",PSI.threshold,"\ncanonical SJs min = ",min.canonical.control.SJs )
