#############
## CONTROL ##
#############

library("dplyr")
library("data.table")
library("stringr")

STAR_control_SJ_list <- c("ctr1.SJ.out.tab","ctr2.SJ.out.tab","ctr3.SJ.out.tab") 

STAR_control_total_SJ_counts <- merge.SJ.files(STAR_control_SJ_list)
head(STAR_control_total_SJ_counts)
write.table(STAR_control_total_SJ_counts, "STAR_tdp43_cnp_scn_SJs_control.tab", quote=F,sep="\t",row.names=F)
STAR_dexseq.res <- "24may_scn_SignificantExons_STAR.csv"
STAR_crypt.res <- as.data.frame(fread(STAR_dexseq.res))
nrow(STAR_crypt.res)

# obtaining cryptic exons with FDR < 0.05
STAR_crypt.res <- subset(STAR_crypt.res,grepl("i",exonID) & FDR < 0.05)
nrow(STAR_crypt.res) #492 rows

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
# need to find the counts of canonical junctions that were assigned by canonical_junction_detector in the case dataset.
#This function assumes that the results file has been appended with the canonical start and end coordinates at positions ??? and ??? respectively.
STAR_upstream_results_control <- bridging_junction_finder(SJ.summary = STAR_control_total_SJ_counts,
                                                     results.df = STAR_crypt.res,
                                                     query.type = "upstream")
head(STAR_upstream_results_control)
nrow(STAR_upstream_results_control) #158 rows

head(STAR_crypt.res)
nrow(STAR_crypt.res) # 486 rows

STAR_crypt.res$upstream.control.canonical.start <- STAR_upstream_results_control$canonical.start[match(rownames(STAR_crypt.res),rownames(STAR_upstream_results_control))]
STAR_crypt.res$upstream.control.cryptic.5prime <- STAR_upstream_results_control$cryptic.5prime[match(rownames(STAR_crypt.res),rownames(STAR_upstream_results_control))]
STAR_crypt.res$upstream.control.strand <- STAR_upstream_results_control$upstream.strand[match(rownames(STAR_crypt.res),rownames(STAR_upstream_results_control))]
STAR_crypt.res$upstream.control.intron.motif <- STAR_upstream_results_control$intron.motif[match(rownames(STAR_crypt.res),rownames(STAR_upstream_results_control))]
STAR_crypt.res$upstream.control.mean.SJ <- STAR_upstream_results_control$upstream.unique.count[match(rownames(STAR_crypt.res),rownames(STAR_upstream_results_control))]
STAR_crypt.res$upstream.control.mean.SJ <- STAR_crypt.res$upstream.control.mean.SJ / length(STAR_control_SJ_list)


## COUNTING DOWNSTREAM JUNCTIONS IN CONTROLS
STAR_downstream_results_control <- bridging_junction_finder(SJ.summary = STAR_control_total_SJ_counts, 
                                                       results.df = STAR_crypt.res, 
                                                       query.type = "downstream")
head(STAR_downstream_results_control)
nrow(STAR_downstream_results_control) #159 rows

head(STAR_crypt.res)
nrow(STAR_crypt.res) # 486 rows

STAR_crypt.res$downstream.control.cryptic.3prime <- STAR_downstream_results_control$cryptic.3prime[match(rownames(STAR_crypt.res),rownames(STAR_downstream_results_control))]
STAR_crypt.res$downstream.control.canonical.end <- STAR_downstream_results_control$canonical.end[match(rownames(STAR_crypt.res),rownames(STAR_downstream_results_control))]
STAR_crypt.res$downstream.control.strand <- STAR_downstream_results_control$downstream.strand[match(rownames(STAR_crypt.res),rownames(STAR_downstream_results_control))]
STAR_crypt.res$downstream.control.intron.motif <- STAR_downstream_results_control$intron.motif[match(rownames(STAR_crypt.res),rownames(STAR_downstream_results_control))]
STAR_crypt.res$downstream.control.mean.SJ <- STAR_downstream_results_control$downstream.unique.count[match(rownames(STAR_crypt.res),rownames(STAR_downstream_results_control))]
STAR_crypt.res$downstream.control.mean.SJ <- STAR_crypt.res$downstream.control.mean.SJ / length(STAR_control_SJ_list)

##########
## CASE ##
##########

STAR_case_SJ_list <- c("cKO1.SJ.out.tab","cKO2.SJ.out.tab","cKO3.SJ.out.tab") 
STAR_case_total_SJ_counts <- merge.SJ.files(STAR_case_SJ_list)

head(STAR_case_total_SJ_counts)
write.table(STAR_case_total_SJ_counts, "STAR_tdp43_cnp_scn_SJs_case.tab", quote=F,sep="\t",row.names=F)

# This function looks for junctions in the case dataset that exactly match those discovered in the control dataset.
# STAR_crypt.res is now augmented with the coordinates discovered previously for the canonical sites
head(STAR_case_total_SJ_counts)
head(STAR_crypt.res[10])
head(STAR_crypt.res[11])
head(STAR_crypt.res[12])
STAR_canonical_results_case <- canonical_junction_detector(SJ.summary = STAR_case_total_SJ_counts, 
                                                      results.df = STAR_crypt.res, 
                                                      mode = "replication")
head(STAR_canonical_results_case)
# Add numbers of canonical splicing in cases to results
STAR_crypt.res$canonical.case.mean.SJ <- STAR_canonical_results_case$canonical.unique.count[match(rownames(STAR_crypt.res),rownames(STAR_canonical_results_case))]
STAR_crypt.res$canonical.case.mean.SJ <- STAR_crypt.res$canonical.case.mean.SJ / length(STAR_control_SJ_list)

# upstream junctions in cases
STAR_upstream_results_case <- bridging_junction_finder(SJ.summary = STAR_case_total_SJ_counts, 
                                                  results.df = STAR_crypt.res, 
                                                  query.type = "upstream")
head(STAR_upstream_results_case)
head(upstream_results_case)

STAR_crypt.res$upstream.case.canonical.start <- STAR_upstream_results_case$canonical.start[match(rownames(STAR_crypt.res),rownames(STAR_upstream_results_case))]
STAR_crypt.res$upstream.case.cryptic.5prime <- STAR_upstream_results_case$cryptic.5prime[match(rownames(STAR_crypt.res),rownames(STAR_upstream_results_case))]
STAR_crypt.res$upstream.case.strand <- STAR_upstream_results_case$upstream.strand[match(rownames(STAR_crypt.res),rownames(STAR_upstream_results_case))]
STAR_crypt.res$upstream.case.intron.motif <- STAR_upstream_results_case$intron.motif[match(rownames(STAR_crypt.res),rownames(STAR_upstream_results_case))]
STAR_crypt.res <- STAR_crypt.res[,1:35] # had to subset this and re-run since I forgot to run certain columns previously
head(STAR_crypt.res)
STAR_crypt.res$upstream.case.mean.SJ <- STAR_upstream_results_case$upstream.unique.count[match(rownames(STAR_crypt.res),rownames(STAR_upstream_results_case))]
STAR_crypt.res$upstream.case.mean.SJ <- STAR_crypt.res$upstream.case.mean.SJ / length(STAR_case_SJ_list)

# downstream junctions in cases
# need to find the counts of canonical junctions that were assigned by canonical_junction_detector in the case dataset.
#This function assumes that the results file has been appended with the canonical start and end coordinates at positions ??? and ??? respectively.

STAR_downstream_results_case <- bridging_junction_finder(SJ.summary = STAR_case_total_SJ_counts, 
                                                    results.df = STAR_crypt.res, 
                                                    query.type = "downstream")
head(STAR_downstream_results_case)
head(downstream_results_case)
head(STAR_crypt.res)

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
STAR_fixed.gene.names <- fix.gene.names(STAR_crypt.res,"mouse_annotation_infos_biomart.txt")
STAR_crypt.res$fix.gene.names <- as.character(STAR_fixed.gene.names$fixed.gene.id[match(row.names(STAR_crypt.res),row.names(STAR_fixed.gene.names))])
STAR_crypt.res$fix.gene.names <- ifelse(test = is.na(STAR_crypt.res$fix.gene.names),yes = STAR_crypt.res$external_gene_id,no = STAR_crypt.res$fix.gene.names)

STAR_crypt.res$fix.strand <- as.character(STAR_fixed.gene.names$fixed.strand[match(row.names(STAR_crypt.res),row.names(STAR_fixed.gene.names))])
STAR_crypt.res$fix.strand <- ifelse(test=is.na(STAR_crypt.res$fix.strand),yes=STAR_crypt.res$strand,no=as.character(STAR_crypt.res$fix.strand))
STAR_crypt.res <- subset(STAR_crypt.res, !is.na(STAR_crypt.res$fix.gene.names))

STAR_crypt.res$upstream_delta_psi <- ifelse(STAR_crypt.res$fix.strand == "+",
                                       yes = (STAR_crypt.res$upstream.case.mean.SJ / (STAR_crypt.res$upstream.case.mean.SJ + STAR_crypt.res$canonical.case.mean.SJ) )  - (STAR_crypt.res$upstream.control.mean.SJ / (STAR_crypt.res$upstream.control.mean.SJ + STAR_crypt.res$canonical.control.mean.SJ) ),
                                       no = (STAR_crypt.res$downstream.case.mean.SJ / (STAR_crypt.res$downstream.case.mean.SJ + STAR_crypt.res$canonical.case.mean.SJ) ) - (STAR_crypt.res$downstream.control.mean.SJ / (STAR_crypt.res$downstream.control.mean.SJ + STAR_crypt.res$canonical.control.mean.SJ) ) )

head(STAR_crypt.res)

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
