## Splice junction rediscovery 
Input: Cryptic tag

Output: Exact splice junctions of cryptic exon

![image](https://user-images.githubusercontent.com/68455070/124528073-ddc35400-de39-11eb-95c8-43949a50863b.png)

### 1. Merge SJ.out.tab files for controls
```r
STAR_control_SJ_list <- c("ctr1.SJ.out.tab","ctr2.SJ.out.tab","ctr3.SJ.out.tab") 
STAR_control_total_SJ_counts <- merge.SJ.files(STAR_control_SJ_list)
```
`STAR_control_SJ_list`:

![image](https://user-images.githubusercontent.com/68455070/124410094-4815bf00-dd7c-11eb-9aff-b4cda71f0cc8.png)

`STAR_control_total_SJ_counts`:

![image](https://user-images.githubusercontent.com/68455070/124410070-3502ef00-dd7c-11eb-9949-7922b61618c8.png)

### 2. Load in dexseq results containing cryptic tag information
```r
STAR_dexseq.res <- "24may_scn_SignificantExons_STAR.csv"
STAR_crypt.res <- as.data.frame(fread(STAR_dexseq.res))
# obtaining cryptic exons with FDR < 0.05
STAR_crypt.res <- subset(STAR_crypt.res,grepl("i",exonID) & FDR < 0.05)
```
`STAR_crypt.res`:

![image](https://user-images.githubusercontent.com/68455070/124410229-99be4980-dd7c-11eb-94e9-3a4827593edc.png)

### 3. Finding canonical junctions flanking each cryptic tag
```r
STAR_canonical_results_control <- canonical_junction_detector(SJ.summary = STAR_control_total_SJ_counts, 
                                                         results.df = STAR_crypt.res,
                                                         mode = "discovery")
```

Functions involved:

`canonical_junction_detector` function: 
```r
canonical_junction_detector <- function(SJ.summary,results.df,mode="discovery"){
	GRanges_object <-  makeGRangesFromDataFrame(SJ.summary,keep.extra.columns=T) # convert the combined SJ file into a GRanges obj
	if(mode == "discovery"){
		junctions.list <- apply(results.df, MAR=1,FUN=function(x) canonical_junction_query(x[10],x[11],x[12], GRanges_object)) # this is for control
		} # apply canonical_junction_query to each cryptic tag
	if(mode == "replication"){
		junctions.list <- apply(results.df, MAR=1,FUN=function(x) canonical_junction_replication(x[10],x[15],x[16], GRanges_object)) # this is for case
	}# apply canonical_junction_replication to each cryptic tag
	#output is a list of GRange objects - unuseable.
	junctions.list <- unlist(GRangesList(junctions.list)) 
	#convert into a dataframe, extracting the relevent information from the GRanges object.
	canonical.df <- data.frame(row.names=names(junctions.list),
			canonical.chr=seqnames(junctions.list),
			canonical.start=start(junctions.list),
			canonical.end=end(junctions.list),
			canonical.unique.count = score(junctions.list),
			canonical.strand = strand(junctions.list),
			intron.motif = mcols(junctions.list)[3])
	return(canonical.df)
}
```
`canonical_junction_query`:
```r
# this function below allows you to specify which GRange object to query
# this function outputs the canonical junction flanking the cryptic tag (since there could be multiple)
canonical_junction_query <- function(CE.chr,CE.start,CE.end, SJ.GRange){ # iterate each SJ through the GRanges of SJs
        junction <- SJ.GRange[seqnames(SJ.GRange)==CE.chr & start(SJ.GRange) <= as.numeric(CE.start)+1 & end(SJ.GRange) >= as.numeric(CE.end)-1] # obtain that particular SJ
        junction <- head(junction[order(score(junction),decreasing=T)],1) # we only want the most abundantly supported SJ
        return(junction)
}
```

Output `STAR_canonical_results_control`:

![image](https://user-images.githubusercontent.com/68455070/124410451-0f2a1a00-dd7d-11eb-89fa-76b0cfcc0a5e.png)

Step-by-step function breakdown: 

`STAR_crypt.res[10:12]`:

![image](https://user-images.githubusercontent.com/68455070/124412342-dbe98a00-dd80-11eb-8d84-3663c0f5e552.png)

`GRanges_object <-  makeGRangesFromDataFrame(SJ.summary,keep.extra.columns=T)`:

![image](https://user-images.githubusercontent.com/68455070/124412538-55817800-dd81-11eb-9c14-bcf2548c7835.png)

Testing on 1 cryptic tag:

`junction <- SJ.GRange[seqnames(SJ.GRange)==CE.chr & start(SJ.GRange) <= as.numeric(CE.start)+1 & end(SJ.GRange) >= as.numeric(CE.end)-1]`:

![image](https://user-images.githubusercontent.com/68455070/124412897-1d2e6980-dd82-11eb-878a-5fe680469224.png)

`junction <- head(junction[order(score(junction),decreasing=T)],1)`:

![image](https://user-images.githubusercontent.com/68455070/124412992-4c44db00-dd82-11eb-99c7-80b51e159f58.png)

Testing on several cryptic tags: 

`junctions.list <- apply(head(STAR_crypt.res), MAR=1,FUN=function(x) canonical_junction_query(x[10],x[11],x[12], GRanges_object))`:

![image](https://user-images.githubusercontent.com/68455070/124413958-52d45200-dd84-11eb-9fe0-7cde5bcd588b.png)

`junctions.list <- unlist(GRangesList(junctions.list))`: (collapse a GRangesList into a single GRanges object)
![image](https://user-images.githubusercontent.com/68455070/124414094-8fa04900-dd84-11eb-8ad7-8fc801872514.png)

```r
canonical.df <- data.frame(row.names=names(junctions.list),
			canonical.chr=seqnames(junctions.list),
			canonical.start=start(junctions.list),
			canonical.end=end(junctions.list),
			canonical.unique.count = score(junctions.list),
			canonical.strand = strand(junctions.list),
			intron.motif = mcols(junctions.list)[3])
```

Generating a dataframe from the GRanges result above:

![image](https://user-images.githubusercontent.com/68455070/124414227-e1e16a00-dd84-11eb-943c-0cfc4592295a.png)

### 4. Add the CANONICAL splice site coordinates to the cryptic tag results
```r
STAR_crypt.res$canonical.chr <- STAR_canonical_results_control$canonical.chr[match(rownames(STAR_crypt.res),rownames(STAR_canonical_results_control))]
STAR_crypt.res$canonical.start <- STAR_canonical_results_control$canonical.start[match(rownames(STAR_crypt.res),rownames(STAR_canonical_results_control))]
STAR_crypt.res$canonical.end <- STAR_canonical_results_control$canonical.end[match(rownames(STAR_crypt.res),rownames(STAR_canonical_results_control))]
STAR_crypt.res$canonical.strand <- STAR_canonical_results_control$canonical.strand[match(rownames(STAR_crypt.res),rownames(STAR_canonical_results_control))]
STAR_crypt.res$canonical.intron.motif <- STAR_canonical_results_control$intron.motif[match(rownames(STAR_crypt.res),rownames(STAR_canonical_results_control))]
STAR_crypt.res$canonical.control.mean.SJ <- STAR_canonical_results_control$canonical.unique.count[match(rownames(STAR_crypt.res),rownames(STAR_canonical_results_control))]
STAR_crypt.res$canonical.control.mean.SJ <- STAR_crypt.res$canonical.control.mean.SJ / length(STAR_control_SJ_list)
```
`STAR_crypt.res`:

![image](https://user-images.githubusercontent.com/68455070/124410559-56180f80-dd7d-11eb-85f2-c5474aced706.png)

### 5. Discover (5') upstream cryptic junctions in controls
```r
STAR_upstream_results_control <- bridging_junction_finder(SJ.summary = STAR_control_total_SJ_counts,
                                                     results.df = STAR_crypt.res,
                                                     query.type = "upstream")
```
Function breakdown:
```r
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
	#convert into a dataframe, extracting the relevant information from the GRanges object.
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
```
`GRanges_object`:

![image](https://user-images.githubusercontent.com/68455070/124685141-40345700-df03-11eb-8a12-3c289be26dc5.png)

`upstream_junction_query` helper function: 
```r
upstream_junction_query <- function(CE.chr,CE.start,CE.end,canonical.start,SJ.GRange){
	junction <- SJ.GRange[seqnames(SJ.GRange) == CE.chr & start(SJ.GRange) == as.numeric(canonical.start) & end(SJ.GRange) >= as.numeric(CE.start) - 1 & end(SJ.GRange) < as.numeric(CE.end)]
	junction <- head(junction[order(score(junction),decreasing=T)],1)
    return(junction)
}
```
`junction <- SJ.GRange[seqnames(SJ.GRange) == CE.chr & start(SJ.GRange) == as.numeric(canonical.start) & end(SJ.GRange) >= as.numeric(CE.start) - 1 & end(SJ.GRange) < as.numeric(CE.end)]` (discovering 5' cryptic junctions; there may be more than one possible): 

![image](https://user-images.githubusercontent.com/68455070/124686007-e9c81800-df04-11eb-996f-e73bb3985428.png)

`junction <- head(junction[order(score(junction),decreasing=T)],1)` (in the case that there are more than 1 cryptic junctions found, get the most abundant supported cryptic junction. For this example, it should be the junction with a score of 57):

![image](https://user-images.githubusercontent.com/68455070/124686173-33186780-df05-11eb-9460-30f9887249d8.png)

Apply the above function to every row in the cryptic tag dataframe (i.e. for every cryptic tag): 

![image](https://user-images.githubusercontent.com/68455070/124686554-f0a35a80-df05-11eb-871c-74e5bf1ad262.png)

Convert result into a dataframe `STAR_upstream_results_control`:

![image](https://user-images.githubusercontent.com/68455070/124410585-64fec200-dd7d-11eb-8cb5-54b877b4c346.png)

### 6. Add the (5') upstream cryptic junctions coordinates to cryptic tag dataframe
```r
STAR_crypt.res$upstream.control.canonical.start <- STAR_upstream_results_control$canonical.start[match(rownames(STAR_crypt.res),rownames(STAR_upstream_results_control))]
STAR_crypt.res$upstream.control.cryptic.5prime <- STAR_upstream_results_control$cryptic.5prime[match(rownames(STAR_crypt.res),rownames(STAR_upstream_results_control))]
STAR_crypt.res$upstream.control.strand <- STAR_upstream_results_control$upstream.strand[match(rownames(STAR_crypt.res),rownames(STAR_upstream_results_control))]
STAR_crypt.res$upstream.control.intron.motif <- STAR_upstream_results_control$intron.motif[match(rownames(STAR_crypt.res),rownames(STAR_upstream_results_control))]
STAR_crypt.res$upstream.control.mean.SJ <- STAR_upstream_results_control$upstream.unique.count[match(rownames(STAR_crypt.res),rownames(STAR_upstream_results_control))]
STAR_crypt.res$upstream.control.mean.SJ <- STAR_crypt.res$upstream.control.mean.SJ / length(STAR_control_SJ_list)
```
`STAR_crypt.res`:

![image](https://user-images.githubusercontent.com/68455070/124410724-a8f1c700-dd7d-11eb-9dd0-5bd1d76a7d6e.png)

### 7. Discover (3') downstream cryptic junctions in controls
```r
STAR_downstream_results_control <- bridging_junction_finder(SJ.summary = STAR_control_total_SJ_counts, 
                                                       results.df = STAR_crypt.res, 
                                                       query.type = "downstream")
```
`STAR_downstream_results_control`: 

![image](https://user-images.githubusercontent.com/68455070/124410753-b73fe300-dd7d-11eb-9550-00c2009ec9cb.png)

### 8. Add the (3') downstream cryptic junctions coordinates to cryptic tag dataframe
```r
STAR_crypt.res$downstream.control.cryptic.3prime <- STAR_downstream_results_control$cryptic.3prime[match(rownames(STAR_crypt.res),rownames(STAR_downstream_results_control))]
STAR_crypt.res$downstream.control.canonical.end <- STAR_downstream_results_control$canonical.end[match(rownames(STAR_crypt.res),rownames(STAR_downstream_results_control))]
STAR_crypt.res$downstream.control.strand <- STAR_downstream_results_control$downstream.strand[match(rownames(STAR_crypt.res),rownames(STAR_downstream_results_control))]
STAR_crypt.res$downstream.control.intron.motif <- STAR_downstream_results_control$intron.motif[match(rownames(STAR_crypt.res),rownames(STAR_downstream_results_control))]
STAR_crypt.res$downstream.control.mean.SJ <- STAR_downstream_results_control$downstream.unique.count[match(rownames(STAR_crypt.res),rownames(STAR_downstream_results_control))]
STAR_crypt.res$downstream.control.mean.SJ <- STAR_crypt.res$downstream.control.mean.SJ / length(STAR_control_SJ_list)
```
`STAR_crypt.res`:

![image](https://user-images.githubusercontent.com/68455070/124410840-d9396580-dd7d-11eb-867a-7cd5d2e410d0.png)

Once Controls are done, we move on to cKO.

### 1. Merge SJ.out.tab files for cKO
```r
STAR_case_SJ_list <- c("cKO1.SJ.out.tab","cKO2.SJ.out.tab","cKO3.SJ.out.tab") 
STAR_case_total_SJ_counts <- merge.SJ.files(STAR_case_SJ_list)
```
`STAR_case_total_SJ_counts`:

![image](https://user-images.githubusercontent.com/68455070/124410859-e5252780-dd7d-11eb-8299-d0e5c0112a32.png)

STAR_crypt.res is now augmented with the coordinates discovered previously for the canonical sites. The subsequent function looks for canonical junctions in the case dataset that exactly match those discovered in the control dataset. Why do we need to do this? --> Because we want to find the SJ splicing from the canonical splice sites into the cryptic tag region.

### 2. Match canonical junctions between case and controls
```r
# We want to find junctions in the merged SJ summary file that matches exactly the canonical splice junctions discovered earlier in Controls analysis
STAR_canonical_results_case <- canonical_junction_detector(SJ.summary = STAR_case_total_SJ_counts, 
                                                      results.df = STAR_crypt.res, 
                                                      mode = "replication")
```

We see that in the `canonical_junction_detector` function, the option of `mode = "replication"`. In the actual `canonical_junction_detector` function: 
```r
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
```

The only difference between `discovery` mode and `replication` mode is that `discovery` mode is used for controls, when canonical splice junctions are being discovered. For `replication` mode, we simply want to match the canonical splice junctions already found in the controls, and recover them from the total splice junction files for cKO. The columns `x[10],x[15],x[16]` under the `replication` mode refer to `canonical_chr`, `canonical_start` and `canonical_end` respectively found earlier in controls. The code for `canonical_junction_replication` is similar to the `canonical_junction_query` function as explained above. 

`STAR_canonical_results_case`:

![image](https://user-images.githubusercontent.com/68455070/124410893-f2421680-dd7d-11eb-99ef-78c159894f01.png)

### 3. Repeat the same upstream and downstream cryptic splice site discovery for cKO (as performed on controls)

## Resultant splice junction rediscovery output file
`splicing_analysis.xlsx`:

![image](https://user-images.githubusercontent.com/68455070/124410951-1271d580-dd7e-11eb-912b-90010f8a326b.png)

![image](https://user-images.githubusercontent.com/68455070/124410983-24ec0f00-dd7e-11eb-8405-cce99e4711f1.png)

![image](https://user-images.githubusercontent.com/68455070/124411030-3e8d5680-dd7e-11eb-8432-2cf93726bda5.png)

![image](https://user-images.githubusercontent.com/68455070/124411815-d3448400-dd7f-11eb-9890-e51e8663234c.png)

![image](https://user-images.githubusercontent.com/68455070/124411848-e5262700-dd7f-11eb-8b8c-cf35424d8548.png)

![image](https://user-images.githubusercontent.com/68455070/124411870-f2431600-dd7f-11eb-9691-5677c08d8238.png)

![image](https://user-images.githubusercontent.com/68455070/124411897-ff600500-dd7f-11eb-8ee5-2b93ff7c1efa.png)

![image](https://user-images.githubusercontent.com/68455070/124411924-0b4bc700-dd80-11eb-881a-15fa48f35322.png)

