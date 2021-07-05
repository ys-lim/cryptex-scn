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
`24may_scn_SignificantExons_STAR.csv`:
`STAR_crypt.res`:
### 3. Finding canonical junctions flanking each cryptic tag
```r
STAR_canonical_results_control <- canonical_junction_detector(SJ.summary = STAR_control_total_SJ_counts, 
                                                         results.df = STAR_crypt.res,
                                                         mode = "discovery")
```
`STAR_canonical_results_control`:
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
### 5. Counting (5') upstream cryptic junctions in controls
```r
STAR_upstream_results_control <- bridging_junction_finder(SJ.summary = STAR_control_total_SJ_counts,
                                                     results.df = STAR_crypt.res,
                                                     query.type = "upstream")
```
`STAR_upstream_results_control`:
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
### 7. Counting (3') downstream cryptic junctions in cKO
```r
STAR_downstream_results_control <- bridging_junction_finder(SJ.summary = STAR_control_total_SJ_counts, 
                                                       results.df = STAR_crypt.res, 
                                                       query.type = "downstream")
```
`STAR_downstream_results_control`: 
`STAR_crypt.res`:
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

Once Controls are done, we move on to cKO.

### 1. Merge SJ.out.tab files for cKO
```r
STAR_case_SJ_list <- c("cKO1.SJ.out.tab","cKO2.SJ.out.tab","cKO3.SJ.out.tab") 
STAR_case_total_SJ_counts <- merge.SJ.files(STAR_case_SJ_list)
```
`STAR_case_total_SJ_counts`:

STAR_crypt.res is now augmented with the coordinates discovered previously for the canonical sites. The subsequent function looks for canonical junctions in the case dataset that exactly match those discovered in the control dataset. Why do we need to do this? --> Because we want to find the SJ splicing from the canonical splice sites into the cryptic tag region.

### 2. Match canonical junctions between case and controls
```r
# We want to find junctions in the merged SJ summary file that matches exactly the canonical splice junctions discovered earlier in Controls analysis
STAR_canonical_results_case <- canonical_junction_detector(SJ.summary = STAR_case_total_SJ_counts, 
                                                      results.df = STAR_crypt.res, 
                                                      mode = "replication")
```
`STAR_canonical_results_case`:
### 3. Repeat the same upstream and downstream cryptic splice site discovery for cKO (as performed on controls)

## Resultant SJ rediscovery output file
`splicing_analysis.xlsx`:
