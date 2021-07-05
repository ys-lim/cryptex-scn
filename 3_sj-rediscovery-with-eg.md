### 1. Merge SJ.out.tab files
```r
STAR_control_SJ_list <- c("ctr1.SJ.out.tab","ctr2.SJ.out.tab","ctr3.SJ.out.tab") 
STAR_control_total_SJ_counts <- merge.SJ.files(STAR_control_SJ_list)
head(STAR_control_total_SJ_counts)
```
### 2. Load in dexseq results containing cryptic tag information
```r
STAR_dexseq.res <- "24may_scn_SignificantExons_STAR.csv"
STAR_crypt.res <- as.data.frame(fread(STAR_dexseq.res))
# obtaining cryptic exons with FDR < 0.05
STAR_crypt.res <- subset(STAR_crypt.res,grepl("i",exonID) & FDR < 0.05)
```

### 3. Finding canonical junctions flanking each cryptic tag
```r
STAR_canonical_results_control <- canonical_junction_detector(SJ.summary = STAR_control_total_SJ_counts, 
                                                         results.df = STAR_crypt.res, 
```                                                         mode = "discovery") 
