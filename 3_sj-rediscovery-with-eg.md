### 1. Merge SJ.out.tab files
```r
STAR_control_SJ_list <- c("ctr1.SJ.out.tab","ctr2.SJ.out.tab","ctr3.SJ.out.tab") 
STAR_control_total_SJ_counts <- merge.SJ.files(STAR_control_SJ_list)
head(STAR_control_total_SJ_counts)
```
