```r
> gtf_test <- import.gff2("gencode.vM27.primary_assembly.annotation.dexseq.chr3.gtf")
> head(gtf_test)
GRanges object with 6 ranges and 7 metadata columns:
      seqnames          ranges strand |                       source           type     score     phase              gene_id
         <Rle>       <IRanges>  <Rle> |                     <factor>       <factor> <numeric> <integer>          <character>
  [1]     chr3 3069070-3069169      - | dexseq_prepare_annotation.py aggregate_gene      <NA>      <NA> ENSMUSG00002075012.1
  [2]     chr3 3069070-3069169      - | dexseq_prepare_annotation.py    exonic_part      <NA>      <NA> ENSMUSG00002075012.1
  [3]     chr3 3092718-3092817      - | dexseq_prepare_annotation.py aggregate_gene      <NA>      <NA> ENSMUSG00002075999.1
  [4]     chr3 3092718-3092817      - | dexseq_prepare_annotation.py    exonic_part      <NA>      <NA> ENSMUSG00002075999.1
  [5]     chr3 3253093-3253219      + | dexseq_prepare_annotation.py aggregate_gene      <NA>      <NA> ENSMUSG00000103416.2
  [6]     chr3 3253093-3253219      + | dexseq_prepare_annotation.py    exonic_part      <NA>      <NA> ENSMUSG00000103416.2
               transcripts exonic_part_number
               <character>        <character>
  [1]                 <NA>               <NA>
  [2] ENSMUST00020182553.1                001
  [3]                 <NA>               <NA>
  [4] ENSMUST00020181724.1                001
  [5]                 <NA>               <NA>
  [6] ENSMUST00000191986.2                001
  -------
  seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
