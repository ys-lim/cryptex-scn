## Identification of cryptic tags from RNA-seq BAM files (subset to Nfasc gene)

### 1. Prepare flattened exon and intron annotation file
```bash
python dexseq_prepare_annotation.py Mus_musculus.GRCm38.94.gtf Mus_musculus.GRCm38.94.dexseq.gtf
```
<p align="center">
    Mus_musculus.GRCm38.94.gtf:
</p>

![image](https://user-images.githubusercontent.com/68455070/123911249-20a2a900-d9ae-11eb-920f-62bea91b9be8.png)

<p align="center">
    Mus_musculus.GRCm38.94.dexseq.gtf:
</p>

![image](https://user-images.githubusercontent.com/68455070/123911122-f3ee9180-d9ad-11eb-9bf4-a5635b0e532f.png)

<p align="center">
    Mus_musculus.GRCm38.94.dexseq.gtf:
</p>
### 1. Extract spliced reads (reads with N in CIGAR)
```bash
samtools view -h -F 256 ../nfascTrue/ctrNfascReads.bam | awk '$1~/@/ || $6~/N/' | samtools view -bh > ctrNfascReads_spliced.bam
```
### 2. Intersect spliced reads with exon gtf to extract spliced exonic reads
```bash
bedtools intersect -a ctrNfascReads_spliced.bam -b Mus_musculus.GRCm38.94.dexseq.exons.only.gtf > ctrNfasc_spliced_exons.bam
```
<p align="center">
    ctrNfasc_spliced_exons.bam:
</p>

![image](https://user-images.githubusercontent.com/68455070/123912024-1c2ac000-d9af-11eb-8f0f-863f3e1ce091.png)

### 3. Convert bam to bed file
```bash
bedtools bamtobed -i ctrNfasc_spliced_exons.bam -split | sort -k1,1 -k2,2n > ctrNfasc_spliced_exons.bed
```
<p align="center">
    ctrNfasc_spliced_exons.bed:
</p>

![image](https://user-images.githubusercontent.com/68455070/123912145-3f556f80-d9af-11eb-89a2-ab3e9ce69dfd.png)

## 4. Inverse intersect spliced reads with exon gtf to extract spliced intronic reads
```bash
bedtools intersect -a ctrNfasc_spliced_exons.bed -b Mus_musculus.GRCm38.94.dexseq.exons.only.gtf -v > ctrNfasc_spliced_introns.bed
```
<p align="center">
    ctrNfasc_spliced_introns.bed:
</p>

![image](https://user-images.githubusercontent.com/68455070/123913513-e25ab900-d9b0-11eb-94e4-abb608a704ca.png)

## 5. Concatenate all spliced intronic reads and sort them
```bash
cat ctrNfasc_spliced_introns.bed ckoNfasc_spliced_introns.bed | sort -S 50% -k1,1 -k2,2n > grouped_spliced_introns_sorted.bed
```

<p align="center">
    grouped_spliced_introns_sorted.bed:
</p>

![image](https://user-images.githubusercontent.com/68455070/123913700-19c96580-d9b1-11eb-8158-b5bbc7d04c3b.png)

## 6. Merge intronic reads within 500 bp of each other to produce a pre-cryptic tag
```bash
bedtools merge -i grouped_spliced_introns_sorted.bed -d 500 -c 1 -o count > spliced_introns_merged.bed
```

![image](https://user-images.githubusercontent.com/68455070/123913868-4ed5b800-d9b1-11eb-842f-df577c0485d8.png)

## 7. Intersect pre-cryptic tag with intron gtf

<p align="center">
    Mus_musculus.GRCm38.94.dexseq.introns.only.bed:
</p>

![image](https://user-images.githubusercontent.com/68455070/123914012-7fb5ed00-d9b1-11eb-8bec-4b1259afaf3b.png)
