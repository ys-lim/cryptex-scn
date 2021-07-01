## Identification of cryptic tags from RNA-seq BAM files (using Nfasc gene as example)

### 1. Prepare flattened exon and intron annotation file
```bash
python dexseq_prepare_annotation.py Mus_musculus.GRCm38.94.gtf Mus_musculus.GRCm38.94.dexseq.gtf
```
`Mus_musculus.GRCm38.94.gtf`:

![image](https://user-images.githubusercontent.com/68455070/123911249-20a2a900-d9ae-11eb-920f-62bea91b9be8.png)


`Mus_musculus.GRCm38.94.dexseq.gtf`:

![image](https://user-images.githubusercontent.com/68455070/123911122-f3ee9180-d9ad-11eb-9bf4-a5635b0e532f.png)

### 2. Extract spliced reads (reads with N in CIGAR)
```bash
samtools view -h -F 256 ../nfascTrue/ctrNfascReads.bam | awk '$1~/@/ || $6~/N/' | samtools view -bh > ctrNfascReads_spliced.bam
```
### 3. Intersect spliced reads with flattened exon annotation gtf to extract spliced exonic reads
```bash
bedtools intersect -a ctrNfascReads_spliced.bam -b Mus_musculus.GRCm38.94.dexseq.exons.only.gtf > ctrNfasc_spliced_exons.bam
```

`ctrNfasc_spliced_exons.bam`:

![image](https://user-images.githubusercontent.com/68455070/123912024-1c2ac000-d9af-11eb-8f0f-863f3e1ce091.png)

### 4. Convert bam to bed file
```bash
bedtools bamtobed -i ctrNfasc_spliced_exons.bam -split | sort -k1,1 -k2,2n > ctrNfasc_spliced_exons.bed
```
`-split` allows RNA-seq reads to be split based on N in CIGAR, that is, one read spanning an intronic region will be split into two distinct bed intervals.

`ctrNfasc_spliced_exons.bed`:

![image](https://user-images.githubusercontent.com/68455070/123912145-3f556f80-d9af-11eb-89a2-ab3e9ce69dfd.png)

## 5. Inverse intersect spliced reads with exon gtf to extract spliced intronic reads
```bash
bedtools intersect -a ctrNfasc_spliced_exons.bed -b Mus_musculus.GRCm38.94.dexseq.exons.only.gtf -v > ctrNfasc_spliced_introns.bed
```
This removes the exonic mapping region of a spliced read, leaving only the intronic mapping region. 

`ctrNfasc_spliced_introns.bed`:

![image](https://user-images.githubusercontent.com/68455070/123913513-e25ab900-d9b0-11eb-94e4-abb608a704ca.png)

## 6. Concatenate all spliced intronic reads and sort them
```bash
cat ctrNfasc_spliced_introns.bed ckoNfasc_spliced_introns.bed | sort -S 50% -k1,1 -k2,2n > grouped_spliced_introns_sorted.bed
```

`grouped_spliced_introns_sorted.bed`:


![image](https://user-images.githubusercontent.com/68455070/123913700-19c96580-d9b1-11eb-8158-b5bbc7d04c3b.png)

## 7. Merge intronic reads within 500 bp of each other to produce a pre-cryptic tag
```bash
bedtools merge -i grouped_spliced_introns_sorted.bed -d 500 -c 1 -o count > spliced_introns_merged.bed
```

`spliced_introns_merged.bed`:

![image](https://user-images.githubusercontent.com/68455070/123913868-4ed5b800-d9b1-11eb-842f-df577c0485d8.png)

## 8. Intersect pre-cryptic tag with intron gtf

```bash
bedtools intersect -a spliced_introns_merged.bed -b ../Mus_musculus.GRCm38.94.dexseq.introns.only.bed -wb | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$8,$9}' | sort -k1,1V -k5,5n > cryptics_merged.bed
```

This prevents cryptic tags that span across an exon of less than 500bp to be wrongly merged together. The intron gtf contains intervals of introns between exonic regions.


`Mus_musculus.GRCm38.94.dexseq.introns.only.bed`:


![image](https://user-images.githubusercontent.com/68455070/123914012-7fb5ed00-d9b1-11eb-8bec-4b1259afaf3b.png)


`bedtools intersect -a spliced_introns_merged.bed -b ../Mus_musculus.GRCm38.94.dexseq.introns.only.bed -wb`:


![image](https://user-images.githubusercontent.com/68455070/124051609-606d9d00-da4f-11eb-8328-c33841f6f867.png)


`cryptics_merged.bed`:

![image](https://user-images.githubusercontent.com/68455070/124051495-2e5c3b00-da4f-11eb-9190-7e628ff2fd4c.png)

### 9. Grab the cryptic tag ID and filter for unique ones

```bash
cat cryptics_merged.bed | awk '{print $6}' | sort -V | uniq > unique_gene_introns.tab
```

`unique_gene_introns.tab`:

![image](https://user-images.githubusercontent.com/68455070/124051714-9874e000-da4f-11eb-9592-b5efca253f2a.png)

### 10. Create a unique list of cryptic exons annotated with user-created unique ID
```bash
N=8
for entry in `cat unique_gene_introns.tab`; do
((i=i%N)); ((i++==0)) && wait # this line prevents overloading for no. of running forks
grep $entry cryptics_merged.bed | awk 'BEGIN{s=1}{print $0"i"s;s+=1}' >> cryptics_merged_annotated.bed
done
```

`cryptics_merged_annotated.bed`:

![image](https://user-images.githubusercontent.com/68455070/124051815-c8bc7e80-da4f-11eb-977f-97e1453804e4.png)

### 11. Convert cryptic tag bed file into cryptic tag gff for merging with annotated gtf

```bash
cat cryptics_merged_annotated.bed | awk 'BEGIN{OFS="\t"}{split($6,a,"_");print $1, "Mus_musculus.GRCm38.94.gtf", "exonic_part", $2, $3, ".", $5, ".", "transcripts \"cryptic_exon\"; exonic_part_number \""a[2]"\"; gene_id \""a[1]"\"" }' | sort -k1,1 -k2,2n > cryptic_exons.gff
```

`cryptic_exons.gff`:

![image](https://user-images.githubusercontent.com/68455070/124051920-fb667700-da4f-11eb-8a86-b2a77f13d531.png)

### 12. Concatenate cryptic exons within annotated exon gff & sort - all_exons_and_cryptics.gff forms the annotation for subsequent cryptic tag hunting

```bash
cat cryptic_exons.gff Mus_musculus.GRCm38.94.dexseq.gtf | sort -k1,1V -k4,4n -k5,5n | awk '$14 ~ /ENS/' > all_exons_and_cryptics.gff
```

`all_exons_and_cryptics.gff`:

![image](https://user-images.githubusercontent.com/68455070/124052087-50a28880-da50-11eb-8d88-fd67641b4fc3.png)

### 13. Read counting for samples (STAR) using actual all_exons_and_cryptics.gff and dexseq

```bash
python dexseq_count.py -s reverse -p yes -r name -f bam all_exons_and_cryptics.gff ctr/ctr1_sorted.bam ctr/ctr1_pycount.txt &> ctr/ctr1_pycount_report.txt
```

`ctr1_pycount.txt`:

![image](https://user-images.githubusercontent.com/68455070/124052655-711f1280-da51-11eb-841b-c0f15a5cfa0a.png)
