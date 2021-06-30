## Identification of cryptic tags from RNA-seq BAM files (subset to Nfasc gene)

### 1. Prepare flattened exon and intron annotation file
```bash
python dexseq_prepare_annotation.py Mus_musculus.GRCm38.94.gtf Mus_musculus.GRCm38.94.dexseq.gtf
```
Mus_musculus.GRCm38.94.gtf:
![image](https://user-images.githubusercontent.com/68455070/123911249-20a2a900-d9ae-11eb-920f-62bea91b9be8.png)

Mus_musculus.GRCm38.94.dexseq.gtf:
![image](https://user-images.githubusercontent.com/68455070/123911122-f3ee9180-d9ad-11eb-9bf4-a5635b0e532f.png)
### 1. Extract spliced reads
```bash
samtools view -h -F 256 ../nfascTrue/ctrNfascReads.bam | awk '$1~/@/ || $6~/N/' | samtools view -bh > ctrNfascReads_spliced.bam
```
