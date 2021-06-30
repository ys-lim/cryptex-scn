# Identification of cryptic tags from RNA-seq BAM files (subset to Nfasc gene)

1. Extract spliced reads
```bash
samtools view -h -F 256 ../nfascTrue/ctrNfascReads.bam | awk '$1~/@/ || $6~/N/' | samtools view -bh > ctrNfascReads_spliced.bam
```
