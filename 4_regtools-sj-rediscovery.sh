### TEST
# 1. Merge all bam files for each treatment (ctrl and cko)
# 2. For each cryptic tag, identify the exon before CE & the exon after CE; i.e the exons flanking the cryptic tag
# - Obtain genomic region between before.exon.end and after.exon.start (bedtools with all_cryptics_and_exons?)
# 3. Use regtools to discover all the SJs within this genomic region (i.e. the region that the CE sits) from the Nfasc bam files (splice junction rediscovery)
# 4. Filter SJs:
# - SJ counts < 10
# - Categorise: upstream SJ = those btween before.exon.end and middle of CE, downstream SJ = those btween middle of CE and after.exon.start
# (do i need to match the start and end splice junctions?)
# 5. Calculate PSI for each SJ somehow

###############
#### NFASC TEST
###############

# Nfasc bam alignments in mm10 (GRCm38): NC_000067.6: chr1: 132546974 - 132759564 https://www.ncbi.nlm.nih.gov/genome/gdv/browser/gene/?id=269116
samtools view /mnt/gtklab01/ira/tdp43_cnp_scn/star/ctr.bam 1:132546974-132759564 -o ctr-nfasc-star.bam

samtools view /mnt/gtklab01/ira/tdp43_cnp_scn/star/cKO.bam 1:132546974-132759564 -o cko-nfasc-star.bam

# Go to 28may_sciatic_nerve_cryptic_exon_report --> filter for nfasc cryptic exons
# there are 4 cryptic tags E025, E008, E024, E009, E008
Nfasc_E025
Nfasc_E008
Nfasc_E024
Nfasc_E009
Nfasc_E008

# from all_exons_and_cryptics: all exons (cryptic included) in gff file
awk '{if($1=="1" && $4>132546975 && $5<132759565) print $0}' /mnt/gtklab01/yongshan/sciatic_nerve/star/all_exons_and_cryptics.gff > test.gff

# regtools to extract junctions within the exonic region that flanks the cryptic tag
# splice junction rediscovery step
# block: junction anchor overhang on either ends of the junction
# determind strand settings: https://rnabio.org/module-09-appendix/0009/12/01/StrandSettings/
# tool to check strandedness: https://github.com/betsig/how_are_we_stranded_here
# A gene can live on a DNA strand in one of two orientations. The gene is said to have a coding strand (also known as its sense strand), and a template strand (also known as its antisense strand). For 50% of genes, its coding strand will correspond to the chromosome's forward strand, and for the other 50% it will correspond to the reverse strand.
# Annotations such as Ensembl and UCSC are concerned with the coding sequences of genes, so when they say a gene is on the forward strand, it means the gene's coding sequence is on the forward strand. To follow through again, that means that during transcription of this forward-strand gene, the gene's template sequence is read from the reverse strand, producing an mRNA that matches the sequence on the forward strand.

## Extracting splice junctions btween the exon before and after cryptic tag
./regtools junctions extract -r 1:132607047-132608383 -s 1 -o ~/cryptex-junctions-test/cko-nfasc-junctions.bed ~/cryptex-junctions-test/cko-nfasc-star.bam

./regtools junctions extract -r 1:132607047-132608383 -s 1 -o ~/cryptex-junctions-test/ctr-nfasc-junctions.bed ~/cryptex-junctions-test/ctr-nfasc-star.bam
# Check this during next meeting: having the same repeated junctions = reads mapping to both strands! Qns: Can I ignore positive strand arbitrarily just because there are more reads on negative strand? Is there a reason why positive strand produces much less reads supporting the same junction compared to negative strand? Does this relate to being first stranded?
# ANS: choose the strand that Nfasc gene is on --> in this case, the negative strand


# PSI calculation:
# "calculated as the ratio of inclusion reads (reads spliced to exon of interest) to the sum of inclusion and exclusion reads (reads spliced to other exons)"

## WRITE A SCRIPT TO AUTOMATE THIS

samtools view -H ctr-nfasc-star.bam |\
   sed -e 's/SN:1/SN:chr1/' | sed -e 's/SN:2/SN:chr2/' | \
   sed -e 's/SN:3/SN:chr3/' | sed -e 's/SN:4/SN:chr4/' | \
   sed -e 's/SN:5/SN:chr5/' | sed -e 's/SN:6/SN:chr6/' | \
   sed -e 's/SN:7/SN:chr7/' | sed -e 's/SN:8/SN:chr8/' | \
   sed -e 's/SN:9/SN:chr9/' | sed -e 's/SN:10/SN:chr10/' | \
   sed -e 's/SN:11/SN:chr11/' | sed -e 's/SN:12/SN:chr12/' | \
   sed -e 's/SN:13/SN:chr13/' | sed -e 's/SN:14/SN:chr14/' | \
   sed -e 's/SN:15/SN:chr15/' | sed -e 's/SN:16/SN:chr16/' | \
   sed -e 's/SN:17/SN:chr17/' | sed -e 's/SN:18/SN:chr18/' | \
   sed -e 's/SN:19/SN:chr19/' | sed -e 's/SN:20/SN:chr20/' | \
   sed -e 's/SN:21/SN:chr21/' | sed -e 's/SN:22/SN:chr22/' | \
   sed -e 's/SN:X/SN:chrX/' | sed -e 's/SN:Y/SN:chrY/' | \
   sed -e 's/SN:MT/SN:chrM/' | samtools reheader - ctr-nfasc-star.bam > ctr-nfasc-star_chr.bam
   
   
samtools view -H cko-nfasc-star.bam |\
   sed -e 's/SN:1/SN:chr1/' | sed -e 's/SN:2/SN:chr2/' | \
   sed -e 's/SN:3/SN:chr3/' | sed -e 's/SN:4/SN:chr4/' | \
   sed -e 's/SN:5/SN:chr5/' | sed -e 's/SN:6/SN:chr6/' | \
   sed -e 's/SN:7/SN:chr7/' | sed -e 's/SN:8/SN:chr8/' | \
   sed -e 's/SN:9/SN:chr9/' | sed -e 's/SN:10/SN:chr10/' | \
   sed -e 's/SN:11/SN:chr11/' | sed -e 's/SN:12/SN:chr12/' | \
   sed -e 's/SN:13/SN:chr13/' | sed -e 's/SN:14/SN:chr14/' | \
   sed -e 's/SN:15/SN:chr15/' | sed -e 's/SN:16/SN:chr16/' | \
   sed -e 's/SN:17/SN:chr17/' | sed -e 's/SN:18/SN:chr18/' | \
   sed -e 's/SN:19/SN:chr19/' | sed -e 's/SN:20/SN:chr20/' | \
   sed -e 's/SN:21/SN:chr21/' | sed -e 's/SN:22/SN:chr22/' | \
   sed -e 's/SN:X/SN:chrX/' | sed -e 's/SN:Y/SN:chrY/' | \
   sed -e 's/SN:MT/SN:chrM/' | samtools reheader - cko-nfasc-star.bam > cko-nfasc-star_chr.bam
