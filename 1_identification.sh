# Annotation file: Mus_musculus.GRCm38.94.gtf 

############################
## PART 0: ANNOTATION PREP
############################
# Run dexseq prepare annotation.py on gtf annotation file
# C:\Users\yongs\R-3.6.0\library\DEXSeq\python_scripts\dexseq_prepare_annotation.py
# C:\Users\yongs\dexseq-test\alignments\Mus_musculus.GRCm38.94.gtf.fai
# Ran in atlas
python dexseq_prepare_annotation.py Mus_musculus.GRCm38.94.gtf Mus_musculus.GRCm38.94.dexseq.gtf

# Run dpryan R script to generate intronic intervals (intron GFF)
# http://seqanswers.com/forums/archive/index.php/t-42420.html
# Full code in dpryan.r file

# Get introns only
grep intronic Mus_musculus.GRCm38.94.dexseq.introns.gtf | awk 'BEGIN{OFS="\t"}{split($12,a,"\"");split($NF,b,"\"");print $1,$4,$5,$7,b[2]"_"a[2]}' > Mus_musculus.GRCm38.94.dexseq.introns.only.bed

grep exonic_part Mus_musculus.GRCm38.94.dexseq.gtf > Mus_musculus.GRCm38.94.dexseq.exons.only.gtf

############################
## PART 1: SPLICE EXTRACTION
############################
# samtools view -F 256 should keep out secondary giving primary aligned only
# awk gets the spliced reads "N" in cigar
## PURPOSE: Extract spliced reads from bam alignment files
samtools view -h -F 256 /mnt/gtklab01/ira/tdp43_cnp_scn/star/ctr1/ctr1.bam | awk '$1~/@/ || $6~/N/' | samtools view -bh > ctr/ctr1_spliced_reads.bam

samtools view -h -F 256 /mnt/gtklab01/ira/tdp43_cnp_scn/star/ctr2/ctr2.bam | awk '$1~/@/ || $6~/N/' | samtools view -bh > ctr/ctr2_spliced_reads.bam

samtools view -h -F 256 /mnt/gtklab01/ira/tdp43_cnp_scn/star/ctr3/ctr3.bam | awk '$1~/@/ || $6~/N/' | samtools view -bh > ctr/ctr3_spliced_reads.bam

samtools view -h -F 256 /mnt/gtklab01/ira/tdp43_cnp_scn/star/cKO1/cKO1.bam | awk '$1~/@/ || $6~/N/' | samtools view -bh > cko/cko1_spliced_reads.bam

samtools view -h -F 256 /mnt/gtklab01/ira/tdp43_cnp_scn/star/cKO2/cKO2.bam | awk '$1~/@/ || $6~/N/' | samtools view -bh > cko/cko2_spliced_reads.bam

samtools view -h -F 256 /mnt/gtklab01/ira/tdp43_cnp_scn/star/cKO3/cko3.bam | awk '$1~/@/ || $6~/N/' | samtools view -bh > cko/cko3_spliced_reads.bam

## PURPOSE: Intersect spliced reads with exon gtf file
bedtools intersect -a ctr1_spliced_reads.bam -b ../../Mus_musculus.GRCm38.94.dexseq.exons.only.gtf > ctr1_spliced_exons.bam

bedtools intersect -a ctr2_spliced_reads.bam -b ../../Mus_musculus.GRCm38.94.dexseq.exons.only.gtf > ctr2_spliced_exons.bam

bedtools intersect -a ctr3_spliced_reads.bam -b ../../Mus_musculus.GRCm38.94.dexseq.exons.only.gtf > ctr3_spliced_exons.bam

bedtools intersect -a cko1_spliced_reads.bam -b ../../Mus_musculus.GRCm38.94.dexseq.exons.only.gtf > cko1_spliced_exons.bam

bedtools intersect -a cko2_spliced_reads.bam -b ../../Mus_musculus.GRCm38.94.dexseq.exons.only.gtf > cko2_spliced_exons.bam

bedtools intersect -a cko3_spliced_reads.bam -b ../../Mus_musculus.GRCm38.94.dexseq.exons.only.gtf > cko3_spliced_exons.bam

## PURPOSE: Convert bam to bed file
bedtools bamtobed -i ctr1_spliced_exons.bam -split | sort -k1,1 -k2,2n > ctr1_spliced_exons.bed

bedtools bamtobed -i ctr2_spliced_exons.bam -split | sort -k1,1 -k2,2n > ctr2_spliced_exons.bed

bedtools bamtobed -i ctr3_spliced_exons.bam -split | sort -k1,1 -k2,2n > ctr3_spliced_exons.bed

bedtools bamtobed -i cko1_spliced_exons.bam -split | sort -k1,1 -k2,2n > cko1_spliced_exons.bed

bedtools bamtobed -i cko2_spliced_exons.bam -split | sort -k1,1 -k2,2n > cko2_spliced_exons.bed

bedtools bamtobed -i cko3_spliced_exons.bam -split | sort -k1,1 -k2,2n > cko3_spliced_exons.bed

## PURPOSE: Inverse intersect spliced exons bed with same exon gff file
bedtools intersect -a ctr1_spliced_exons.bed -b ../../Mus_musculus.GRCm38.94.dexseq.exons.only.gtf -v > ctr1_spliced_introns.bed

bedtools intersect -a ctr2_spliced_exons.bed -b ../../Mus_musculus.GRCm38.94.dexseq.exons.only.gtf -v > ctr2_spliced_introns.bed

bedtools intersect -a ctr3_spliced_exons.bed -b ../../Mus_musculus.GRCm38.94.dexseq.exons.only.gtf -v > ctr3_spliced_introns.bed

bedtools intersect -a cko1_spliced_exons.bed -b ../../Mus_musculus.GRCm38.94.dexseq.exons.only.gtf -v > cko1_spliced_introns.bed

bedtools intersect -a cko2_spliced_exons.bed -b ../../Mus_musculus.GRCm38.94.dexseq.exons.only.gtf -v > cko2_spliced_introns.bed

bedtools intersect -a cko3_spliced_exons.bed -b ../../Mus_musculus.GRCm38.94.dexseq.exons.only.gtf -v > cko3_spliced_introns.bed

#################################
## PART 2: MERGING & GFF CREATION
#################################

## PURPOSE: Group all intronic mapping reads bed together (ko and control)
cat cko1_spliced_introns.bed cko2_spliced_introns.bed cko3_spliced_introns.bed ../ctr/ctr1_spliced_introns.bed ../ctr/ctr2_spliced_introns.bed ../ctr/ctr3_spliced_introns.bed | sort -S 50% -k1,1 -k2,2n > grouped_spliced_introns_sorted.bed

## PURPOSE: Merge spliced introns within 500bp of each other
bedtools merge -i grouped_spliced_introns_sorted.bed -d 500 -c 1 -o count > spliced_introns_merged.bed

## PURPOSE: Intersect merged spliced introns file with introns BED file 
bedtools intersect -a spliced_introns_merged.bed -b ../Mus_musculus.GRCm38.94.dexseq.introns.only.bed -wb | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$8,$9}' | sort -k1,1V -k5,5n > cryptics_merged.bed

## PURPOSE: Filter unique cryptic tags
cat cryptics_merged.bed | awk '{print $6}' | sort -V | uniq > unique_gene_introns.tab

## PURPOSE: create .sh script to create unique list of cryptic exons annotated with unique id
N=8
for entry in `cat unique_gene_introns.tab`; do
((i=i%N)); ((i++==0)) && wait # this line prevents overloading for no. of running forks
grep $entry cryptics_merged.bed | awk 'BEGIN{s=1}{print $0"i"s;s+=1}' >> cryptics_merged_annotated.bed
done

## PURPOSE: Convert cryptic exon bed file into gff file for merging
cat cryptics_merged_annotated.bed | awk 'BEGIN{OFS="\t"}{split($6,a,"_");print $1, "Mus_musculus.GRCm38.94.gtf", "exonic_part", $2, $3, ".", $5, ".", "transcripts \"cryptic_exon\"; exonic_part_number \""a[2]"\"; gene_id \""a[1]"\"" }' | sort -k1,1 -k2,2n > cryptic_exons.gff

## PURPOSE: Concatenate cryptic exons within total exon gff & sort - all_exons_and_cryptics.gff forms the annotation for subsequent cryptic tag hunting
cat cryptic_exons.gff Mus_musculus.GRCm38.94.dexseq.gtf | sort -k1,1V -k4,4n -k5,5n | awk '$14 ~ /ENS/' > all_exons_and_cryptics.gff

#########################
## PART 3: READ COUNTING
#########################

samtools sort -n /mnt/gtklab01/ira/tdp43_cnp_scn/star/ctr1/ctr1.bam -o ctr1_sorted.bam

samtools sort -n /mnt/gtklab01/ira/tdp43_cnp_scn/star/ctr1/ctr2.bam -o ctr2_sorted.bam 

samtools sort -n /mnt/gtklab01/ira/tdp43_cnp_scn/star/ctr1/ctr3.bam -o ctr3_sorted.bam

samtools sort -n /mnt/gtklab01/ira/tdp43_cnp_scn/star/cKO3/cko1.bam -o cko1_sorted.bam

samtools sort /mnt/gtklab01/ira/tdp43_cnp_scn/star/cKO3/cko2.bam -o cko2_sorted.bam

samtools sort /mnt/gtklab01/ira/tdp43_cnp_scn/star/cKO3/cko3.bam -o cko3_sorted.bam


# read counting for scn samples (star) using actual all_exons_and_cryptics.gff 
python dexseq_count.py -s reverse -p yes -r name -f bam all_exons_and_cryptics.gff ctr/ctr1_sorted.bam ctr/ctr1_pycount.txt &> ctr/ctr1_pycount_report.txt &

python dexseq_count.py -s reverse -p yes -r name -f bam all_exons_and_cryptics.gff ctr/ctr2_sorted.bam ctr/ctr2_pycount.txt &> ctr/ctr2_pycount_report.txt &

python dexseq_count.py -s reverse -p yes -r name -f bam all_exons_and_cryptics.gff ctr/ctr3_sorted.bam ctr/ctr3_pycount.txt &> ctr/ctr3_pycount_report.txt &

python dexseq_count.py -s reverse -p yes -r name -f bam all_exons_and_cryptics.gff cko/cko1_sorted.bam cko/cko1_pycount.txt &> cko/cko1_pycount_report.txt &

python dexseq_count.py -s reverse -p yes -r name -f bam all_exons_and_cryptics.gff cko/cko2_sorted.bam cko/cko2_pycount.txt &> cko/cko2_pycount_report.txt &

python dexseq_count.py -s reverse -p yes -r name -f bam all_exons_and_cryptics.gff cko/cko3_sorted.bam cko/cko3_pycount.txt &> cko/cko3_pycount_report.txt &
