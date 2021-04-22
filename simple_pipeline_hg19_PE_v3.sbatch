#!/bin/bash -l

#SBATCH -A snic2020-15-9
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 8:00:00


bedpath="/proj/snic2020-6-3/bed_hg19/"
indexpath="/proj/snic2020-6-3/bowtie2-index"

filebase=$1

module load bioinfo-tools
module load bowtie2
module load samtools
module load IGVtools
module load FastQC
module load picard
module load deepTools
module load BEDTools/2.29.2

#filebase="pipeline/$2/$2"

#mkdir "pipeline"
#mkdir "pipeline/$2"

R1=$filebase"_R1.fastq.gz"
R2=$filebase"_R2.fastq.gz"

#if [ ! -f $filebase".hg19.all.bam" ]
#then
bowtie2 -p 8 -x $indexpath/hg19 -q -1 $R1 -2 $R2 --fast 2> $filebase".hg19.all.bowtie.log" | samtools view -bS -F 4  - | samtools sort -o $filebase".hg19.all.bam" -
samtools index $filebase".hg19.all.bam"
samtools flagstat $filebase".hg19.all.bam" > $filebase".hg19.all.flagstat.txt"
samtools idxstats $filebase".hg19.all.bam" > $filebase".hg19.all.idxstats.txt"
#fi

#if [ ! -f $filebase".hg19.fltd.bam" ]
#then
java -jar $PICARD_HOME/picard.jar MarkDuplicates REMOVE_DUPLICATES=TRUE I=$filebase".hg19.all.bam" O=$filebase".hg19.dedup.bam" M=$filebase".hg19.all.MarkDuplicates.txt"
samtools index $filebase".hg19.dedup.bam"
bedtools intersect -v -abam $filebase".hg19.dedup.bam" -b $bedpath"/blacklist.hg19.bed" > $filebase".hg19.fltd.bam"
samtools index $filebase".hg19.fltd.bam"

samtools flagstat $filebase".hg19.fltd.bam" > $filebase".hg19.fltd.flagstat.txt"
samtools idxstats $filebase".hg19.fltd.bam" > $filebase".hg19.fltd.idxstats.txt"
#fi

#if [ ! -f $filebase".hg19.bw" ]
#then
#	java -jar $PICARD_HOME/picard.jar CollectInsertSizeMetrics I=$filebase".hg19.fltd.bam" O=$filebase".hg19.insert_sizes.txt" H=$filebase".hg19.insert_sizes.pdf" M=0.5 &
#	igvtools count -e 60 $filebase".hg19.fltd.bam" $filebase".hg19.tdf" hg19.chrom.sizes &
bamCoverage --bam $filebase".hg19.fltd.bam" --binSize 5 --normalizeUsing RPGC --effectiveGenomeSize 2864785220 -o $filebase".hg19.bw"
bamCoverage -e --binSize 25 --smoothLength 75 --bam $filebase".hg19.fltd.bam" --normalizeUsing RPGC --effectiveGenomeSize 2864785220 -o $filebase".hg19.smooth.bw"
bamCoverage -e --minFragmentLength 30 --maxFragmentLength 100 --bam $filebase".hg19.fltd.bam" --normalizeUsing RPGC --effectiveGenomeSize 2864785220 -o $filebase".hg19.small.bw"
bamCoverage -e --minFragmentLength 150 --maxFragmentLength 300  --bam $filebase".hg19.fltd.bam" --normalizeUsing RPGC --effectiveGenomeSize 2864785220 -o $filebase".hg19.nucl.bw"
#fi


wait

module load ngsplot

rm -f "$filebase.config.txt"
for class in "RefGene_uniqueTSS.lt10kb.hg19.bed"
do
	echo -e $filebase".hg19.fltd.bam\t"$class"\t"$class >> "$filebase.config.txt"
done
ngs.plot.r -G mm9 -R bed -T "Genes $name" -C "$filebase.config.txt" -O $filebase".genes" -P 0 -SC global -L 4000 -IN 1 -MQ 20 -RB 0.01 -FL 150 -GO max
rm -rf $filebase".genes"
rm -f $filebase".genes.zip"
rm -f "$filebase.config.txt"
