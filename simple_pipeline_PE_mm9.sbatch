#!/bin/bash -l

#SBATCH -A snic2020-15-9
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 10:00:00
#SBATCH -M snowy


bedpath="/proj/snic2020-6-3/bed_mm9/"
indexpath="/proj/snic2020-6-3/bowtie2-index"


filebase=$1


module load bioinfo-tools
module load bowtie2
module load samtools
module load IGVtools
module load FastQC
module load picard
module load deepTools
module load BEDTools



R1=$filebase"_R1.fastq.gz"
R2=$filebase"_R2.fastq.gz"

#unique to Cut&Tag: remove 19bp mosaic adaptor sequence
bowtie2 -p 8 -x $indexpath/mm9 -q -1 $R1 -2 $R2 --fast-local 2> $filebase".mm9.all.bowtie.log" | samtools view -bS -F 4  - | samtools sort -o $filebase".mm9.all.bam" -
samtools index $filebase".mm9.all.bam"
samtools flagstat $filebase".mm9.all.bam" > $filebase".mm9.all.flagstat.txt"
samtools idxstats $filebase".mm9.all.bam" > $filebase".mm9.all.idxstats.txt"
java -jar $PICARD_HOME/picard.jar MarkDuplicates REMOVE_DUPLICATES=TRUE I=$filebase".mm9.all.bam" O=$filebase".mm9.dedup.bam" M=$filebase".mm9.all.MarkDuplicates.txt"

samtools index $filebase".mm9.dedup.bam"
bedtools intersect -v -abam $filebase".mm9.dedup.bam" -b $bedpath"/Blacklist.bed" > $filebase".mm9.fltd.bam"

rm -f $filebase".mm9.dedup.bam"
rm -f $filebase".mm9.dedup.bam.bai"

samtools index $filebase".mm9.fltd.bam"
samtools flagstat $filebase".mm9.fltd.bam" > $filebase".mm9.fltd.flagstat.txt"
samtools idxstats $filebase".mm9.fltd.bam" > $filebase".mm9.fltd.idxstats.txt"

samtools view -bh -q 20 $filebase".mm9.fltd.bam" > $filebase".mm9.fltd.uniq.bam"

samtools index $filebase".mm9.fltd.uniq.bam"
samtools flagstat $filebase".mm9.fltd.uniq.bam" > $filebase".mm9.fltd.uniq.flagstat.txt"
samtools idxstats $filebase".mm9.fltd.uniq.bam" > $filebase".mm9.fltd.uniq.idxstats.txt"

samtools view -c -L $bedpath"/refGene.bed" $filebase".mm9.fltd.bam" > $filebase".mm9.genes.bedstats.txt"
samtools view -c -L $bedpath"/GSE17051_ciLAD_regions.bed" $filebase".mm9.fltd.bam" > $filebase".mm9.ciLAD.bedstats.txt"
samtools view -c -L $bedpath"/GSE17051_cLAD_regions.bed" $filebase".mm9.fltd.bam" > $filebase".mm9.cLAD.bedstats.txt"
samtools view -c -L $bedpath"/ESC_Enhancers.bed" $filebase".mm9.fltd.bam" > $filebase".mm9.enhancers.bedstats.txt"
samtools view -c -L $bedpath"/rmsk_IAP.bed" $filebase".mm9.fltd.bam" > $filebase".mm9.IAP.bedstats.txt"
samtools view -c -L $bedpath"/IAPEz.bed" $filebase".mm9.fltd.bam" > $filebase".mm9.IAPEz.bedstats.txt"

java -jar $PICARD_HOME/picard.jar CollectInsertSizeMetrics I=$filebase".mm9.fltd.bam" O=$filebase".mm9.insert_sizes.txt" H=$filebase".mm9.insert_sizes.pdf" M=0.5 &

igvtools count --pairs $filebase".mm9.fltd.bam" $filebase".mm9.fltd.tdf" mm9.chrom.sizes &
wait

bamCoverage -p 2 --normalizeUsing RPGC --effectiveGenomeSize 2150570000 -bs 5 -e --minMappingQuality 20 -b $filebase".mm9.fltd.bam" -o $filebase".mm9.pe.uniq.bw" &
bamCoverage -p 2 --normalizeUsing RPGC --effectiveGenomeSize 2150570000 -bs 5 -e --minMappingQuality 20 --minFragmentLength 50 --maxFragmentLength 300 -b $filebase".mm9.fltd.bam" -o $filebase".mm9.pe.50-300.uniq.bw" &
bamCoverage -p 2 --normalizeUsing RPGC --effectiveGenomeSize 2150570000 -bs 5 -b $filebase".mm9.fltd.bam" -o $filebase".mm9.bw" &

plotFingerprint -b $filebase".mm9.fltd.bam" --minMappingQuality 30 --skipZeros --region 19 --numberOfSamples 50000 --plotFile $filebase".mm9.fingerprint.png"

bowtie2 -p 8 -x $indexpath/rodrep -q -1 $R1 -2 $R2 --fast-local 2> $filebase".rodrep.bowtie.log" | samtools view -bS -F 4 - | samtools sort -o $filebase".rep.bam" -

samtools index $filebase".rep.bam"
samtools flagstat $filebase".rep.bam" > $filebase".rep.flagstat.txt"
samtools idxstats $filebase".rep.bam" > $filebase".rep.idxstats.txt"


#----- write out statistics


total=$(( $(zcat $R1 | wc -l) / 2 ))
mapped=$(grep '+ 0 mapped' $filebase".mm9.all.flagstat.txt" | sed -r 's/ \+ 0 mapped.+//g' -)
rep=$(grep '+ 0 mapped' $filebase".rep.flagstat.txt" | sed -r 's/ \+ 0 mapped.+//g' -)
telomere=$(grep 'Rep_Telomere' $filebase".rep.idxstats.txt" | sed -r 's/Rep_Telomere\t[0-9]+\t//g' - | sed -r 's/\t0//g' -)
gsat=$(grep 'Rep_GSAT_MM_SAT_Mus_musculus' $filebase".rep.idxstats.txt" | sed -r 's/Rep_GSAT_MM_SAT_Mus_musculus\t[0-9]+\t//g' - | sed -r 's/\t0//g' -)
iapi=$(grep 'Rep_IAPLTR1a_I_MM_ERV2_Mus_musculus' $filebase".rep.idxstats.txt" | sed -r 's/Rep_IAPLTR1a_I_MM_ERV2_Mus_musculus\t[0-9]+\t//g' - | sed -r 's/\t0//g' -)

echo -e $filebase"\t"$total"\t"$mapped"\t"$rep"\t"$gsat"\t"$telomere"\t"$iapi >> umifilterstats.txt

mapped=$(grep '+ 0 mapped' $filebase".mm9.all.flagstat.txt" | sed -r 's/ \+ 0 mapped.+//g' -)
fltd=$(grep '+ 0 mapped' $filebase".mm9.fltd.flagstat.txt" | sed -r 's/ \+ 0 mapped.+//g' -)
rep=$(grep '+ 0 mapped' $filebase".rep.flagstat.txt" | sed -r 's/ \+ 0 mapped.+//g' -)
telomere=$(grep 'Rep_Telomere' $filebase".rep.idxstats.txt" | sed -r 's/Rep_Telomere\t[0-9]+\t//g' - | sed -r 's/\t0//g' -)
gsat=$(grep 'Rep_GSAT_MM_SAT_Mus_musculus' $filebase".rep.idxstats.txt" | sed -r 's/Rep_GSAT_MM_SAT_Mus_musculus\t[0-9]+\t//g' - | sed -r 's/\t0//g' -)
iapi=$(grep 'Rep_IAPLTR1a_I_MM_ERV2_Mus_musculus' $filebase".rep.idxstats.txt" | sed -r 's/Rep_IAPLTR1a_I_MM_ERV2_Mus_musculus\t[0-9]+\t//g' - | sed -r 's/\t0//g' -)
genes=$(cat $filebase".mm9.genes.bedstats.txt")
ciLAD=$(cat $filebase".mm9.ciLAD.bedstats.txt")
cLAD=$(cat $filebase".mm9.cLAD.bedstats.txt")
enh=$(cat $filebase".mm9.enhancers.bedstats.txt")
IAP=$(cat $filebase".mm9.IAP.bedstats.txt")
IAPEz=$(cat $filebase".mm9.IAPEz.bedstats.txt")

echo -e $filebase"\t"$mapped"\t"$fltd"\t"$genes"\t"$enh"\t"$IAP"\t"$IAPEz"\t"$ciLAD"\t"$cLAD"\t"$rep"\t"$gsat"\t"$telomere"\t"$iapi >> bamstats.txt

module load MACS

mkdir macs2
macs2 callpeak -t $filebase".mm9.fltd.uniq.bam" -n $filebase -g mm

wait

module load ngsplot
rm -f "$1.config.txt"
for class in "genes_hi_lt10kb" "genes_mid_lt10kb" "genes_lo_lt10kb" "K27_bivalent_genes.mm9"
do
	echo -e $1".mm9.fltd.bam\t"$bedpath$class".bed\t"$class >> "$1.config.txt"
done
ngs.plot.r -G mm9 -R bed -T "Bivalent Genes $name" -C "$1.config.txt" -O $1".bivalent_genes" -P 0 -SC global -L 4000 -IN 1 -MQ 20 -RB 0.01 -FL 150 -GO max
rm -rf $1".bivalent_genes"
rm -f $1".bivalent_genes.zip"

rm -f "$1.config.txt"

rm -f $name.mm9.bam.cnt
rm -f $name.mm9.fltd.bam.cnt
