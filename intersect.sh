
#return all intersections
bedtools intersect -a $1 -b $2 -wao > $1.anno.temp

#for each original interval, retain only one intersection with the biggest overlap in bp
cat $1.anno.temp | awk '{ if ($8 > max[$1$2$3]) i[$1$2$3]=$1"\t"$2"\t"$3"\t"$7; max[$1$2$3]=$8} END {for (x in i) {print i[x]}}' > $1.anno.bed

#count
cut -f4 $1.anno.bed | sort | uniq -c | sort -nr > $1.anno.counts.txt
cut -f4 $1.anno.bed | sed -E 's/[0-9]+_//' | sort | uniq -c | sort -nr > $1.anno.simple.txt

echo "results:"
cat $1.anno.counts.txt

echo "simplified:"
cat $1.anno.simple.txt
