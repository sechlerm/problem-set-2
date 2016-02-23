#! usr/bin/env bash

## Question 1

#Use BEDtools intersect to identify the size of the largest overlap between
#CTCF and H3K4me3 locations.
#print last column in awk- so you don't have to count columns {print $NF}
data='/Users/marybethsechler/Desktop/GAW/data-sets'

overlap=$(gzcat $data/bed/peaks.chr22.bed.gz\
 | awk '$4 == "CTCF" '\
 | bedtools intersect -a - -b $data/bed/encode.h3k4me3.hela.chr22.bed.gz -wo\
 | sort -k15n | tail -n1 | awk '{print $NF}')

echo 'answer-1 ' $overlap


## Question 2

#Use BEDtools to calculate the GC content of nucleotides 19,000,000 to
#19,000,500 on chr22 of `hg19` genome build. Report the GC content
#as a fraction (e.g., 0.50).


echo -e "chr22\t19000000\t19000500" > interval.bed
gc_content=$(bedtools nuc -fi $data/fasta/hg19.chr22.fa -bed interval.bed | grep -v "#" | cut -f5)
echo 'answer-2 ' $gc_content 



## Question 3

#Use BEDtools to identify the length of the CTCF ChIP-seq peak (i.e.,
#interval) that has the largest mean signal in `ctcf.hela.chr22.bg.gz`.
gzcat $data/bed/peaks.chr22.bed.gz\
 | awk '$4 == "CTCF"' > CTCFpeaks.bed

sort -i CTCFpeaks.bed > temp.bed
mv temp.bed CTCFpeaks.bed

largesig=$(bedtools map -a CTCFpeaks.bed -b $data/bedtools/ctcf.hela.chr22.bg.gz -c 4 -o mean | awk '$5 != "."' | sort -k5n | tail -n1 |  awk '{print $3 - $2} ' )

echo 'answer-3 ' $largesig


## Question 4

#Use BEDtools to identify the gene promoter (defined as 1000 bp upstream of
#a TSS) with the highest median signal in `ctcf.hela.chr22.bg.gz`.  Report
#the gene name (e.g., 'ABC123')

gzcat $data/bed/genes.hg19.bed.gz | awk 'BEGIN {OFS="\t"} ($6 == "+") {print $1, $2, $2 + 1, $4, $5, $6 }' > tss.bed

gzcat $data/bed/genes.hg19.bed.gz | awk 'BEGIN {OFS="\t"} ($6 == "-") {print $1, $3, $3 + 1, $4, $5, $6 }' >> tss.bed

bedtools sort -i tss.bed > tmp.bed
mv tmp.bed tss.bed


promoter=$(bedtools flank -i tss.bed -g $data/genome/hg19.genome -s -b 1000\
 |bedtools sort -i -\
 |  bedtools  map -a - -b $data/bedtools/ctcf.hela.chr22.bg.gz -c 4 -o median | awk '$7 != "."' | sort -k7n | tail -n1 | cut -f4)

echo "answer-4 " $promoter 


## Question 5

#Use BEDtools to identify the longest interval on `chr22` that is not
#covered by `genes.hg19.bed.gz`. Report the interval like `chr1:100-500`.

#gzcat $data/bed/genes.hg19.bed.gz | awk '$1 == "chr22"' > chr22q6.bed
#sort -i chr22q6.bed > tmp.bed
#mv tmp.bed chr22q6.bed

#bedtools complement -i chr22q6.bed -g $data/genome/hg19.genome\
# | awk ' $1 == "chr22"'\
# | awk '{print $1, $2, $3, $3 - $2}'\
# | sort - k4n
 
notcov=$(bedtools complement -i $data/bed/genes.hg19.bed.gz -g $data/genome/hg19.genome\
 | awk 'BEGIN {OFS = "\t"} {print $1, $2, $3, $3-$2}'\
 | sort -k4n\
 | tail -n1\
 | awk ' { print $1":"$2"-"$3}') 

echo "answer-5 " $notcov

## Question 6 (extra credit)
# I will use Jaccard to show the statistic representing the overlap of two sets, what would be found using intersect as in question 1

overlapstat=$(gzcat $data/bed/peaks.chr22.bed.gz\
 | awk '$4 == "CTCF" '\
 | bedtools jaccard -a - -b $data/bed/encode.h3k4me3.hela.chr22.bed.gz | grep -v "intersection" | cut -f3)

echo "answer-6  the overlap statistic from question 1 using Jaccard : " $overlapstat 




#Use one or more BEDtools that we haven't covered in class. Be creative.



