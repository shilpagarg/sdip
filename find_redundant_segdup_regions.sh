#Convert tsv file to bed file
awk 'OFS="\t" {print $2, $3, $4, $1}' segdups.similar.tsv > segdups.similar.bed

#Compute a table of good overlaps between different segdup regions
#Segdup region in column 8 is contained in segdup region in column 4
bedtools intersect -a segdups.similar.bed  -b segdups.similar.bed -wa -wb -F 0.9 | awk '$4!=$8' | sort -k 8,8n | less
