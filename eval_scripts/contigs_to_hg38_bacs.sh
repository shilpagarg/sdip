#Number and sequence length of contigs that align end to end to hg38 with identity >= 99.8%
cat /project/heller_d-data/shilpa/NA12878/regions/eval/t5.b5.d2/tables/polished.diploid.to.hg38.tbl | awk '{if ($2 == 0 && $3 == $4 && $11 >= 99.8) {print $0} }' | wc -l
cat /project/heller_d-data/shilpa/NA12878/regions/eval/t5.b5.d2/tables/polished.diploid.to.hg38.tbl | awk '{if ($2 == 0 && $3 == $4 && $11 >= 99.8) {SUM += $4} } END {print SUM}'

#Number and sequence length of contigs that align end to end to hg38 with identity >= 90%
cat /project/heller_d-data/shilpa/NA12878/regions/eval/t5.b5.d2/tables/polished.diploid.to.hg38.tbl | awk '{if ($2 == 0 && $3 == $4 && $11 >= 90) {print $0} }' | wc -l
cat /project/heller_d-data/shilpa/NA12878/regions/eval/t5.b5.d2/tables/polished.diploid.to.hg38.tbl | awk '{if ($2 == 0 && $3 == $4 && $11 >= 90) {SUM += $4} } END {print SUM}'

#Get list of unaligned contigs
join -v 2 -1 1 -2 1 <(cat /project/heller_d-data/shilpa/NA12878/regions/eval/t5.b5.d2/tables/polished.diploid.to.hg38.tbl | awk '{if ($2 == 0 && $3 == $4 && $11 >= 99.8) {print $0} }' | sort -k 1,1) <(sort -k 1,1 /project/heller_d-data/shilpa/NA12878/regions/contigs/pooled.t5.b5.d2.polished.diploid.fa.fai) | cut -f 1 -d' ' > unaligned.hg38.998.txt
join -v 2 -1 1 -2 1 <(cat /project/heller_d-data/shilpa/NA12878/regions/eval/t5.b5.d2/tables/polished.diploid.to.hg38.tbl | awk '{if ($2 == 0 && $3 == $4 && $11 >= 90) {print $0} }' | sort -k 1,1) <(sort -k 1,1 /project/heller_d-data/shilpa/NA12878/regions/contigs/pooled.t5.b5.d2.polished.diploid.fa.fai) | cut -f 1 -d' ' > unaligned.hg38.90.txt

#Number of unaligned contigs that align end to end to BACs with identity >= 99.8%
join -1 1 -2 1 -t$'\t' <(sort -k 1,1 /project/heller_d-data/shilpa/NA12878/regions/eval/t5.b5.d2/tables/polished.diploid.to.bacs.tbl) unaligned.hg38.998.txt | awk '{if ($2 == 0 && $3 == $4 && $11 >= 99.8) {print $0} }' | wc -l
#Number of unaligned contigs that align end to end to BACs with identity >= 90%
join -1 1 -2 1 -t$'\t' <(sort -k 1,1 /project/heller_d-data/shilpa/NA12878/regions/eval/t5.b5.d2/tables/polished.diploid.to.bacs.tbl) unaligned.hg38.90.txt | awk '{if ($2 == 0 && $3 == $4 && $11 >= 90) {print $0} }' | wc -l