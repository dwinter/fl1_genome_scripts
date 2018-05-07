#assign genes to fragments for HiC comparisons
#note some genes will be on two fragments
# frags.txt comes from HiCDatR digest fo the genome (which is GUI only)
awk -v OFS='\t' '{print $2,$3,$4,$1'} frags.txt  | tail -n +2 > digest.bed
bedtools intersect -wb -a ../../annotation/M3_nuc.bed -b digest.bed  > gene_by_fragment.tsv
