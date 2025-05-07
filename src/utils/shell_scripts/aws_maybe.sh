## if I want to run some long-running jobs on AWS, I was getting started on this script to do so.
## it still needs a fair amount of work.

for i in {1..4} x y m; do
    grep -w $i data/gtfs/human_genome/gencode.v19.annotation.gtf > data/gtfs/human_genome/gencode.v19.annotation.chr$i.gtf
done

bash src/utils/shell_scripts/index_all_chr.sh

for i in {1..4} x y m; do
data/human_genome/
    samtools faidx data/human_genome/chr$i/chr$i.fna
done

for i in {1..4} x y m; do
    bash src/utils/shell_scripts/bwa_align_and_sort_by_chr.sh chr$i SRR413984
done

for i in {1..4} x y m; do
samtools index data/bams/chr$i/SRR413984_chr$i.sorted.bam
done

for i in {1..4} x y m; do
bash src/utils/shell_scripts/querysort_bam_by_chr.sh chr$i SRR413984
done

for i in {1..4} x y m; do
Circle-Map ReadExtractor -i data/bams/chr1/SRR413984_chr1.querysorted.bam -o data/bams/chr1/SRR413984_chr1.candidates.bam
 done

for i in {1..4} x y m; do 
bash src/utils/shell_scripts/sort_and_index_candidates_by_chr.sh chr$i SRR413984
 done

for i in {1..4} x y m; do
Circle-Map Realign \
  -i data/bams/chr$i/SRR413984_chr$i.candidates.sorted.bam \
  -qbam data/bams/chr$i/SRR413984_chr$i.querysorted.bam \
  -sbam data/bams/chr$i/SRR413984_chr$i.sorted.bam \
  -fasta data/human_genome/chr$i/chr$i.fna \
  -o data/beds/chr$i/SRR413984_chr$i.eccdna.bed \
  -t 12
done


for i in {1..4} x y m; do
if [ "$i" = "1" ]; then
    sed "s/^NC_000001.10/chr1/" data/beds/chr1/SRR413984_chr1.eccdna.cleaned.bed > data/beds/chr1/SRR413984_chr1.eccdna.cleaned.ucsc.bed
elif [ "$i" = "2" ]; then
    sed "s/^NC_000002.11/chr2/" data/beds/chr2/SRR413984_chr2.eccdna.cleaned.bed > data/beds/chr2/SRR413984_chr2.eccdna.cleaned.ucsc.bed
elif [ "$i" = "3" ]; then
    sed "s/^NC_000003.11/chr3/" data/beds/chr3/SRR413984_chr3.eccdna.cleaned.bed > data/beds/chr3/SRR413984_chr3.eccdna.cleaned.ucsc.bed
elif [ "$i" = "4" ]; then
    sed "s/^NC_000004.11/chr4/" data/beds/chr4/SRR413984_chr4.eccdna.cleaned.bed > data/beds/chr4/SRR413984_chr4.eccdna.cleaned.ucsc.bed
elif [ "$i" = "x" ]; then
    sed "s/^NC_000023.10/chrx/" data/beds/chrx/SRR413984_chrx.eccdna.cleaned.bed > data/beds/chrx/SRR413984_chrx.eccdna.cleaned.ucsc.bed
elif [ "$i" = "y" ]; then
    sed "s/^NC_000024.9/chry/" data/beds/chry/SRR413984_chry.eccdna.cleaned.bed > data/beds/chry/SRR413984_chry.eccdna.cleaned.ucsc.bed
elif [ "$i" = "m" ]; then
    sed "s/^NC_012920.1/chrm/" data/beds/chrm/SRR413984_chrm.eccdna.cleaned.bed > data/beds/chrm/SRR413984_chrm.eccdna.cleaned.ucsc.bed
fi
done



