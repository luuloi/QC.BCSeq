#!/bin/bash

module load gi/bedtools/2.22.0
module load gi/samtools/1.2
module load gi/picard-tools/1.138
module load phuluu/R/3.1.2
module load phuluu/python/2.7.8
picar=/home/phuluu/bin/picard/dist/picard.jar
hg19=/home/phuluu/methods/darlo/annotations/hg19/      #hg19.fa #chrom.sizes
annotation=/share/Temp/phuluu/captureSeq/BisCaptureSeqQC/annotation/
mem="Xmx64g"

script=/home/phuluu/methods/001NGS/CaptureSeq/script/
path=/share/ClusterScratch/phuluu/captureSeq/15092015/ #15092015/ #18082015/
#path=/share/Temp/phuluu/captureSeq/18082015/
bam=${path}bams/
out=${path}tmp/
off=${out}bedtools/off_target/
on_target=${path}target/150416_HG19_LeanObese_SN_EPI_primary_targets.bed
off_target=${off}150416_HG19_LeanObese_SN_EPI_primary_targets.off.bed
on_target_capture=${path}target/150416_HG19_LeanObese_SN_EPI_capture_targets.bed
off_target_capture=${off}150416_HG19_LeanObese_SN_EPI_capture_targets.off.bed
DMR=${path}target/DMRs/
# ${DMR}primary_lean_obese.bed

# make directory
mkdir -p $out
mkdir -p ${out}bedtools/
mkdir -p ${out}bedtools/overlap_target/
mkdir -p ${out}bedtools/CpG_target/
mkdir -p $off

# make a off target bed
bedtools complement -i ${on_target}.sorted -g ${hg19}chrom.sizes > $off_target
bedtools complement -i ${on_target_capture}.sorted -g ${hg19}chrom.sizes > $off_target_capture

# Prepare all CpG coordinates
bigTable=/share/ClusterScratch/phuluu/captureSeq/15092015/bigTable/bigTable.tsv
awk 'NR>1{print $1"\t"$2"\t"$2}' $bigTable > ${annotation}CpGs.bed

# 0. Targets
# 0.1 Total number of targets
sort -k1,1 -k2,2n $on_target| uniq -u > ${on_target}.sorted
# 0.2 Total length of targets
awk '{print $3-$2}' ${on_target}.sorted > ${out}bedtools/overlap_target/length.primary

# 0.1 Total number of targets
sort -k1,1 -k2,2n $on_target_capture| uniq -u > ${on_target_capture}.sorted
# 0.2 Total length of targets
awk '{print $3-$2}' ${on_target_capture}.sorted > ${out}bedtools/overlap_target/length.capture

# 0.3 count the number of overlapping intervals in itself
oi=$(awk '{print $0"\t"1}' ${on_target}.sorted| bedtools merge -i - -c 5 -o sum| awk 'BEGIN{io=0}{if($4>1){io=io+1}}END{print io}')
oi_capture=$(awk '{print $0"\t"1}' ${on_target_capture}.sorted| bedtools merge -i - -c 5 -o sum| awk 'BEGIN{io=0}{if($4>1){io=io+1}}END{print io}')

# 0.4 Number of CpGs 
nCpG=$(bedtools intersect -a ${on_target}.sorted -b ${annotation}CpGs.bed -c| tee ${out}bedtools/CpG_target/count.tsv | awk 'BEGIN{sum=0}{sum+=$5}END{print sum}')
nCpG_capture=$(bedtools intersect -a ${on_target_capture}.sorted -b ${annotation}CpGs.bed -c| tee ${out}bedtools/CpG_target/count_capture.tsv | awk 'BEGIN{sum=0}{sum+=$5}END{print sum}')

# plot 
Rscript ${script}00TargetOverview/plot.nCpG.target.r ${out}bedtools/CpG_target/count.tsv

# plot metrics table and distribution of target widths
Rscript ${script}00TargetOverview/plot.metrics.table.width.target.r ${out}bedtools/overlap_target/ $nCpG $nCpG_capture $oi $oi_capture

# overlaping between primary and capture
bedtools intersect -a ${on_target}.sorted -b ${on_target_capture}.sorted -c| awk '{print $5}'| sort -k1,1n | uniq -c > ${out}bedtools/overlap_target/count.tsv
# plot 
Rscript ${script}00TargetOverview/plot.overlap.target.r ${out}bedtools/overlap_target/count.tsv

