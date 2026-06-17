#!/bin/bash
基础配置
REF_DIR=/data2/chenh/bsa/REF/PP64_Cgla_REF
REF_FILE=${REF_DIR}/PP64.genome.fa
fq_data_DIR=/data2/chenh/bsa/F2_3WORKSPACE/fq_data/F2_3_hunchi
bam_DIR=/data2/chenh/bsa/F2_3WORKSPACE/bam/F2_3_hunchi/PP64_Cgla_REF_bam
gvcf_DIR=/data2/chenh/bsa/F2_3WORKSPACE/gvcf/F2_3_hunchi/PP64_REF_gvcf
scripts_DIR=/data2/chenh/bsa/F2_3WORKSPACE/scripts
log_DIR=/data2/chenh/bsa/F2_3WORKSPACE/log  
sample_list_FILE=/data2/chenh/bsa/F2_3WORKSPACE/sample_list/sample.list
# 使用 Samtools 建立 BWA 索引
# bwa index 默认使用 bwa index [OPTIONS] <fastafile>
# Samtools 可以直接使用 BWA 索引进行比对
# 对于 BWA 索引，不需要额外 Samtools 索引
# Samtools 用于后续的排序、变异检测等
# 建立索引
echo "Building BWA index..."
bwa index $REF_FILE
echo "Building Samtools index..."
samtools index $REF_FILE
gatk CreateSequenceDictionary -R $REF_FILE

#用法,注意不要遗忘sample文件夹
#bash call.test.sh /data2/chenh/bsa/F2_3WORKSPACE/sample_list/sample.list
#fq2gvcf主流程
for i in $(cat $1);
do
    cmd="docker run --gpus all --rm --volume /data2/chenh/bsa/F2_3WORKSPACE/fq_data/F2_3_hunchi/:/1 --volume /data2/chenh/bsa/REF/PP64_Cgla_REF/:/2 --volume /data2/chenh/bsa/F2_3WORKSPACE/bam/F2_3_hunchi/PP64_Cgla_REF_bam/:/3 nvcr.io/nvidia/clara/clara-parabricks:4.1.1-1 pbrun fq2bam  --ref /2/PP64.genome.fa --in-fq /1/${i}_1_clean.fq.gz /1/${i}_2_clean.fq.gz '@RG\tID:foo\tLB:lib1\tPL:bar\tSM:${i}\tPU:unit1' --out-bam /3/${i}.bam" 
    eval $cmd
    cmd="docker run --gpus all --rm --volume /data2/chenh/bsa/F2_3WORKSPACE/fq_data/F2_3_hunchi/:/1 --volume /data2/chenh/bsa/REF/PP64_Cgla_REF/:/2 --volume /data2/chenh/bsa/F2_3WORKSPACE/gvcf/F2_3_hunchi/PP64_REF_gvcf/:/3 nvcr.io/nvidia/clara/clara-parabricks:4.1.1-1 pbrun haplotypecaller --ref /2/PP64.genome.fa --in-bam /3/${i}.bam --gvcf --out-variant /3/${i}.g.vcf.gz"
    eval $cmd
done
#生成gvcf.list

cd $gvcf_DIR
for f in *.g.vcf.gz; do
    sample=$(basename "$f" .g.vcf.gz)
    echo -e "${sample}\t$(realpath "$f")"
  done > gvcf_sample.list

#注意不要忘记亲本的fq2gvcf，在最开始的sample.list就写好

