Referencefile=./BWA_index/mouse_ensembl_ref.fa
Vcf_file=./mus_musculus.vcf
oridic=./WES_out/filtered_data/
temdic=./WES_out/temp_data/
report_folder=./WES_out/data_QC/
rawdic=./WES_data/
finaldic=./ready_bam/
interval=./mm10_exon.interval_list
sh_folder=./process_sh/
Picard=/home/xudeshu/picard/picard.jar
GATK3=/home/xudeshu/gatk3.8/GenomeAnalysisTK.jar
cat wes_configure_1.txt | while read id
do
arr=(${id})
sample_Name=${arr[0]}
Lane=${arr[1]}
node=${arr[2]}
Read_1=$(ls ${rawdic} | grep ${sample_Name} | grep 'R1')
Read_2=$(ls ${rawdic} | grep ${sample_Name} | grep 'R2')
Read_3=${sample_Name}_R1.fq.gz
Read_4=${sample_Name}_R2.fq.gz
echo "#This pipline is writen for transforing fastq data into ready_sort_BAM
#xudeshu
#2020/11/7
#v7
#software: GATK3.8 /bwa 0.7.17-r1188 / samtools 1.7 /Picard
#environment: conda python 3.6
source /home/xudeshu/.bashrc
/home/xudeshu/fastp/fastp \
-i ${rawdic}${Read_1} \
-o ${oridic}${Read_3} \
-I ${rawdic}${Read_2} \
-O ${oridic}${Read_4} \
-j ${report_folder}${sample_Name}_fastp.json \
-w 4 \
-q 19 \
-u 50 \
-l 36 \
-x
conda activate py3.6
##################### Map to Reference###################################
bwa mem -t 6 -R \"@RG\tID:${sample_Name}\tPL:Illumina\tLB:${Lane}\tPU:${sample_Name}\tSM:${sample_Name}\" \
 ${Referencefile} \
 ${oridic}${Read_3} \
 ${oridic}${Read_4} \
 | samtools view -S-b > ${temdic}${sample_Name}.bam
##################### Mark Duplication###################################
samtools sort -@ 4 \
-m 4G \
-O bam \
-o ${temdic}${sample_Name}.sorted.bam \
${temdic}${sample_Name}.bam
samtools index ${temdic}${sample_Name}.sorted.bam
java -Xmx20g -Djava.io.tmpdir=${temdic} \
-jar ${Picard} MarkDuplicates \
  I=${temdic}${sample_Name}.sorted.bam \
  O=${temdic}${sample_Name}.sorted.markdup.bam \
  M=${temdic}${sample_Name}.markdup_metrics.txt
samtools index ${temdic}${sample_Name}.sorted.markdup.bam
#####################       BQSR      ###################################
java -Xmx20g -Djava.io.tmpdir=${temdic} -jar ${GATK3} \
    -T BaseRecalibrator  \
    -R ${Referencefile} \
    -I ${temdic}${sample_Name}.sorted.markdup.bam \
    -L ${interval} \
    -knownSites ${Vcf_file} \
    -o ${finaldic}${sample_Name}.recal_data1.table
java -Xmx20g -Djava.io.tmpdir=${temdic} -jar ${GATK3} \
    -T BaseRecalibrator  \
    -R ${Referencefile} \
    -I ${temdic}${sample_Name}.sorted.markdup.bam \
    -L ${interval}\
    -knownSites ${Vcf_file} \
    -o ${finaldic}${sample_Name}.recal_data2.table \
    -BQSR ${finaldic}${sample_Name}.recal_data1.table
java -Xmx20g -Djava.io.tmpdir=${temdic} -jar ${GATK3} \
    -T PrintReads \
    -R ${Referencefile} \
    -I ${temdic}${sample_Name}.sorted.markdup.bam \
    -BQSR ${finaldic}${sample_Name}.recal_data1.table \
    -o ${finaldic}${sample_Name}.sorted.markdup.realign.BQSR.bam
" > ${sh_folder}${sample_Name}.sh
echo "qsub -cwd -l vf=20g,h=${node} ${sh_folder}${sample_Name}.sh" >> ${sh_folder}submit.sh
done
sh ${sh_folder}submit.sh
