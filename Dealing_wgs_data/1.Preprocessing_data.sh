source /home/xudeshu/.bashrc
raw_dic=./WGS_data/
oridic=./WGS_output/filtered_data/
temdic=./WGS_output/temp_data/
finaldic=./WGS_output/ready_bam/
Referencefile=./mouse_ensembl_ref.fa
Vcf_file=./vcf_NCBI_dpsnp/00-All.vcf
gatk_temp=./WGS_output/gatk_temp/
report_folder=./WGS_output/data_QC/
sh_folder=./WGS_output/sh_folder/
cat configure.txt | while read id
do
arr=(${id})
sample_Name=${arr[0]}
Lane=${arr[1]}
node=${arr[2]}
Read_1=${sample_Name}_R1.fq.gz
Read_2=${sample_Name}_R2.fq.gz
echo "source /home/xudeshu/.bashrc
/home/xudeshu/fastp/fastp \
-i ${raw_dic}${Read_1} \
-o ${oridic}${Read_1} \
-I ${raw_dic}${Read_2} \
-O ${oridic}${Read_2} \
-j ${report_folder}${sample_Name}_fastp.json \
-w 4 \
-q 19 \
-u 50 \
-l 36 \
-x
conda activate py3.6
##################### Map to Reference###################################
bwa mem -t 6 -R \"@RG\tID:${Lane}\tPL:Illumina\tLB:${sample_Name}\tSM:${sample_Name}\" \
 ${Referencefile} \
 ${oridic}/${Read_1} \
 ${oridic}/${Read_2} \
 | samtools view -S-b > ${temdic}${sample_Name}.bam
##################### Mark Duplication###################################
gatk --java-options \"-Xmx12G\" MarkDuplicatesSpark \
 -I ${temdic}${sample_Name}.bam \
 -O ${temdic}${sample_Name}__Marked_duplicate.bam \
 -M ${temdic}${sample_Name}__Marked_dup_metrics.txt \
 --conf 'spark.executor.cores=5' \
 --tmp-dir ${gatk_temp}
#####################       BQSR      ###################################
gatk --java-options \"-Xmx12G\" BaseRecalibrator \
 -I ${temdic}${sample_Name}__Marked_duplicate.bam \
 -R ${Referencefile} \
 --known-sites ${Vcf_file} \
 -O ${finaldic}${sample_Name}_recal_data.table \
 --tmp-dir ${gatk_temp}
gatk --java-options \"-Xmx12G\" ApplyBQSR \
 -R ${Referencefile} \
 -I ${temdic}${sample_Name}__Marked_duplicate.bam \
 --bqsr-recal-file ${finaldic}${sample_Name}_recal_data.table \
 -O ${finaldic}${sample_Name}_ready.bam \
 --tmp-dir ${gatk_temp}" > ${sh_folder}${sample_Name}.sh
 echo "qsub -cwd -l vf=15g,h=${node} ${sh_folder}${sample_Name}.sh" >> ${sh_folder}submid.sh
done
sh ${sh_folder}submid.sh
