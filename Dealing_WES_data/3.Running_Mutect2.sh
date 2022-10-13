Referencefile=./mouse_ensembl_ref.fa
Vcf_file=./mus_musculus.vcf
ready_bam=./ready_bam/
GATK3=/home/xudeshu/gatk3.8/GenomeAnalysisTK.jar
interval=./mm10_exon.interval_list
panel_dic=./build_normal_panel/
result_dic=./mutect_vcf/
temp_data=./temp_data/
sh_folder=./vcf_sh/
cat wes_configure2.txt | while read id
do
arr=(${id})
tumor_name=${arr[0]}
normal_name=${arr[1]}
even_name=${arr[2]}
node_name=${arr[3]}
echo "source /home/xudeshu/.bashrc
conda activate py3.6
java -Xmx10g -Djava.io.tmpdir=${temp_data} \
-jar ${GATK3} -T MuTect2 \
-I:tumor ${ready_bam}${tumor_name}.sorted.markdup.realign.BQSR.bam \
-I:normal ${ready_bam}${normal_name}.sorted.markdup.realign.BQSR.bam \
--dbsnp ${Vcf_file} \
--normal_panel ${panel_dic}3_pon_siteonly.vcf.gz \
--output_mode EMIT_VARIANTS_ONLY \
-o ${result_dic}${even_name}.vcf.gz \
-L ${interval} \
-R ${Referencefile}" > ${sh_folder}${even_name}.sh
echo "qsub -cwd -l vf=10g,h=${node_name} ${sh_folder}${even_name}.sh" >> ${sh_folder}submit.sh
done
sh ${sh_folder}submit.sh
