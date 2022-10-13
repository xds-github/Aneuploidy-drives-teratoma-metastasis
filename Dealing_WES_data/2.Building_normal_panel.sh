Referencefile=./mouse_ensembl_ref.fa
Vcf_file=./mus_musculus.vcf
interval=./mm10_exon.interval_list
GATK3=./GenomeAnalysisTK.jar
bam_file=./ready_bam/
panel_dic=./build_normal_panel/
fn=$(ls ${bam_file} | grep "bam$" | grep "N.*sorted" | cut -d "." -f 1)
for i in ${fn}
do
echo "#This pipline is writen for building normal panel
#xudeshu
#2020/11/17
#v12
#software: GATK3.8 /bwa 0.7.17-r1188 / samtools 1.7 /Picard
#environment: conda python 3.6
source /home/xudeshu/.bashrc
conda activate py3.6
java -Xmx10g -Djava.io.tmpdir=/temp_data -jar ${GATK3} \
 -T MuTect2 \
 -I:tumor ${bam_file}${i}.sorted.markdup.realign.BQSR.bam \
 --dbsnp ${Vcf_file} \
 --artifact_detection_mode \
 -o ${panel_dic}${i}.vcf.gz \
 -L ${interval} \
 -R ${Referencefile}" > ${panel_dic}${i}.sh
echo "qsub -cwd -l vf=10g,h=node2 ${panel_dic}${i}.sh" >> ${panel_dic}submid.sh
done
sh ${panel_dic}submid.sh
source /home/xudeshu/.bashrc
conda activate py3.6
fn=$(ls ${panel_dic} | grep "vcf.gz$" | cut -d "." -f 1)
# Create input.list
for i in ${fn}
do
echo "-V ${panel_dic}${i}.vcf.gz \ " >> ${panel_dic}inputs.list
done
java -Xmx15g -Djava.io.tmpdir=./temp_data -jar ${GATK3} \
-T CombineVariants \
--arg_file ${panel_dic}inputs.list \
-minN 2 \
--setKey "null" \
--filteredrecordsmergetype KEEP_IF_ANY_UNFILTERED \
-o ${panel_dic}2_pon_combinevariants.vcf.gz \
-R ${Referencefile}
# Finally, generate a sites-only VCF with Picard's MakeSitesOnlyVcf
java -Xmx15g -jar ${Picard} MakeSitesOnlyVcf \
I=${panel_dic}2_pon_combinevariants.vcf.gz \
O=${panel_dic}3_pon_siteonly.vcf.gz
