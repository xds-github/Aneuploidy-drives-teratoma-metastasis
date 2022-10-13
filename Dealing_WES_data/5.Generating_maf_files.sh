# output vcf files from mutect2 will annotated by VEP and the result will be converted to maf files
# Xudeshu 2020/9/13
# v6
source /home/xudeshu/.bashrc
#conda activate WGS
input_dic=./filtered_vcf/
output_dic=./vep_annotation/
sh_folder=./vep_sh/
fn=$(ls ${input_dic} | grep 'filtered.vcf$' | cut -d '_' -f 1,2)
for i in ${fn}
do
echo "source /home/xudeshu/.bashrc
conda activate WGS
perl ./vcf2maf.pl \
--input-vcf ${input_dic}${i}_filtered.vcf \
--output-maf ${output_dic}${i}.vep.maf \
--tumor-id ${i} \
--normal-id temp \
--ncbi-build GRCm38 \
--cache-version 101 \
--vep-data ./mm10_vep/ \
--vep-path /home/xudeshu/anaconda3/pkgs/ensembl-vep-101.0-pl526hecda079_0/bin \
--species mus_musculus \
--filter-vcf ./mus_musculus.vcf.gz \
--ref-fasta ./Mus_musculus.GRCm38.dna.toplevel.fa.gz" > ${sh_folder}${i}.sh
echo "sh ${sh_folder}${i}.sh" >> ${sh_folder}submit.sh
done
sh ${sh_folder}submit.sh
cat ${output_dic}WES_Ts11-1-M2P.vep.maf | head -n 2 > ${output_dic}total.maf
fn=$( ls ${output_dic} | grep 'maf$')
for i in ${fn}
do
    cat ${output_dic}${i} | grep -v '#' | grep -v 'Hugo_Symbol' >> ${output_dic}total.maf
done
