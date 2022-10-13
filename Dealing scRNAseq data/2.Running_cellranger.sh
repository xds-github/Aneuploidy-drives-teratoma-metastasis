fastq_dic=./scRNA_data/
result_dic=./scRNA_output/
reference_dic=./Mouse_GFP_Lu_reference/
sh_folder=./sh_file/
cat sc_configure.txt | while read id
do
arr=(${id})
fasq_id=${arr[0]}
sample_id=${arr[1]}
echo "source /home/xudeshu/.bashrc
cd ${result_dic}
./cellranger-3.1.0/cellranger count --id=${sample_id} --fastqs=${fastq_dic} --sample=${fasq_id} --transcriptome=${reference_dic} --localmem=128 --localcores=32" > ${sh_folder}${sample_id}.sh
echo "qsub -cwd -l vf=128g  ${sh_folder}${sample_id}.sh" >> ${sh_folder}submid.sh
done
sh ${sh_folder}submid.sh
