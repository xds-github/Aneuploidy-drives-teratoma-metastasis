# 只适用于小鼠非链特异性建库的转录组数据
# 直接在命令行运行下述程序，会自动生成任务投递脚本并执行
# 最终的结果以及数据质量报告都储存在result_foler里面
# 许德澍 2021/10/16 v2
################################ 直接命令行运行 ################################
your_folder=xudeshu/ #本人所在的根目录
input_folder=data/ # 原始数据目录
temp_folder=temp/ # 临时文件夹
result_foler=result/ #结果文件夹
mkdir -p ${result_foler}
mkdir -p ${temp_folder}
mkdir -p ${temp_folder}filter_dic/
mkdir -p ${temp_folder}report_dic/
filtered_data=${temp_folder}filter_dic/
report_folder=${temp_folder}report_dic/
STAR_index=/share/soft/mouse_star_index/
gtf_file=Mus_musculus.GRCm39.104.gtf
cd ${temp_folder}
ls ${input_folder} | grep '_R1.fq.gz$' > 1
ls ${input_folder} | grep '_R2.fq.gz$' > 2
ls ${input_folder} | grep '_R1.fq.gz$' | cut -d '_' -f 1 > 0
paste 0 1 2 > configure
cat configure | while read id
do
arr=(${id})
Read_1=${arr[1]}
Read_2=${arr[2]}
sample_Name=${arr[0]}
cp /share/soft/slurm_template.txt ${temp_folder}${sample_Name}.sh
echo "#SBATCH --job-name=rnaseq_${sample_Name}
#SBATCH --partition compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=6
docker run --rm --cpus=6 -v ${your_folder}:${your_folder} -v /share/soft:/share/soft  -v /data:/data f0651ca08acd /root/miniconda3/bin/fastp -i ${input_folder}${Read_1} -o ${filtered_data}${Read_1} -I ${input_folder}${Read_2} -O ${filtered_data}${Read_2} -j ${report_folder}${sample_Name}_fastp.json -w 6 -q 19 -u 50 -l 36 -x 
docker run --rm --cpus=6 -v ${your_folder}:${your_folder} -v /share/soft:/share/soft -v /data:/data f0651ca08acd /root/miniconda3/bin/STAR --runThreadN 6 --genomeDir ${STAR_index} --readFilesIn ${filtered_data}${Read_1} ${filtered_data}${Read_2} --readFilesCommand gunzip -c --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --outFileNamePrefix ${temp_folder}${sample_Name}_
cp ${temp_folder}${sample_Name}_ReadsPerGene.out.tab ${result_foler}${sample_Name}_ReadsPerGene.out.tab" >> ${temp_folder}${sample_Name}.sh
echo "sbatch ${temp_folder}${sample_Name}.sh" >> ${temp_folder}submit.sh
done
sh ${temp_folder}submit.sh
docker run --rm --cpus=1 -v ${your_folder}:${your_folder} -v /share/soft:/share/soft f0651ca08acd /root/miniconda3/bin/multiqc ${temp_folder}. -o ${result_foler}
