source /home/xudeshu/.bashrc
refer_fa=./mouse_ensembl_ref.fa
wig_gz=./WGS_output/seq_cnv//mm10.gc200Base.wig.gz
data_folder=./WGS_output/ready_bam/
output_dic=./WGS_output/seq_cnv/
cat ${output_dic}cnv_list.txt | while read id
#sequenza-utils gc_wiggle -w 200 --fasta ./mouse_ensembl_ref.fa -o ./WGS_output/seq_cnv//mm10.gc200Base.wig.gz
do
echo "source /home/xudeshu/.bashrc
sequenza-utils bam2seqz -n ./WGS_output/ready_bam/WT_cellline_ready.bam -t ${data_folder}${id}_ready.bam --fasta ${refer_fa} \
-gc ${wig_gz} -o ${output_dic}${id}.seqz.gz -f illumina
sequenza-utils seqz_binning --seqz ${output_dic}${id}.seqz.gz -w 200 -o ${output_dic}${id}_small.seqz.gz" > ${output_dic}${id}.sh
echo "library(sequenza)
c_list <- c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','X')
test <- sequenza.extract(\"${output_dic}${id}_small.seqz.gz\", verbose = FALSE, chromosome.list = c_list, parallel = 5)
CP <- sequenza.fit(test, ratio.priority = T, chromosome.list = c_list,ploidy = seq(1, 2, 0.1))
sequenza.results(sequenza.extract = test,
    cp.table = CP, sample.id = \"${id}\",
    out.dir=\"${output_dic}${id}\", CNt.max = 4, ratio.priority = T,chromosome.list = c_list)" > ${output_dic}${id}.R
echo "conda activate sc
Rscript ${output_dic}${id}.R" >> ${output_dic}${id}.sh
echo "qsub -cwd -l vf=5g,h=node6 ${output_dic}${id}.sh" >> ${output_dic}submit.sh
done
sh ${output_dic}submit.sh
