./cellranger-3.1.0/cellranger mkgtf \
Mus_musculus.GRCm38.102.Gfp.Lu.gtf \
Mus_musculus.GRCm38.102.Gfp.Lu.filetered.gtf \
--attribute=gene_biotype:protein_coding

./cellranger-3.1.0/cellranger mkref --genome=Mouse_GFP_Lu \
--fasta=Mus_musculus.GRCm38.dna.primary_assembly.Gfp.Lu.fa \
--genes=Mus_musculus.GRCm38.102.Gfp.Lu.filetered.gtf
