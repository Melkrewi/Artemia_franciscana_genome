# Fst analysis
```
module load samtools
module load STAR
srun STAR --runThreadN 20 --runMode genomeGenerate --genomeDir asm_fran --genomeFastaFiles asm_np_female_mkf02_01_09_2023_renamed_final_fin.fa
STAR --runThreadN 40 --genomeDir ./asm_fran --readFilesIn /nfs/scistore18/vicosgrp/melkrewi/Artemia_franciscana_genome_V1/5.FST/male/male_1.fastq /nfs/scistore18/vicosgrp/melkrewi/Artemia_franciscana_genome_V1/5.FST/male/male_2.fastq  --outFileNamePrefix afran_male.
samtools view -q 20 -bS afran_male.Aligned.out.sam | samtools sort -o afran_male_q20.sam.sort
STAR --runThreadN 40 --genomeDir ./asm_fran --readFilesIn /nfs/scistore18/vicosgrp/melkrewi/Artemia_franciscana_genome_V1/5.FST/female/female_1.fastq /nfs/scistore18/vicosgrp/melkrewi/Artemia_franciscana_genome_V1/5.FST/female/female_2.fastq  --outFileNamePrefix afran_female.
samtools view -q 20 -bS afran_female.Aligned.out.sam | samtools sort -o afran_female_q20.sam.sort
```
my_bams.fofn contains the following two lines 
```
afran_male_q20.sam.sort
afran_female_q20.sam.sort
```
generate the pileup file:
```
module load samtools
module load parallel
samtools faidx asm_np_female_mkf02_01_09_2023_renamed_final_fin.fa
samtools index afran_male_q20.sam.sort
samtools index afran_female_q20.sam.sort
parallel --jobs 40 --keep-order --colsep '\t' samtools mpileup -b my_bams.fofn -r {1} :::: asm_np_female_mkf02_01_09_2023_renamed_final_fin.fa.fai > my.pileup
```
generate sync file using grenedalf:
```
export PATH=/nfs/scistore18/vicosgrp/melkrewi/Artemia_franciscana_genome_V1/5.FST/grenedalf/bin/:$PATH
export PATH=/nfs/scistore18/vicosgrp/melkrewi/Artemia_franciscana_genome_V1/7.poopolation/popoolation2_1201/:$PATH
grenedalf sync --pileup-path my.pileup --pileup-min-base-qual 20 --threads 100 --file-prefix afran --reference-genome-fasta-file asm_np_female_mkf02_01_09_2023_renamed_final_fin.fa --allow-file-overwriting
```
Use gtf to get FST per gene:
```
export PATH=/nfs/scistore18/vicosgrp/melkrewi/Artemia_franciscana_genome_V1/5.FST/grenedalf/bin/:$PATH
export PATH=/nfs/scistore18/vicosgrp/melkrewi/Artemia_franciscana_genome_V1/7.poopolation/popoolation2_1201/:$PATH
grep exons braker_filconagat_isoORF3_sup100compgen.gtf > braker_filconagat_isoORF3_sup100compgen_exons.gtf
perl /nfs/scistore18/vicosgrp/melkrewi/Artemia_franciscana_genome_V1/7.poopolation/popoolation2_1201/create-genewise-sync.pl --input afransync.sync --gtf braker_filconagat_isoORF3_sup100compgen_exons.gtf --output afransync_exons.sync
perl /nfs/scistore18/vicosgrp/melkrewi/Artemia_franciscana_genome_V1/7.poopolation/popoolation2_1201/fst-sliding.pl --min-count 5 --min-coverage 30 --max-coverage 2000 --pool-size 10 --min-covered-fraction 0.0 --window-size 1000000 --step-size 1000000 --input afransync_exons.sync --suppress-noninformative --output afransync_exons.fst
```

