# blobtools

### blast against nt database:

```
module load blast
blastn \
-task megablast \
-query asmnp_female_mkf02_renamed.fa \
-db /nfs/scistore18/vicosgrp/melkrewi/panorpa_assembly_v2/blobtools/ncbi-db/nt \
-outfmt '6 qseqid staxids bitscore std' \
-num_threads 60 \
-max_target_seqs 1 \
-max_hsps 1 \
-evalue 1e-25 \
-out assembly.vs.nt.megablast.out
```
### sort sam file from coverage analysis
```
samtools view -bS CC2U_6.sam | samtools sort /dev/stdin -o CC2U_6.sorted.bam

samtools index CC2U_6.sorted.bam
```
### run blobtools
```
module load anaconda3/2022.05
source /mnt/nfs/clustersw/Debian/bullseye/anaconda3/2022.05/activate_anaconda3_2022.05.txt
conda activate blobtools
/nfs/scistore18/vicosgrp/melkrewi/panorpa_project/blobology/blobtools/blobtools create -i asmnp_female_mkf02_renamed.fa -b CC2U_6.sorted.bam -t assembly.vs.nt.megablast.out -o test && /nfs/scistore18/vicosgrp/melkrewi/panorpa_project/blobology/blobtools/blobtools view -i test.blobDB.json && /nfs/scistore18/vicosgrp/melkrewi/panorpa_project/blobology/blobtools/blobtools plot -i test.blobDB.json
conda deactivate
```
### snail plot
```
module load anaconda3/2022.05
source /mnt/nfs/clustersw/Debian/bullseye/anaconda3/2022.05/activate_anaconda3_2022.05.txt
conda activate btk
blobtools create --fasta asm_np_female_mkf02_01_09_2023_renamed_final_fin.fa Afran_assembly
blobtools add \
    --busco /nfs/scistore18/vicosgrp/melkrewi/Artemia_franciscana_genome_V3/9.disable_post_join/round2/round_3/purge_haplotigs/merge_arcs_rascaf/redundans/redundans_limit/chromonomer/Genome_assembly_with_diff/scaffold_whole_genome_round2/blobtools_final/busco/busco/run_arthropoda_odb10/full_table.tsv \
    Afran_assembly
blobtools view --view snail Afran_assembly
conda deactivate
```

