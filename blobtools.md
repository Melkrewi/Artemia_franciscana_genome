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
