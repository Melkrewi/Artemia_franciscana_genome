# Artemia_franciscana_genome

### Generate consensus from raw bam file:

The first step is to generate the fastq file of the consensus reads from the raw reads using the CCS tool (here we can try --min-passes=3 and --min-passes=8):

Make sure you set the TMP directory to the directory you're working in (adjust the first line)
```
export TMPDIR=/nfs/scistore18/vicosgrp/melkrewi/panorpa_paper/15.CCS/ccs_8_passes/
module load anaconda3/2022.05
source /mnt/nfs/clustersw/Debian/bullseye/anaconda3/2022.05/activate_anaconda3_2022.05.txt
conda activate pbccs
ccs -j 100 --all /nfs/scistore18/vicosgrp/melkrewi/artemia_franciscana_genome_data/r64046_20230525_135433_C02/m64046_230528_074955.subreads.bam ccs_all.fastq.gz
conda deactivate
```
### Assemble the reads using hifiasm:
```
module load hifiasm
hifiasm  -t 100 -o artemia_franciscana.asm ccs.fastq.gz 
```
Get fasta file from gfa. There are three assemblies (the primary, first haplotype and second haplotype, lets run three commands, 1 for each). Make sure you change the names accordingly:
```
awk '/^S/{print ">"$2;print $3}' artemia_franciscana.asm.bp.p_ctg.gfa > artemia_franciscana.asm.bp.p_ctg.fasta
awk '/^S/{print ">"$2;print $3}' BRLI.asm.hic.hap1.p_ctg.gfa > BRLI.asm.hic.hap1.p_ctg.fasta
awk '/^S/{print ">"$2;print $3}' BRLI.asm.hic.hap2.p_ctg.gfa > BRLI.asm.hic.hap2.p_ctg.fasta
```
Let's look at the assembly stats and busco at this stage
```
module load assembly-stats/20170224
assembly-stats artemia_franciscana.asm.bp.p_ctg.fasta
```
We can filter the Pacbio reads using the female short reads:
```
/nfs/scistore18/vicosgrp/melkrewi/panorpa_paper/11.filter_reads/Filtlong/bin/filtlong -1 CC2U_7_1.fastq.gz -2 CC2U_7_2.fastq.gz --min_length 1000 --trim --split 500 ccs.fastq.gz | gzip > output_short.fastq.gz
hifiasm  -t 100 -o artemia_franciscana_filtshortlong.asm output_short.fastq.gz
```
The stats when assembling using the filtered reads:
```
assembly-stats artemia_franciscana_filtshortlong.asm.bp.p_ctg.fasta
stats for artemia_franciscana_filtshortlong.asm.bp.p_ctg.fasta
sum = 1591112320, n = 19113, ave = 83247.65, largest = 1378558
N50 = 122757, n = 3898
N60 = 98595, n = 5344
N70 = 76150, n = 7177
N80 = 54827, n = 9639
N90 = 36805, n = 13184
N100 = 3799, n = 19113
N_count = 0
Gaps = 0
```
The basic statistics for minimum passes 3 and filtering using female short reads:
```
stats for artemia_franciscana_short_minpasses3.asm.bp.p_ctg.fasta
sum = 1462319637, n = 11286, ave = 129569.35, largest = 1953135
N50 = 278035, n = 1473
N60 = 203051, n = 2088
N70 = 142127, n = 2947
N80 = 88917, n = 4248
N90 = 51132, n = 6431
N100 = 3351, n = 11286
N_count = 0
Gaps = 0
```
Code to run BUSCO:
```
module load anaconda3/2022.05

source /mnt/nfs/clustersw/Debian/bullseye/anaconda3/2022.05/activate_anaconda3_2022.05.txt

conda activate busco2

busco -f -i artemia_franciscana.asm.bp.p_ctg.fasta -o busco -l arthropoda_odb10 -m geno -c 40

conda deactivate
```
Add BUSCO score here:
```
BUSCO score
```
### Purge haplotigs 
Purging duplicates with the short reads:
```
module load java
module load bwa
module load samtools
module load python/3.7
export PATH=/nfs/scistore18/vicosgrp/melkrewi/Improved_genome_project/minimap/minimap2/minimap2:$PATH
bwa index A_franciscana.genome.fasta
bwa mem -t 60 A_franciscana.genome.fasta CC2U_7_1_paired.fastq CC2U_7_2_paired.fastq | samtools view -b -o - > CC2U_7.bam
/nfs/scistore18/vicosgrp/melkrewi/Project_confirm_genome_assembly/purge/purge_dups/bin/ngscstat CC2U_7.bam # The program will generate two/three outputs, TX.stat and TX.base.cov which functions the same way as PB.stat and PB.base.cov respectively.  
/nfs/scistore18/vicosgrp/melkrewi/Project_confirm_genome_assembly/purge/purge_dups/bin/calcuts TX.stat > cutoffs 2>calcults.log

/nfs/scistore18/vicosgrp/melkrewi/Project_confirm_genome_assembly/purge/purge_dups/bin/split_fa A_franciscana.genome.fasta  > A_franciscana.genome.fasta.split
/nfs/scistore18/vicosgrp/melkrewi/Improved_genome_project/minimap/minimap2/minimap2 -t 40 -xasm5 -DP A_franciscana.genome.fasta.split A_franciscana.genome.fasta.split | gzip -c - > A_franciscana.genome.fasta.split.self.paf.gz

/nfs/scistore18/vicosgrp/melkrewi/Project_confirm_genome_assembly/purge/purge_dups/bin/purge_dups -2 -T cutoffs -c TX.base.cov A_franciscana.genome.fasta.split.self.paf.gz > dups.bed 2> purge_dups.log

/nfs/scistore18/vicosgrp/melkrewi/Project_confirm_genome_assembly/purge/purge_dups/bin/get_seqs -e dups.bed A_franciscana.genome.fasta
```
Let's look at the assembly stats and busco at this stage too.
Add stats here:
```
Stats
```
Add BUSCO score here.
```
BUSCO score
```
### Scaffolding using ARCS
```
module load samtools
module load minimap2
module load abyss
export PATH=/nfs/scistore18/vicosgrp/melkrewi/panorpa_assembly_v2/12.ARKS/arcs-1.2.5/Arcs/:$PATH
export PATH=/nfs/scistore18/vicosgrp/melkrewi/panorpa_assembly_v2/links/bin/:$PATH
export PATH=/nfs/scistore18/vicosgrp/melkrewi/panorpa_assembly_v2/12.ARKS/arcs-1.2.5/::$PATH
/nfs/scistore18/vicosgrp/melkrewi/panorpa_assembly_v2/12.ARKS/arcs-1.2.5/bin/arcs-make arcs-long draft=purged reads=output_short_minpasses3 t=100 m=8-10000 s=70 c=2 l=2 a=0.7
```

### Scaffold using the franciscana linkage map
If I remember correctly, we would need to remove Ns from the start and ends of sequences
```
/nfs/scistore18/vicosgrp/melkrewi/Improved_genome_project/quickmerge/KPI_and_bmc_2/purge/round2/chromonomer/seqkit -is replace -p "^n+|n+$" -r "" purged.fa > purged_clean.fa
```
Then we generate an agp file:
```
python /nfs/scistore18/vicosgrp/melkrewi/project_save_the_genome_project/chromonomer/chromonomer-1.13/scripts/fasta2agp.py --fasta purged_clean.fa > test.agp
```
Map the markers to the assembly (the marker files are in this repository):
```
module load bwa

module load samtools

bwa index purged_clean.fa

bwa mem -M -t 30 purged_clean.fa markers_first100.fa markers_last100_rc.fa > aligned_paired.sam
```
Then scaffold using chromonomer:
Initially, I think we should try without rescaffold. The linkage_map.tsv file is in this repository. Note: You need to make the output folder before running. (output_paired)
```
module load bwa

module load samtools

mkdir output_paired
/nfs/scistore18/vicosgrp/melkrewi/project_save_the_genome_project/chromonomer/chromonomer-1.13/chromonomer -p linkage_map.tsv --out_path output_paired --alns aligned_paired.sam -a test.agp --fasta purged_clean.fa 
#--rescaffold
```
