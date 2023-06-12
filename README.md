# Artemia_franciscana_genome

### Generate consensus from raw bam file (if we get raw bam file):

The first step is to generate the fastq file of the consensus reads from the raw reads using the CCS tool (here we can try --min-passes=3 and --min-passes=8):
```
module load anaconda3/2022.05
source ~/anaconda3/2022.05/activate_anaconda3_2022.05.txt
conda activate pbccs
ccs --minLength=100 ~/r64046_20230118_140133_C02/m64046_230121_101921.subreads.bam ccs.fastq.gz -j 40
conda deactivate
```
### Assemble the reads using hifiasm:
```
module load hifiasm
hifiasm  -t 100 -o artemia_franciscana.asm ccs.fastq.gz 
```
Let's look at the assembly stats and busco at this stage
### Purge haplotigs 
We can use the primary assembly for now I guess, and purge it to remove the duplicates and haplotigs.
```
module load python/3.7

export PATH=~/minimap2/minimap2:$PATH

~/minimap2/minimap2 -t 50 -xmap-pb artemia_franciscana.asm.bp.p_ctg.fasta ccs.fastq.gz | gzip -c - > artemia_franciscana.asm.bp.p_ctg.paf.gz
~/purge_dups/bin/pbcstat *.paf.gz
~/purge_dups/bin/calcuts PB.stat > cutoffs 2>calcults.log
~/purge_dups/bin/split_fa artemia_franciscana.asm.bp.p_ctg.fasta  > artemia_franciscana.asm.bp.p_ctg.fasta.split
~/minimap2/minimap2 -t 50 -xasm5 -DP artemia_franciscana.asm.bp.p_ctg.fasta.split artemia_franciscana.asm.bp.p_ctg.fasta.split | gzip -c - > artemia_franciscana.asm.bp.p_ctg.fasta.split.self.paf.gz
~/purge_dups/bin/purge_dups -2 -T cutoffs -c PB.base.cov artemia_franciscana.asm.bp.p_ctg.fasta.split.self.paf.gz > dups.bed 2> purge_dups.log
~/purge_dups/bin/get_seqs -e dups.bed artemia_franciscana.asm.bp.p_ctg.fasta
```
### Scaffold using the franciscana linkage map
Map the markers to the assembly:
```
module load bwa

module load samtools

bwa index purged.fa

bwa mem -M -t 30 purged.fa markers_first100.fa markers_last100_rc.fa > aligned_paired.sam
```
The scaffold using chromonomer:
```
module load bwa

module load samtools

/nfs/scistore03/vicosgrp/melkrewi/project_save_the_genome_project/chromonomer/chromonomer-1.13/chromonomer -p linkage_map.tsv --out_path output_paired --alns aligned_paired.sam -a test.agp --fasta purged_clean.fa 
#--rescaffold
```


