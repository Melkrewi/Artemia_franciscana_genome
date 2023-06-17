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
Get fasta file from gfa. There are three assemblies (the primary, first haplotype and second haplotype, lets run three commands, 1 for each). Make sure you change the names accordingly:
```
awk '/^S/{print ">"$2;print $3}' artemia_franciscana.asm.bp.p_ctg.gfa > artemia_franciscana.asm.bp.p_ctg.fasta
awk '/^S/{print ">"$2;print $3}' BRLI.asm.hic.hap1.p_ctg.gfa > BRLI.asm.hic.hap1.p_ctg.fasta
awk '/^S/{print ">"$2;print $3}' BRLI.asm.hic.hap2.p_ctg.gfa > BRLI.asm.hic.hap2.p_ctg.fasta
```
Let's look at the assembly stats and busco at this stage
### Purge haplotigs 
We can use the primary assembly for now I guess, and purge it to remove the duplicates and haplotigs.
```
module load python/3.7
module load minimap2

minimap2 -t 50 -xmap-pb artemia_franciscana.asm.bp.p_ctg.fasta ccs.fastq.gz | gzip -c - > artemia_franciscana.asm.bp.p_ctg.paf.gz
/nfs/scistore18/vicosgrp/melkrewi/Project_confirm_genome_assembly/purge/purge_dups/bin/pbcstat *.paf.gz
/nfs/scistore18/vicosgrp/melkrewi/Project_confirm_genome_assembly/purge/purge_dups/bin/calcuts PB.stat > cutoffs 2>calcults.log
/nfs/scistore18/vicosgrp/melkrewi/Project_confirm_genome_assembly/purge/purge_dups/bin/split_fa artemia_franciscana.asm.bp.p_ctg.fasta  > artemia_franciscana.asm.bp.p_ctg.fasta.split
minimap2 -t 50 -xasm5 -DP artemia_franciscana.asm.bp.p_ctg.fasta.split artemia_franciscana.asm.bp.p_ctg.fasta.split | gzip -c - > artemia_franciscana.asm.bp.p_ctg.fasta.split.self.paf.gz
/nfs/scistore18/vicosgrp/melkrewi/Project_confirm_genome_assembly/purge/purge_dups/bin/purge_dups -2 -T cutoffs -c PB.base.cov artemia_franciscana.asm.bp.p_ctg.fasta.split.self.paf.gz > dups.bed 2> purge_dups.log
/nfs/scistore18/vicosgrp/melkrewi/Project_confirm_genome_assembly/purge/purge_dups/bin/get_seqs -e dups.bed artemia_franciscana.asm.bp.p_ctg.fasta
```
### Scaffold using the franciscana linkage map
If I remember correctly, we would need to remove Ns from the start and ends of sequences
```
./seqkit -is replace -p "^n+|n+$" -r "" purged.fa > purged_clean.fa
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
Initially, I think we should try without rescaffold. The linkage_map.tsv file is in this repository
```
module load bwa

module load samtools

/nfs/scistore18/vicosgrp/melkrewi/project_save_the_genome_project/chromonomer/chromonomer-1.13/chromonomer -p linkage_map.tsv --out_path output_paired --alns aligned_paired.sam -a test.agp --fasta purged_clean.fa 
#--rescaffold
```


