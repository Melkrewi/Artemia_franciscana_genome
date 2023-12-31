# Artemia_franciscana_genome

# Currently testing a different tool to call consensensus: deepconsensus
https://github.com/google/deepconsensus/blob/r1.2/docs/quick_start.md

# Goals:
1. Clean coverage in the S0 (high importance)
2. N50 >= 1MB before scaffolding (medium improtance, but we can survive without it)
3. BUSCO score complete >= 90% and a reasonable number of duplicated buscos (medium importance)
4. Some of the W-chromosome (the icing on the cake).
# Current Imaginary Plan:
1. Use SRY to identify the W/female-specific reads.
2. Remove those reads from the pacbio fastq, and assemble using Hifiasm and Nextdenovo.
3. Purge using long-reads/short reads, depending on what works better.
4. Scaffold using Arcs.
5. Scaffold using Chromonomer with the modified linkage map, but hopefully without the xb-1 (if we manage to remove the W specific reads).
6. Check stats, busco, and coverage.
7. Assemble the W reads using the same pipeline to see if we get anything, and then add it to the assembly (hopefully we can assemble xb-1).
8. Run blobtools to check for contamination, and remove any if found.
9. Do final stats, coverage and busco.
10. Annotate the genome.
11. Make final figures and write the paper/genome report.

# Things that we tried/could still potentially try:
1. I tested using the ultra long option with reads > 50000 bp, and it seems promising (N50 improved a bit for 8 passes), but let's see if we can get good results without it. The ultralong reads option is for a specific ONT technology, but I've seen other people use the option with the longest pacbio subreads.
2. I tried to self-correct the reads using Fec, but the assembly was really bad.
3. I tried to correct the reads with Lordec (short reads), but it took forever, and I cancelled the job after three days. (We can start it again at some point, but I think it might take 10 days or so to run).
4. We can try other ways to align the markers for chromonomer (low priority).
5. I checked for adapter contamination, but it was only in around 200 reads out of 5 million, so not a problem.
6. I tried hifiasm with smaller kmer size, but it did not improve the assembly.
7. I could not find how the -x and -y parameters in hifiasm work, but we can try changing them at some point (very low priority).
8. There are many CLR assemblers, which we can try (Miniasm/Raven/Flye). I tried flye but with hifi option, not CLR.

# Initial pipeline

### Generate consensus from raw bam file (if we get raw bam file):

The first step is to generate the fastq file of the consensus reads from the raw reads using the CCS tool (here we can try --min-passes=3 and --min-passes=8):

Make sure you set the TMP directory to the directory you're working in (adjust the first line)
```
export TMPDIR=/nfs/scistore18/vicosgrp/melkrewi/panorpa_paper/15.CCS/ccs_8_passes/
module load anaconda3/2022.05
source /mnt/nfs/clustersw/Debian/bullseye/anaconda3/2022.05/activate_anaconda3_2022.05.txt
conda activate pbccs
ccs -j 100 --min-passes 8 ~/r64046_20230118_140133_C02/m64046_230121_101921.subreads.bam ccs.fastq.gz
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
Let's look at the assembly stats and busco at this stage too.
Add stats here:
```
Stats
```
Add BUSCO score here.
```
BUSCO score
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
### Additional Stuff (Filter the pacbio reads for quality/using short reads) with filtlong:
We can filter for length > 10000:
```
/nfs/scistore18/vicosgrp/melkrewi/panorpa_paper/11.filter_reads/Filtlong/bin/filtlong ccs.fastq --min_length 10000 | gzip > output_filtlong.fastq.gz
```
We can use the female franciscana short reads here (change the file names)
```
export TMPDIR=/nfs/scistore18/vicosgrp/melkrewi/panorpa_paper/11.filter_reads/

/nfs/scistore18/vicosgrp/melkrewi/panorpa_paper/11.filter_reads/Filtlong/bin/filtlong -1 199208_S1_L003_1_paired.fastq -2 199208_S1_L003_2_paired.fastq --min_length 1000 --trim --split 500 ccs.fastq.gz | gzip > output.fastq.gz
```
### Flye (hifi option)
```
export TMPDIR=/nfs/scistore18/vicosgrp/melkrewi/panorpa_paper/9.random_sample_on_steroids/
module load python
python /nfs/scistore18/vicosgrp/melkrewi/panorpa_paper/5.Flye/Flye/bin/flye --pacbio-hifi ccs.fastq.gz --genome-size 1g --threads 100 --out-dir ./assembly --scaffold
```
### hifiasm (auto coverage)
```
export PATH=/nfs/scistore18/vicosgrp/melkrewi/panorpa_assembly_v2/63.hifiasm_updated/hifiasm-0.19.4/:$PATH
hifiasm -t 100 --hom-cov 20 --n-hap 4 -o artemia_franciscana.asm ccs.fastq.gz 
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
### Correcting the reads using Hicanu:
```
module load perl

module load java
module load gnuplot

export TMPDIR=/nfs/scistore18/vicosgrp/melkrewi/artemia_franciscana_genome_data/hicanu/
module load python
export PATH=/nfs/scistore18/vicosgrp/melkrewi/artemia_franciscana_genome_data/hicanu/canu-2.2/bin/:$PATH
canu -correct \
 -p asm -d artemia_hifi_correct \
 genomeSize=1g \
 -pacbio output_short_minpasses3.fastq.gz maxThreads=30 maxMemory=200G useGrid=false
```
### Something else definetly worth testing
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
### Most recent idea
Currently thinking of how to phase the Z and W haplotypes. I tried trio-binning and it seems to kinda work, but not so well. I also tried filtering with male reads to remove W and assemble just the Z, but that did not work.
I am testing SRY now to get the W specific reads, then the idea would be to remove them from the reads file to get the Z seperately and then attempt to assemble the W on its own. 
# Promising Assembler Nextdenovo
I tested it with the 8 passes and filtered only, and it did slightly better than hifiasm, but let's see. I'll keep you updated.
Assembly with no scaffolding:
```
sum = 947684741, n = 6815, ave = 139058.66, largest = 1503388
N50 = 187352, n = 1293
N60 = 147642, n = 1867
N70 = 115369, n = 2589
N80 = 88576, n = 3528
N90 = 64194, n = 4781
N100 = 12472, n = 6815
N_count = 0
Gaps = 0
```
Busco score with no scaffolding same assembly:
```
C:77.6%[S:72.7%,D:4.9%],F:11.2%,M:11.2%,n:1013
```
