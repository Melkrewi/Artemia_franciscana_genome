# Artemia_franciscana_genome

### Generate consensus from raw bam file:

The first step is to generate the fastq file of the consensus reads from the raw reads using the CCS tool:

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
export PATH=/nfs/scistore18/vicosgrp/melkrewi/artemia_franciscana_genome_data/hifiasm_updated/hifiasm-0.19.5/:$PATH
hifiasm -t 100 --hg-size 1g --n-hap 2 -s 0.75 -o artemia_franciscana.asm ccs_all.fastq.gz
awk '/^S/{print ">"$2;print $3}' artemia_franciscana.asm.bp.p_ctg.gfa > artemia_franciscana.asm.bp.p_ctg.fasta
```
### Find Female specific scaffolds:
We find the female specific kmers from the short reads, and use them to remove putative scaffolds from the assembly:
```
module load java

module load bbmap

kmercountexact.sh k=31 in1=CC2U_7_1.fastq.gz in2=CC2U_7_2.fastq.gz out=sfemale_mer_fran.fa mincount=2

bbduk.sh k=31 in=sfemale_mer_fran.fa out=female_specific_mers_31.fasta ref=CC2U_6_1.fastq.gz,CC2U_6_2.fastq.gz -Xmx300g
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
### Scaffolding using ntLink+ARCS
```
export TMPDIR=/nfs/scistore18/vicosgrp/melkrewi/Artemia_franciscana_genome_assembly/7.longstitch/
module load samtools
module load minimap2
module load abyss
export PATH=/nfs/scistore18/vicosgrp/melkrewi/panorpa_assembly_v2/12.ARKS/arcs-1.2.5/Arcs/:$PATH
export PATH=/nfs/scistore18/vicosgrp/melkrewi/panorpa_assembly_v2/links/bin/:$PATH
export PATH=/nfs/scistore18/vicosgrp/melkrewi/panorpa_assembly_v2/12.ARKS/arcs-1.2.5/::$PATH
module load anaconda3/2023.04
source /mnt/nfs/clustersw/Debian/bullseye/anaconda3/2023.04/activate_anaconda3_2023.04.txt
conda activate longstitch_new
/nfs/scistore18/vicosgrp/melkrewi/Artemia_franciscana_genome_assembly/7.longstitch/longstitch-1.0.4/longstitch ntLink-arks draft=purged reads=ccs_all_without_W G=1g t=100 k_arks=20 j=0.05 c=2 l=2 a=0.8
```

### Scaffold using the franciscana linkage map
If I remember correctly, we would need to remove Ns from the start and ends of sequences
```
/nfs/scistore18/vicosgrp/melkrewi/Improved_genome_project/quickmerge/KPI_and_bmc_2/purge/round2/chromonomer/seqkit -is replace -p "^n+|n+$" -r "" purged.fa > purged_clean.fa
```
Then we generate an agp file:
```
module load python/2.7
python /nfs/scistore18/vicosgrp/melkrewi/project_save_the_genome_project/chromonomer/chromonomer-1.13/scripts/fasta2agp.py --fasta purged_clean.fa > test.agp
```
Map the markers to the assembly:

I wrote some commands to clean the excel file with the markers to make life easier. I used csvkit to convert the xlsx file to csv (pip install csvkit) :
```
module load python
wget https://ars.els-cdn.com/content/image/1-s2.0-S0044848621003549-mmc3.xlsx
in2csv 1-s2.0-S0044848621003549-mmc3.xlsx > slaf_markers.csv
sed -i '1,1d' slaf_markers.csv
tr -d ',' < slaf_markers.csv | sed -n '/>Marker/,+1p' > slaf_markers_clean.csv
sed -E '/>/! s/^(.{100}).*/\1/' slaf_markers_clean.csv > markers_1.fa
sed -Ee 's/^.*(.{101})$/\1/' slaf_markers_clean.csv > markers_2.fa
```
mapping:
```
module load bwa

module load samtools

bwa index purged_clean.fa

bwa mem -M -t 30 purged_clean.fa markers_first100.fa markers_last100_rc.fa > aligned_paired.sam
```
Then scaffold using chromonomer:
Initially, I think we should try without rescaffold. The modified linkage_map is in the folder below. Note: You need to make the output folder before running. (output_paired)
```
cp /nfs/scistore18/vicosgrp/melkrewi/Artemia_franciscana_genome_assembly/7.longstitch/ntlinks+arks/chromonomer/linkage_map_modified_no_XB1.tsv .
```
```
module load bwa

module load samtools

mkdir output_paired
/nfs/scistore18/vicosgrp/melkrewi/project_save_the_genome_project/chromonomer/chromonomer-1.13/chromonomer -p linkage_map_modified_no_XB1.tsv --out_path output_paired --alns aligned_paired.sam -a test.agp --fasta purged.fa.k32.w100.z1000.ntLink.scaffolds_c2_m8-10000_cut250_k20_r0.05_e30000_z1000_l2_a0.8.scaffolds_clean.fa
```
### Gap filling using TGSgapcloser 
```
export TMPDIR=/nfs/scistore18/vicosgrp/melkrewi/Artemia_franciscana_genome_V1/3.TGSgapfiller/mkf0.2_female/
#./smrtlink_11.1.0.166339.run
module load anaconda3/2022.05
source /mnt/nfs/clustersw/Debian/bullseye/anaconda3/2022.05/activate_anaconda3_2022.05.txt
conda activate TGScloser
sed '/^>/ s/,.*//' CHRR_integrated.fa > CHRR_integrated_clean.fa
tgsgapcloser --scaff CHRR_integrated_clean.fa --reads ccs_all_without_W_0.2.fastq.gz --output CHRR_integrated_TGS --thread 80 --ne --tgstype pb
conda deactivate
```
### Polish with racon + merfin
```
export TMPDIR=/nfs/scistore18/vicosgrp/melkrewi/Artemia_franciscana_genome_V1/3.TGSgapfiller/mkf0.2_female/T2T/
module load samtools
module load bcftools
module load minimap2
export PATH=/nfs/scistore18/vicosgrp/melkrewi/Artemia_franciscana_genome_assembly/12.T2T/racon/build/bin/:$PATH
export PATH=/nfs/scistore18/vicosgrp/melkrewi/Artemia_franciscana_genome_assembly/12.T2T/merfin/build/bin/:$PATH
export PATH=/nfs/scistore18/vicosgrp/melkrewi/Artemia_franciscana_genome_assembly/12.T2T/Winnowmap/bin/:$PATH
export PATH=/nfs/scistore18/vicosgrp/melkrewi/Artemia_franciscana_genome_assembly/12.T2T/T2T-Polish-1.0/automated_polishing/:$PATH
module load anaconda3/2022.05
source /mnt/nfs/clustersw/Debian/bullseye/anaconda3/2022.05/activate_anaconda3_2022.05.txt
conda activate T2Tpolishing
k=21
thread=80
# Construct k-mer db
cat CC2U_6_1.fastq.gz CC2U_6_2.fastq.gz > CC2U_6.fastq.gz
meryl count k=${k} threads=${thread} CC2U_6.fastq.gz output reads.meryl
# Collect histogram for GenomeScope
meryl histogram reads.meryl > reads.hist
# Exclude frequency = 1 k-mers
meryl greater-than 1 reads.meryl output reads.gt1.meryl
automated-polishing.sh ${thread} 1 CHRR_integrated_TGS.scaff_seqs ccs_all_without_W_0.2.fastq.gz reads.gt1.meryl T2T_polished
conda deactivate
```
### Polish using Nextpolish2:
```
export TMPDIR=/nfs/scistore18/vicosgrp/melkrewi/Artemia_franciscana_genome_assembly/9.nextpolish2/
module load samtools
module load minimap2
export PATH=/nfs/scistore18/vicosgrp/melkrewi/Artemia_franciscana_genome_assembly/9.nextpolish2/NextPolish2/target/release/:$PATH
export PATH=/nfs/scistore18/vicosgrp/melkrewi/artemia_franciscana_genome_data/hifiasm_updated/yak/:$PATH
minimap2 -ax map-hifi -t 40 T2T_polished.iter_1.consensus.fasta ccs_all_without_W_0.2.fastq.gz | samtools sort -o hifi.map.sort.bam -
samtools index hifi.map.sort.bam

yak count -o k21.yak -k 21 -b 37 -t 40 <(zcat CC2U_6_*.fastq.gz) <(zcat CC2U_6_*.fastq.gz)
yak count -o k31.yak -k 31 -b 37 -t 40 <(zcat CC2U_6_*.fastq.gz) <(zcat CC2U_6_*.fastq.gz)

nextPolish2 -r -t 5 hifi.map.sort.bam T2T_polished.iter_1.consensus.fasta k21.yak k31.yak > asm.np_female_mkf0.2.fa
```
Let's look at the assembly stats and busco at this stage
```
module load assembly-stats/20170224
assembly-stats artemia_franciscana_all_no_W.asm.bp.p_ctg.fasta
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
