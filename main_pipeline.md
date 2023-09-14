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
### Find Female specific kmers:
We find the female specific kmers from the short reads, and use them to remove putative W-reads from the fastq file:
```
#mkf=0.2
module load bbmap

kmercountexact.sh k=21 in1=CC2U_7_1.fastq in2=CC2U_7_2.fastq out=sfemale_mer_fran.fa mincount=2

bbduk.sh k=21 in=sfemale_mer_fran.fa out=female_specific_mers.fasta ref=CC2U_6_1.fastq,CC2U_6_2.fastq -Xmx300g

bbduk.sh k=21 in=ccs_all.fastq outm=ccs_female_specific_0.2.fastq ref=female_specific_mers.fasta mkf=0.2 -Xmx300g

cat ccs_female_specific_0.2.fastq | awk 'NR%4==1' | sed 's/@//' > ccs_female_specific_0.2.fastq.readsID
cat ccs_all.fastq | awk 'NR%4==1' | sed 's/@//' > ccs_all.fastq.readsID
grep -f ccs_female_specific_0.2.fastq.readsID ccs_all.fastq.readsID -v > remaining.list
module load seqtk
seqtk subseq ccs_all.fastq remaining.list | gzip - > ccs_all_without_W_0.2.fastq.gz
```

### Assemble the reads using hifiasm:
```
export PATH=/nfs/scistore18/vicosgrp/melkrewi/panorpa_assembly_v2/63.hifiasm_updated/hifiasm-0.19.4/:$PATH
hifiasm -t 100 --hg-size 1g --n-hap 4 -r 5 -s 0.7 -N 150 -o artemia_franciscana_allnoWmkf02.asm ccs_all_without_W_0.2.fastq.gz
awk '/^S/{print ">"$2;print $3}' artemia_franciscana_allnoWmkf02.asm.bp.p_ctg.gfa > artemia_franciscana_allnoWmkf02.asm.bp.p_ctg.fasta
```
### Purge haplotigs 
Purging duplicates with the short reads:
```
module load java
module load bwa
module load samtools
module load python/3.7
export PATH=/nfs/scistore18/vicosgrp/melkrewi/Improved_genome_project/minimap/minimap2/minimap2:$PATH
bwa index artemia_franciscana_allnoWmkf02.asm.bp.p_ctg.fasta
bwa mem -t 80 artemia_franciscana_allnoWmkf02.asm.bp.p_ctg.fasta CC2U_7_1.fastq.gz CC2U_7_2.fastq.gz | samtools view -b -o - > CC2U_7.bam
/nfs/scistore18/vicosgrp/melkrewi/Project_confirm_genome_assembly/purge/purge_dups/bin/ngscstat CC2U_7.bam
/nfs/scistore18/vicosgrp/melkrewi/Project_confirm_genome_assembly/purge/purge_dups/bin/calcuts TX.stat > cutoffs 2>calcults.log
/nfs/scistore18/vicosgrp/melkrewi/Project_confirm_genome_assembly/purge/purge_dups/bin/split_fa artemia_franciscana_allnoWmkf02.asm.bp.p_ctg.fasta > artemia_franciscana_allnoWmkf02.asm.bp.p_ctg.fasta.split
/nfs/scistore18/vicosgrp/melkrewi/Improved_genome_project/minimap/minimap2/minimap2 -t 80 -xasm5 -DP artemia_franciscana_allnoWmkf02.asm.bp.p_ctg.fasta.split artemia_franciscana_allnoWmkf02.asm.bp.p_ctg.fasta.split | gzip -c - > artemia_franciscana_allnoWmkf02.asm.bp.p_ctg.fasta.split.self.paf.gz
/nfs/scistore18/vicosgrp/melkrewi/Project_confirm_genome_assembly/purge/purge_dups/bin/purge_dups -2 -T cutoffs -c TX.base.cov  artemia_franciscana_allnoWmkf02.asm.bp.p_ctg.fasta.split.self.paf.gz > dups.bed 2> purge_dups.log
/nfs/scistore18/vicosgrp/melkrewi/Project_confirm_genome_assembly/purge/purge_dups/bin/get_seqs -e dups.bed artemia_franciscana_allnoWmkf02.asm.bp.p_ctg.fasta
```
### Scaffolding using Rascaf
First map the male RNA reads to the genome:
```
module load samtools
module load hisat2
export TMPDIR=/nfs/scistore18/vicosgrp/melkrewi/Artemia_franciscana_genome_V3/9.disable_post_join/round2/round_3/purge_haplotigs/rascaf/
mkdir hisat2
hisat2-build purged.fa genome_index
for i in *_1.fastq
do
    prefix=$(basename $i _1.fastq)
    hisat2 --phred33 -p 50 --novel-splicesite-outfile hisat2/${prefix}_splicesite.txt -S hisat2/${prefix}_accepted_hits.sam -x genome_index -1 ${prefix}_1.fastq -2 ${prefix}_2.fastq --rna-strandness RF
    samtools view -@ 25 -bS -o hisat2/${prefix}_accepted_hits.bam hisat2/${prefix}_accepted_hits.sam
    samtools sort -@ 25 -o hisat2/${prefix}_accepted_hits.sorted.bam hisat2/${prefix}_accepted_hits.bam
done
```
Then scaffold the assembly as follows:
```
export TMPDIR=/nfs/scistore18/vicosgrp/melkrewi/Artemia_franciscana_genome_V3/9.disable_post_join/round2/round_3/purge_haplotigs/rascaf/

export PATH=/nfs/scistore18/vicosgrp/melkrewi/project_w/december_2021/scaffod_assembly/rescaf/rascaf/:$PATH

rascaf -b ./hisat2/60541_accepted_hits.sorted.bam -o 60541_rascaf -f purged.fa

rascaf -b ./hisat2/60542_accepted_hits.sorted.bam -o 60542_rascaf -f purged.fa

rascaf -b ./hisat2/60543_accepted_hits.sorted.bam -o 60543_rascaf -f purged.fa

rascaf -b ./hisat2/60544_accepted_hits.sorted.bam -o 60544_rascaf -f purged.fa

rascaf-join -r 60541_rascaf.out -r 60542_rascaf.out -r 60543_rascaf.out -r 60544_rascaf.out -o assembly_scaffold_males
```
### Scaffolding using longstitch
```
export TMPDIR=/nfs/scistore18/vicosgrp/vbett/Artemia_hififinal/genome_mkf0.2
module load samtools
module load minimap2
module load abyss
export PATH=/nfs/scistore18/vicosgrp/melkrewi/panorpa_assembly_v2/12.ARKS/arcs-1.2.5/Arcs/:$PATH
export PATH=/nfs/scistore18/vicosgrp/melkrewi/panorpa_assembly_v2/links/bin/:$PATH
export PATH=/nfs/scistore18/vicosgrp/melkrewi/panorpa_assembly_v2/12.ARKS/arcs-1.2.5/::$PATH
module load anaconda3/2023.04
source /mnt/nfs/clustersw/Debian/bullseye/anaconda3/2023.04/activate_anaconda3_2023.04.txt
source activate /nfs/scistore18/vicosgrp/melkrewi/.conda/envs/longstitch_new
/nfs/scistore18/vicosgrp/melkrewi/Artemia_franciscana_genome_assembly/7.longstitch/longstitch-1.0.4/longstitch ntLink-arks draft=purged reads=ccs_all_without_W_0.2 G=1g t=100 k_arks=20 j=0.05 c=2 l=2 a=0.8
```
### Identify and keep the merges that appear in both approaches:
First map the two scaffolded assemblies to the purged assembly using minimap2:
```
export PATH=/nfs/scistore18/vicosgrp/melkrewi/Improved_genome_project/minimap/minimap2/minimap2:$PATH
/nfs/scistore18/vicosgrp/melkrewi/Improved_genome_project/minimap/minimap2/minimap2 -x asm5 -t 60 assembly_scaffold_males.fa purged.fa > rascaf_vs_purged.paf
/nfs/scistore18/vicosgrp/melkrewi/Improved_genome_project/minimap/minimap2/minimap2 -x asm5 -t 60 purged.fa.k32.w100.z1000.ntLink.scaffolds_c2_m8-10000_cut250_k20_r0.05_e30000_z1000_l2_a0.8.scaffolds.fa purged.fa > longstitch_vs_purged.paf
```
Then get the consensus of the merges using this python code:
```
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy
import seaborn as sns
import os
rascaf=pd.read_csv("rascaf_vs_purged.paf",sep="\s+",header=None)
rascaf=rascaf.sort_values(10, ascending=False).drop_duplicates([0])
rascaf=rascaf[rascaf[10]>=0.3*rascaf[1]]
longstitch=pd.read_csv("longstitch_vs_purged.paf",sep="\s+",header=None)
longstitch=longstitch.sort_values(10, ascending=False).drop_duplicates([0])
longstitch=longstitch[longstitch[10]>=0.3*longstitch[1]]
longstitch_vs_rascaf=pd.merge(rascaf,longstitch,on=0)
longstitch_vs_rascaf=pd.merge(rascaf,longstitch,on=0)[[0,'1_x','5_x','6_x','10_x','5_y','6_y','10_y']].drop_duplicates([0,'1_x','5_x','6_x','5_y','6_y'])
longstitch_vs_rascaf_results=longstitch_vs_rascaf[longstitch_vs_rascaf.duplicated(['5_x','5_y'], keep=False)].sort_values('5_y')
d = []
for i in np.arange(0,len(longstitch_vs_rascaf_results[0])):
    ref=longstitch_vs_rascaf_results[longstitch_vs_rascaf_results[0]==longstitch_vs_rascaf_results.iloc[i][0]]
    d.append(
        {
            'contigs': ','.join(longstitch_vs_rascaf_results[(longstitch_vs_rascaf_results['5_x']==ref['5_x'].iloc[0])&(longstitch_vs_rascaf_results['5_y']==ref['5_y'].iloc[0])][0].to_list()),
            'rascaf_scaffold': longstitch_vs_rascaf_results[(longstitch_vs_rascaf_results['5_x']==ref['5_x'].iloc[0])&(longstitch_vs_rascaf_results['5_y']==ref['5_y'].iloc[0])]['5_x'].iloc[0],
            'longstitch_scaffold':  longstitch_vs_rascaf_results[(longstitch_vs_rascaf_results['5_x']==ref['5_x'].iloc[0])&(longstitch_vs_rascaf_results['5_y']==ref['5_y'].iloc[0])]['5_y'].iloc[0]
        }
    )
final_results=pd.DataFrame(d).drop_duplicates()
final_results.to_csv('final_results_longstitch_vs_rascaf.txt',sep='\t',index=False)
```
The resulting scaffolds were generated using ragtag with the following commands:
```
module load samtools
module load anaconda3/2022.05
source /mnt/nfs/clustersw/Debian/bullseye/anaconda3/2022.05/activate_anaconda3_2022.05.txt
conda activate ragtag
for i in {2..689}
do
 cat final_results_longstitch_vs_rascaf.txt | sed "${i}q;d" | cut -f1 | sed 's/,/\n/g' > $i.contig.txt
 cat final_results_longstitch_vs_rascaf.txt | sed "${i}q;d" | cut -f2 > $i.reference.txt
 perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' $i.contig.txt purged.fa > $i.contig.fasta
 perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' $i.reference.txt assembly_scaffold_males.fa > $i.reference.fasta
 ragtag.py scaffold $i.reference.fasta $i.contig.fasta -t 40 -o ./$i.ragtag
 cat $i.contig.fasta >> pre_scaffolded_sequences.fasta
 cat ./$i.ragtag/ragtag.scaffold.fasta >> scaffolded_sequences.fasta
 rm $i.contig.txt
 rm $i.reference.txt
 rm $i.contig.fasta
 rm $i.reference.fasta
 rm -r $i.ragtag
 rm $i.contig.fasta.fai
 echo $i
done
perl fasta_len.pl > pre_scaffolded_sequences.fasta.len
cat pre_scaffolded_sequences.fasta.len | cut -f1 > remove.list
perl -ne 'if(/^>(\S+)/){$c=!$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' remove.list purged.fa > purged_without_scaffolded.fasta
cat purged_without_scaffolded.fasta scaffolded_sequences.fasta > purged_with_scaffolded.fasta
awk '/^>/{print ">scaffold_" ++i; next}{print}' < purged_with_scaffolded.fasta > purged_with_scaffolded_renamed.fasta
```
### Scaffolding using redundans and mate pairs:
```
export TMPDIR=/nfs/scistore18/vicosgrp/melkrewi/Artemia_franciscana_genome_V3/9.disable_post_join/round2/round_3/purge_haplotigs/merge_arcs_rascaf/redundans/redundans_limit
module load anaconda3/2023.04
source /mnt/nfs/clustersw/Debian/bullseye/anaconda3/2023.04/activate_anaconda3_2023.04.txt
conda activate redundans
redundans.py -v -i SRR6980924_1_fixed.fastq SRR6980924_2_fixed.fastq -f purged_with_scaffolded_renamed.fasta -o run_short-scaffolding -t 60 --limit 1 --noreduction --nogapclosing -j 10 --nocleaning
```
## Coverage analysis to identify S0
```
module load bowtie2/2.4.4
module load soap/coverage
bowtie2-build scaffolds.fa Afran_new_genome
bowtie2 -x Afran_new_genome -1 CC2U_6_1.fastq.gz -2 CC2U_6_2.fastq.gz --end-to-end --sensitive -p 100 -S CC2U_6.sam
bowtie2 -x Afran_new_genome -1 CC2U_7_1.fastq.gz -2 CC2U_7_2.fastq.gz --end-to-end --sensitive -p 100 -S CC2U_7.sam
grep -vw "XS:i" CC2U_6.sam > CC2U_6_unique.sam
grep -vw "XS:i" CC2U_7.sam > CC2U_7_unique.sam
soap.coverage -sam -cvg -i CC2U_6_unique.sam -onlyuniq -p 100 -refsingle scaffolds.fa -window CC2U_6_window 10000
soap.coverage -sam -cvg -i CC2U_7_unique.sam -onlyuniq -p 100 -refsingle scaffolds.fa -window CC2U_7_window 10000
```
Identify S0 scaffolds from female/male coverage with the following python code:
```
#Coverage scaffolds

artemia_male_scaf=pd.read_csv("CC2U_6_window_redundans_all_coverage_S0",sep="\s+",header=None)
artemia_female_scaf=pd.read_csv("CC2U_7_window_redundans_all_coverage_S0",sep="\s+",header=None)
artemia_male_scaf['log2fsm']=np.log2((artemia_female_scaf[3]+0.001)/(artemia_male_scaf[3]+0.001))#artemia_male[3]/artemia_female[3]
artemia_scaf=artemia_male_scaf
artemia_scaf=artemia_scaf.groupby([0]).median().reset_index()
artemia_scaf[artemia_scaf['log2fsm']<np.median(artemia_scaf['log2fsm'])-0.5][0].to_csv('artemia_Z_scaf_rar_all.txt',index=False)
```
Scaffold the S0 region first using the linkage map:
```
grep LG6 linkage_map_modified_no_XB1.tsv > LG6_map.txt
cat LG6_map.txt | cut -f2 > LG6_markers.txt
perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' LG6_markers.txt markers_1.fa > markers_LG6_1.fa
perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' LG6_markers.txt markers_2.fa > markers_LG6_2.fa
perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' artemia_Z_scaf_rar.txt scaffolds.fa  > Artemia_Z_scaf_rar.fasta
perl -ne 'if(/^>(\S+)/){$c=!$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' artemia_Z_scaf_rar.txt scaffolds.fa > purged_without_differentiated_region.fasta
```
As some of the S0 scaffolds had regions with higher coverage that appear to be misassembled, we break the S0 scaffolds into contigs:
```
perl /nfs/scistore18/vicosgrp/melkrewi/Artemia_franciscana_genome_V3/9.disable_post_join/round2/round_3/purge_haplotigs/merge_arcs_rascaf/redundans/redundans_limit/chromonomer/Genome_assembly_with_diff/GenoToolBox/SeqTools/BreakScaffolds -f Artemia_Z_scaf_rar.fasta -o broken.fasta -n 100
```
Map LG6 markers to the contigs:
```
module load bwa
bwa index broken.fasta.contigs.fasta
bwa mem -M -t 50 broken.fasta.contigs.fasta markers_LG6_1.fa markers_LG6_2.fa > aligned_paired_broken_LG6_markers.sam
```
generate the agp file:
```
module load python/2.7
python /nfs/scistore18/vicosgrp/melkrewi/project_save_the_genome_project/chromonomer/chromonomer-1.13/scripts/fasta2agp.py --fasta broken.fasta.contigs.fasta > LG6_broken.agp
```
Scaffold using Chromonomer:
```
mkdir output_paired_diff_LG6_markers_broken
/nfs/scistore18/vicosgrp/melkrewi/project_save_the_genome_project/chromonomer/chromonomer-1.13/chromonomer -p LG6_map.txt --out_path output_paired_diff_LG6_markers_broken --alns aligned_paired_broken_LG6_markers.sam -a LG6_broken.agp --fasta broken.fasta.contigs.fasta
```
We add the resulting scaffolds to the genome:
```
cat ../output_paired_diff_LG6_markers_broken/CHRR_integrated.fa > Z_diff.fasta
cat ../purged_without_differentiated_region.fasta Z_diff.fasta > purged_with_differentiated_region.fasta
```
Then we do two rounds of scaffolding using Chromonomer, the first to identify the location of the differentiated region without splitting, and then we modify the agp accordingly and run the second round with splittling and rescaffolding to allow for fixing missasemblies:
```
module load bwa
module load samtools
bwa index purged_with_differentiated_region.fasta
bwa mem -M -t 50 purged_with_differentiated_region.fasta markers_1.fa markers_2.fa > aligned_paired.sam
```
```
module load python/2.7
python /nfs/scistore18/vicosgrp/melkrewi/project_save_the_genome_project/chromonomer/chromonomer-1.13/scripts/fasta2agp.py --fasta purged_with_differentiated_region.fasta > whole_genome.agp
```
Run chromonomer:
```
mkdir output_paired_no_split
/nfs/scistore18/vicosgrp/melkrewi/project_save_the_genome_project/chromonomer/chromonomer-1.13/chromonomer -p linkage_map_modified_no_XB1.tsv --out_path output_paired_no_split --alns aligned_paired.sam -a whole_genome.agp --fasta purged_with_differentiated_region.fasta --disable_splitting
```
Modify linkage map to prevent the splitting of S0 region:
```
cat aligned_paired.sam | cut -f1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16 > aligned_paired_no_XA.sam
grep LG6 aligned_paired_no_XA.sam | cut -f1 > diff_region_markers.txt
grep -wvf diff_region_markers.txt linkage_map_modified_no_XB1.tsv > new_linkage_map.tsv
#find the location and markers of S0 using the following command and add it to the new linkage map manually
awk '$4 == "LG6"' ./output_paired_no_split/CHRR_linkage_map.tsv | cut -f1,2,3 > S0_markers.txt
```
Run chromonomer with splitting and rescaffolding:
```
mkdir output_paired_split
/nfs/scistore18/vicosgrp/melkrewi/project_save_the_genome_project/chromonomer/chromonomer-1.13/chromonomer -p new_linkage_map.tsv --out_path output_paired_split --alns aligned_paired.sam -a whole_genome.agp --fasta purged_with_differentiated_region.fasta --rescaffold 
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

### getting the markers used in the analysis:
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
