# Expression analysis in *A. franciscana*
We use star to align and quantify transcript abundance

First we load necessary tools
```
module load STAR
module load samtools
```

Then we generate genome indexes

```
STAR --runThreadN 50 --runMode genomeGenerate --genomeDir augenome --genomeFastaFiles asm_np_female_mkf02_01_09_2023_renamed_final_fin.fa --sjdbGTFfile braker_filconagat_isoORF3_sup100compgen.gtf --sjdbOverhang 124 --genomeSAindexNbases 13
```

then we ran the alignment of RNA reads 

We included the options `--quantMode TranscriptomeSAM` to generate alignments that have translated transcript coordinates 

`--quantTranscriptomeBan IndelSoftclipSingleend` - this will ensure no indels or soft clips are included in the alignments 

`--quantMode GeneCounts` - this generate reads count per gene and can be used to determine the RNA library strandedness

```
for base in $(ls /nfs/scistore18/vicosgrp/vbett/Artemia_franhifiNewGenomeNew/Expression/*_1.fastq | sed -r 's/_1.fastq//' | uniq)
do
STAR --runThreadN 50 --genomeDir augenome --readFilesIn "${base}_1.fastq" "${base}_2.fastq" --outFileNamePrefix "${base}_brakmasked" --twopassMode Basic --sjdbGTFfile braker_filconagat_isoORF3_sup100compgen.gtf --sjdbScore 2 --sjdbOverhang 124 --limitSjdbInsertNsj 1000000 --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.025 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMunmapped Within --outFilterType BySJout --outSAMattributes NH HI AS NM MD --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --quantTranscriptomeBan IndelSoftclipSingleend
done
rm -R *_STARgenome *_STARpass1 *_STARtmp
```
We then check whether we will use strandedness (no, yes or reverse) during RSEM quantification

The first column of `ReadsPerGene.out.tab` represent gene name
The second column indicate RNA-seq reads count for unstranded library (`--stranded no`)
The third column represent reads count for first strand (`--stranded yes`)
The fourth column represent reads count for reverse strand (`--stranded reverse`)

This is detailed here (https://ycl6.gitbook.io/guide-to-rna-seq-analysis/raw-read-processing/mapping/alignment-based-method)

```
head(60544_brakmaskedReadsPerGene.out.tab)
N_unmapped      5789789 5789789 5789789
N_multimapping  13087023        13087023        13087023
N_noFeature     12424675        37432194        12575987
N_ambiguous     337703  1676    319465
jg591   3       0       3
jg592   3       0       3
jg593   3       0       3
jg594   1       0       1
jg595   0       0       0
jg596   1       0       1
jg597   0       0       0
jg598   0       0       0
jg599   2       2       0
jg600   0       0       0
jg601   1       1       0
jg603   2       2       0
jg604   0       0       0
jg605   3       0       3
jg606   3       1       2
jg607   1       0       1
jg608   3       0       3
jg609   25      3       22
jg610   545     5       540
jg611   1374    13      1361
jg612   245     4       241
jg613   182     23      159
jg614   335     3       332
jg615   10      0       10
jg616   173     1       172
jg617   1754    9       1745
```
