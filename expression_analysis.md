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
for base in $(ls /Expression/*_1.fastq | sed -r 's/_1.fastq//' | uniq)
do
STAR --runThreadN 50 --genomeDir augenome --readFilesIn "${base}_1.fastq" "${base}_2.fastq" --outFileNamePrefix "${base}_brakmasked" --twopassMode Basic --sjdbGTFfile braker_filconagat_isoORF3_sup100compgen.gtf --sjdbScore 2 --sjdbOverhang 124 --limitSjdbInsertNsj 1000000 --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.025 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMunmapped Within --outFilterType BySJout --outSAMattributes NH HI AS NM MD --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --quantTranscriptomeBan IndelSoftclipSingleend
done
rm -R *_STARgenome *_STARpass1 *_STARtmp
```

### We then check whether we will use strandedness (no, yes or reverse) during RSEM quantification

The first column of `ReadsPerGene.out.tab` represent gene name

The second column indicate RNA-seq reads count for unstranded library (`--stranded no`)

The third column represent reads count for first strand (`--stranded yes`)

The fourth column represent reads count for reverse strand (`--stranded reverse`)

```
head 60544_brakmaskedReadsPerGene.out.tab
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
# RSEM Quantification

Let's load necessary tools

```
module load anaconda3/2022.05 
source /mnt/nfs/clustersw/Debian/bullseye/anaconda3/2022.05/activate_anaconda3_2022.05.txt
conda activate rsem
module load bowtie2/2.4.5
```

### We build RSEM references
```
rsem-prepare-reference --gtf Rsem_brak/braker_filconagat_isoORF3_sup100compgen.gtf --bowtie2 Rsem_brak/asm_np_female_mkf02_01_09_2023_renamed_final_fin.fa Rsem_brak/asmnp_female_mkf02 
```

### We calculate Expression Values

```
for base in $(ls /Expression/*_1.fastq | sed -r 's/_1.fastq//' | uniq)
do
rsem-calculate-expression -p 40 --paired-end --alignments --estimate-rspd --ci-memory 100000 --no-bam-output --strandedness reverse "${base}_brakmaskedAligned.toTranscriptome.out.bam" Rsem_brak/asmnp_female_mkf02 "${base}_brakmasked" 
done
```

head 60541_brakmasked.genes.results
```
gene_id transcript_id(s)        length  effective_length        expected_count  TPM     FPKM
jg10    jg10.t1 2796.00 2407.31 0.00    0.00    0.00
jg10000 jg10000.t1      480.00  110.77  3.00    0.79    1.49
jg10001 jg10001.t1      615.00  230.50  16.00   2.03    3.82
```

We then normalize the expression for heads and gonads separately
```
setwd('Expression/NormalizedExpression')
mh1<-read.table("60541_brakmasked.genes.results", head=T, sep="\t")
mh2<-read.table("60542_brakmasked.genes.results", head=T, sep="\t")
fh1<-read.table("60545_brakmasked.genes.results", head=T, sep="\t")
fh2<-read.table("60546_brakmasked.genes.results", head=T, sep="\t")
all <- merge(mh1[, c('gene_id', 'TPM')],mh2[, c('gene_id', 'TPM')],by.x="gene_id", by.y="gene_id")
all <- merge(all,fh1[, c('gene_id', 'TPM')],by.x="gene_id", by.y="gene_id")
all <- merge(all,fh2[, c('gene_id', 'TPM')],by.x="gene_id", by.y="gene_id")
write.csv(all, file = "Expression_afranfhmh.txt",row.names = FALSE)
library(dplyr)
exp<-read.table("Expression_afranfhmh.txt", head=T, sep=",")
colnames(exp) <- c("gene_id","mh1","mh2","fh1","fh2")
expf<-exp[,-1]
head(expf)
rownames(expf)<-exp[,1]
par(mfrow=c(2,1))
par(mar=c(3,3,0,0))
pdf('Expression_visualizationfhmh_beforeNorm.pdf',width=7,height=7)
boxplot(log2(expf+1))
dev.off()
bolFMat<-as.matrix(expf, nrow = nrow(expf), ncol = ncol(expf))
library(NormalyzerDE)
expf2<-performQuantileNormalization(bolFMat, noLogTransform = T)
rownames(expf2)<-rownames(expf)
colnames(expf2)<-colnames(expf)
expf2<-as.data.frame(expf2)
pdf('Expression_visualizationNormfhmh.pdf',width=7,height=7)
boxplot(log2(expf2+1))
dev.off()
write.table(expf2, file = "Expression_afranciscana_headsnormalized.txt", quote=F)
```

References:
(https://ycl6.gitbook.io/guide-to-rna-seq-analysis/raw-read-processing/mapping/alignment-based-method)
