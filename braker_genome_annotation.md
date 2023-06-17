# Alignment of RNA reads to the genome with hisat

### Let's load necessary modules

*module load STAR*

*module load hisat2*

*module load BRAKER/2.1.6*

*module load samtools*

## Actual alignment of RNA reads to the genome

### We include option --dta in order to report alignments tailored for transcript assemblers

`srun hisat2-build genome_clean_sorted.fasta genome_clean_sorted_hisat2`

```
for base in $(ls /nfs/scistore18/vicosgrp/vbett/Artemia_franEMdata/EMReads_analysis/Expression_dir/*_1.fastq | sed -r 's/_1.fastq//' | uniq)

do

srun hisat2 -q -x genome_clean_sorted_hisat2 -1 "${base}_1.fastq" -2 "${base}_2.fastq" -p 30 --dta --qc-filter -S "${base}_hisatr_alignment.sam"

samtools view -bS "${base}_hisatr_alignment.sam" > "${base}_hisatr_alignment.bam"

srun samtools sort -o "${base}_hisatr_alignment_sorted.bam" "${base}_hisatr_alignment.bam"

samtools index "${base}_hisatr_alignment_sorted.bam"

done
```


# We need to softmasked the repetitive elements in the genome

## First is to identify repeats de novo from your reference genome using RepeatModeler

*module load RepeatModeler*

*module load RepeatMasker*

```
srun BuildDatabase -name Afranciscana_genomeclean -engine ncbi genome_clean_sorted.fasta

srun RepeatModeler -threads 60 -database Afranciscana_genomeclean -engine ncbi

srun RepeatMasker -pa 60 -e ncbi -lib Afranciscana_genomeclean-families.fa -xsmall -dir Afran_maskallrepeats genome_clean_sorted.fasta
```

# We then run genome annotation using BRAKER 

This is detailed in this link (https://github.com/Gaius-Augustus/BRAKER)

*module load BRAKER/2.1.6*

`srun braker.pl --genome=genome_clean_sorted.fasta.masked --bam=60541_hisatr_alignment_sorted.bam,60542_hisatr_alignment_sorted.bam,60543_hisatr_alignment_sorted.bam,60544_hisatr_alignment_sorted.bam,60545_hisatr_alignment_sorted.bam,60546_hisatr_alignment_sorted.bam,60547_hisatr_alignment_sorted.bam,60548_hisatr_alignment_sorted.bam -gff3 --useexisting --species=Artemia_francisc --cores=30 --min_contig=5000 --softmasking --workingdir=/nfs/scistore18/vicosgrp/vbett/Artemia_franEMdata/EMReads_analysis/Expression_dir/braker_output`

# We can also run genome annotation with RNA bam and protein sequences and compare the results

## First we download all arthropoda protein sequences 

`wget https://v100.orthodb.org/download/odb10_arthropoda_fasta.tar.gz`

`tar xvfz odb10_arthropoda_fasta.tar.gz`

`cat arthropoda/Rawdata/* > proteins.fasta`

### We then generate a protein hint that will be used in BRAKER

Let's load necessary modules


```
module load java/8 bamtools/2.5.2 lpsolve/5.5.2.11 perl/5.34.0 genemark
module load GSL/2.7 htslib/1.13 augustus/3.4.0 python/3.9.7 BRAKER/2.1.6
module load ProtHint/2.6.0
```

`prothint.py genome_clean_sorted_masked.fasta proteins.fasta --workdir ProhintDir --threads 50`

This produces 3 output files 

- prothint.gff Gff file with all reported hints (introns, starts and stops)
- evidence.gff High confidence subset of prothint.gff

An output which is ready to be used in BRAKER and AUGUSTUS is also generated:
- prothint_augustus.gff

This is detailed in this link (https://github.com/gatech-genemark/ProtHint#protein-database-preparation)

## We then run braker using both RNA and protein hint

`srun braker.pl --genome=genome_clean_sorted.fasta.masked --bam=60541_hisatr_alignment_sorted.bam,60542_hisatr_alignment_sorted.bam,60543_hisatr_alignment_sorted.bam,60544_hisatr_alignment_sorted.bam,60545_hisatr_alignment_sorted.bam,60546_hisatr_alignment_sorted.bam,60547_hisatr_alignment_sorted.bam,60548_hisatr_alignment_sorted.bam -gff3 --useexisting --species=Artemia_francisca --etpmode --cores=30 --min_contig=5000 --softmasking --hints=prothint_augustus.gff --workingdir=/nfs/scistore18/vicosgrp/vbett/Artemia_franEMdata/EMReads_analysis/Expression_all/braker_outputetp`







