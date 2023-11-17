# Alignment of RNA reads to the genome with STAR

### Let's load necessary modules

*module load STAR*

*module load BRAKER/2.1.6*

*module load samtools*

## Actual alignment of RNA reads to the genome

`STAR --runMode genomeGenerate --genomeDir STARgenome --runThreadN 50 --genomeFastaFiles asm_np_female_mkf02_01_09_2023_renamed_final_fin.fa`

```
for base in $(ls /Annotation/*_1.fastq | sed -r 's/_1.fastq//' | uniq)

do
srun STAR --twopassMode Basic --genomeDir STARgenome --runThreadN 50 --readFilesIn "${base}_1.fastq" "${base}_2.fastq" --outFileNamePrefix "${base}_RNAalignment_" --outSAMtype BAM Unsorted

samtools sort -@ 50 -T ./ -o "${base}_RNAalignment_Aligned.sorted.bam" "${base}_RNAalignment_Aligned.out.bam"

samtools index -@ 50 "${base}_RNAalignment_Aligned.sorted.bam"

done
```


# We need to softmasked the repetitive elements in the genome

## First is to identify repeats de novo from your reference genome using RepeatModeler

*module load RepeatModeler*

*module load RepeatMasker*

```
BuildDatabase -name asm_np_female_mkf02_01_09_2023_renamed_final -engine ncbi asm_np_female_mkf02_01_09_2023_renamed_final_fin.fa

RepeatModeler -threads 100 -database asm_np_female_mkf02_01_09_2023_renamed_final -engine ncbi

```
## We then annotate the repeats using the library databases created above
```
RepeatMasker -pa 60 -lib asm_np_female_mkf02_01_09_2023_renamed_final-families.fa -xsmall asm_np_female_mkf02_01_09_2023_renamed_final_fin.fa
```
# We can also run genome annotation with RNA bam and protein sequences and compare the results

## First we download all arthropoda protein sequences 

`wget https://v100.orthodb.org/download/odb10_arthropoda_fasta.tar.gz`

`tar xvfz odb10_arthropoda_fasta.tar.gz`

`cat arthropoda/Rawdata/* > proteins.fasta`

### We then generate a protein hint that will be used in BRAKER

Let's load necessary modules


```
module load java bamtools/2.5.2 lpsolve perl genemark
module load GSL htslib augustus/3.4.0 python/3.9.7 BRAKER/2.1.6
module load ProtHint
```

`prothint.py asm_np_female_mkf02_01_09_2023_renamed_final_fin.fa.masked proteins.fasta --workdir ProhintDirmasked --threads 50`

This produces 3 output files 

- prothint.gff Gff file with all reported hints (introns, starts and stops)
- evidence.gff High confidence subset of prothint.gff

An output which is ready to be used in BRAKER and AUGUSTUS is also generated:
- prothint_augustus.gff

This is detailed in this link (https://github.com/gatech-genemark/ProtHint#protein-database-preparation)

## We then run braker using both RNA bam files and protein hint

`braker.pl --species=artfrannewhifiremasked --genome=asm_np_female_mkf02_01_09_2023_renamed_final_fin.fa.masked --bam=60541_RNAalignment_Aligned.sorted.bam,60542_RNAalignment_Aligned.sorted.bam,60543_RNAalignment_Aligned.sorted.bam,60544_RNAalignment_Aligned.sorted.bam,60545_RNAalignment_Aligned.sorted.bam,60546_RNAalignment_Aligned.sorted.bam,60547_RNAalignment_Aligned.sorted.bam,60548_RNAalignment_Aligned.sorted.bam -gff3 --useexisting --etpmode --hints=/nfs/scistore18/vicosgrp/vbett/Artemia_franhifiNewGenomeNew/Annotation/ProhintDirmasked/prothint_augustus.gff --cores=30 --softmasking --workingdir=/nfs/scistore18/vicosgrp/vbett/Artemia_franhifiNewGenomeNew/Annotation/braker_staralign`

### It is recommended to run braker at least two times in order to reduce possibility of false positives. Therefore we will run above script but we will include the hint output of the braker

`braker.pl --species=artfranewhifiremaskedrn2 --genome=asm_np_female_mkf02_01_09_2023_renamed_final_fin.fa.masked --bam=60541_RNAalignment_Aligned.sorted.bam,60542_RNAalignment_Aligned.sorted.bam,60543_RNAalignment_Aligned.sorted.bam,60544_RNAalignment_Aligned.sorted.bam,60545_RNAalignment_Aligned.sorted.bam,60546_RNAalignment_Aligned.sorted.bam,60547_RNAalignment_Aligned.sorted.bam,60548_RNAalignment_Aligned.sorted.bam -gff3 --useexisting --etpmode --hints=/nfs/scistore18/vicosgrp/vbett/Artemia_franhifiNewGenomeNew/Annotation/ProhintDirmasked/prothint_augustus.gff --cores=30 --softmasking --hints=/nfs/scistore18/vicosgrp/vbett/Artemia_franhifiNewGenomeNew/Annotation/braker_staralign/hintsfile.gff --workingdir=/nfs/scistore18/vicosgrp/vbett/Artemia_franhifiNewGenomeNew/Annotation/braker_staralignrnd2`

Braker produces many output files: 
- augustus.hints.aa
- bam_header.map
- errors
- genome_header.map
- augustus.hints.codingseq
- braker.gff3
- GeneMark-ETP
- hintsfile.gff
- augustus.hints.gff3       
- braker.gtf
- genemark_evidence.gff
- species
- augustus.hints.gtf
- braker.log
- genemark_hintsfile.gff
-  what-to-cite.txt

```
sed 's/file_1_file_1_//g' braker.gff3 > braker_fil.gff3
sed -i 's/file_1_file_2_/gm/g' braker_fil.gff3
sed -i 's/file_2_/gm/g' braker_fil.gff3
agat_convert_sp_gxf2gxf.pl -g braker_fil.gff3 -o braker_filconagat.gff3
agat_sp_keep_longest_isoform.pl --gff braker_filconagat.gff3 -o braker_filconagat_iso.gff3
agat_sp_filter_by_ORF_size.pl --gff braker_filconagat_iso.gff3 -o braker_filconagat_isoORF.gff3

agat_sp_filter_incomplete_gene_coding_models.pl --gff braker_filconagat_isoORF3_sup100.gff --fasta ../asm_np_female_mkf02_01_09_2023_renamed_final_fin.fa.masked -o braker_filconagat_isoORF3_sup100compgen.gff
```
