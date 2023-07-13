# Coverage Analysis

```
module load java

module load bowtie2/2.4.4

bowtie2-build asm.np_female_mkf0.2.fa Afran_new_genome

bowtie2 -x Afran_new_genome -1 CC2U_6_1.fastq.gz -2 CC2U_6_2.fastq.gz --end-to-end --sensitive -p 100 -S CC2U_6.sam

bowtie2 -x Afran_new_genome -1 CC2U_7_1.fastq.gz -2 CC2U_7_2.fastq.gz --end-to-end --sensitive -p 40 -S CC2U_7.sam

grep -vw "XS:i" CC2U_6.sam > CC2U_6_unique.sam

grep -vw "XS:i" CC2U_7.sam > CC2U_7_unique.sam

module load soap/coverage

soap.coverage -sam -cvg -i CC2U_6_unique.sam -onlyuniq -p 100 -refsingle asm.np_female_mkf0.2.fa -window CC2U_6_window 10000

soap.coverage -sam -cvg -i CC2U_7_unique.sam -onlyuniq -p 8 -refsingle CHRR_integrated.fa -window CC2U_7_window 10000

```
