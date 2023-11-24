# De novo identification of repeat families using RepeatModeler  

`module load RepeatModeler`

`BuildDatabase -name asm_np_female_mkf02_01_09_2023_renamed_finalwscaff -engine ncbi asm_np_female_mkf02_01_09_2023_renamed_final_fin_wscaff.fa`

`RepeatModeler -threads 80 -database asm_np_female_mkf02_01_09_2023_renamed_finalwscaff -engine ncbi`

### Classification of unknown TEs from RepeatModeler using DeepTE

Retaining first names from sequences fasta families

`awk '/^>/ {$0=$1} 1' asm_np_female_mkf02_01_09_2023_renamed_finalwscaff-families.fa > asm_np_female_mkf02_01_09_2023_renamed_finalwscaffren-families.fa`

Getting unknown fasta headers

`grep -e "Unknown" asm_np_female_mkf02_01_09_2023_renamed_finalwscaffren-families.fa > asmfinalwscaff_unknown-families.txt`

getting fasta sequences for unknown TEs

sed -i 's/>//g' asmfinalwscaff_unknown-families.txt

seqtk subseq asm_np_female_mkf02_01_09_2023_renamed_finalwscaffren-families.fa asmfinalwscaff_unknown-families.txt > asm_np_female_mkf02_01_09_2023_renamed_finalwscaff_unknown-families.fa

module load python/3.7
module load hmmer

export TMPDIR=/nfs/scistore18/vicosgrp/vbett/Tools/DeepTE

srun python /nfs/scistore18/vicosgrp/vbett/Tools/DeepTE/DeepTE.py -i asm_np_female_mkf02_01_09_2023_renamed_finalwscaff_unknown-families.fa -sp M -m M -o /nfs/scistore18/vicosgrp/vbett/Artemia_franhifiNewGenomeNewWscaff/TEAnnotation/DeepTE_output -d /nfs/scistore18/vicosgrp/vbett/Artemia_franhifiNewGenomeNewWscaff/TEAnnotation/DeepTE_working

grep ">"  asm_np_female_mkf02_01_09_2023_renamed_finalwscaffren-families.fa > asmfinalwscaff_all-families.txt

grep -v "Unknown" asmfinalwscaff_all-families.txt > asmfinalwscaff_known-families.txt

sed -i 's/>//g' asmfinalwscaff_known-families.txt

seqtk subseq asm_np_female_mkf02_01_09_2023_renamed_finalwscaffren-families.fa asmfinalwscaff_known-families.txt > asm_np_female_mkf02_01_09_2023_renamed_finalwscaff_known-families.fa

We then combine known sequences with just annotated TEs (annotated with DeepTE) that were unknown 

cat DeepTE_output/opt_DeepTE.fasta asm_np_female_mkf02_01_09_2023_renamed_finalwscaff_known-families.fa > asm_np_female_mkf02_01_09_2023_renamed_finalwscaff_all-families.fa 


sed 's/ClassII_//g' asm_np_female_mkf02_01_09_2023_renamed_finalwscaff_all-families.fa > asm_np_female_mkf02_01_09_2023_renamed_finalwscaff_all-familiesfil.fa

sed -i 's/ClassIII_//g' asm_np_female_mkf02_01_09_2023_renamed_finalwscaff_all-familiesfil.fa

sed -i 's/ClassI_//g' asm_np_female_mkf02_01_09_2023_renamed_finalwscaff_all-familiesfil.fa

sed -i 's/ /\//g' asm_np_female_mkf02_01_09_2023_renamed_finalwscaff_all-familiesfil.fa

sed -i 's/DNA_/DNA /g' asm_np_female_mkf02_01_09_2023_renamed_finalwscaff_all-familiesfil.fa

sed -i 's/LTR_/LTR /g' asm_np_female_mkf02_01_09_2023_renamed_finalwscaff_all-familiesfil.fa

sed -i 's/nLTR_/nLTR /g' asm_np_female_mkf02_01_09_2023_renamed_finalwscaff_all-familiesfil.fa

sed -i 's/ /\//g' asm_np_female_mkf02_01_09_2023_renamed_finalwscaff_all-familiesfil.fa
