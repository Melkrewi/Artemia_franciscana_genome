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

export TMPDIR=/home/DeepTE

srun python DeepTE.py -i asm_np_female_mkf02_01_09_2023_renamed_finalwscaff_unknown-families.fa -sp M -m M -o DeepTE_output -d DeepTE_working

Get DNA transposons sequences

seqtk subseq asm_np_female_mkf02_01_09_2023_renamed_finalwscaff_known-families.fa DNA_transposons.txt > DNA_transposons.fa

Get Class II from output of DeepTE

grep -e "ClassII" opt_DeepTE.fasta > opt_DeepTE_classII.txt

seqtk subseq opt_DeepTE.fasta opt_DeepTE_classII.txt > opt_DeepTE_classII.fa

Combine with DNA transposons

cat DNA_transposons.fa opt_DeepTE_classII.fa > DNA_transposons_DeepTE_classIIrep.fa

Run DeepTE

export TMPDIR=/home/DeepTE
srun python DeepTE.py -i DNA_transposons_DeepTE_classIIrep.fa -sp M -m M  -fam ClassII -o DeepTE_transposons 

cat opt_DeepTE_classIIfilrep.fasta asm_np_femaleNoDNA_transposons.fa opt_DeepTE_classI.fasta > asm_np_female_TEclass_super-families.fa
