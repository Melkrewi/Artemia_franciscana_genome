# dNdS

### finding the homologs for the W-genes:

```
#extract W-genes from codingseq file
perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' Artemia_fran_W_names.txt braker_filconagatwgenes_isoORF3_sup100compgen.codingseq > afran_W.fasta

#removing the w-genes from the codingseq file
perl -ne 'if(/^>(\S+)/){$c=!$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' Artemia_fran_W_names.txt braker_filconagatwgenes_isoORF3_sup100compgen.codingseq > afran_cds_without_w.fasta

#Reciprocal besthits
module load blast
blastn -query afran_cds_without_w.fasta -subject afran_W.fasta -out Z_vs_W.out -outfmt 6
blastn -query afran_W.fasta -subject afran_cds_without_w.fasta -out W_vs_Z.out -outfmt 6
perl Reciprocal_BLAST.pl --blast_output_1=W_vs_Z.out --blast_output_2=Z_vs_W.out --output_file=reciprocal_besthits_ZW_afran.txt
sort reciprocal_besthits_ZW_afran.txt | uniq > reciprocal_besthits_ZW_afran_uniq.txt
cat reciprocal_besthits_ZW_afran_uniq.txt | cut -f1 | uniq > Artemia_fran_Z_names.txt

perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' Artemia_fran_Z_names.txt braker_filconagatwgenes_isoORF3_sup100compgen.codingseq > afran_Z.fasta
```
# dNdS analysis
```
module load wise/2.4.1

module load blat
module load seqtk/20210915
mkdir ORF2
for (( i = 1; i <= 300; i++ ))

do

sed -n "$i"p  reciprocal_besthits_ZW_afran_uniq.txt | cut -f 1 | sort | uniq | seqtk subseq afran_Z.fasta /dev/stdin > ORF2/"$i"_z.fa

sed -n "$i"p  reciprocal_besthits_ZW_afran_uniq.txt | cut -f 2 | sort | uniq | seqtk subseq afran_W.fasta /dev/stdin > ORF2/"$i"_w.fa

cat ORF2/"$i"_w.fa ORF2/"$i"_z.fa > ORF2/"$i"_wz.fa

done

module load muscle

export PATH=~/Gblocks_0.91b/:$PATH

for file in ./ORF2/*wz.fa

do

perl ~/translatorx_vLocal.pl -i ${file} -o ${file}.fa.tx -t T -g "-s"

done

for file in ./ORF2/*.fa.tx.nt_cleanali.fasta

do

perl ~/parseFastaIntoAXT.pl ${file}
~/KaKs_Calculator2.0/src/KaKs_Calculator -i ${file}.axt -o ${file}.kaks -m NG -m YN
done
cd ORF2
for (( i = 1; i <= 300; i++ ))
do
cat "$i"_*.kaks | cut -f 1,2,3,4,5,7 >> tableofKaKs_wz.txt
done
```
