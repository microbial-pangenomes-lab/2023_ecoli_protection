cat /dev/null > out/phylogroup.txt
mkdir -p out/phylogroups
for i in $(cat out/mash_input.txt);
do
  sid=$(basename $i .fasta);
  mkdir -p out/phylogroups/$sid;
  makeblastdb -in $i -input_type fasta -out out/phylogroups/$sid/db -dbtype nucl;
  blastn -query ClermonTyping/data/primers.fasta -perc_identity 90 -task blastn -outfmt 5 -db out/phylogroups/$sid/db -out out/phylogroups/$sid/blast.xml;
  pg=$(python3 ClermonTyping/bin/clermont.py -x out/phylogroups/$sid/blast.xml | awk -F '\t' '{print $NF}');
  echo -e $sid"\t"$pg;
  echo -e $sid"\t"$pg >> out/phylogroup.txt;
done
sed -i 's/Non Escherichia/Non_Escherichia/g' out/phylogroup.txt

rm -rf out/phylogroups
