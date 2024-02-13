for i in $(ls data/gffs); do name=$(basename $i .gff); echo $name;python workflow/scripts/gff2faa.py data/gffs/$i data/fastas/$name.fasta > data/faas/$name.faa; done
python workflow/scripts/extract_pks.py data/pks.gb > out/pks.faa
for i in $(ls data/faas); do name=$(basename $i .faa); blastp -query out/pks.faa -subject data/faas/$i -evalue 1E-4 -outfmt "6 qseqid sseqid pident qlen slen length nident evalue" > out/pks/$name.tsv; echo $name; done
