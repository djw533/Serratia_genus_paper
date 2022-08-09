#mkdir gbk_files

mkdir only_fasta
mkdir only_annotation

for file in *.gff
do
  	strain=$( echo ${file} | rev | cut -d '.' -f 2- | rev )
        head -n $( expr $(grep -n '##FASTA' ${file} | cut -d ':' -f 1) - 1) ${file} >  only_annotation/${strain}.gff
        tail -n +$( expr $(grep -n '##FASTA' ${file} | cut -d ':' -f 1) + 1) ${file} >  only_fasta/${strain}.fasta
done
