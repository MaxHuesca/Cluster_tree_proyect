#!/bin/bash
set -eou pipefail

#this program has the objective of create a tabular file after obtaining the local alignments usign blatp program 
#Arguments  
#$1 = file with the fasta sequences 
#$2 = number of sequences you want to use from the fasta 
#$3 (optional) = out dir 

#asign the variables 
input_fasta=$1
number_seqs=$2 
if [[ $# < 3 ]]; then 
    out_dir="../results"
    if [[ ! -d "$out_dir" ]]; then 
        mkdir "$out_dir"
    fi
else 
    out_dir=$3
fi

out_baseline=${input_fasta##*/} 
out_baseline=${out_baseline%.*}
data_baseline=${input_fasta%/*}

#we clean the headers of the fasta if it is neccesary 
fasta_clean="${data_baseline}/${out_baseline}_clean.faa"
awk '/^>/{sub(/\|.*/,"",$0)}1' "$input_fasta" > "$fasta_clean"

#now we are going to cut the fasta in order to obtain onyu 150 sequences 
fasta_cuted="${data_baseline}/${out_baseline}_cuted.faa"
awk -v num="$number_seqs" 'BEGIN{count_seq=0} /^>/{count_seq++}count_seq<=num{print $0}' "$fasta_clean" > "$fasta_cuted"

#one time we have filtered the data now we can make the blast procesing 
makeblastdb -in "$fasta_cuted" -dbtype prot -parse_seqids 
#now we can do the blast all for all  
out_blast="${out_dir}/${out_baseline}_blastp.out"
blastp -query "$fasta_cuted" -db "$fasta_cuted" -outfmt 7 -max_hsps 1 -use_sw_tback > "$out_blast"

# make the tsv 
out_tsv="${out_dir}/${out_baseline}_aling.tsv"
awk 'BEGIN{FS=OFS="\t"} !/^#/ {print $1, $2, $12}' "$out_blast" > "$out_tsv" 

#finally we delete all the blast db files 
rm -rf "${data_baseline}/*.*.*"
