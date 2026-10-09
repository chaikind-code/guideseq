#!/bin/bash
SRR=$1; TOT=$2; NC=$3
U="https://sra-pub-run-odp.s3.amazonaws.com/sra/$SRR/$SRR"
PER=$((TOT/NC))
mkdir -p L_$SRR && rm -f L_$SRR/*.cnt
work(){
  i=$1; SRR=$2; U=$3; PER=$4
  s=$((i*PER+1)); e=$(((i+1)*PER))
  fastq-dump -Z --split-spot -I -N $s -X $e "$U" 2>/dev/null \
   | /home/user/work/demux id - \
   | bwa mem -p -t 1 -v 1 /home/user/work/loci.fa - 2>/dev/null \
   | samtools view -F 2308 -q 30 - \
   | awk -F'\t' '{split($1,a,"_"); st = (int($2/16)%2) ? "-" : "+"; print a[1]"\t"$3"\t"$4"\t"st"\t"$6}' \
   > L_$SRR/$i.cnt
}
export -f work
seq 0 $((NC-1)) | xargs -P 4 -I{} bash -c 'work {} '"$SRR $U $PER"
cat L_$SRR/*.cnt | wc -l
