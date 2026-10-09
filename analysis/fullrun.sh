#!/bin/bash
# Full-depth demultiplex of selected libraries from one SRA run, streamed from S3.
# usage: fullrun.sh <SRR> <f|r> <barcodes.tsv> <outdir> <nspots> <nchunks>
set -euo pipefail
SRR=$1; LAY=$2; BC=$3; OUT=$4; TOT=$5; NC=$6
U="https://sra-pub-run-odp.s3.amazonaws.com/sra/$SRR/$SRR"
PER=$(( (TOT + NC - 1) / NC ))
mkdir -p "$OUT"
work(){
  i=$1; U=$2; PER=$3; OUT=$4; BC=$5; LAY=$6
  s=$((i*PER+1)); e=$(((i+1)*PER))
  d="$OUT/chunk_$i"; mkdir -p "$d"
  fastq-dump -Z --split-spot -I -N $s -X $e "$U" 2>/dev/null \
    | /home/user/work/demux split "$BC" "$d" 1 "$LAY" 2> "$d/demux.log"
}
export -f work
seq 0 $((NC-1)) | xargs -P 6 -I{} bash -c 'work {} '"\"$U\" $PER \"$OUT\" \"$BC\" $LAY"
# concatenate per-chunk gzip members (valid gzip concatenation)
while read -r name b1 b2; do
  [ -z "${name:-}" ] && continue
  for rd in r1 r2; do
    cat "$OUT"/chunk_*/"$name.$rd.fastq.gz" > "$OUT/$name.$rd.fastq.gz"
  done
done < "$BC"
rm -rf "$OUT"/chunk_*
echo "=== read pairs per library ==="
for f in "$OUT"/*.r1.fastq.gz; do
  n=$(( $(zcat "$f" | wc -l) / 4 ))
  echo -e "$(basename "$f" .r1.fastq.gz)\t$n"
done
