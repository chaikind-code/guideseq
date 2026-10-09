#!/bin/bash
# usage: drive.sh <manifest.yaml> <sample> <shared-control-sample-name>
set -euo pipefail
MF=$1; S=$2; SRC=$3
cd /home/user/work/an
mkdir -p aligned
# A shared untreated control only needs aligning once; point this sample's
# expected control alignment at the one already produced.
if [ -s "aligned/Control_$SRC.dedup.bam" ] && [ ! -e "aligned/Control_$S.dedup.bam" ]; then
  ln -s "Control_$SRC.dedup.bam" "aligned/Control_$S.dedup.bam"
  [ -s "aligned/Control_$SRC.dedup.bam.bai" ] && ln -s "Control_$SRC.dedup.bam.bai" "aligned/Control_$S.dedup.bam.bai" || true
fi
PYTHONPATH=/home/user/guideseq/guideseq /home/user/venv/bin/python \
  /home/user/guideseq/guideseq/guideseq.py main -m "$MF" --sample "$S" --step align+identify
