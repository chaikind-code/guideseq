import sys, os
# usage: mk_manifest.py <out.yaml> <fastqdir> <name:target:control> ...
out, fqdir = sys.argv[1], os.path.abspath(sys.argv[2])
hdr = f"""reference_genome: /home/user/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa
blacklist: /home/user/work/an/empty_blacklist.bed
genome: hg38
PAM: NGG
analysis_folder: {os.path.dirname(os.path.abspath(out))}
bwa: bwa
samtools: samtools
bedtools: bedtools
umi_tools: /home/user/venv/bin/umi_tools
window_size: 25
max_score: 7
mapq_threshold: 50
njobs: 4
save_pickle: False
samples:
"""
body=[]
for spec in sys.argv[3:]:
    name, target, ctrl = spec.split(':')
    body.append(f"""    {name}:
        target: {target}
        read1: {fqdir}/{name}.r1.fastq.gz
        read2: {fqdir}/{name}.r2.fastq.gz
        controlread1: {fqdir}/{ctrl}.r1.fastq.gz
        controlread2: {fqdir}/{ctrl}.r2.fastq.gz
        description: {name}
""")
open(out,'w').write(hdr+"".join(body))
print("wrote", out)
