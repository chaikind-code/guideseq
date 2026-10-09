"""Independent GUIDE-seq quantification straight from a deduplicated BAM.

For the tag-bearing mate (read 2, whose 5' end marks the dsODN integration
point) collect 5' positions, cluster them, and for every cluster report
bidirectional read support plus the best gapless match of the guide protospacer
in the surrounding window. Deliberately simple and independent of the pipeline
so the two can be compared.
"""
import collections, subprocess, sys, re, argparse

cig = re.compile(r'(\d+)([MIDNSHP=X])')
def reflen(c):
    return sum(int(n) for n, o in cig.findall(c) if o in 'MDN=X')

COMP = str.maketrans("ACGTN", "TGCAN")
def rc(s): return s.translate(COMP)[::-1]

def mismatches(guide, window):
    """Best (n_mismatch, position, strand) for a gapless guide+PAM match."""
    best = (99, None, None)
    L = len(guide)
    for strand in '+-':
        g = guide if strand == '+' else rc(guide)
        for i in range(0, len(window) - L + 1):
            sub = window[i:i+L]
            if strand == '+':
                if sub[-2:] != 'GG':
                    continue
            else:
                if sub[:2] != 'CC':
                    continue
            n = 0
            for a, b in zip(g, sub):
                if a == 'N' or b == 'N':
                    continue
                if a != b:
                    n += 1
                    if n >= best[0]:
                        break
            if n < best[0]:
                best = (n, i, strand)
    return best

ap = argparse.ArgumentParser()
ap.add_argument('bam'); ap.add_argument('guide')   # guide incl. NGG, e.g. GAAG...NGG
ap.add_argument('--ref', default='/home/user/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa')
ap.add_argument('--mapq', type=int, default=50)
ap.add_argument('--window', type=int, default=25)
ap.add_argument('--out', default=None)
a = ap.parse_args()

pos = collections.defaultdict(lambda: [0, 0])     # (chrom,pos) -> [plus, minus]
p = subprocess.Popen(['samtools', 'view', '-q', str(a.mapq), '-f', '128', '-F', '2308', a.bam],
                     stdout=subprocess.PIPE, universal_newlines=True)
n = 0
for line in p.stdout:
    f = line.split('\t', 9)
    flag = int(f[1]); chrom = f[2]
    if chrom == '*':
        continue
    rev = flag & 16
    five = int(f[3]) + reflen(f[5]) - 1 if rev else int(f[3])
    pos[(chrom, five)][1 if rev else 0] += 1
    n += 1
p.stdout.close()
if p.wait() != 0:
    raise SystemExit('samtools view failed')
sys.stderr.write(f"tag reads used: {n}, distinct 5' positions: {len(pos)}\n")

# cluster positions within `window` bp on the same chromosome
clusters = []
for (chrom, q) in sorted(pos):
    pl, mi = pos[(chrom, q)]
    if clusters and clusters[-1][0] == chrom and q - clusters[-1][2] <= a.window:
        c = clusters[-1]
        c[2] = q; c[3] += pl; c[4] += mi
        if pl + mi > c[6]:
            c[5] = q; c[6] = pl + mi
    else:
        clusters.append([chrom, q, q, pl, mi, q, pl + mi])
sys.stderr.write(f"clusters: {len(clusters)}\n")

import pyfaidx
fa = pyfaidx.Fasta(a.ref, as_raw=True, sequence_always_upper=True)
guide = a.guide
rows = []
for chrom, lo, hi, pl, mi, peak, peakn in clusters:
    if pl + mi < 2:
        continue
    s = max(0, peak - a.window - len(guide))
    e = peak + a.window + len(guide)
    try:
        win = str(fa[chrom][s:e])
    except Exception:
        continue
    nm, off, strand = mismatches(guide, win)
    rows.append((pl + mi, pl, mi, chrom, peak, nm, strand, (pl*mi) ** 0.5))
rows.sort(reverse=True)
out = open(a.out, 'w') if a.out else sys.stdout
print("reads\tplus\tminus\tchrom\tpeak_pos\tmismatches_vs_guide\tstrand\tbidirectional_geomean", file=out)
for r in rows:
    print("\t".join(str(x) for x in r), file=out)
if a.out:
    out.close()
