"""Control-subtract and summarise a GUIDE-seq library.

Takes quantify.py output for a treated library and its untreated (dsODN-only)
control, pairs clusters by position, and reports the on-target plus the
off-target sites that survive: guide match within MAXMM mismatches, reads on
both strands, and enrichment over the control.
"""
import sys, argparse, collections

ap = argparse.ArgumentParser()
ap.add_argument('sample'); ap.add_argument('control')
ap.add_argument('--on-chrom', required=True); ap.add_argument('--on-pos', type=int, required=True)
ap.add_argument('--maxmm', type=int, default=6)
ap.add_argument('--min-reads', type=int, default=10)
ap.add_argument('--min-enrich', type=float, default=5.0)
ap.add_argument('--tol', type=int, default=25)
ap.add_argument('--label', default='')
a = ap.parse_args()

def load(p):
    rows = []
    with open(p) as f:
        next(f)
        for line in f:
            r = line.rstrip('\n').split('\t')
            rows.append((r[3], int(r[4]), int(r[0]), int(r[1]), int(r[2]), int(r[5])))
    return rows

smp = load(a.sample)
ctl = load(a.control)
cidx = collections.defaultdict(list)
for ch, pos, n, pl, mi, mm in ctl:
    cidx[(ch, pos // 1000)].append((pos, n))
def ctl_reads(ch, pos):
    tot = 0
    for b in (pos // 1000 - 1, pos // 1000, pos // 1000 + 1):
        for p, n in cidx.get((ch, b), ()):
            if abs(p - pos) <= a.tol:
                tot += n
    return tot

smp_total = sum(r[2] for r in smp)
on = None; offs = []
for ch, pos, n, pl, mi, mm in smp:
    c = ctl_reads(ch, pos)
    enr = n / max(c, 1)
    is_on = (ch == a.on_chrom and abs(pos - a.on_pos) <= a.tol)
    rec = (n, pl, mi, ch, pos, mm, c, enr)
    if is_on:
        if on is None or n > on[0]:
            on = rec
    elif mm <= a.maxmm and pl > 0 and mi > 0 and n >= a.min_reads and enr >= a.min_enrich:
        offs.append(rec)
offs.sort(reverse=True)
L = a.label or a.sample
onr = on[0] if on else 0
offtot = sum(o[0] for o in offs)
print(f"===== {L} =====")
print(f"tag reads in library            : {smp_total}")
if on:
    print(f"ON-TARGET {on[3]}:{on[4]}  reads={on[0]} (+{on[1]}/-{on[2]})  mm={on[5]}  control={on[6]}  enrichment={on[7]:.0f}x")
else:
    print("ON-TARGET: NOT DETECTED")
print(f"off-target sites (mm<={a.maxmm}, bidirectional, >={a.min_reads} reads, >={a.min_enrich}x over control): {len(offs)}")
print(f"off-target reads total         : {offtot}")
if onr:
    print(f"on-target / off-target reads   : {onr/max(offtot,1):.3f}")
    print(f"on-target as % of on+off reads : {100.0*onr/(onr+offtot):.1f}%")
print("top off-targets:")
print("  reads\tplus\tminus\tchrom\tpos\tmm\tctrl\tenrich")
for o in offs[:12]:
    print(f"  {o[0]}\t{o[1]}\t{o[2]}\t{o[3]}\t{o[4]}\t{o[5]}\t{o[6]}\t{o[7]:.0f}x")
print()
