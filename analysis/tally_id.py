import collections, glob, re, sys
OFF={'AAVS1reg':55105000,'CXCR4reg':136105000}
cuts={'AAVS1_site_10':('AAVS1reg',55116279),'AAVS1_site_3':('AAVS1reg',55115828),
      'AAVS1_site_13':('AAVS1reg',55115774),'AAVS1_site_14':('AAVS1reg',55115755),
      'CXCR4_site_3':('CXCR4reg',136115797)}
LAYOUT={'SRR11667141':'first','SRR11667144':'first','SRR11667142':'last','SRR11667143':'last'}
def rc(s): return s.translate(str.maketrans("ACGTN","TGCAN"))[::-1]
cig=re.compile(r'(\d+)([MIDNSHP=X])')
def reflen(c):
    return sum(int(n) for n,o in cig.findall(c) if o in 'MDN=X')

srr=sys.argv[1]
which=LAYOUT[srr]
tot=collections.Counter(); hit=collections.defaultdict(collections.Counter)
pos_detail=collections.defaultdict(collections.Counter)
for f in glob.glob(f'id_{srr}/*.cnt'):
    for ln in open(f):
        p=ln.rstrip('\n').split('\t')
        if len(p)<5: continue
        bcfull,ref,pos,st,c=p
        if len(bcfull)<24: continue
        i7=bcfull[:8]; i2=bcfull[8:24]
        bc5 = rc(i2[8:16]) if which=='last' else i2[0:8]
        bc=f"{i7}+{bc5}"
        pos=int(pos)
        five = pos if st=='+' else pos+reflen(c)-1
        g=OFF[ref]+five
        tot[bc]+=1
        for gn,(r,cu) in cuts.items():
            if r==ref and abs(g-cu)<=25:
                hit[bc][gn]+=1
                pos_detail[gn][g]+=1

print(f"=== {srr}: barcodes with on-target integration signal (±25bp of cut site) ===")
rows=[]
for bc,h in hit.items():
    rows.append((sum(h.values()),bc,h,tot[bc]))
rows.sort(reverse=True)
for s,bc,h,t in rows[:25]:
    print(f"{bc}  locusreads={t:7d}  " + "  ".join(f"{k}={v}" for k,v in h.most_common()))
print()
for gn in cuts:
    print(f"--- {gn}: top integration positions (all barcodes): " + ", ".join(f"{p}:{c}" for p,c in pos_detail[gn].most_common(6)))
