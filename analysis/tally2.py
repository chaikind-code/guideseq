import collections, glob, re, sys
OFF={'AAVS1reg':55105000,'CXCR4reg':136105000}
REGCHR={'AAVS1reg':'chr19','CXCR4reg':'chr2'}
sites=[]
for ln in open('allguide_sites.tsv'):
    n,st,c,s0,cut=ln.split()
    cut=int(cut)
    for reg,ch in REGCHR.items():
        if c==ch and OFF[reg] < cut < OFF[reg]+25000:
            sites.append((n,reg,cut,st))
LAYOUT={'SRR11667141':'first','SRR11667144':'first','SRR11667142':'last','SRR11667143':'last'}
def rc(s): return s.translate(str.maketrans("ACGTN","TGCAN"))[::-1]
cigre=re.compile(r'(\d+)([MIDNSHP=X])')
def reflen(c): return sum(int(n) for n,o in cigre.findall(c) if o in 'MDN=X')
TOL=3
srr=sys.argv[1]; which=LAYOUT[srr]
tot=collections.Counter(); hit=collections.defaultdict(collections.Counter)
for f in glob.glob(f'id_{srr}/*.cnt'):
    for ln in open(f):
        p=ln.rstrip('\n').split('\t')
        if len(p)<5: continue
        bcfull,ref,pos,st,c=p
        if len(bcfull)<24: continue
        bc5 = rc(bcfull[16:24]) if which=='last' else bcfull[8:16]
        bc=f"{bcfull[:8]}+{bc5}"
        five=int(pos) if st=='+' else int(pos)+reflen(c)-1
        g=OFF[ref]+five
        tot[bc]+=1
        for n,reg,cut,gst in sites:
            if reg==ref and abs(g-cut)<=TOL: hit[bc][n]+=1
print(f"=== {srr}  (±{TOL}bp of exact cut sites; i5 window={which}) ===")
rows=sorted(((sum(h.values()),bc,h) for bc,h in hit.items()), reverse=True)
for s,bc,h in rows:
    if s<40: break
    print(f"{bc}  regionreads={tot[bc]:7d}  total_at_cuts={s:6d}  " + "  ".join(f"{k}={v}" for k,v in h.most_common(5)))
