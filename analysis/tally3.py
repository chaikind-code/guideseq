import collections, glob, re, sys
cuts=collections.defaultdict(list)
for ln in open('regcuts.txt'):
    r,g,o=ln.split(); cuts[r].append((g,int(o)))
LAYOUT={'SRR11667141':'first','SRR11667144':'first','SRR11667142':'last','SRR11667143':'last'}
def rc(s): return s.translate(str.maketrans("ACGTN","TGCAN"))[::-1]
cigre=re.compile(r'(\d+)([MIDNSHP=X])')
def reflen(c): return sum(int(n) for n,o in cigre.findall(c) if o in 'MDN=X')
TOL=3
srr=sys.argv[1]; which=LAYOUT[srr]
tot=collections.Counter(); hit=collections.defaultdict(collections.Counter)
nspot=collections.Counter()
for f in glob.glob(f'L_{srr}/*.cnt'):
    for ln in open(f):
        p=ln.rstrip('\n').split('\t')
        if len(p)<5: continue
        bcfull,ref,pos,st,c=p
        if len(bcfull)<24: continue
        bc5 = rc(bcfull[16:24]) if which=='last' else bcfull[8:16]
        bc=f"{bcfull[:8]}+{bc5}"
        five=int(pos) if st=='+' else int(pos)+reflen(c)-1
        tot[bc]+=1
        best=None;bd=99
        for g,o in cuts.get(ref,[]):
            d=abs(five-o)
            if d<=TOL and d<bd: bd=d; best=g
        if best: hit[bc][best]+=1
print(f"=== {srr}: sample libraries identified by dominant on-target cut (±{TOL}bp) ===")
rows=sorted(((sum(h.values()),bc,h) for bc,h in hit.items()), reverse=True)
for s,bc,h in rows:
    if s<30: break
    top=h.most_common(1)[0]
    frac=100.0*top[1]/s
    print(f"{bc}  locireads={tot[bc]:7d}  cuts={s:6d}  -> {top[0]:16s} ({top[1]}, {frac:.0f}% of its cut reads)   others: " + ", ".join(f"{k}={v}" for k,v in h.most_common(4)[1:]))
