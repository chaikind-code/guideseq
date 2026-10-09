import collections, glob, re, sys
cuts=collections.defaultdict(list)
for ln in open('regcuts.txt'):
    r,g,o=ln.split(); cuts[r].append((g,int(o)))
LAYOUT={'SRR11667141':'first','SRR11667144':'first','SRR11667142':'last','SRR11667143':'last'}
def rc(s): return s.translate(str.maketrans("ACGTN","TGCAN"))[::-1]
cigre=re.compile(r'(\d+)([MIDNSHP=X])')
def reflen(c): return sum(int(n) for n,o in cigre.findall(c) if o in 'MDN=X')
srr=sys.argv[1]; which=LAYOUT[srr]; want=set(sys.argv[2:])
tot=collections.Counter(); hit=collections.defaultdict(collections.Counter)
for f in glob.glob(f'L_{srr}/*.cnt'):
    for ln in open(f):
        p=ln.rstrip('\n').split('\t')
        if len(p)<5: continue
        bcfull,ref,pos,st,c=p
        if len(bcfull)<24: continue
        bc5 = rc(bcfull[16:24]) if which=='last' else bcfull[8:16]
        bc=f"{bcfull[:8]}+{bc5}"
        if want and bc not in want: continue
        five=int(pos) if st=='+' else int(pos)+reflen(c)-1
        tot[bc]+=1
        best=None;bd=99
        for g,o in cuts.get(ref,[]):
            d=abs(five-o)
            if d<=3 and d<bd: bd=d; best=g
        if best: hit[bc][best]+=1
for bc in sys.argv[2:]:
    h=hit[bc]
    print(f"{bc}  locireads={tot[bc]:6d}  cutreads={sum(h.values()):5d}  " + ", ".join(f"{k}={v}" for k,v in h.most_common(5)))
