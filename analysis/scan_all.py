FA="/home/user/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa"
guides={}
for ln in open('allguides.tsv'):
    n,s=ln.split(); guides[n]=s
comp=str.maketrans("ACGTNacgtn","TGCANtgcan")
def rc(s): return s.translate(comp)[::-1]
pats=[]
for n,s in guides.items():
    pats.append((n,"+",s)); pats.append((n,"-",rc(s)))
out=[]
name=None; chunks=[]
def process(name, seq):
    if name is None or not seq: return
    S=seq.upper()
    for n,strand,p in pats:
        i=S.find(p)
        while i!=-1:
            if strand=="+":
                if S[i+21:i+23]=="GG":
                    # 0-based site start; cut between protospacer pos17/18 -> 1-based cut coord
                    out.append((n,strand,name,i, i+17))
            else:
                if S[i-3:i-1]=="CC":
                    out.append((n,strand,name,i-3, i-3+6+1))
            i=S.find(p,i+1)
with open(FA) as f:
    for line in f:
        if line[0]=='>':
            process(name,"".join(chunks)); name=line[1:].split()[0]; chunks=[]
        else: chunks.append(line.strip())
process(name,"".join(chunks))
for n,st,c,s0,cut in sorted(out):
    print(f"{n}\t{st}\t{c}\t{s0}\t{cut}")
