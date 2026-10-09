import sys
FA="/home/user/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa"
guides={
 "AAVS1_site_10":"GGGAACCCAGCGAGTGAAGA",
 "AAVS1_site_13":"GTCCCCTCCACCCCACAGTG",
 "AAVS1_site_14":"GGGGCCACTAGGGACAGGAT",
 "AAVS1_site_3" :"GAGCCACATTAACCGGCCCT",
 "CXCR4_site_3" :"GAAGATGATGGAGTAGATGG",
}
comp=str.maketrans("ACGTNacgtn","TGCANtgcan")
def rc(s): return s.translate(comp)[::-1]
pats=[]
for n,s in guides.items():
    pats.append((n,"+",s))
    pats.append((n,"-",rc(s)))

hits=[]
name=None; chunks=[]
def process(name, seq):
    if name is None: return
    S=seq.upper()
    for n,strand,p in pats:
        i=S.find(p)
        while i!=-1:
            if strand=="+":
                pam=S[i+21:i+23]
                ok = pam=="GG"
                pstart,pend=i,i+23
            else:
                pam=S[i-3:i-1]
                ok = pam=="CC"
                pstart,pend=i-3,i+20
            if ok:
                hits.append((n,strand,name,pstart,pend,S[max(0,pstart-40):pend+40]))
            i=S.find(p,i+1)

with open(FA) as f:
    for line in f:
        if line[0]=='>':
            process(name,"".join(chunks))
            name=line[1:].split()[0]; chunks=[]
        else:
            chunks.append(line.strip())
process(name,"".join(chunks))

for h in hits:
    print("\t".join(map(str,h)))
