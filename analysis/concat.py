"""Concatenate per-chunk gzip members into one file per library, deleting each
chunk member as soon as it is appended so peak disk stays ~1x the output."""
import os, sys, shutil, glob, re
base = os.path.realpath(sys.argv[1])          # fq dir
bcfile = sys.argv[2]
assert base.startswith('/home/user/work/'), base
names = [l.split()[0] for l in open(bcfile) if l.strip() and not l.startswith('#')]
chunks = sorted(glob.glob(os.path.join(base, 'chunk_*')),
                key=lambda p: int(re.search(r'chunk_(\d+)$', p).group(1)))
print(f"{len(chunks)} chunks, {len(names)} libraries")
for name in names:
    for rd in ('r1', 'r2'):
        out = os.path.join(base, f"{name}.{rd}.fastq.gz")
        with open(out, 'wb') as o:
            for c in chunks:
                f = os.path.join(c, f"{name}.{rd}.fastq.gz")
                rp = os.path.realpath(f)
                if not rp.startswith(base + os.sep) or not rp.endswith('.fastq.gz'):
                    raise SystemExit(f"refusing to touch {rp}")
                if os.path.isfile(rp):
                    if os.path.getsize(rp):
                        with open(rp, 'rb') as i:
                            shutil.copyfileobj(i, o, 1 << 22)
                    os.remove(rp)
        print(f"built {out} {os.path.getsize(out)/1e9:.2f} GB", flush=True)
for c in chunks:
    rp = os.path.realpath(c)
    if rp.startswith(base + os.sep) and re.search(r'/chunk_\d+$', rp):
        leftover = os.listdir(rp)
        if leftover:
            print("leftover in", rp, leftover[:3])
        shutil.rmtree(rp)
print("done")
