/* GUIDE-seq demultiplexer: reads `fastq-dump -Z --split-spot -I` stream
   (4 records per spot: I1(8), I2(16), R1, R2) from stdin.
   Mode "id":   emit interleaved FASTQ to stdout, read name = BC1+BC2 (for alignment-based barcode ID)
   Mode "split": write per-sample gzipped R1/R2 FASTQ, read name = <spotname>_<umi10>  (guideseq convention)
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#define MAXS 64
#define LB 4096

static char names[MAXS][128], b1[MAXS][16], b2[MAXS][16];
static gzFile o1[MAXS], o2[MAXS];
static long long cnt[MAXS];
static int ns = 0, mism = 1;
static int rev_i5 = 0;   /* 1 => sample i5 is revcomp(I2[8:16]) and UMI is I2[0:8] */

static inline int ham8(const char *a, const char *b){
    int d = 0;
    for (int i = 0; i < 8; i++) if (a[i] != b[i]) { if (++d > mism) return d; }
    return d;
}

static inline int assign(const char *s1, const char *s2){
    for (int i = 0; i < ns; i++)
        if (ham8(s1, b1[i]) <= mism && ham8(s2, b2[i]) <= mism) return i;
    return -1;
}

/* strip trailing newline; return length */
static inline int chomp(char *s){
    int l = strlen(s);
    while (l > 0 && (s[l-1] == '\n' || s[l-1] == '\r')) s[--l] = 0;
    return l;
}

int main(int argc, char **argv){
    if (argc < 3){ fprintf(stderr,"usage: demux <id|split> <barcodes.tsv|-> [outdir] [mismatch] [f|r]\n"); return 1; }
    int idmode = strcmp(argv[1], "id") == 0;
    const char *bcfile = argv[2];
    const char *pref = argc > 3 ? argv[3] : "out";
    if (argc > 4) mism = atoi(argv[4]);
    if (argc > 5) rev_i5 = (argv[5][0] == 'r');

    if (!idmode){
        FILE *f = fopen(bcfile, "r");
        if (!f){ perror("barcodes"); return 1; }
        char ln[512];
        while (fgets(ln, sizeof ln, f)){
            chomp(ln);
            if (!ln[0] || ln[0]=='#') continue;
            char nm[128], x1[32], x2[32];
            if (sscanf(ln, "%127s %31s %31s", nm, x1, x2) != 3) continue;
            strcpy(names[ns], nm); strncpy(b1[ns], x1, 15); strncpy(b2[ns], x2, 15);
            char p[512];
            snprintf(p, sizeof p, "%s/%s.r1.fastq.gz", pref, nm); o1[ns] = gzopen(p, "wb1");
            snprintf(p, sizeof p, "%s/%s.r2.fastq.gz", pref, nm); o2[ns] = gzopen(p, "wb1");
            if (!o1[ns] || !o2[ns]){ fprintf(stderr,"cannot open outputs for %s\n", nm); return 1; }
            ns++;
            if (ns >= MAXS){ fprintf(stderr,"too many samples\n"); return 1; }
        }
        fclose(f);
        fprintf(stderr, "loaded %d samples\n", ns);
    }

    char l[16][LB];
    long long spots = 0, kept = 0;
    char ob[1<<16];
    while (1){
        int got = 1;
        for (int i = 0; i < 16; i++) if (!fgets(l[i], LB, stdin)){ got = 0; break; }
        if (!got) break;
        spots++;
        char *i1s = l[1],  *i2s = l[5];
        char *r1s = l[9],  *r1q = l[11];
        char *r2s = l[13], *r2q = l[15];
        chomp(i1s); int l2 = chomp(i2s);
        chomp(r1s); chomp(r1q); chomp(r2s); chomp(r2q);
        if (strlen(i1s) < 8 || l2 < 8) continue;

        /* Index-read layout differs between the runs in this study:
           "f": I2 = i5 barcode (8) + molecular UMI (8)            [HiSeq-style runs]
           "r": I2 = molecular UMI (8) + revcomp(i5 barcode) (8)   [NextSeq-style runs]
           bc2 holds the forward-orientation i5 used for matching; umi_p the UMI. */
        char bc2[16];
        const char *umi_p;
        if (rev_i5 && l2 >= 16){
            static const char cmpl[128] = {['A']='T',['C']='G',['G']='C',['T']='A',['N']='N'};
            for (int k = 0; k < 8; k++){
                unsigned char ch = (unsigned char) i2s[15 - k];
                bc2[k] = (ch < 128 && cmpl[ch]) ? cmpl[ch] : 'N';
            }
            bc2[8] = 0;
            umi_p = i2s;            /* first 8 bases */
        } else {
            memcpy(bc2, i2s, 8); bc2[8] = 0;
            umi_p = i2s + (l2 - 8 > 0 ? l2 - 8 : 0);
        }

        /* spot name: "@SRRxxxx.N.1 INSTRUMENT:... length=8" -> take 2nd token (instrument id) */
        char *nm = l[0] + 1;
        char *sp = strchr(nm, ' ');
        char *inst = nm;
        if (sp){ *sp = 0; inst = sp + 1; char *sp2 = strchr(inst, ' '); if (sp2) *sp2 = 0; }
        chomp(inst);

        if (idmode){
            int n = snprintf(ob, sizeof ob,
                "@%.8s%.16s_%lld/1\n%s\n+\n%s\n@%.8s%.16s_%lld/2\n%s\n+\n%s\n",
                i1s, i2s, spots, r1s, r1q, i1s, i2s, spots, r2s, r2q);
            fwrite(ob, 1, n, stdout);
            kept++;
        } else {
            int s = assign(i1s, bc2);
            if (s < 0) continue;
            const char *umi = umi_p;
            cnt[s]++; kept++;
            int n = snprintf(ob, sizeof ob, "@%s_%.8s\n%s\n+\n%s\n", inst, umi, r1s, r1q);
            gzwrite(o1[s], ob, n);
            n = snprintf(ob, sizeof ob, "@%s_%.8s\n%s\n+\n%s\n", inst, umi, r2s, r2q);
            gzwrite(o2[s], ob, n);
        }
        if (spots % 20000000 == 0) fprintf(stderr, "..%lld spots, %lld kept\n", spots, kept);
    }
    for (int i = 0; i < ns; i++){ gzclose(o1[i]); gzclose(o2[i]); fprintf(stderr, "%s\t%lld\n", names[i], cnt[i]); }
    fprintf(stderr, "TOTAL_SPOTS\t%lld\nKEPT\t%lld\n", spots, kept);
    return 0;
}
