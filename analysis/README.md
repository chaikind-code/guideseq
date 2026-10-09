# GUIDE-seq re-analysis of SRX8228066–SRX8228069

Goal: rank a panel of candidate guide RNAs by on-target dsODN integration versus
off-target activity, using the four GUIDE-seq sequencing runs deposited under
BioProject PRJNA625995 (study SRP258158).

## Data access

NCBI web/eutils/FTP hosts are not reachable from this environment, so run
accessions were recovered from the SRA metadata snapshot published in the AWS
Open Data bucket `s3://sra-pub-metadata-us-east-1/sra/metadata_json/`, and the
runs themselves are read **in place** over HTTPS from
`s3://sra-pub-run-odp/sra/<SRR>/<SRR>` (199 GB total, far more than local disk).
`sra-tools` streams them with range requests once remote accession resolution is
disabled:

    /repository/remote/disabled = "true"

| Experiment  | Run         | Size   | Spots |
|-------------|-------------|--------|-------|
| SRX8228066  | SRR11667144 | 64 GB  | 349,788,426 |
| SRX8228067  | SRR11667143 | 33 GB  | 273,282,111 |
| SRX8228068  | SRR11667142 | 51 GB  | 389,762,591 |
| SRX8228069  | SRR11667141 | 51 GB  | 263,193,390 |

Each run is an **undemultiplexed** NextSeq/HiSeq run with four reads per spot:
I1 (8 nt, i7), I2 (16 nt), R1 and R2 (147–151 nt).

## The missing index→guide mapping, and the work-around

The deposited metadata lists which guides were pooled in a run, but not which
index pair belongs to which guide, so `guideseq demultiplex` cannot be driven
from a manifest. The sample barcodes are nevertheless present in the data, so
they are recovered and assigned empirically:

1. Tally I1/I2 across each run to find the real sample barcodes. The index-read
   layout is **not** the same in every run:
   * SRR11667141, SRR11667144: `I2 = i5(8) + UMI(8)`
   * SRR11667142, SRR11667143: `I2 = UMI(8) + revcomp(i5)(8)`
   In each run the barcode window resolves to TruSeq D7xx/D5xx indices with a
   sharp abundance break separating real samples from sequencing-error barcodes.
2. Locate every guide in the panel in GRCh38 by exact protospacer+NGG match
   (`scan_all.py`) and compute its Cas9 blunt cut coordinate, 3 bp from the PAM.
3. Align a subsample of each run to a 24 kb mini-reference of those loci
   (`idrun2.sh`) and, per barcode, attribute each dsODN integration position to
   the nearest cut site within 3 bp (`tally3.py`). A library's guide is the cut
   site that dominates it.

Tight tolerances matter: the AAVS1 guides are packed into ~170 bp
(cut sites 55115755 / 55115761 / 55115773 / 55115774 / 55115782 / 55115828),
so a loose window misattributes one guide's on-target reads to its neighbours.
AAVS1_site_2 and AAVS1_site_13 cut 1 bp apart and are not separable by position
alone.

## Files

| File | Purpose |
|------|---------|
| `demux.c` | Streaming demultiplexer for `fastq-dump --split-spot` output; handles both index layouts, appends the UMI to the read name as guideseq expects |
| `scan_all.py` | Exact genome-wide protospacer+NGG search; emits site start and cut coordinate |
| `allguides.tsv` | Guide names and 20 nt spacers from the supplementary table |
| `allguide_sites.tsv` | Resulting on-target coordinates and strands |
| `idrun.sh`, `idrun2.sh` | Parallel range-streamed alignment of a run against a mini-reference |
| `tally3.py` | Per-barcode assignment of libraries to guides |

## Reference

GRCh38 plus the pre-built BWA index from
`s3://1000genomes/technical/reference/GRCh38_reference_genome/`. The
supplementary-table coordinates are GRCh38, so they are directly comparable.

## Caveats

* The ENCODE hg38 blacklist is not reachable from this environment, so blacklist
  filtering is skipped. It is skipped identically for every guide, and
  sequence-match and control filtering are retained.

## Which mate carries the dsODN tag

`identifyOfftargetSites` inspects only the mate with `flag & 128` (read 2) and
treats its 5' end as the dsODN integration point. That is correct for these
runs, verified by aligning each mate separately and histogramming 5' ends for
CXCR4_site_3 (predicted cut site chr2:136,115,797):

| mate | top 5' positions |
|------|------------------|
| R1 | 136115634 (18), 136115914 (15), 136128216 (12) — diffuse |
| R2 | **136115798 (641), 136115797 (428)**, 136115796 (56) — sharp |

A FASTQ grep is misleading here: the dsODN primer appears at the *start* of 29%
of R1 reads but those are adapter-dimer / ODN-only fragments that do not align
and are dropped at MAPQ >= 50.

## Primer annotation caveat

Tag reads at the cut site have CIGAR `22S129M` with the soft-clipped prefix
`ACATATGACAACTCAATTAAAC`, which is `dsODN_primer_revcomp[12:]`: these runs used
a shorter nested ODN primer than `default.yaml` describes, so the read begins 12
nt inside the configured primer. `match_dsODN` anchors at the read start, so it
reports `nomatch` for most reads and the primer columns in the output are
uninformative.

This does not affect site calling. `addPositionBarcode` counts every read into
the strand totals irrespective of primer class, and a window is retained when
`barcode_geometric_mean > 0`, i.e. it has reads on both strands — the standard
GUIDE-seq bidirectional criterion. The behaviour is identical for every guide
compared here.

These libraries also lack the chr2:99,357,763-99,357,801 spike-in amplicon that
`control_primer_coord` points at, so `control_counts` falls back to its -1
sentinel and the `normlization_ratio` column is meaningless. It is not used for
any conclusion drawn here.

## Independent cross-check

`quantify.py` recomputes the same quantities straight from the deduplicated BAM
— tag-read 5' ends, clustered within the window, with per-cluster plus/minus
support and the best gapless protospacer+PAM match in the surrounding sequence.
It shares no code with the pipeline, so agreement between the two is meaningful.
