#!/usr/bin/env python3
"""
27_identity_matrices_nucmer.py
==============================
Recomputes Data Files 2-4 (chromosome-by-chromosome coverage and identity of
MRL8996 against Fol4287, of Fo47 against MRL8996, and of MRL8996 against
itself) from NUCmer alignments.

Why this replaces the earlier matrices
--------------------------------------
The earlier matrices were computed from minimap2 PAF files as
sum(column 10) / sum(column 11). Those PAF files were written without
base-level alignment (no cg: tag), so column 10 is only the number of bases
covered by chained minimizer anchors, not the number of identical bases, and
column 11 includes every gap inside the chain. The ratio is therefore an
under-estimate of nucleotide identity that depends on repeat and indel
density (values of 5-80 % for chromosome pairs that are > 98 % identical).
The earlier script also renamed the MRL8996 scaffolds by sorting their
PGA_scaffold names, which swaps Chr03 and Chr04 relative to the official
naming (Chr03 = PGA_scaffold_4, Chr04 = PGA_scaffold_3).

Method
------
* Sequences are renamed to the official chromosome names first
  (MRL8996: mapping_cromosomi_NCBI.tsv; Fol4287: NC_030986-NC_031000 =
  chromosomes 1-15; Fo47: NC_072840-NC_072851 = chromosomes I-XII).
* NUCmer (MUMmer 3.1, default parameters), one query chromosome at a time.
  For the self matrix each chromosome is aligned against the other 15, so the
  trivial self-alignment is excluded.
* Identity = length-weighted mean %IDY of the alignments of the pair, the
  AvgIdentity statistic of dnadiff. It is computed on the 1-to-1 alignment set
  (delta-filter -1), as in dnadiff's "1-to-1" report.
* Coverage = fraction of the query chromosome covered by the merged query
  intervals of the many-to-many alignment set (delta-filter -m), i.e. every
  homologous copy is counted, as in the earlier coverage matrix.
* Identity is left empty where the 1-to-1 set has less than MIN_ALN bp for
  the pair, because a weighted mean over a few short repeat hits is not a
  chromosome-level identity.

out: Revision_2026/tables/identity_nucmer/Matrix_*.csv, per-pair long tables,
     genome_wide_avg_identity.tsv
"""
import glob
import os
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

REV = Path("/home/usuario/nvme_data/home/usuario/adoddi/MRL8996_genome_paper/Revision_2026")
WORK = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("nuc")   # holds out/*.delta
OUT = REV / "tables" / "identity_nucmer"
MIN_ALN = 10_000


def coords(delta, mode):
    """delta-filter mode ('-1' or '-m') | show-coords; returns a DataFrame."""
    filt = subprocess.run(["delta-filter", mode, str(delta)], capture_output=True,
                          text=True, check=True).stdout
    tmp = str(delta) + mode.replace("-", ".f")
    with open(tmp, "w") as fh:
        fh.write(filt)
    sc = subprocess.run(["show-coords", "-rclTH", tmp], capture_output=True,
                        text=True, check=True).stdout
    os.remove(tmp)
    rows = [l.split("\t") for l in sc.strip().split("\n") if l]
    cols = ["S1", "E1", "S2", "E2", "LEN1", "LEN2", "IDY", "LENR", "LENQ",
            "COVR", "COVQ", "REF", "QRY"]
    df = pd.DataFrame(rows, columns=cols)
    for c in cols[:11]:
        df[c] = pd.to_numeric(df[c])
    return df


def merged_len(iv):
    iv = sorted((min(a, b), max(a, b)) for a, b in iv)
    tot, cs, ce = 0, None, None
    for s, e in iv:
        if cs is None or s > ce:
            if cs is not None:
                tot += ce - cs + 1
            cs, ce = s, e
        else:
            ce = max(ce, e)
    if cs is not None:
        tot += ce - cs + 1
    return tot


def pairs(prefix):
    recs = []
    for d in sorted(glob.glob(str(WORK / "out" / (prefix + "_*.delta")))):
        one, many = coords(d, "-1"), coords(d, "-m")
        if many.empty:
            continue
        qlen = int(many.LENQ.iloc[0])
        for (r, q), g in many.groupby(["REF", "QRY"]):
            g1 = one[(one.REF == r) & (one.QRY == q)]
            aln1 = int(g1.LEN2.sum())
            idy = (g1.IDY * g1.LEN2).sum() / aln1 if aln1 else float("nan")
            recs.append(dict(query=q.split("_", 1)[1], target=r.split("_", 1)[1],
                             query_len=qlen,
                             coverage_pct=100 * merged_len(zip(g.S2, g.E2)) / qlen,
                             aligned_bp_1to1=aln1, n_aln_1to1=len(g1),
                             identity_pct=idy))
    return pd.DataFrame(recs)


def matrices(df, qnames, tnames, self_diag=False):
    cov = df.pivot(index="query", columns="target", values="coverage_pct")
    idt = df.assign(identity_pct=df.identity_pct.where(df.aligned_bp_1to1 >= MIN_ALN)) \
            .pivot(index="query", columns="target", values="identity_pct")
    cov = cov.reindex(index=qnames, columns=tnames).fillna(0.0).round(2)
    idt = idt.reindex(index=qnames, columns=tnames).round(2)
    if self_diag:
        for c in qnames:
            cov.loc[c, c] = 100.0
            idt.loc[c, c] = 100.0
    return cov, idt


def genome_wide(prefix):
    tot_l = tot_w = 0
    for d in glob.glob(str(WORK / "out" / (prefix + "_*.delta"))):
        g = coords(d, "-1")
        tot_l += g.LEN2.sum()
        tot_w += (g.IDY * g.LEN2).sum()
    return tot_w / tot_l, tot_l


def write(cov, idt, path, title):
    with open(path, "w") as fh:
        fh.write(title + "\n\n--- COVERAGE MATRIX (%) ---\n")
        cov.to_csv(fh)
        fh.write("\n\n--- IDENTITY MATRIX (%) ---\n")
        idt.to_csv(fh)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    mrl = ["Chr%02d" % i for i in range(1, 17)]
    fol = ["Chr%02d" % i for i in range(1, 16)]
    fo47 = ["Chr%02d" % i for i in range(1, 13)]
    gw = []
    for prefix, q, t, fn, title, diag in [
        ("MRLvsFOL", mrl, fol, "Matrix_MRL8996_vs_Fol4287_nucmer.csv",
         "--- QUERY MRL8996 (rows) vs TARGET Fol4287 (columns); NUCmer ---", False),
        ("FO47vsMRL", fo47, mrl, "Matrix_Fo47_vs_MRL8996_nucmer.csv",
         "--- QUERY Fo47 (rows) vs TARGET MRL8996 (columns); NUCmer ---", False),
        ("SELF", mrl, mrl, "Matrix_MRL8996_self_nucmer.csv",
         "--- MRL8996 vs MRL8996, each chromosome against the other 15; NUCmer ---", True)]:
        df = pairs(prefix)
        df.to_csv(OUT / fn.replace(".csv", "_pairs.tsv"), sep="\t", index=False)
        cov, idt = matrices(df, q, t, diag)
        write(cov, idt, OUT / fn, title)
        a, l = genome_wide(prefix)
        gw.append(dict(comparison=prefix, avg_identity_1to1_pct=round(a, 2),
                       aligned_bp_1to1=int(l)))
        print(prefix, "pairs", len(df), "genome-wide 1-to-1 AvgIdentity %.2f %% over %d bp" % (a, l))
    pd.DataFrame(gw).to_csv(OUT / "genome_wide_avg_identity.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
