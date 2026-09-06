#!/usr/bin/env python3
"""
08_FOSC_phylogeny.py
====================
Multi-locus ML phylogeny placing MRL8996 within the Fusarium oxysporum species
complex, using the local genome collection (no downloads) and the FOSC probe
sequences already assembled for the F012 project.

Panel: MRL8996 plus one representative genome per taxon listed in PANEL below -
Fol4287, the endophyte Fo47, the clinical/human-associated isolate NRRL 32931,
a spread of plant-pathogenic formae speciales, and F. verticillioides as
outgroup. Representatives are resolved by organism name from
fusarium_genomes/genomes_scaffold_up/genomeID_mapping.tsv, so the accession of
every tip is recorded in the metadata table rather than hard-coded here.

Pipeline: blastn of each locus probe set against each genome -> best hit region
extracted with flanks -> per-locus MAFFT alignment -> trimAl gappyout ->
concatenation with a partition file -> IQ-TREE partitioned ML with 1000
ultrafast bootstrap replicates -> figure.

Requires blastn, mafft, trimal and iqtree. All four are present in the
chr10-phylo conda environment on this machine:

    conda run -n chr10-phylo python3 08_FOSC_phylogeny.py --threads 12

Runtime: BLAST databases dominate (~30-60 s per genome, cached in
Revision_2026/results/phylogeny/blastdb), so allow ~20-30 min for a first run
and a couple of minutes for re-runs.

Outputs (all under Revision_2026/):
    results/phylogeny/loci/<locus>.aln.fasta      per-locus alignments
    results/phylogeny/concatenated.fasta          concatenated supermatrix
    results/phylogeny/partitions.nex              locus partitions
    results/phylogeny/FOSC_ML.treefile            IQ-TREE ML tree
    tables/Table_FOSC_phylogeny_taxa.tsv          tip -> accession/organism
    figures/Figure_FOSC_phylogeny.pdf / .png      600 dpi

Usage:
    conda run -n chr10-phylo python3 08_FOSC_phylogeny.py --threads 12
    conda run -n chr10-phylo python3 08_FOSC_phylogeny.py --plot-only
"""

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

ADODDI = Path("/home/usuario/nvme_data/home/usuario/adoddi")
REV = ADODDI / "MRL8996_genome_paper" / "Revision_2026"
GENOMES = ADODDI / "fusarium_genomes" / "formatted_genomes"
MAPPING = ADODDI / "fusarium_genomes" / "genomes_scaffold_up" / "genomeID_mapping.tsv"
PROBES = ADODDI / "F012_HiC_data" / "final_results" / "Phylogeny_files"
MRL_FASTA = ADODDI / "MRL8996_HiC_data" / "MRL8996_chromosome_level.fasta"

OUT = REV / "results" / "phylogeny"
LOCI = ["tef1", "rpb1", "rpb2", "tub2", "cam", "his3", "acl1", "its"]

# tip label -> exact organism string to match in the mapping table
PANEL = {
    "Fol4287": "Fusarium oxysporum f. sp. lycopersici 4287",
    "Fol_MN25": "Fusarium oxysporum f. sp. lycopersici MN25",
    "Fo47_endophyte": "Fusarium oxysporum Fo47",
    "NRRL32931_clinical": "Fusarium oxysporum NRRL 32931",
    "Fo_melonis": "Fusarium oxysporum f. sp. melonis",
    "Fo_fragariae": "Fusarium oxysporum f. sp. fragariae",
    "Fo_niveum": "Fusarium oxysporum f. sp. niveum",
    "Fo_cubense": "Fusarium oxysporum f. sp. cubense",
    "Fo_cucumerinum": "Fusarium oxysporum f. sp. cucumerinum",
    "Fo_lini": "Fusarium oxysporum f. sp. lini",
    "Fo_vasinfectum": "Fusarium oxysporum f. sp. vasinfectum",
    "Fo_conglutinans": "Fusarium oxysporum f. sp. conglutinans",
    "Fo_radicis_lycopersici": "Fusarium oxysporum f. sp. radicis-lycopersici",
    "Fo_cepae": "Fusarium oxysporum f. sp. cepae",
    "Fo_albedinis": "Fusarium oxysporum f. sp. albedinis",
    "Fo_raphani": "Fusarium oxysporum f. sp. raphani 54005",
    "F_verticillioides_outgroup": "Fusarium cf. verticillioides",
}
FLANK = 0
MIN_ALN_FRAC = 0.5      # discard blast hits shorter than half the probe


def need(tool):
    p = shutil.which(tool)
    if p is None:
        sys.exit("%s not on PATH - run inside the chr10-phylo environment "
                 "(conda run -n chr10-phylo ...)" % tool)
    return p


def read_fasta(path):
    name, seq, out = None, [], {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    out[name] = "".join(seq)
                name, seq = line[1:].split()[0], []
            elif line:
                seq.append(line)
    if name:
        out[name] = "".join(seq)
    return out


def resolve_panel():
    m = pd.read_csv(MAPPING, sep="\t")
    rows = [{"tip": "MRL8996", "genomeID": "MRL8996", "accession": "this study",
             "organism": "Fusarium oxysporum MRL8996 (keratitis isolate)",
             "fasta": str(MRL_FASTA)}]
    for tip, organism in PANEL.items():
        sub = m[m["organism"] == organism]
        chosen = None
        for _, r in sub.iterrows():
            f = GENOMES / ("%s_genomic.fna" % r["genomeID"])
            if f.exists():
                chosen = (r, f)
                break
        if chosen is None:
            print("WARNING: no local genome for %s (%s)" % (tip, organism),
                  file=sys.stderr)
            continue
        r, f = chosen
        rows.append({"tip": tip, "genomeID": r["genomeID"],
                     "accession": r["accession"], "organism": r["organism"],
                     "fasta": str(f)})
    taxa = pd.DataFrame(rows)
    (REV / "tables").mkdir(parents=True, exist_ok=True)
    taxa.drop(columns=["fasta"]).to_csv(
        REV / "tables" / "Table_FOSC_phylogeny_taxa.tsv", sep="\t", index=False)
    print("panel: %d tips" % len(taxa))
    return taxa


def blastdb(fasta, tip):
    db = OUT / "blastdb" / tip
    db.parent.mkdir(parents=True, exist_ok=True)
    if not (db.with_suffix(".nsq").exists() or Path(str(db) + ".nsq").exists()):
        subprocess.run(["makeblastdb", "-in", str(fasta), "-dbtype", "nucl",
                        "-out", str(db)], check=True,
                       stdout=subprocess.DEVNULL)
    return db


def extract_locus(db, probe_file, genome_seqs):
    """Best blastn hit of any probe; returns the subject region as a string."""
    r = subprocess.run(
        ["blastn", "-query", str(probe_file), "-db", str(db),
         "-outfmt", "6 qseqid sseqid pident length sstart send bitscore qlen",
         "-max_target_seqs", "5", "-evalue", "1e-20", "-num_threads", "2"],
        capture_output=True, text=True, check=True)
    best = None
    for line in r.stdout.splitlines():
        f = line.split("\t")
        length, qlen, bits = int(f[3]), int(f[7]), float(f[6])
        if length < MIN_ALN_FRAC * qlen:
            continue
        if best is None or bits > best[0]:
            best = (bits, f[1], int(f[4]), int(f[5]))
    if best is None:
        return None
    _, sid, s0, s1 = best
    seq = genome_seqs.get(sid)
    if seq is None:
        return None
    lo, hi = sorted((s0, s1))
    frag = seq[max(0, lo - 1 - FLANK):hi + FLANK]
    if s1 < s0:                      # minus strand
        comp = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
        frag = frag.translate(comp)[::-1]
    return frag


def build_alignments(taxa, threads):
    (OUT / "loci").mkdir(parents=True, exist_ok=True)
    genome_cache = {}
    for locus in LOCI:
        probe = PROBES / ("probes_%s.fasta" % locus)
        if not probe.exists():
            print("skipping %s (no probe file)" % locus, file=sys.stderr)
            continue
        raw = OUT / "loci" / ("%s.fasta" % locus)
        with open(raw, "w") as out:
            for _, t in taxa.iterrows():
                if t["tip"] not in genome_cache:
                    genome_cache[t["tip"]] = read_fasta(t["fasta"])
                db = blastdb(t["fasta"], t["tip"])
                frag = extract_locus(db, probe, genome_cache[t["tip"]])
                if frag is None:
                    print("  %s: no %s hit" % (t["tip"], locus), file=sys.stderr)
                    continue
                out.write(">%s\n%s\n" % (t["tip"], frag))
        aln = OUT / "loci" / ("%s.aln.fasta" % locus)
        with open(aln, "w") as fh:
            subprocess.run(["mafft", "--auto", "--thread", str(threads), str(raw)],
                           stdout=fh, stderr=subprocess.DEVNULL, check=True)
        trimmed = OUT / "loci" / ("%s.trim.fasta" % locus)
        subprocess.run(["trimal", "-in", str(aln), "-out", str(trimmed),
                        "-gappyout"], check=True, stdout=subprocess.DEVNULL)
        print("  %s aligned" % locus)
        genome_cache = {k: v for k, v in genome_cache.items()}   # keep cache


def concatenate(taxa):
    blocks, parts, tips = {}, [], list(taxa["tip"])
    pos = 1
    for locus in LOCI:
        f = OUT / "loci" / ("%s.trim.fasta" % locus)
        if not f.exists():
            continue
        seqs = read_fasta(f)
        L = len(next(iter(seqs.values())))
        for tip in tips:
            blocks.setdefault(tip, []).append(seqs.get(tip, "-" * L))
        parts.append((locus, pos, pos + L - 1))
        pos += L
    cat = OUT / "concatenated.fasta"
    with open(cat, "w") as fh:
        for tip in tips:
            fh.write(">%s\n%s\n" % (tip, "".join(blocks[tip])))
    nex = OUT / "partitions.nex"
    with open(nex, "w") as fh:
        fh.write("#nexus\nbegin sets;\n")
        for name, a, b in parts:
            fh.write("  charset %s = %d-%d;\n" % (name, a, b))
        fh.write("end;\n")
    print("supermatrix: %d tips x %d sites, %d partitions"
          % (len(tips), pos - 1, len(parts)))
    return cat, nex


def run_iqtree(cat, nex, threads):
    exe = next((shutil.which(x) for x in ("iqtree2", "iqtree3", "iqtree")
                if shutil.which(x)), None)
    if exe is None:
        sys.exit("iqtree not on PATH - run inside chr10-phylo")
    subprocess.run([exe, "-s", str(cat), "-p", str(nex), "-m", "MFP",
                    "-B", "1000", "-T", str(threads), "--prefix",
                    str(OUT / "FOSC_ML"), "-o", "F_verticillioides_outgroup",
                    "-redo"], check=True)


def plot_tree():
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from Bio import Phylo

    tree = Phylo.read(str(OUT / "FOSC_ML.treefile"), "newick")
    tree.root_with_outgroup({"name": "F_verticillioides_outgroup"})
    tree.ladderize()

    mpl.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                         "pdf.fonttype": 42, "ps.fonttype": 42})
    fig, ax = plt.subplots(figsize=(6.8, 6.2))
    highlight = {"MRL8996": "#b2182b", "Fol4287": "#1f5c99",
                 "Fo47_endophyte": "#2e7d32",
                 "NRRL32931_clinical": "#b2182b"}

    def label(clade):
        if not clade.name:
            return None
        return clade.name.replace("_", " ")

    Phylo.draw(tree, axes=ax, do_show=False, label_func=label,
               show_confidence=False,
               label_colors=lambda n: highlight.get(n.replace(" ", "_"), "black")
               if n else "black")
    # UFBoot support (IQ-TREE writes it as the internal-node label) printed for
    # every internal node resolved at >= 70%; weaker nodes are left unlabelled.
    depths = tree.depths()
    ypos = {}
    for i, t in enumerate(tree.get_terminals(), 1):
        ypos[t] = i

    def y_of(clade):
        if clade in ypos:
            return ypos[clade]
        v = sum(y_of(c) for c in clade.clades) / len(clade.clades)
        ypos[clade] = v
        return v

    n_sup = 0
    for clade in tree.get_nonterminals():
        sup = clade.confidence
        if sup is None:
            continue
        if sup <= 1:
            sup *= 100
        if sup >= 70:
            ax.text(depths.get(clade, 0), y_of(clade) - 0.22, "%d" % round(sup),
                    fontsize=6, ha="right", va="center", color="0.35")
            n_sup += 1
    print("support labels drawn (UFBoot >= 70): %d" % n_sup)
    ax.set_xlabel("substitutions per site")
    ax.set_ylabel("")
    ax.set_title("Multi-locus ML placement of MRL8996 in the F. oxysporum "
                 "species complex\n(%d loci, IQ-TREE partitioned model, 1000 "
                 "ultrafast bootstrap replicates)" % len(LOCI),
                 fontsize=9, loc="left")
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.get_yaxis().set_visible(False)
    fig.tight_layout()
    (REV / "figures").mkdir(parents=True, exist_ok=True)
    fig.savefig(REV / "figures" / "Figure_FOSC_phylogeny.pdf", bbox_inches="tight")
    fig.savefig(REV / "figures" / "Figure_FOSC_phylogeny.png", dpi=600,
                bbox_inches="tight")
    print("figure written")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--plot-only", action="store_true")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    if not a.plot_only:
        for t in ("makeblastdb", "blastn", "mafft", "trimal"):
            need(t)
        taxa = resolve_panel()
        build_alignments(taxa, a.threads)
        cat, nex = concatenate(taxa)
        run_iqtree(cat, nex, a.threads)
    plot_tree()


if __name__ == "__main__":
    main()
