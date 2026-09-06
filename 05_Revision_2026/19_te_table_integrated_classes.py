#!/usr/bin/env python3
"""
Table 2, rebuilt on the integrated EDTA + DeepTE classification (R2.20).

The previous version of the table used the EDTA/RepeatMasker label alone,
which left 2.33 % of the core and 15.13 % of the accessory compartment as
"unknown", and split single superfamilies across several near-duplicate
labels (DNA/DTA vs TIR/hAT vs MITE/DTA; unknown vs Unknown; SINE vs SINE?).

This script reproduces the cross-reference the DeepTE run applied
(09_TE_annotation/export_final_integrated_stats.py): every RepeatMasker
interval labelled "unknown" is relabelled with the DeepTE Fungi-model
class of the DeepTE interval it overlaps.  The resulting labels from both
tools are then mapped onto one hierarchy (class -> superfamily), intervals
are merged per compartment so nothing is double-counted, and the table is
written with the compartment sizes and the totals in the table itself.

Inputs (all pre-existing, nothing is re-run):
    09_TE_annotation/.../MRL8996_clean_for_EDTA.fasta.mod.out   (RepeatMasker)
    09_TE_annotation/DeepTE/DeepTE-master/output_dir/opt_DeepTE.txt
    08_circos_features/chrom.sizes

Output:
    Revision_2026/tables/Table_TE_classes_core_vs_accessory.tsv
"""
from collections import defaultdict
from pathlib import Path

import pandas as pd

ADODDI = Path("/home/usuario/nvme_data/home/usuario/adoddi")
HIC = ADODDI / "MRL8996_HiC_data"
REV = ADODDI / "MRL8996_genome_paper" / "Revision_2026"
TE_DIR = HIC / "09_TE_annotation"
RM_OUT = (TE_DIR / "MRL8996_clean_for_EDTA.fasta.mod.EDTA.anno"
          / "MRL8996_clean_for_EDTA.fasta.mod.out")
DEEPTE = TE_DIR / "DeepTE/DeepTE-master/output_dir/opt_DeepTE.txt"
SIZES = HIC / "08_circos_features" / "chrom.sizes"
ACCESSORY = {"Chr12", "Chr13", "Chr14", "Chr15", "Chr16"}

# label (EDTA/RepeatMasker or DeepTE) -> (class, superfamily)
RETRO = "Class I (retrotransposons)"
DNA = "Class II (DNA transposons)"
HEL = "Class III (Helitrons)"
OTH = "Other / unclassified"
MAP = {
    # --- Class I, LTR
    "LTR/Gypsy/Ty3": (RETRO, "LTR / Gypsy"),
    "LTR/Gypsy": (RETRO, "LTR / Gypsy"),
    "ClassI_LTR_Gypsy": (RETRO, "LTR / Gypsy"),
    "LINE/Gypsy/Ty3": (RETRO, "LTR / Gypsy"),
    "LTR/Copia": (RETRO, "LTR / Copia"),
    "ClassI_LTR_Copia": (RETRO, "LTR / Copia"),
    "LTR/unknown": (RETRO, "LTR, superfamily undet."),
    "ClassI_LTR": (RETRO, "LTR, superfamily undet."),
    # --- Class I, non-LTR
    "LINE/Tad1": (RETRO, "LINE / Tad1"),
    "LINE/CRE": (RETRO, "LINE / CRE"),
    "ClassI_nLTR_LINE_L1": (RETRO, "LINE / L1"),
    "ClassI_nLTR_DIRS": (RETRO, "DIRS"),
    "SINE/Foxy": (RETRO, "SINE / Foxy"),
    "ClassI_nLTR_SINE_tRNA": (RETRO, "SINE / tRNA-derived"),
    "ClassI_nLTR_SINE_7SL": (RETRO, "SINE / 7SL-derived"),
    "ClassI_nLTR_SINE": (RETRO, "SINE, family undet."),
    "SINE?": (RETRO, "SINE, family undet."),
    "ClassI_nLTR": (RETRO, "Non-LTR, order undet."),
    "ClassI": (RETRO, "Retro, order undet."),
    # --- Class II, TIR superfamilies
    "TIR/hAT": (DNA, "TIR / hAT"),
    "DNA/DTA": (DNA, "TIR / hAT"),
    "MITE/DTA": (DNA, "TIR / hAT"),
    "DNA/hAT-Restless": (DNA, "TIR / hAT"),
    "ClassII_DNA_hAT_nMITE": (DNA, "TIR / hAT"),
    "ClassII_DNA_hAT_MITE": (DNA, "TIR / hAT"),
    "ClassII_DNA_hAT_unknown": (DNA, "TIR / hAT"),
    "DNA/DTC": (DNA, "TIR / CACTA"),
    "MITE/DTC": (DNA, "TIR / CACTA"),
    "ClassII_DNA_CACTA_nMITE": (DNA, "TIR / CACTA"),
    "ClassII_DNA_CACTA_MITE": (DNA, "TIR / CACTA"),
    "DNA/DTH": (DNA, "TIR / PIF-Harbinger"),
    "MITE/DTH": (DNA, "TIR / PIF-Harbinger"),
    "ClassII_DNA_Harbinger_nMITE": (DNA, "TIR / PIF-Harbinger"),
    "ClassII_DNA_Harbinger_MITE": (DNA, "TIR / PIF-Harbinger"),
    "DNA/DTM": (DNA, "TIR / Mutator"),
    "MITE/DTM": (DNA, "TIR / Mutator"),
    "TIR/Mule": (DNA, "TIR / Mutator"),
    "ClassII_DNA_Mutator_nMITE": (DNA, "TIR / Mutator"),
    "ClassII_DNA_Mutator_MITE": (DNA, "TIR / Mutator"),
    "DNA/DTT": (DNA, "TIR / Tc1-Mariner"),
    "MITE/DTT": (DNA, "TIR / Tc1-Mariner"),
    "TIR/Tc1/mariner": (DNA, "TIR / Tc1-Mariner"),
    "ClassII_DNA_TcMar_nMITE": (DNA, "TIR / Tc1-Mariner"),
    "ClassII_DNA_TcMar_MITE": (DNA, "TIR / Tc1-Mariner"),
    "ClassII_DNA_TcMar_unknown": (DNA, "TIR / Tc1-Mariner"),
    "DNA/PiggyBac": (DNA, "TIR / PiggyBac"),
    "ClassII_MITE": (DNA, "MITE, superfamily undet."),
    "ClassII_nMITE": (DNA, "TIR, superfamily undet."),
    # --- Class III
    "DNA/Helitron": (HEL, "Helitron"),
    "Helitron/Helitron": (HEL, "Helitron"),
    "ClassIII_Helitron": (HEL, "Helitron"),
    # --- other
    "SDR/Crypton": (OTH, "Crypton"),
    "Tac/Tac": (OTH, "Tac (Fusarium-specific)"),
    "Marsu/Marsu": (OTH, "Marsu (Fusarium-specific)"),
    "tRNA": (OTH, "tRNA-related repeat"),
    "unknown": (OTH, "Unclassified"),
    "Unknown/Unknown": (OTH, "Unclassified"),
}
CLASS_ORDER = [RETRO, DNA, HEL, OTH]


def merged_bp(iv):
    iv = sorted(iv)
    tot = 0
    cs, ce = iv[0]
    for s, e in iv[1:]:
        if s <= ce + 1:
            ce = max(ce, e)
        else:
            tot += ce - cs + 1
            cs, ce = s, e
    return tot + ce - cs + 1


def main():
    sizes = {}
    for line in open(SIZES):
        c, l = line.split()[:2]
        sizes[c] = int(l)
    core = {c: l for c, l in sizes.items() if c not in ACCESSORY}
    acc = {c: l for c, l in sizes.items() if c in ACCESSORY}

    deepte = defaultdict(list)
    for line in open(DEEPTE):
        if ":" not in line:
            continue
        loc, cls = line.split()[:2]
        chrom, rng = loc.split(":")
        s, e = rng.split("-")
        deepte[chrom].append((int(s), int(e), cls))

    iv = defaultdict(list)          # (compartment, class, superfamily) -> intervals
    tot_iv = defaultdict(list)      # compartment -> all repeat intervals
    unmapped = defaultdict(int)
    with open(RM_OUT) as fh:
        for line in fh:
            f = line.split()
            if len(f) < 11 or not f[4].startswith("Chr"):
                continue
            chrom, s, e, lab = f[4], int(f[5]), int(f[6]), f[10]
            if chrom not in sizes:
                continue
            if "Simple_repeat" in lab or "Low_complexity" in lab:
                continue
            if lab == "unknown":
                for ds, de, dc in deepte.get(chrom, ()):
                    if ds <= s <= de or ds <= e <= de:
                        lab = dc
                        break
            comp = "accessory" if chrom in ACCESSORY else "core"
            key = MAP.get(lab)
            if key is None:
                unmapped[lab] += e - s + 1
                key = (OTH, "Unclassified")
            iv[(comp, key[0], key[1])].append((chrom, s, e))
            tot_iv[comp].append((chrom, s, e))

    def bp(entries):
        per = defaultdict(list)
        for c, s, e in entries:
            per[c].append((s, e))
        return sum(merged_bp(v) for v in per.values())

    keys = sorted({(k[1], k[2]) for k in iv},
                  key=lambda k: (CLASS_ORDER.index(k[0]),
                                 -bp(iv[("accessory", k[0], k[1])] or
                                     iv[("core", k[0], k[1])])))
    core_len, acc_len = sum(core.values()), sum(acc.values())
    gen_len = core_len + acc_len
    rows = []
    for cls, sf in keys:
        cb = bp(iv[("core", cls, sf)])
        ab = bp(iv[("accessory", cls, sf)])
        cp, ap = 100 * cb / core_len, 100 * ab / acc_len
        gb = bp(iv[("core", cls, sf)] + iv[("accessory", cls, sf)])
        rows.append({"class": cls, "superfamily": sf,
                     "core_bp": cb, "core_pct": cp,
                     "acc_bp": ab, "acc_pct": ap,
                     "gen_bp": gb, "gen_pct": 100 * gb / gen_len,
                     "ratio": (ap / cp) if cp else float("nan")})
    tab = pd.DataFrame(rows)

    # class subtotals and the all-repeat total, from merged intervals
    sub = []
    for cls in CLASS_ORDER:
        cb = bp([x for k, v in iv.items() if k[0] == "core" and k[1] == cls for x in v])
        ab = bp([x for k, v in iv.items() if k[0] == "accessory" and k[1] == cls for x in v])
        cp, ap = 100 * cb / core_len, 100 * ab / acc_len
        gb = bp([x for k, v in iv.items() if k[1] == cls for x in v])
        sub.append({"class": cls, "superfamily": "all %s" % cls.split(" (")[0],
                    "core_bp": cb, "core_pct": cp, "acc_bp": ab, "acc_pct": ap,
                    "gen_bp": gb, "gen_pct": 100 * gb / gen_len,
                    "ratio": ap / cp if cp else float("nan")})
    tcb, tab_ = bp(tot_iv["core"]), bp(tot_iv["accessory"])
    tgb = bp(tot_iv["core"] + tot_iv["accessory"])
    total = {"class": "All repeats", "superfamily": "total (non-redundant)",
             "core_bp": tcb, "core_pct": 100 * tcb / core_len,
             "acc_bp": tab_, "acc_pct": 100 * tab_ / acc_len,
             "gen_bp": tgb, "gen_pct": 100 * tgb / gen_len,
             "ratio": (100 * tab_ / acc_len) / (100 * tcb / core_len)}

    out = []
    for cls in CLASS_ORDER:
        for r in tab[tab["class"] == cls].to_dict("records"):
            r = dict(r); r["label"] = "   " + r["superfamily"]
            out.append(r)
        for r in sub:
            if r["class"] == cls:
                r = dict(r); r["label"] = cls.replace(" (retrotransposons)", " retrotransposons").replace(" (DNA transposons)", " DNA transposons").replace(" (Helitrons)", " Helitrons") + ", all"
                out.append(r)
    total = dict(total); total["label"] = "All repeats, total (non-redundant)"
    out.append(total)
    fin = pd.DataFrame(out)

    kb = lambda v: "0" if v == 0 else ("<1" if v < 500 else "{:,}".format(int(round(v / 1000))))
    fmt = pd.DataFrame({
        "Class and superfamily": fin["label"],
        "Core (kb)": fin.core_bp.map(kb),
        "Core (%)": fin.core_pct.map("{:.2f}".format),
        "Acc. (kb)": fin.acc_bp.map(kb),
        "Acc. (%)": fin.acc_pct.map("{:.2f}".format),
        "Genome (kb)": fin.gen_bp.map(kb),
        "Genome (%)": fin.gen_pct.map("{:.2f}".format),
        "Acc. : core": fin.ratio.map(lambda v: "-" if v != v else
                                     ("%.1f" % v if v >= 1 else "%.2f" % v)),
    })
    p = REV / "tables" / "Table_TE_classes_core_vs_accessory.tsv"
    fmt.to_csv(p, sep="\t", index=False)
    print("core %.2f Mb, accessory %.2f Mb, genome %.2f Mb; repeats %.2f%% / %.2f%% / %.2f%% of genome; rows %d"
          % (core_len / 1e6, acc_len / 1e6, gen_len / 1e6, total["core_pct"],
             total["acc_pct"], total["gen_pct"], len(fmt)))
    print("unclassified after DeepTE cross-ref: core %.2f%%, accessory %.2f%%"
          % (float(tab[tab.superfamily == "Unclassified"].core_pct.iloc[0]),
             float(tab[tab.superfamily == "Unclassified"].acc_pct.iloc[0])))
    if unmapped:
        print("labels not in MAP:", dict(sorted(unmapped.items(), key=lambda x: -x[1])[:10]))


if __name__ == "__main__":
    main()
