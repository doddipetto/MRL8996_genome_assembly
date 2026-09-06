#!/usr/bin/env python3
"""
01_pan_effectorome_cdhit.py
===========================
Three-strain pan-effectorome clustering for the MRL8996 genome paper revision.

Builds a pan-effectorome from the final (SignalP + EffectorP + DeepTMHMM-filtered)
effector sets of:
    - Fusarium oxysporum MRL8996  (clinical keratitis isolate, this study)
    - F. oxysporum f. sp. lycopersici 4287 (Fol4287, tomato pathogen)
    - F. oxysporum Fo47 (endophytic biocontrol strain)

and clusters it with CD-HIT using exactly the parameters used for the companion
FO12 data descriptor (-c 0.5 -n 3 -d 0; parameters back-inferred from
F012_HiC_data/final_results/Effectorome_CDHIT/clustered_effectors.clstr, whose
minimum within-cluster identity is 50.0% and whose minimum within-cluster length
ratio is 0.199, i.e. no -s/-aS length constraint was applied).

Outputs (all written under Revision_2026/, nothing outside it is modified):
    results/pan_effectorome_MRL8996.fasta      strain-tagged concatenated input
    results/clustered_effectors_MRL8996        CD-HIT representatives
    results/clustered_effectors_MRL8996.clstr  CD-HIT clusters
    tables/cdhit_cluster_matrix.tsv            cluster x strain presence/absence
    tables/venn_counts.tsv                     counts for the 7 Venn regions
    results/lists/1_core_effectors_all_three.txt
    results/lists/2_MRL8996_and_Fol4287.txt
    results/lists/3_MRL8996_and_Fo47.txt
    results/lists/4_MRL8996_unique.txt

Usage:  python3 01_pan_effectorome_cdhit.py
"""

import os
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

# ----------------------------------------------------------------------------
# Paths
# ----------------------------------------------------------------------------
ADODDI = Path("/home/usuario/nvme_data/home/usuario/adoddi")
REV = ADODDI / "MRL8996_genome_paper" / "Revision_2026"

MRL_PROT_DIR = (ADODDI / "MRL8996_genome_paper" / "JGI_reference" / "FoxMRL8996"
                / "Mycocosm" / "Annotation" / "Filtered_Models___best__" / "Proteins")
MRL_IDS = MRL_PROT_DIR / "FoxMRL8996_Final_Effectorome_IDs.txt"
MRL_FASTA = MRL_PROT_DIR / "SignalP_EffectorP_intersect.fasta"

F012_FINAL = ADODDI / "F012_HiC_data" / "final_results" / "Effectorome_Final_Results"
FOL_FASTA = F012_FINAL / "Fol4287_final_effectors.fasta"
FO47_FASTA = F012_FINAL / "Fo47_final_effectors.fasta"

OUT_RES = REV / "results"
OUT_TAB = REV / "tables"
OUT_LISTS = OUT_RES / "lists"

CDHIT_PARAMS = ["-c", "0.5", "-n", "3", "-d", "0"]

STRAINS = ["MRL8996", "Fol4287", "Fo47"]


# ----------------------------------------------------------------------------
# FASTA helpers
# ----------------------------------------------------------------------------
def read_fasta(path):
    """Yield (header, sequence) tuples."""
    header, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq)
                header, seq = line[1:], []
            else:
                seq.append(line.strip())
    if header is not None:
        yield header, "".join(seq)


def mrl_tag(header):
    """jgi|FoxMRL8996|10088|CE10087_16778 ... -> MRL8996|10088_CE10087_16778"""
    first = header.split()[0]
    parts = first.split("|")
    if len(parts) >= 4 and parts[1] == "FoxMRL8996":
        return "MRL8996|%s_%s" % (parts[2], parts[3])
    raise ValueError("unexpected MRL8996 header: %s" % header)


# ----------------------------------------------------------------------------
# Build the pan-effectorome
# ----------------------------------------------------------------------------
def build_pan(pan_path):
    counts = {}
    written = set()
    with open(pan_path, "w") as out:

        # --- MRL8996: restrict the 535-entry SignalP/EffectorP intersect to the
        #     513 IDs that survived DeepTMHMM filtering (the final effectorome).
        wanted = set()
        with open(MRL_IDS) as fh:
            for line in fh:
                line = line.strip()
                if line:
                    wanted.add(line.split()[0])
        n = 0
        for header, seq in read_fasta(MRL_FASTA):
            if header.split()[0] not in wanted:
                continue
            tag = mrl_tag(header)
            if tag in written:
                continue
            written.add(tag)
            out.write(">%s\n%s\n" % (tag, seq))
            n += 1
        counts["MRL8996"] = n
        if n != len(wanted):
            print("[WARN] %d of %d MRL8996 final effector IDs recovered from %s"
                  % (n, len(wanted), MRL_FASTA.name), file=sys.stderr)

        # --- Fol4287 and Fo47: RefSeq protein accessions
        for strain, fasta in (("Fol4287", FOL_FASTA), ("Fo47", FO47_FASTA)):
            n = 0
            for header, seq in read_fasta(fasta):
                acc = header.split()[0]
                tag = "%s|%s" % (strain, acc)
                if tag in written:
                    continue
                written.add(tag)
                out.write(">%s\n%s\n" % (tag, seq))
                n += 1
            counts[strain] = n

    return counts


# ----------------------------------------------------------------------------
# CD-HIT
# ----------------------------------------------------------------------------
def run_cdhit(pan_path, out_prefix, log_path):
    cmd = ["cd-hit", "-i", str(pan_path), "-o", str(out_prefix)] + CDHIT_PARAMS
    print("[cd-hit] " + " ".join(cmd))
    with open(log_path, "w") as log:
        log.write(" ".join(cmd) + "\n")
        log.flush()
        subprocess.run(cmd, check=True, stdout=log, stderr=subprocess.STDOUT)
    return cmd


def parse_clstr(clstr_path):
    """-> {cluster_id: [member_tag, ...]}, {cluster_id: representative_tag}"""
    members = defaultdict(list)
    reps = {}
    cid = None
    with open(clstr_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">Cluster"):
                cid = int(line.split()[1])
                continue
            m = re.search(r">(\S+?)\.\.\.", line)
            if not m:
                continue
            tag = m.group(1)
            members[cid].append(tag)
            if line.rstrip().endswith("*"):
                reps[cid] = tag
    return dict(members), reps


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def main():
    for d in (OUT_RES, OUT_TAB, OUT_LISTS, REV / "logs"):
        d.mkdir(parents=True, exist_ok=True)

    pan = OUT_RES / "pan_effectorome_MRL8996.fasta"
    counts = build_pan(pan)
    print("Input effectors: " + ", ".join("%s=%d" % (s, counts[s]) for s in STRAINS))

    prefix = OUT_RES / "clustered_effectors_MRL8996"
    run_cdhit(pan, prefix, REV / "logs" / "cdhit.log")

    members, reps = parse_clstr(Path(str(prefix) + ".clstr"))
    print("CD-HIT clusters: %d" % len(members))

    # --- cluster x strain matrix
    matrix_path = OUT_TAB / "cdhit_cluster_matrix.tsv"
    region_of = {}
    with open(matrix_path, "w") as out:
        out.write("Cluster_ID\tRepresentative\tCluster_size\t"
                  + "\t".join("n_%s" % s for s in STRAINS)
                  + "\tPattern\tRegion\t"
                  + "\t".join("%s_members" % s for s in STRAINS) + "\n")
        for cid in sorted(members):
            by_strain = defaultdict(list)
            for tag in members[cid]:
                strain, ident = tag.split("|", 1)
                by_strain[strain].append(ident)
            present = tuple(1 if by_strain.get(s) else 0 for s in STRAINS)
            pattern = "".join(str(x) for x in present)
            region = {
                (1, 1, 1): "core_all_three",
                (1, 1, 0): "MRL8996_Fol4287",
                (1, 0, 1): "MRL8996_Fo47",
                (1, 0, 0): "MRL8996_unique",
                (0, 1, 1): "Fol4287_Fo47",
                (0, 1, 0): "Fol4287_unique",
                (0, 0, 1): "Fo47_unique",
            }[present]
            region_of[cid] = region
            out.write("\t".join([
                str(cid), reps.get(cid, "NA"), str(len(members[cid])),
                *[str(len(by_strain.get(s, []))) for s in STRAINS],
                pattern, region,
                *[",".join(by_strain.get(s, [])) or "-" for s in STRAINS],
            ]) + "\n")

    # --- Venn region counts (clusters, and protein counts per strain)
    region_clusters = defaultdict(int)
    region_prots = defaultdict(lambda: defaultdict(int))
    for cid, region in region_of.items():
        region_clusters[region] += 1
        for tag in members[cid]:
            region_prots[region][tag.split("|", 1)[0]] += 1

    order = ["core_all_three", "MRL8996_Fol4287", "MRL8996_Fo47", "Fol4287_Fo47",
             "MRL8996_unique", "Fol4287_unique", "Fo47_unique"]
    with open(OUT_TAB / "venn_counts.tsv", "w") as out:
        out.write("Region\tn_clusters\t" + "\t".join("n_%s_proteins" % s for s in STRAINS) + "\n")
        for region in order:
            out.write("\t".join([region, str(region_clusters[region])]
                                + [str(region_prots[region][s]) for s in STRAINS]) + "\n")
        out.write("TOTAL\t%d\t%s\n" % (len(members),
                                       "\t".join(str(counts[s]) for s in STRAINS)))

    # --- ID lists (MRL8996-centred, as requested by the reviewer response)
    list_files = {
        "core_all_three": OUT_LISTS / "1_core_effectors_all_three.txt",
        "MRL8996_Fol4287": OUT_LISTS / "2_MRL8996_and_Fol4287.txt",
        "MRL8996_Fo47": OUT_LISTS / "3_MRL8996_and_Fo47.txt",
        "MRL8996_unique": OUT_LISTS / "4_MRL8996_unique.txt",
    }
    for region, path in list_files.items():
        with open(path, "w") as out:
            out.write("MRL8996_effector\tCluster_ID\tCluster_size\tFol4287_members\tFo47_members\n")
            for cid in sorted(members):
                if region_of[cid] != region:
                    continue
                by_strain = defaultdict(list)
                for tag in members[cid]:
                    s, ident = tag.split("|", 1)
                    by_strain[s].append(ident)
                for ident in by_strain["MRL8996"]:
                    out.write("\t".join([
                        ident, str(cid), str(len(members[cid])),
                        ",".join(by_strain.get("Fol4287", [])) or "-",
                        ",".join(by_strain.get("Fo47", [])) or "-",
                    ]) + "\n")

    print("\nVenn regions (clusters):")
    for region in order:
        print("  %-18s %4d clusters  (%s)" % (
            region, region_clusters[region],
            ", ".join("%s:%d" % (s, region_prots[region][s]) for s in STRAINS
                      if region_prots[region][s])))


if __name__ == "__main__":
    main()
