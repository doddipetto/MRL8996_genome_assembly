#!/usr/bin/env python3
"""Patch the author version 2.0 of the manuscript with the corrections that are
still missing after the reviewers' comments, marking every insertion or change in
orange.  The document is edited in place (runs spliced, paragraphs inserted) so
that Word citation fields, table formatting and styles survive untouched.

in :  MRL8996_genome_paper/MRL8996_2026_v2.0.docx
out:  Revision_2026/manuscript/MRL8996_2026_v2.1_highlighted.docx
"""
import copy
import shutil
from pathlib import Path

import docx
from docx.oxml.ns import qn
from docx.oxml import OxmlElement
from docx.shared import RGBColor

ROOT = Path("/home/usuario/nvme_data/home/usuario/adoddi/MRL8996_genome_paper")
SRC = ROOT / "MRL8996_2026_v2.0.docx"
OUT = ROOT / "Revision_2026" / "manuscript" / "MRL8996_2026_v2.1_highlighted.docx"

ORANGE = RGBColor(0xC0, 0x50, 0x00)
SHADE = "FFE0B2"


def mark(run):
    """Orange text on a light-orange character shading."""
    run.font.color.rgb = ORANGE
    rPr = run._element.get_or_add_rPr()
    shd = OxmlElement("w:shd")
    shd.set(qn("w:val"), "clear")
    shd.set(qn("w:color"), "auto")
    shd.set(qn("w:fill"), SHADE)
    rPr.append(shd)
    return run


def norm(s):
    """Length-preserving normalisation so that ASCII needles match Word's
    non-breaking spaces, curly quotes and dashes."""
    for a, b in (("\xa0", " "), ("\u2019", "'"), ("\u2018", "'"),
                 ("\u201c", '"'), ("\u201d", '"'),
                 ("\u2013", "-"), ("\u2014", "-"), ("\u2212", "-")):
        s = s.replace(a, b)
    return s


def para(doc, anchor, nth=0):
    hits = [p for p in doc.paragraphs if norm(anchor) in norm(p.text)]
    if not hits:
        raise SystemExit("anchor not found: %r" % anchor[:60])
    return hits[nth]


def splice(p, old, new):
    """Replace `old` by `new` inside paragraph p, touching only the runs that
    overlap the match.  `new` is inserted as an orange run."""
    runs = p.runs
    text = norm("".join(r.text for r in runs))
    i = text.find(norm(old))
    if i < 0:
        raise SystemExit("text not found: %r" % old[:70])
    j = i + len(old)
    pos, first, tail_run = 0, None, None
    for r in runs:
        a, b = pos, pos + len(r.text)
        pos = b
        if b <= i or a >= j:
            continue
        head = r.text[: max(0, i - a)]
        tail = r.text[max(0, j - a):] if b > j else ""
        if first is None:
            first = r
            r.text = head
            newrun = copy.deepcopy(r._element)
            r._element.addnext(newrun)
            nr = docx.text.run.Run(newrun, p)
            nr.text = new
            mark(nr)
            if tail:
                trun = copy.deepcopy(first._element)
                newrun.addnext(trun)
                tail_run = docx.text.run.Run(trun, p)
                tail_run.text = tail
                tail_run.font.color.rgb = None
        else:
            r.text = tail
    return p


def renumber(doc, old, new):
    """Replace every occurrence of `old` in the document by `new`."""
    n = 0
    for p in doc.paragraphs:
        while norm(old) in norm("".join(r.text for r in p.runs)):
            splice(p, old, new)
            n += 1
    if n == 0:
        raise SystemExit("renumber found nothing: %r" % old)
    return n


def append(p, text):
    r = p.add_run(text)
    mark(r)
    return p


def insert_after(p, text, style=None, bold=False):
    new_p = copy.deepcopy(p._p)
    for child in list(new_p):
        if child.tag == qn("w:r") or child.tag == qn("w:hyperlink") or child.tag == qn("w:fldSimple"):
            new_p.remove(child)
    p._p.addnext(new_p)
    np = docx.text.paragraph.Paragraph(new_p, p._parent)
    if style is not None:
        np.style = style
    r = np.add_run(text)
    r.bold = bold
    mark(r)
    return np


def drop(p):
    p._p.getparent().remove(p._p)


def main():
    doc = docx.Document(str(SRC))

    # ------------------------------------------------------------------ E7
    # The twelve large tables leave the manuscript files and are cited as data
    # files of the repository deposition; the three small new tables added in
    # this revision take the Supplementary Table numbering.
    for old, new in [("Supplementary Table 11", "Data File 11"),
                     ("Supplementary Table 12", "Data File 12"),
                     ("Supplementary Table 1", "Data File 1"),
                     ("Supplementary Tables 2 and 3", "Data Files 2 and 3"),
                     ("Supplementary Table 4", "Data File 4"),
                     ("Supplementary Table 5", "Data File 5"),
                     ("Supplementary Table 6", "Data File 6"),
                     ("Supplementary Table 7", "Data File 7"),
                     ("Supplementary Table 8-9", "Data Files 8-10")]:
        renumber(doc, old, new)


    # ---------------------------------------------------------------- R2.6
    splice(para(doc, "strain MRL8996 (NRRL 47514), originally isolated"),
           "strain MRL8996 (NRRL 47514)",
           "strain MRL8996 (deposited as NRRL 47514 in the Agricultural Research "
           "Service Culture Collection, Peoria, IL, USA)")

    p18 = para(doc, "The currently available draft genome of MRL8996")
    splice(p18, "The currently available draft genome of MRL8996(PRJNA554890) is fragmented",
           "The currently available draft genome of MRL8996 (v1.0; GenBank GCA_009746015.1, "
           "BioProject PRJNA554890), generated as part of the first comparative-genomic "
           "analysis of this isolate [CIT: Zhang et al. 2020], is fragmented")

    # ---------------------------------------------------------------- R1.1
    append(p18,
           " The phylogenetic position of MRL8996 within the FOSC has already been "
           "established and is not re-examined here: in a maximum-likelihood phylogeny of "
           "55 conserved single-copy orthologous genes of the genus Fusarium, rooted on "
           "F. verticillioides, MRL8996 (NRRL 47514) falls in the same clade as the "
           "clinical isolate NRRL 32931, within a subclade supported by 100% bootstrap "
           "that also contains the two reference strains used for comparison here, "
           "Fol4287 (NRRL 34936) and Fo47 (NRRL 54002) [CIT: Zhang et al. 2020].")

    # ------------------------------------------------------- E2, section order
    p20 = para(doc, "The improved MRL8996 chromosome-level genome assembly was generated")
    old_head = para(doc, "DATA OVERVIEW")
    head_style, head_p = old_head.style, old_head._p
    new_head = copy.deepcopy(head_p)
    for child in list(new_head):
        if child.tag in (qn("w:r"), qn("w:hyperlink"), qn("w:fldSimple")):
            new_head.remove(child)
    p20._p.addprevious(new_head)
    nh = docx.text.paragraph.Paragraph(new_head, p20._parent)
    nh.style = head_style
    mark(nh.add_run("DATA OVERVIEW"))
    drop(old_head)

    # ------------------------------------------------------- R2.19, R2.7, R2.8
    p23 = para(doc, "The hybrid assembly pipeline produced 16 chromosome-scale scaffolds")
    splice(p23, "In addition, StainedGlass was used to confirm the detection of centromeric regions.",
           "Candidate centromeric regions were called from the intersection of GC-content "
           "minima, the Hi-C cis-interaction maximum of each chromosome and bands of "
           "84-93% self-identity in a StainedGlass map, and are listed in Supplementary "
           "Table 1 (see Methods).")
    splice(p23, "with a scaffold N50 of 4.82 Mb, an L50 of 5,",
           "with a scaffold N50 of 4.82 Mb, an auN of 4.15 Mb, an L50 of 5,")
    splice(p23, "(Figure 1B), where the core chromosomes shared across the F. oxysporum "
                "species complex are clearly conserved.",
           "(Figure 1B). This overview does not itself demonstrate conservation of the "
           "core compartment; that is shown by the synteny analysis below (Figure 3A).")

    # ---------------------------------------------------------- R2.9 renumbering
    p25 = para(doc, "A direct comparison with the previously available draft genome")
    for _ in range(4):
        splice(p25, "Table 2", "Table 1")
    append(p25,
           " The area under the Nx curve, which summarises contiguity across the whole "
           "size distribution instead of at a single quantile, increases from 1.87 Mb in "
           "v1.0 to 4.15 Mb in v2.0 (Table 1).")
    splice(para(doc, "Whole-genome alignment of the two assemblies using D-GENIES"),
           "Table 2", "Table 1")

    # ---------------------------------------------------------------- R2.22
    p30 = para(doc, "Overall, these findings define 11 stable core chromosomes")
    insert_after(p30,
        "Gene and repeat content differ systematically between the two compartments. "
        "Chr01-Chr11 carry 309-366 genes per Mb, 45-50% coding sequence and 3.5-10.5% "
        "repeat-masked sequence, and contain all 56 predicted biosynthetic gene clusters, "
        "whereas Chr12-Chr16 carry 217-277 genes per Mb, 21-28% coding sequence, "
        "23.7-50.3% repeat-masked sequence and no biosynthetic gene cluster "
        "(Supplementary Table 2). The five accessory chromosomes together span 6.73 Mb, "
        "12.9% of the assembly, and account for 1,688 of the 16,612 annotated genes. "
        "Repeats occupy 6.0% of the core and 34.5% of the accessory compartment, and every "
        "repeat class detected is proportionally more abundant in the accessory "
        "compartment; the largest contributions there are unclassified elements (15.1% of "
        "accessory sequence), LTR elements (5.5%), LINEs (4.2%) and DNA transposons (3.5%) "
        "(Table 2). Predicted secreted proteins are correspondingly less dense on the "
        "accessory chromosomes (37 candidates, 5.5 per Mb) than on the core chromosomes "
        "(476 candidates, 10.4 per Mb).")

    # ------------------------------------------------------------ R2.1 rationale
    p32 = para(doc, "Beyond chromosome-level architecture")
    splice(p32,
           "Beyond chromosome-level architecture, the MRL8996 genome encodes a large "
           "repertoire of predicted secreted proteins that may be relevant for its lifestyle.",
           "The dataset also includes a catalogue of predicted secreted proteins. It is "
           "included for two reasons. Small secreted proteins are the gene class most often "
           "implicated in the interaction of F. oxysporum with its hosts and in adaptation "
           "to new ones, so a catalogue for a clinical isolate is the part of the annotation "
           "most likely to be reused; and it is the part least accessible from a fragmented "
           "assembly, because assignment of a candidate to the core or the accessory "
           "compartment requires chromosome-level coordinates. The catalogue is therefore "
           "described here as a component of the resource, without functional inference.")
    splice(p32, "(Supplementary Image 2)", "(Figure 4A)")

    # -------------------------------------------------- R2.1 three-strain result
    p36 = para(doc, "Growing evidence suggests that fungal effectors may have originated")
    insert_after(p36,
        "To let users separate the shared from the isolate-specific part of the catalogue, "
        "the 513 MRL8996 candidates were clustered at 70% amino-acid identity together with "
        "the equivalent catalogues of Fol4287 (500 candidates) and Fo47 (393 candidates). "
        "The 1,406 sequences collapse into 625 clusters, of which 252 contain members of all "
        "three strains and 131 contain only MRL8996 proteins (135 proteins, 26.3% of the "
        "MRL8996 catalogue; Figure 4B, Supplementary Table 3). Of these 135 sequence-unique "
        "candidates, 124 lie on core and 11 on accessory chromosomes, and 61 are predicted "
        "antimicrobial; in both respects they do not differ in proportion from the shared "
        "candidates. Sequence-level uniqueness does not imply structural novelty: of the 122 "
        "sequence-unique candidates with a confident model, 90 have a structural counterpart "
        "elsewhere in the MRL8996 catalogue at a TM-score of at least 0.5, and 90 have one "
        "among the 500 Fol4287 models, leaving 26 candidates with no structural counterpart "
        "in either set. In the structural network (Figure 5B), 44 of these 122 candidates are "
        "among the 131 nodes that share no edge above the threshold, and the family most "
        "enriched in isolate-specific members is the KP4-like family (Family 4, 13 of 25 "
        "members).")

    # ------------------------------------------------------------------- E3
    append(para(doc, "The MRL8996 strain used in this study was originally isolated"),
           " The isolate was provided by Prof. Li-Jun Ma (University of Massachusetts "
           "Amherst, USA). No human subjects, patient material or cell lines were handled "
           "in this study; the strain is a fungal culture derived from a clinical specimen "
           "collected in 2006.")

    # ------------------------------------------------------------------ R2.13
    p44 = para(doc, "For genomic DNA isolation, mycelium was collected by filtration")
    splice(p44,
           "However, an additional cleanup step was introduced to obtain "
           "ultra-high-molecular-weight (UHMW) DNA for our sequencing analyses. Following "
           "DNA precipitation, the samples were treated with the Short Read Eliminator Kit "
           "(PacBio, California, USA) according to the manufacturer's instructions, "
           "effectively depleting DNA fragments <25 kb.",
           "Two DNA preparations were made from the same culture protocol. The first "
           "followed the published protocol without modification and provided the DNA for "
           "the sequencing run deposited under SRR39391754, from which the assembly "
           "reported here was built. In the second preparation a size-selection step was "
           "added after DNA precipitation, in which the sample was treated with the Short "
           "Read Eliminator Kit (PacBio, California, USA) according to the manufacturer's "
           "instructions to deplete fragments below 25 kb; this preparation was used only "
           "to test whether the read-length distribution could be shifted further and "
           "contributed no reads to the assembly. The size selection removes short "
           "fragments and does not extract molecules longer than those the published "
           "protocol yields, so the two preparations differ in read-length distribution "
           "rather than in the length of the molecules extracted.")

    p45 = para(doc, "Library preparation with the Ligation Sequencing Kit")
    splice(p45,
           " To this end, after DNA precipitation, the Short Reads Eliminator Kit (PacBio, "
           "California, USA) was used following the manufacturer's protocol to select "
           "DNA fragments >25 kb.", "")

    # ------------------------------------------------------------------ R2.14
    p49 = para(doc, "To obtain a first draft assembly, the reads from Nanopore sequencing")
    splice(p49,
           " Hi-C library construction of F. oxysporum strain MRL8996 was performed "
           "according to the standard protocol and sequenced on the Illumina NovaSeq 6000 "
           "platform (Phase Genomics Inc., Seattle, WA, USA).", "")
    append(para(doc, "For Hi-C sequencing, an aliquot of MRL8996 spores"),
           " The Hi-C library was then constructed according to the standard Phase Genomics "
           "protocol and sequenced on an Illumina NovaSeq 6000 platform (Phase Genomics "
           "Inc., Seattle, WA, USA).")

    # ------------------------------------------------------------- R2.7 Methods
    p_syn = para(doc, "Comparative genome alignment and synteny")
    body = para(doc, "Synteny between MRL8996 and the reference assemblies")
    h = insert_after(body, "Candidate centromeric regions", style=p_syn.style, bold=False)
    insert_after(h,
        "F. oxysporum centromeres are not marked by a satellite repeat; they are AT-rich, "
        "gene-poor regions of a few tens of kb that behave as a single interaction focus in "
        "Hi-C maps. Candidates were therefore called from the intersection of three signals: "
        "minima of GC content in 20-kb windows, at least three percentage points below the "
        "chromosome median; the position of the strongest cis-interaction focus in the Hi-C "
        "contact matrix; and bands of 84-93% self-identity in a StainedGlass map computed "
        "from an all-versus-all alignment of the assembly in 2-kb windows. Calls supported "
        "by at least two of the three signals are reported in Supplementary Table 1, and "
        "the per-chromosome profiles are shown in Supplementary Image 4.")

    # ---------------------------------------------- Table and image renumbering
    p52 = para(doc, "Gene models for the chromosome-level assembly")
    splice(p52, "(Table 3, Data File 7)", "(Table 2, Data File 7)")
    splice(p52, "(Supplementary Image 3)", "(Supplementary Image 2)")
    splice(para(doc, "The accuracy and completeness of the assembly were evaluated"),
           "(Supplementary Image 4)", "(Supplementary Image 3)")

    # -------------------------------------------------- network: keep singletons
    splice(para(doc, "A structural-similarity network (Figure 4B) was constructed independently"),
           "To filter out background noise, nodes with a degree <3 were removed.",
           "All 482 nodes were retained, including the 131 that share no edge above the "
           "threshold, so that candidates without a structural counterpart remain visible "
           "in the layout.")

    # ------------------------------------------------ E5, E6, E7, R2.2 Data Records
    p65 = para(doc, "The raw sequencing data of Nanopore and HiC have been deposited")
    append(p65,
           " Each of these datasets is cited in the reference list as a data citation: the "
           "assembly GCA_009746015.2 [CIT: data citation], the Nanopore run SRR39391754 "
           "[CIT: data citation] and the Hi-C run SRR39391753 [CIT: data citation]. The "
           "annotation and all derived data are deposited as a single record at Figshare "
           "under DOI [AUTHOR ACTION: insert DOI] [CIT: data citation]. That record holds "
           "the gene annotation in GFF3 format lifted over to the chromosome-level "
           "coordinates and the predicted proteome; the transposable-element annotation "
           "with the non-redundant TE library and its DeepTE classification; the effector "
           "catalogue with the per-candidate SignalP, EffectorP, DeepTMHMM, pLDDT, "
           "structural-family and AMAPEC assignments, the all-versus-all TM-score matrix "
           "and the structural-similarity network; the CD-HIT clustering of the MRL8996, "
           "Fol4287 and Fo47 catalogues with the isolate-specific candidates and their "
           "structural comparisons; the chromosome-by-chromosome alignment matrices and "
           "the minimap2 alignments behind them; the candidate centromeric regions with "
           "their supporting signals; the per-chromosome content table; the biosynthetic "
           "gene cluster coordinates; the read-depth windows, BUSCO, Merqury and QUAST "
           "summaries and the Hi-C contact matrix; and the alignment and trees of the "
           "FOSC phylogeny. Every file is listed with its checksum and a one-line "
           "description in the manifest of the deposition. The twelve large data tables "
           "that were previously supplied as manuscript files are part of this record, "
           "where they are named Data File 1 to Data File 12 and are cited in the text "
           "under those names; they are no longer distributed with the manuscript.")

    # ------------------------------------------------------------------ R2.15
    p75 = para(doc, "The genome assembly of F. oxysporum strain MRL8996 has been deposited")
    for child in list(p75._p):
        if child.tag in (qn("w:r"), qn("w:hyperlink"), qn("w:fldSimple"), qn("w:proofErr")):
            p75._p.remove(child)
    mark(p75.add_run(
        "All data listed in Data Records are publicly available without restriction and "
        "under the accessions and the DOI given there. No permissions or material transfer "
        "agreements are required for their reuse."))

    # -------------------------------------------------------- figure legends
    splice(para(doc, "Figure 1. Overview of the chromosome-level genome assembly"),
           "(A) Genome-wide Hi-C contact map showing the interaction matrix among the 16 "
           "assembled chromosomes.",
           "(A) Genome-wide Hi-C contact map showing the interaction matrix among the 16 "
           "assembled chromosomes; the boxes along the diagonal delimit the 16 "
           "pseudomolecules, and the candidate centromeric region called for each "
           "chromosome (Supplementary Table 1) coincides with the strongest "
           "cis-interaction focus inside the corresponding block.")

    splice(para(doc, "Figure 4. Whole-genome distribution of the MRL8996 predicted effector"),
           "while red vertical lines indicate predicted antimicrobial effectors.",
           "while red vertical lines indicate predicted antimicrobial effectors. Black "
           "triangles above the chromosome bars mark the 135 candidates that are specific "
           "to MRL8996 in the three-strain comparison shown in (B).")

    splice(para(doc, "Figure 5. Structure-based classification of the predicted effector"),
           "nodes with a degree < 3 were excluded.",
           "all 482 nodes are shown, including the 131 that share no edge above the "
           "threshold, which are drawn as unconnected nodes so that the position of the "
           "MRL8996-specific candidates within the network can be seen.")

    append(para(doc, "Supplementary Image 1. Telomere identification across the 16"),
           " Bar colour gives the Nanopore read depth in 200-kb windows as a multiple of "
           "the genome-median depth of 47.5x, on the six-level scale of tapestry (0, 0.5x, "
           "1x, 1.5x, 2x, >=2.5x) shown in the key below the plot. Three regions exceed one "
           "genomic copy of depth: the terminal 70 kb of Chr05, which carries eight tandem "
           "rDNA units in the assembly against an estimated ~124 units in the reads; a "
           "4.7-kb mitochondrial insertion at Chr05:4.55 Mb; and the first contig of Chr13 "
           "(0-0.49 Mb), which is present at about 1.8x uniquely-mapping depth and is the "
           "one segment of the assembly that appears to represent two genomic copies. The "
           "first two are collapsed tandem repeats of known composition rather than "
           "misassemblies of unique sequence.")

    last_supp = para(doc, "Supplementary Image 3. Analysis of k-mer spectra using Merqury")
    insert_after(last_supp,
        "Supplementary Image 4. Per-chromosome centromere diagnostics. For each of the 16 "
        "chromosomes, the density of off-diagonal 2-kb window pairs at 84-93% identity in "
        "the StainedGlass self-alignment (red, left axis) and the GC content in 20-kb "
        "windows (blue, right axis); the dotted line marks the chromosome median GC content "
        "minus three percentage points, the depth threshold used to accept a GC trough, and "
        "the shaded band marks the call reported in Supplementary Table 1.")

    # --------------------------------------------------------- Table 3 -> Table 2
    splice(para(doc, "Table 3. Statistics of repeat sequences"), "Table 3.", "Table 2.")


    # ------------------------------------- figure numbering of the structural panels
    for anchor, pairs in [
        ("Three-dimensional models for the 513 candidate effectors were predicted",
         [("hierarchical clustering (Figure 4A)", "hierarchical clustering (Figure 5A)"),
          ("DALI network organisation (Figure 4B)", "DALI network organisation (Figure 5B)")]),
        ("To further validate these groupings, an independent force-directed",
         [("network (Figure 4B)", "network (Figure 5B)"),
          ("(LPMO; Family 5) (Figure 4C)", "(LPMO; Family 5) (Figure 5C)")]),
        ("Growing evidence suggests that fungal effectors may have originated",
         [("annotated on the heatmap (Figure 4A)", "annotated on the heatmap (Figure 5A)")]),
        ("For the structural heatmap", [("structural heatmap (Figure 4A)", "structural heatmap (Figure 5A)")]),
        ("A structural-similarity network", [("network (Figure 4B)", "network (Figure 5B)")]),
        ("network analysis and visualisation were performed",
         [("structural network (Figure 4C)", "structural network (Figure 5C)")]),
        ("Antimicrobial activity was predicted with AMAPEC",
         [("beneath the heatmap (Figure 4A", "beneath the heatmap (Figure 5A")]),
    ]:
        p = para(doc, anchor)
        for old, new in pairs:
            splice(p, old, new)

    # ------------------------------------------------------- R2.7 / R2.10 validation
    splice(para(doc, "The obtained genome assembly was manually corrected"),
           "We utilised StainedGlass and tapestry v1.0.1 to position all centromeres and "
           "almost all telomeres (27/32) (TTAGGG) (Supplementary Image 1).",
           "tapestry v1.0.1 was used to locate the telomeric repeat (TTAGGG) at 27 of the "
           "32 chromosome ends (Supplementary Image 1), and StainedGlass self-identity, GC "
           "content and the Hi-C cis-interaction maxima were combined to call candidate "
           "centromeric regions, which are supported by at least two of the three signals "
           "on 10 of the 16 chromosomes (Supplementary Table 1, Supplementary Image 4).")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    doc.save(str(OUT))
    print("written", OUT)
    print("orange runs:", sum(
        1 for p in doc.paragraphs for r in p.runs
        if r.font.color is not None and r.font.color.rgb == ORANGE))


if __name__ == "__main__":
    main()
