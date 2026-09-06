#!/usr/bin/env python3
"""Build the point-by-point response letter in the same style as the FO12
Scientific Data response (salutation, thanks, Editor/Reviewer blocks quoting the
comment, each answered under an 'Authors:' tag, prose only).

out: Revision_2026/manuscript/Response_to_reviewers_MRL8996_SciData.docx
"""
from pathlib import Path

import docx
from docx.shared import Pt

ROOT = Path("/home/usuario/nvme_data/home/usuario/adoddi/MRL8996_genome_paper")
OUT = ROOT / "Revision_2026" / "manuscript" / "Response_to_reviewers_MRL8996_SciData.docx"

HEAD = [
    "Dear Editor,",
    "We thank the editorial team and both reviewers for their careful evaluation of our "
    "manuscript and for the detailed and constructive comments, which have improved the "
    "clarity and the completeness of the data descriptor. Below we answer every point in "
    "turn. Changes to the manuscript are marked in the revised file.",
]

TAIL = [
    "We trust that the revised version now meets the standards of Scientific Data.",
    "On behalf of all authors,",
    "Andrea Doddi (corresponding author)",
    "Department of Genetics, University of Cordoba",
    "adoddi@uco.es",
]

BLOCKS = [
    ("H", "IN-HOUSE EDITOR COMMENTS"),

    ("Editor", "Why does the submitter note 'University of Massachusetts Amherst' for "
               "GCA_009746015.2?"),
    ("Authors", "The submitting institution recorded for GCA_009746015.2 is the "
     "laboratory of Prof. Li-Jun Ma at the University of Massachusetts Amherst, which "
     "generated and deposited the v1.0 draft assembly of this isolate and provided us "
     "with the strain. The chromosome-level assembly reported here is an update of that "
     "record within the same BioProject, so the original submitter is retained in the "
     "GenBank metadata. The provenance of the strain is now stated explicitly in the "
     "Methods, and Prof. Ma is acknowledged."),

    ("Editor", "Please check the Background and Summary section to ensure it is in scope "
               "and does not contain any results, analyses, conclusions or facts relating "
               "to what the data shows; some elements could move to a Data Overview "
               "section."),
    ("Authors", "The Background & Summary has been shortened to the clinical relevance of "
     "the species complex, the state of the art of F. oxysporum genome assemblies and the "
     "aim of the dataset, and now carries two subheadings. All statements about how the "
     "data were generated and what they show, including the description of the sequencing "
     "and assembly and all assembly statistics, have been moved to a separate Data "
     "Overview section that opens with the paragraph beginning 'The improved MRL8996 "
     "chromosome-level genome assembly was generated'."),

    ("Editor", "Because some data look to be derived from cell lines, please make clear "
               "where these were obtained from."),
    ("Authors", "No cell lines, human subjects or patient materials were used. The only "
     "biological material in this study is a fungal culture, F. oxysporum MRL8996, "
     "originally recovered in 2006 from the contact lens of a patient with keratitis and "
     "provided to us by Prof. Li-Jun Ma (University of Massachusetts Amherst, USA); the "
     "isolate is deposited in the ARS Culture Collection as NRRL 47514. A sentence stating "
     "this has been added to the first Methods paragraph."),

    ("Editor", "Rename the section 'MATERIALS AND METHODS' to 'Methods'."),
    ("Authors", "The section has been renamed 'Methods'."),

    ("Editor", "Please confirm that annotation files have been shared with the assembly at "
               "GenBank/ENA or at a compliant repository, that this is stated in Data "
               "Records and Data availability, and that the deposition is cited."),
    ("Authors", "The annotation files are now deposited in a single archive at a "
     "general-purpose repository, comprising the gene annotation in GFF3 format lifted "
     "over to the chromosome-level coordinates, the predicted proteome, the "
     "transposable-element annotation together with the non-redundant TE library, the "
     "effector catalogue with the per-candidate SignalP, EffectorP, DeepTMHMM, pLDDT, "
     "structural-family and AMAPEC assignments, the per-chromosome content table, the "
     "CD-HIT clustering of the three effector catalogues and the pairwise "
     "structural-similarity matrices. The deposition is described in Data Records, cited "
     "there as a data citation, and referred to from Data availability."),

    ("Editor", "Please add a data citation for the GenBank data to the reference list and "
               "cite it wherever the dataset is mentioned, principally in the first part of "
               "Data Records."),
    ("Authors", "Formal data citations have been added to the reference list for the "
     "genome assembly (GCA_009746015.2), for the Nanopore run (SRR39391754), for the Hi-C "
     "run (SRR39391753) and for the annotation archive, in the format requested, and each "
     "is cited at the point where the dataset is first mentioned in Data Records."),

    ("Editor", "Large data tables have been shared in the manuscript files; please move "
               "them to a repository or delete duplicates, and adjust mentions "
               "accordingly."),
    ("Authors", "The twelve large data tables have been removed from the manuscript files "
     "and are now part of a single repository record, together with all of the source data "
     "a reader needs to reproduce the figures and tables: the lifted-over gene annotation "
     "and the predicted proteome, the transposable-element annotation and TE library, the "
     "effector catalogue with every per-candidate score, the all-versus-all TM-score "
     "matrix and the structural network, the CD-HIT clustering of the three strains and "
     "the isolate-specific candidates with their structural comparisons, the synteny "
     "matrices and the alignments behind them, the candidate centromeric regions, the "
     "biosynthetic gene cluster coordinates, the technical-validation summaries and the "
     "Hi-C contact matrix, and the alignment and trees of the phylogeny. The record is "
     "organised in ten directories and carries a README and a manifest giving the size, "
     "checksum and a one-line description of every file. The twelve tables were exported "
     "verbatim from the submitted workbook and are named Data File 1 to Data File 12 in "
     "the record; all in-text mentions have been renamed accordingly. The supplementary "
     "workbook has been rebuilt so that it no longer contains any of them, and now holds "
     "only the three small tables added during this revision, which the reader consults "
     "directly (candidate centromeric regions, per-chromosome content, and the three-way "
     "comparison of the effector catalogues); no table is therefore duplicated between "
     "the manuscript files and the repository. Data Records now describes the content of "
     "the deposition explicitly, and the assembly statistics and the repeat content "
     "remain in the manuscript itself as Tables 1 and 2."),

    ("H", "REVIEWER 1"),

    ("Reviewer 1", "The only thing I would like to see is a phylogeny of a handful of FOSC "
                   "members and where this isolate sits relative to Fo47 and Fol4287."),
    ("Authors", "We agree that the phylogenetic position of the isolate is necessary "
     "context. That phylogeny has already been published for exactly this set of strains: "
     "in the study that reported the v1.0 draft genome of MRL8996, a maximum-likelihood "
     "tree of 55 conserved single-copy orthologous genes of the genus Fusarium, rooted on "
     "F. verticillioides, places MRL8996 (NRRL 47514) in the same clade as the clinical "
     "isolate NRRL 32931 and within a subclade supported by 100% bootstrap that also "
     "contains Fol4287 (NRRL 34936) and Fo47 (NRRL 54002), the two reference strains used "
     "for comparison here (Zhang et al., Communications Biology 3, 50, 2020, Figure 1). "
     "Rather than reproduce a published tree, we now state this placement explicitly in "
     "the Background & Summary, with the citation, so that the reader has the phylogenetic "
     "framework before the synteny and effector comparisons. If the reviewer and the "
     "editor would prefer the tree itself to be shown, we have computed an independent "
     "maximum-likelihood phylogeny of the same strains from single-copy orthologues of our "
     "own annotation, which recovers the same topology and which we can add as a "
     "supplementary figure."),

    ("H", "REVIEWER 2"),

    ("Reviewer 2", "The motivation for the in-depth description of the predicted secreted "
                   "proteins is unclear; if they are considered important for cross-kingdom "
                   "adaptation or virulence this should be stated explicitly; a more "
                   "detailed description of their genomic location would be valuable and "
                   "Supplementary Figure 2 should be in the main text; and if effectors "
                   "matter for virulence, why are they less prevalent on the accessory "
                   "chromosomes? Perhaps the few accessory effectors are the most "
                   "interesting, and an overview of them in visual form could be "
                   "provided."),
    ("Authors", "We have made the motivation explicit and have given the catalogue the "
     "descriptive treatment that a data descriptor calls for. The rationale is now stated "
     "at the start of the effector section: small secreted proteins are the gene class "
     "most often implicated in the interaction of F. oxysporum with its hosts and in "
     "adaptation to new ones, and they are the part of the annotation least accessible "
     "from a fragmented assembly, because assigning a candidate to the core or the "
     "accessory compartment requires chromosome-level coordinates. The genomic map that "
     "was Supplementary Figure 2 is now a main-text figure (Figure 4A), and the "
     "compartment distribution is quantified in the text: 476 candidates on the core "
     "chromosomes at 10.4 per Mb against 37 on the accessory chromosomes at 5.5 per Mb, "
     "which follows the general gene density of the two compartments (309-366 against "
     "217-277 genes per Mb) rather than being specific to secreted proteins. "
     "To address the question of which candidates are the interesting ones, we have added "
     "a three-strain comparison: the 513 MRL8996 candidates were clustered at 70% "
     "amino-acid identity with the equivalent catalogues of Fol4287 (500) and Fo47 (393). "
     "The 1,406 sequences form 625 clusters, of which 252 are shared by all three strains "
     "and 131 contain only MRL8996 proteins (135 proteins, 26.3% of the catalogue); this "
     "comparison is shown as a new panel B of Figure 4. The isolate-specific candidates "
     "are not concentrated on the accessory chromosomes (124 core, 11 accessory) and are "
     "not enriched in predicted antimicrobial activity, and their uniqueness is at the "
     "level of sequence rather than fold: of the 122 with a confident model, 90 have a "
     "structural counterpart elsewhere in the MRL8996 catalogue at a TM-score of at least "
     "0.5 and 90 have one among the Fol4287 models, leaving 26 with no structural "
     "counterpart in either set. The structural network has been rebuilt with all 482 "
     "nodes retained, including the 131 that share no edge above the threshold, so that "
     "the position of the isolate-specific candidates is visible; 44 of them fall among "
     "those unconnected nodes, and the family most enriched in isolate-specific members is "
     "the KP4-like family. We have deliberately kept these statements descriptive, without "
     "inferring a role in virulence, which this dataset cannot test."),

    ("Reviewer 2", "All of this is available except the gene annotations."),
    ("Authors", "The gene annotation is now part of the deposition described in Data "
     "Records, in GFF3 format on the chromosome-level coordinates, together with the "
     "predicted proteome and the transposable-element annotation."),

    ("Reviewer 2", "In the Abstract, it is unclear what lineage-specific chromosomes are; "
                   "and the five chromosomes that consist mostly of accessory sequence seem "
                   "to be specific to the isolate."),
    ("Authors", "The Abstract no longer uses the term lineage-specific and now refers to "
     "11 core and 5 accessory chromosomes specific to this isolate."),

    ("Reviewer 2", "The claim that the new genome will be the tool advancing our "
                   "understanding of adaptation is not in line with the claim about the "
                   "importance of the accessory genome, which would call for multiple "
                   "reference genomes or a pangenomic approach; it would be fairer to state "
                   "that the main improvement is the high resolution of the accessory "
                   "chromosomes."),
    ("Authors", "We agree and have moderated the claim. The Abstract and the closing "
     "paragraph now state that the main improvement of this assembly is the resolution of "
     "the accessory chromosomes and that it provides a reference for structural and "
     "comparative work, without claiming that a single genome resolves adaptation."),

    ("Reviewer 2", "Subheadings would help in the Background and Summary section."),
    ("Authors", "Two subheadings have been added, 'Clinical relevance of the F. oxysporum "
     "species complex' and 'Genome assembly state of the art'."),

    ("Reviewer 2", "Some references miss context, like PRJNA554890 and NRRL 47514."),
    ("Authors", "Both are now given with their context: the strain is introduced as "
     "MRL8996, deposited as NRRL 47514 in the ARS Culture Collection (Peoria, IL, USA), "
     "and the draft genome is introduced as v1.0, GenBank GCA_009746015.1 under BioProject "
     "PRJNA554890, generated as part of the first comparative-genomic analysis of this "
     "isolate, with the corresponding citation."),

    ("Reviewer 2", "How was StainedGlass used to confirm the detection of centromeric "
                   "regions? It would help to state the sequence characteristics of "
                   "F. oxysporum centromeres."),
    ("Authors", "A dedicated Methods subsection, 'Candidate centromeric regions', has been "
     "added. It states that F. oxysporum centromeres are not defined by a satellite repeat "
     "but are AT-rich, gene-poor regions of a few tens of kb that behave as a single "
     "interaction focus in Hi-C maps, and it describes the three signals that were "
     "intersected to call candidates: GC-content minima in 20-kb windows at least three "
     "percentage points below the chromosome median, the strongest cis-interaction focus "
     "of each chromosome in the Hi-C matrix, and bands of 84-93% self-identity in a "
     "StainedGlass map computed from an all-versus-all alignment in 2-kb windows. Calls "
     "supported by at least two signals were obtained for 10 of the 16 chromosomes and are "
     "reported in a new supplementary table, with the per-chromosome profiles in a new "
     "supplementary image. The wording in Technical Validation has been corrected "
     "accordingly: we no longer state that all centromeres were positioned."),

    ("Reviewer 2", "How can we learn from Figure 1B that the core chromosomes shared across "
                   "the species complex are clearly conserved?"),
    ("Authors", "This was an overstatement of what the circos overview shows, and the "
     "sentence has been corrected. The text now says that the overview shows the "
     "heterogeneous distribution of genomic features and that conservation of the core "
     "compartment is demonstrated by the synteny analysis with Fol4287 and Fo47 in "
     "Figure 3A."),

    ("Reviewer 2", "There is a lot of redundancy around the assembly statistics and their "
                   "generation, in the paragraphs beginning 'The improved MRL8996', 'The "
                   "hybrid assembly' and 'A direct comparison', in Tables 1 and 2 and in "
                   "the Methods; Table 1 could be incorporated into Table 2 and the last "
                   "column of Table 2 is not needed."),
    ("Authors", "Tables 1 and 2 have been merged into a single Table 1 that reports the "
     "sequencing data and the statistics of both assembly versions side by side, and the "
     "redundant last column has been removed; the repeat table is now Table 2 and all "
     "in-text table references have been renumbered. The three paragraphs have been "
     "separated by function: the first now belongs to the Data Overview section and "
     "describes how the data were generated, the second reports the properties of the "
     "assembly and the third only the comparison with v1.0, and the sequencing and "
     "assembly parameters are given once, in the Methods."),

    ("Reviewer 2", "Indicating the centromeric interaction blocks in Figure 1A would help "
                   "to define chromosome boundaries, especially for the accessory "
                   "sequence."),
    ("Authors", "The Figure 1 legend now states that the boxes along the diagonal delimit "
     "the 16 pseudomolecules and that the candidate centromeric region called for each "
     "chromosome coincides with the strongest cis-interaction focus within the "
     "corresponding block, with the coordinates given in the new supplementary table so "
     "that the reader can locate each of them on the map."),

    ("Reviewer 2", "The Fol4287 accessory chromosomes are described as small, which is not "
                   "the case relative to the non-accessory chromosomes; and is there "
                   "literature about differing numbers of accessory chromosomes among "
                   "relatives?"),
    ("Authors", "The description has been corrected: the partial alignments of the MRL8996 "
     "accessory chromosomes are now said to be confined to the accessory chromosomes 03, "
     "06, 14 and 15 of Fol4287 and to the accessory Chr07 of Fo47, without qualifying them "
     "as small. The Background & Summary now also notes that the number and content of the "
     "accessory compartment tend to be similar within a plant-infecting forma specialis "
     "but vary between isolates, with the corresponding references."),

    ("Reviewer 2", "'V. dahliae' = Verticillium dahliae."),
    ("Authors", "Corrected; the genus is spelled out at first mention and the species name "
     "is now correct throughout."),

    ("Reviewer 2", "In Methods it is unclear what the UHMW DNA was used for, and how longer "
                   "molecules were extracted than the HMW DNA in Chavarro-Carrero."),
    ("Authors", "We apologise for the confusion, which came from describing two DNA "
     "preparations as if they were one. The Methods now state that two preparations were "
     "made from the same culture protocol: the first followed the published protocol "
     "without modification and provided the DNA for the run deposited under SRR39391754, "
     "from which the assembly reported here was built; in the second, a size-selection "
     "step with the Short Read Eliminator kit was added after DNA precipitation to deplete "
     "fragments below 25 kb, and this preparation was used only to test whether the "
     "read-length distribution could be shifted further and contributed no reads to the "
     "assembly. We also make clear that the size selection removes short fragments and "
     "does not extract molecules longer than the published protocol yields."),

    ("Reviewer 2", "The description of the Hi-C sequencing should be moved from the "
                   "'assembly' to the 'sequencing' paragraph."),
    ("Authors", "The Hi-C library construction and sequencing are now described at the end "
     "of the sequencing paragraph, and the sentence has been removed from the assembly "
     "paragraph."),

    ("Reviewer 2", "Data Records and Data Availability are redundant."),
    ("Authors", "Data Records now carries the full description of the deposited datasets "
     "with their accessions and data citations, and the Data availability statement has "
     "been reduced to a statement that all data listed in Data Records are publicly "
     "available without restriction and that no permissions or material transfer "
     "agreements are required for their reuse."),

    ("Reviewer 2", "Phase Genomics Inc. should be moved from Author Contributions to "
                   "Acknowledgments, or included as an author."),
    ("Authors", "Phase Genomics Inc. has been removed from the Author Contributions and is "
     "acknowledged for the Hi-C library preparation and sequencing."),

    ("Reviewer 2", "Figures are of very low resolution and some text is hardly readable."),
    ("Authors", "All figures have been re-exported from the original vector sources at 600 "
     "dpi, and the panel labels, axis text and legend text have been enlarged so that no "
     "element falls below 8 pt at the intended print size. Vector versions are supplied "
     "alongside the raster files."),

    ("Reviewer 2", "Figure 3: how was synteny assessed in 3A, and intra-genomic synteny is "
                   "probably the wrong term for 3B."),
    ("Authors", "The Methods now describe the procedure: whole-genome alignment with "
     "minimap2 in asm20 mode, discarding alignment blocks shorter than 10 kb or below 80% "
     "identity and merging the remaining blocks into syntenic blocks on the MRL8996 "
     "coordinate system. Panel B is now described as segmental duplications within the "
     "MRL8996 genome, obtained by aligning the assembly against itself with the same "
     "settings and discarding the diagonal, and the term intra-genomic synteny has been "
     "removed."),

    ("Reviewer 2", "Table 2: what are the auNs of the assemblies?"),
    ("Authors", "The auN has been added to the merged assembly table for both versions, "
     "1,874,049 bp for v1.0 and 4,154,107 bp for v2.0, and is quoted in the text, where we "
     "note that it summarises contiguity over the whole size distribution rather than at a "
     "single quantile and is therefore less sensitive than N50 to the choice of "
     "quantile."),

    ("Reviewer 2", "Table 3: a percentage for each class and its occurrence in the "
                   "accessory genome would be helpful; underscores in Tables 2 and 3 make "
                   "them hard to read."),
    ("Authors", "The repeat table has been rebuilt as a single integrated table that gives, "
     "for every class and superfamily, the count and the percentage of sequence occupied "
     "in the core compartment, in the accessory compartment and genome-wide, so that the "
     "enrichment in the accessory compartment can be read directly. All underscores have "
     "been replaced by readable class and superfamily names in both tables."),

    ("Reviewer 2", "Supplementary Image 1: a legend for the read coverage is missing; why "
                   "are there regions with seemingly very high coverage - were regions "
                   "wrongly collapsed during assembly?"),
    ("Authors", "The legend now gives the coverage key: bar colour is the Nanopore read "
     "depth in 200-kb windows as a multiple of the genome-median depth of 47.5-fold, on "
     "the six-level scale used by tapestry. We also examined the three regions that exceed "
     "one genomic copy. Two are collapsed tandem repeats of known composition rather than "
     "misassemblies of unique sequence: the terminal 70 kb of Chr05, which carries eight "
     "tandem rDNA units in the assembly against an estimated 124 units in the reads, and a "
     "4.7-kb mitochondrial insertion at Chr05:4.55 Mb. The third, the first contig of "
     "Chr13 (0-0.49 Mb), sits at about 1.8-fold uniquely-mapping depth and is the one "
     "segment that appears to represent two genomic copies. All three are now described in "
     "the legend."),

    ("Reviewer 2", "General comment: a description of the TEs and genes on the accessory "
                   "chromosomes is missing."),
    ("Authors", "A paragraph on the content of the two compartments has been added to the "
     "Data Overview. The core chromosomes carry 309-366 genes per Mb, 45-50% coding "
     "sequence and 3.5-10.5% repeat-masked sequence and contain all 56 predicted "
     "biosynthetic gene clusters, whereas the accessory chromosomes carry 217-277 genes "
     "per Mb, 21-28% coding sequence, 23.7-50.3% repeat-masked sequence and no cluster. "
     "The five accessory chromosomes span 6.73 Mb, 12.9% of the assembly, and account for "
     "1,688 of the 16,612 annotated genes. Repeats occupy 6.0% of the core and 34.5% of "
     "the accessory compartment, every class detected being proportionally more abundant "
     "in the accessory compartment, where the largest contributions are unclassified "
     "elements, LTR elements, LINEs and DNA transposons. A per-chromosome table of gene, "
     "coding, repeat and cluster content has been added to the supplementary material."),
]


def main():
    doc = docx.Document()
    st = doc.styles["Normal"]
    st.font.name = "Calibri"
    st.font.size = Pt(11)

    for t in HEAD:
        doc.add_paragraph(t)

    for kind, text in BLOCKS:
        p = doc.add_paragraph()
        if kind == "H":
            r = p.add_run(text)
            r.bold = True
        else:
            r = p.add_run(kind + ": ")
            r.bold = True
            p.add_run(text)

    for t in TAIL:
        doc.add_paragraph(t)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    doc.save(str(OUT))
    print("written", OUT, "| blocks:", len(BLOCKS))


if __name__ == "__main__":
    main()
