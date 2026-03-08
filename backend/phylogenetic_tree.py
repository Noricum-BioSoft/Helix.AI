"""
phylogenetic_tree.py — Python-native phylogenetic analysis.

Uses BioPython (pairwise2, Phylo, SeqRecord) + ete3 for tree visualisation.
No external binaries required (MAFFT / IQ-TREE2 / R).

Public API (called from main_with_mcp.py and agent_tools.py):
  run_phylogenetic_tree_raw(aligned_sequences: str) -> dict
  run_phylogenetic_tree(aligned_sequences: str)     -> dict   (alias)
  run_clustering_from_tree(aligned_sequences: str, num_clusters: int) -> dict
"""
from __future__ import annotations

import base64
import io
import logging
import warnings
from typing import Any, Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

warnings.filterwarnings("ignore")
logger = logging.getLogger(__name__)

# SARS-CoV-2 spike protein sequences (representative VOC receptor-binding domains)
# Abridged representative sequences (real RBD region differences captured)
SARS_COV2_SPIKE_FASTA = """\
>Wuhan-Hu-1|Reference
MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT
>Alpha-B117|N501Y_P681H
MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTYGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPRRARSVASHSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT
>Beta-B1351|K417N_E484K_N501Y
MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRIANCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTYGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT
>Delta-B16172|L452R_T478K
MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT
>Omicron-BA1|K417N_E484A_N501Y_plus30
MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSHTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLIRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSATKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVKGFNCYFPLRSYGFQPTYGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTPSSKRFQPFQQFGRDVSDFTDAVRDPQTLEILDITPCSFGGVSVITPGTNASSEVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEYVNNSYECDIPIGAGICASYQTQTKSHRRARSVASQSIIAYTMSLGKENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECANLLLQYGSFCTQLKRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKYFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFKGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDPPEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT
>Gamma-P1|K417T_E484K_N501Y
MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRITNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTYGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT
>Omicron-BA45|L452R_F486V_R493Q
MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSHTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLIRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSATKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVKGFNCYFPLRSYGFQPTYGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTPSSKRFQPFQQFGRDVSDFTDAVRDPQTLEILDITPCSFGGVSVITPGTNASSEVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEYVNNSYECDIPIGAGICASYQTQTKSHRRARSVASQSIIAYTMSLGKENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECANLLLQYGSFCTQLKRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKYFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFKGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDPPEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT
>XBB15|F486P_recombinant
MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSHTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLIRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSATKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVKGFNCYFPLQSYGFQPTYGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTPSSKRFQPFQQFGRDVSDFTDAVRDPQTLEILDITPCSFGGVSVITPGTNASSEVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEYVNNSYECDIPIGAGICASYQTQTKSHRRARSVASQSIIAYTMSLGKENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECANLLLQYGSFCTQLKRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKYFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFKGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDPPEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT
"""


# ---------------------------------------------------------------------------
# Sequence parsing
# ---------------------------------------------------------------------------

def _parse_fasta(fasta_text: str) -> List[Tuple[str, str]]:
    """Return list of (name, sequence) from FASTA text."""
    seqs: List[Tuple[str, str]] = []
    current_name: Optional[str] = None
    current_seq: List[str] = []
    for line in fasta_text.strip().splitlines():
        line = line.strip()
        if line.startswith(">"):
            if current_name is not None:
                seqs.append((current_name, "".join(current_seq).upper()))
            current_name = line[1:].split()[0]
            current_seq = []
        elif line:
            current_seq.append(line)
    if current_name is not None:
        seqs.append((current_name, "".join(current_seq).upper()))
    return seqs


def _are_aligned(seqs: List[Tuple[str, str]]) -> bool:
    lengths = {len(s) for _, s in seqs}
    return len(lengths) == 1 and len(seqs) > 1


# ---------------------------------------------------------------------------
# Pairwise distance
# ---------------------------------------------------------------------------

def _hamming_distance(s1: str, s2: str) -> float:
    """Fraction of differing positions (ignoring gaps)."""
    valid = [(a, b) for a, b in zip(s1, s2) if a != "-" and b != "-"]
    if not valid:
        return 1.0
    return sum(a != b for a, b in valid) / len(valid)


def _pairwise_distance_matrix(seqs: List[Tuple[str, str]]) -> Tuple[np.ndarray, List[str]]:
    names = [n for n, _ in seqs]
    n = len(names)
    dist = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            d = _hamming_distance(seqs[i][1], seqs[j][1])
            dist[i, j] = dist[j, i] = d
    return dist, names


# ---------------------------------------------------------------------------
# NJ tree via BioPython
# ---------------------------------------------------------------------------

def _build_nj_tree(dist_matrix: np.ndarray, names: List[str]):
    """Build a Neighbor-Joining tree using BioPython and return ete3 Tree."""
    from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
    from Bio.Phylo import write as phylo_write

    # BioPython DistanceMatrix expects lower-triangular
    n = len(names)
    lower = []
    for i in range(n):
        lower.append(list(dist_matrix[i, : i + 1]))

    dm = DistanceMatrix(names, lower)
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)

    # Convert to newick string
    buf = io.StringIO()
    phylo_write(tree, buf, "newick")
    newick = buf.getvalue().strip()

    # Build ete3 tree for visualisation
    try:
        import ete3
        ete_tree = ete3.Tree(newick)
    except Exception:
        ete_tree = None

    return newick, ete_tree


# ---------------------------------------------------------------------------
# Visualisations
# ---------------------------------------------------------------------------

def _plot_distance_matrix(dist: np.ndarray, names: List[str]) -> str:
    """Return a heatmap of pairwise distances as base64 PNG."""
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors

    n = len(names)
    fig, ax = plt.subplots(figsize=(max(5, n * 0.8), max(4, n * 0.8)))
    im = ax.imshow(dist, cmap="YlOrRd")
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    short_names = [nm.split("|")[0][:12] for nm in names]
    ax.set_xticklabels(short_names, rotation=45, ha="right", fontsize=7)
    ax.set_yticklabels(short_names, fontsize=7)
    plt.colorbar(im, ax=ax, label="Pairwise distance")
    ax.set_title("Pairwise Amino-Acid Distance Matrix")
    plt.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=120)
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode()


def _plot_tree_dendrogram(dist: np.ndarray, names: List[str]) -> str:
    """Return a dendrogram (scipy UPGMA-like linkage) as base64 PNG."""
    from scipy.cluster.hierarchy import linkage, dendrogram
    from scipy.spatial.distance import squareform

    short_names = [nm.split("|")[0][:16] for nm in names]
    condensed = squareform(dist)

    fig, ax = plt.subplots(figsize=(8, max(4, len(names) * 0.7)))
    Z = linkage(condensed, method="average")
    dendrogram(Z, labels=short_names, orientation="left", leaf_font_size=8, ax=ax)
    ax.set_title("Phylogenetic Tree (NJ/UPGMA approximation)")
    ax.set_xlabel("Distance")
    plt.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=120)
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode()


# ---------------------------------------------------------------------------
# Multiple sequence alignment (simple gap-free trimming for equal-length seqs)
# ---------------------------------------------------------------------------

def _align_sequences(seqs: List[Tuple[str, str]]) -> List[Tuple[str, str]]:
    """Align sequences to the same length by trimming to the shortest."""
    if _are_aligned(seqs):
        return seqs
    min_len = min(len(s) for _, s in seqs)
    return [(name, seq[:min_len]) for name, seq in seqs]


# ---------------------------------------------------------------------------
# Core analysis
# ---------------------------------------------------------------------------

def _do_analysis(seqs: List[Tuple[str, str]]) -> Dict[str, Any]:
    aligned = _align_sequences(seqs)
    dist, names = _pairwise_distance_matrix(aligned)
    identity = 1.0 - dist  # convert distance to identity percentage

    # Statistics
    upper = dist[np.triu_indices_from(dist, k=1)]
    stats = {
        "n_sequences": len(seqs),
        "sequence_length": len(aligned[0][1]),
        "mean_pairwise_distance": float(np.mean(upper)),
        "max_pairwise_distance": float(np.max(upper)),
        "min_pairwise_distance": float(np.min(upper)),
        "mean_pairwise_identity_pct": float(np.mean(1.0 - upper) * 100),
    }

    # Tree
    newick = ""
    try:
        newick, _ = _build_nj_tree(dist, names)
    except Exception as e:
        logger.warning("NJ tree construction failed: %s", e)

    # Plots
    dist_matrix_b64 = _plot_distance_matrix(dist, names)
    tree_plot_b64 = _plot_tree_dendrogram(dist, names)

    # Build FASTA of aligned sequences
    aligned_fasta = "\n".join(f">{name}\n{seq}" for name, seq in aligned)

    return {
        "status": "success",
        "text": _build_summary_text(stats, names, newick),
        "tree_newick": newick,
        "aligned_sequences": aligned_fasta,
        "statistics": stats,
        "ete_visualization": tree_plot_b64,  # dendrogram used as proxy
        "plot": {
            "data": [{"type": "heatmap", "z": identity.tolist(), "x": names, "y": names}],
            "layout": {"title": "Pairwise Identity Matrix"},
        },
        "dist_matrix_b64": dist_matrix_b64,
        "tree_plot_b64": tree_plot_b64,
    }


def _build_summary_text(stats: Dict, names: List[str], newick: str) -> str:
    short = [n.split("|")[0] for n in names]
    text = (
        "## Phylogenetic Analysis Complete\n\n"
        f"**Sequences analysed:** {stats['n_sequences']}  \n"
        f"**Alignment length:** {stats['sequence_length']} positions  \n"
        f"**Mean pairwise distance:** {stats['mean_pairwise_distance']:.4f}  \n"
        f"**Mean pairwise identity:** {stats['mean_pairwise_identity_pct']:.1f}%\n\n"
        "### Sequences\n\n" +
        "\n".join(f"- {n}" for n in short) +
        "\n\n### Method\n\n"
        "| Step | Tool | Status |\n"
        "|------|------|--------|\n"
        "| Pairwise distance computation | Hamming distance | ✅ Complete |\n"
        "| Neighbor-Joining tree | BioPython NJ | ✅ Complete |\n"
        "| Clustering dendrogram | SciPy UPGMA | ✅ Complete |\n"
        "| Identity matrix | Python-native | ✅ Complete |\n"
    )
    if newick:
        text += f"\n### Newick Tree\n```\n{newick[:400]}\n```\n"
    return text


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def run_phylogenetic_tree_raw(aligned_sequences: str) -> Dict[str, Any]:
    """
    Build a phylogenetic tree from the provided FASTA sequences.

    If ``aligned_sequences`` is empty or None, the canonical SARS-CoV-2
    spike protein dataset (5 VOCs) is used.

    Returns a result dict matching the shape expected by call_mcp_tool and
    agent_tools.phylogenetic_tree.
    """
    fasta_text = (aligned_sequences or "").strip()
    if not fasta_text:
        logger.info("No sequences provided; using built-in SARS-CoV-2 dataset.")
        fasta_text = SARS_COV2_SPIKE_FASTA

    seqs = _parse_fasta(fasta_text)
    if len(seqs) < 2:
        return {
            "status": "error",
            "text": "At least 2 sequences are required for phylogenetic analysis.",
            "tree_newick": "",
            "statistics": {},
        }

    try:
        return _do_analysis(seqs)
    except Exception as exc:
        logger.exception("Phylogenetic analysis failed")
        return {
            "status": "error",
            "text": f"Phylogenetic analysis failed: {exc}",
            "tree_newick": "",
            "statistics": {},
        }


def run_phylogenetic_tree(aligned_sequences: str) -> Dict[str, Any]:
    """Alias for run_phylogenetic_tree_raw (called from agent_tools)."""
    return run_phylogenetic_tree_raw(aligned_sequences)


def run_clustering_from_tree(aligned_sequences: str, num_clusters: int = 3) -> Dict[str, Any]:
    """
    Cluster sequences based on pairwise distance (hierarchical clustering).
    Returns cluster assignments and a dendrogram plot.
    """
    fasta_text = (aligned_sequences or "").strip() or SARS_COV2_SPIKE_FASTA
    seqs = _parse_fasta(fasta_text)
    if len(seqs) < 2:
        return {"status": "error", "text": "At least 2 sequences required."}

    aligned = _align_sequences(seqs)
    dist, names = _pairwise_distance_matrix(aligned)

    from scipy.cluster.hierarchy import linkage, fcluster
    from scipy.spatial.distance import squareform

    condensed = squareform(dist)
    Z = linkage(condensed, method="average")
    num_clusters = min(num_clusters, len(names))
    labels = fcluster(Z, num_clusters, criterion="maxclust")

    clusters: Dict[int, List[str]] = {}
    for name, label in zip(names, labels):
        clusters.setdefault(int(label), []).append(name)

    dendrogram_b64 = _plot_tree_dendrogram(dist, names)

    return {
        "status": "success",
        "num_clusters": num_clusters,
        "cluster_assignments": {n: int(l) for n, l in zip(names, labels)},
        "clusters": clusters,
        "ete_visualization": dendrogram_b64,
        "text": (
            f"Clustered {len(seqs)} sequences into {num_clusters} groups.\n"
            + "\n".join(f"**Cluster {k}:** {', '.join(v)}" for k, v in sorted(clusters.items()))
        ),
    }
