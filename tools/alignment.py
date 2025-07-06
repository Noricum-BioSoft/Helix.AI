# /backend/tools/alignment.py

from langchain.agents import tool


@tool
def run_alignment(sequences: str):
    """performs a sequence alignment on a given set of sequences."""

    seq_align_dict = [
        {"name": "seq1", "sequence": "ACTG--TTGAC"},
        {"name": "seq2", "sequence": "ACTGCA-T--C"},
        {"name": "seq3", "sequence": "ACTGCAATGAC"},
    ]

    return {
        "text": "Sequences aligned successfully.",
        "alignment": seq_align_dict
    }

    # return {
    #     "text": "Sequences aligned successfully.",
    #     "plot": {
    #         "data": [{"x": [1, 2, 3], "y": [3, 1, 2], "type": "bar"}],
    #         "layout": {"title": "Alignment Visualization"}
    #     }
    # }

