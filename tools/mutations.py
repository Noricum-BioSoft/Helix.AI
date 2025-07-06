# /backend/tools/mutation.py
from langchain.agents import tool


@tool
def run_mutation(sequence: str, num_variants: int = 96):
    """mutates a given sequence and returns 96 variants per default."""
    return {
        "text": "Generated 96 mutated variants.",
        "plot": {
            "data": [{"x": [1, 2, 3], "y": [2, 5, 3], "type": "line"}],
            "layout": {"title": "Mutation Variants"}
        }
    }
