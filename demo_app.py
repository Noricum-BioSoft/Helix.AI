import gradio as gr
from PIL import Image, ImageDraw
from typing import Tuple, Optional


def align_sequences(data: str) -> Tuple[str, Image.Image]:
    """
    Simulates sequence alignment and returns a textual and visual representation.

    Parameters
    ----------
    data : str
        User-provided input that may include sequences or an alignment instruction.

    Returns
    -------
    Tuple[str, Image.Image]
        A tuple containing the alignment result as a string and an image visualizing the alignment.
    """
    # Placeholder alignment result
    text_output = (
        "Aligned sequences:\n"
        ">Seq1\nATGCTAGCTA\n"
        ">Seq2\nATG-TAG-TA\n"
        f"Input: {data}"
    )

    # Generate a placeholder alignment image
    img = Image.new("RGB", (400, 100), color="white")
    draw = ImageDraw.Draw(img)
    draw.text((10, 10), "Sequence Alignment Preview", fill="black")

    return text_output, img


def visualize_structure(data: str) -> str:
    """
    Simulates 3D protein structure visualization.

    Parameters
    ----------
    data : str
        User-provided input referring to a protein or gene.

    Returns
    -------
    str
        Placeholder text indicating structure visualization.
    """
    return f"3D structure visualization for: {data}"


def codon_optimize(data: str, organism: str = "E. coli") -> str:
    """
    Simulates codon optimization for a target organism.

    Parameters
    ----------
    data : str
        User-provided gene sequence or gene name.
    organism : str, optional
        Target organism for codon optimization (default is 'E. coli').

    Returns
    -------
    str
        Optimized DNA sequence (placeholder).
    """
    return f"Codon-optimized sequence for {organism}:\n[placeholder output]\nInput: {data}"


def clone_into_vector(data: str, vector: str = "pTet") -> str:
    """
    Simulates cloning of a gene into a given vector.

    Parameters
    ----------
    data : str
        User-provided gene sequence or gene name.
    vector : str, optional
        Cloning vector name (default is 'pTet').

    Returns
    -------
    str
        Cloning plan or output (placeholder).
    """
    return f"Cloned into {vector} vector with annotated restriction sites.\nInput: {data}"


class BioAgent:
    """
    Simple routing agent for handling bioinformatics instructions.
    """

    def __init__(self):
        self.tools = {
            "alignment": align_sequences,
            "structure": visualize_structure,
            "optimization": codon_optimize,
            "cloning": clone_into_vector,
        }

    def route_task(self, prompt: str) -> Tuple[str, Optional[Image.Image]]:
        """
        Determines which tool to use based on user prompt.

        Parameters
        ----------
        prompt : str
            The natural language command from the user.

        Returns
        -------
        Tuple[str, Optional[Image.Image]]
            The result of the operation as text, and optionally an image.
        """
        prompt_lower = prompt.lower()

        if "align" in prompt_lower:
            return self.tools["alignment"](prompt)
        elif "structure" in prompt_lower or "3d" in prompt_lower:
            return self.tools["structure"](prompt), None
        elif "codon optimize" in prompt_lower or "optimize" in prompt_lower:
            return self.tools["optimization"](prompt), None
        elif "clone" in prompt_lower:
            return self.tools["cloning"](prompt), None
        else:
            return "❌ Could not understand the task. Try prompts like 'align', 'optimize', 'clone', or 'structure'.", None


# Instantiate the agent
agent = BioAgent()


def handle_prompt(prompt: str) -> Tuple[str, Optional[Image.Image]]:
    """
    Handles user prompt by routing to the appropriate bioinformatics tool.

    Parameters
    ----------
    prompt : str
        The user's natural language instruction.

    Returns
    -------
    Tuple[str, Optional[Image.Image]]
        The result of the processing as text and optional image.
    """
    return agent.route_task(prompt)


# --- Gradio Interface ---
with gr.Blocks() as demo:
    gr.Markdown("## 🧬 BioPrompt — Natural Language Bioinformatics Platform")

    with gr.Row():
        user_input = gr.Textbox(label="🔍 Prompt", placeholder="e.g., Codon optimize lacZ for E. coli...", lines=2)

    with gr.Row():
        text_output = gr.Textbox(label="📝 Agent Output", lines=10)

    with gr.Row():
        image_output = gr.Image(label="🖼️ Visualization", type="pil")

    user_input.submit(handle_prompt, inputs=user_input, outputs=[text_output, image_output])

# Launch the app
demo.launch()
