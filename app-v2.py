import os
import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import io
import base64
import json

from pathlib import Path
from smolagents import InferenceClientModel

from tools.data_science import DataAnalysisAgent
from tools.data_science import *
from tools.final_answer import FinalAnswerTool
from tools.bio import align_and_visualize_fasta

# ---- Load ENVs ----
from dotenv import load_dotenv
load_dotenv()


# Initialize agent
# model_id = "Qwen2/5-72B-Instruct"
# agent = CodeAgent(
#     name="BioBloomAgent",
#     tools=[align_sequences],
#     model=model_id
# )
# llm = OpenAI(temperature=0.7)

from smolagents import OpenAIServerModel

model = OpenAIServerModel(
    model_id="gpt-4o-mini",
    api_key=os.environ["OPENAI_API_KEY"],
)

# model = InferenceClientModel(
#     provider="hf-inference",
#     api_key=os.environ["HF_TOKEN"]
# )


def load_image(path: Path) -> str:
    with open(path, "rb") as img_file:
        encoded = base64.b64encode(img_file.read()).decode()
    return f"<img src='data:image/png;base64,{encoded}' width='40' style='vertical-align:middle; margin-right:10px;'>"


def main():

    st.set_page_config(page_title="DataBloom.AI - Prompt UI", layout="wide")

    # Initialize session state
    if 'data' not in st.session_state:
        st.session_state['data'] = None
    if 'agent' not in st.session_state:
        st.session_state['agent'] = None

    logo_path = Path(Path(__file__).parent, "images", "logo.png")
    logo_html = load_image(logo_path)
    st.markdown(f"""
    <div style="display:flex; align-items:center; gap:10px;">
        {logo_html}
        <h1 style='display:inline; vertical-align:middle;'>DataBloom.AI — Prompt-Based Data Management</h1>
    </div>
    """, unsafe_allow_html=True)

    # File upload
    uploaded_file = st.file_uploader("Upload data file", type="csv")

    # Prompt input
    user_prompt = st.text_area("Enter your prompt", placeholder="e.g., Visualize protein expression trends over time")

    # Submit button
    st.button("Submit")

    if uploaded_file is not None:
        with st.spinner('Loading and processing your data...'):
            # Load the dataset
            data = pd.read_csv(uploaded_file)
            st.session_state['data'] = data

            # Initialize the agent with the dataset
            # and store it in the session
            st.session_state['agent'] = DataAnalysisAgent(
                dataset=data,
                model=model,
                tools=[analyze_basic_stats, generate_correlation_matrix,
                       analyze_categorical_columns, suggest_features, FinalAnswerTool(),
                       align_and_visualize_fasta],
                additional_authorized_imports=["pandas", "numpy", "numpy.random", "matplotlib", "matplotlib.plt", "seaborn", "io", "io.StringIO", "pypandoc"],
                # stream_outputs=True,
            )

            # st.success(f'Successfully loaded dataset with {data.shape[0]} rows and {data.shape[1]} columns')
            st.subheader("Data Preview")
            st.dataframe(data.head())

    # Submit button
    if st.session_state['data'] is not None and st.session_state['agent'] is not None and user_prompt:
        with st.spinner('Generating and executing code to analyze and visualize this data...'):

            # Load CSV into DataFrame

            df = pd.DataFrame(st.session_state['data'])
            data_string = df.to_csv(index=False)

            # Create prompt for agent
            full_prompt = f"""
        You are given a CSV file containing the following data:
        
        {data_string}
        
        Instruction: {user_prompt}
        
        Only use the return value of the tool. Do not add any additional text or information.
                """

            # Run the agent
            agent = st.session_state['agent']
            results = agent.run(full_prompt)
            print(results)
            st.image(results)

            # json_results = results
            # if isinstance(json_results, str):
            #     json_results = json.loads(results)
            #
            # df = pd.DataFrame(json_results)
            # markdown_table = df.to_markdown(index=False)
            # st.write(markdown_table)

            # # Display generated code (optional for transparency)
            # with st.expander("🔍 Show generated code"):
            #     st.code(code, language="python")
            #
            # try:
            #     # Execute generated code
            #     local_vars = {'pd': pd, 'plt': plt, 'df': df}
            #     exec(code, {}, local_vars)
            #
            #     # Capture the plot
            #     buf = io.BytesIO()
            #     plt.savefig(buf, format='png')
            #     buf.seek(0)
            #     plt.close()
            #
            #     # Display result
            #     st.success("✅ Visualization generated successfully!")
            #     st.image(buf)
            #
            # except Exception as e:
            #     st.error(f"❌ Error: {str(e)}")


if __name__ == "__main__":
    main()