from smolagents import (
    CodeAgent,
    DuckDuckGoSearchTool,
    InferenceClientModel,
    load_tool,
    tool,
)
import datetime
import requests
import pytz
import yaml
from tools.final_answer import FinalAnswerTool
from tools.visit_webpage import VisitWebpageTool
from tools.web_search import DuckDuckGoSearchTool
from tools.file_reader import FileReaderTool

from Gradio_UI import GradioUI

final_answer = FinalAnswerTool()
visit_webpage = VisitWebpageTool()
web_search_tool = DuckDuckGoSearchTool()
file_reader_tool = FileReaderTool()


from smolagents import CodeAgent, InferenceClientModel

model_id = "Qwen2/5-72B-Instruct"

agent = CodeAgent(
    tools=[final_answer, visit_webpage, web_search_tool],
    model=InferenceClientModel(model_id=model_id),
    add_base_tools=True,
)

GradioUI(agent, file_upload_folder=".").launch()
