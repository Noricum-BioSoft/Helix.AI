import PyPDF2

from pathlib import Path
from typing import Any, Optional
# from smolagents.tools import Tool


class FileReaderTool:

    name = "file_reader"
    description = "reads and analyzes the content of a given file."
    inputs = {
        "query": {
            "type": "string",
            "description": "The search query to perform.",
        }
    }
    output_type = "string"

    def __init__(self):
        super().__init__()

    def read_file(self, file_path: Path) -> str:
        with open(file_path, 'rb') as f:
            pdf = PyPDF2.PdfReader(f)
            text = ''
            for page in pdf.pages:
                text += page.extract_text()

            return text

    def forward(self, query: str) -> str:
        print(query)

        return query
        # file_content = self.read_file(file_path)
        # return file_content[:50]

