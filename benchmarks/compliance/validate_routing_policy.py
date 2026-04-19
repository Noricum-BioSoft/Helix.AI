from __future__ import annotations

import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
ROUTER_PATH = REPO_ROOT / "backend" / "command_router.py"


def validate() -> int:
    text = ROUTER_PATH.read_text()

    # Guard against reintroducing the legacy broad fallback:
    # any(... ['align', 'alignment', 'sequences']) -> sequence_alignment
    forbidden_snippets = [
        "if any(phrase in command_lower for phrase in ['align', 'alignment', 'sequences']):",
        'if any(phrase in command_lower for phrase in ["align", "alignment", "sequences"]):',
    ]
    for snippet in forbidden_snippets:
        if snippet in text:
            print("Routing policy validation failed: broad alignment fallback pattern was reintroduced.")
            print(f"Forbidden snippet: {snippet}")
            return 1

    required_snippet = "def _should_sequence_alignment_fallback("
    if required_snippet not in text:
        print("Routing policy validation failed: conservative fallback helper is missing.")
        return 1

    print("Routing policy validation passed.")
    return 0


if __name__ == "__main__":
    sys.exit(validate())
