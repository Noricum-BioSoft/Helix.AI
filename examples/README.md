# Helix examples

Caller-side demos that treat Helix as an **execution plane** (API), not a chat UI.

| File | Purpose |
|------|---------|
| [`labos_client.py`](./labos_client.py) | Thin “LabOS” client: session → intent → approve → runs → bundle |

```bash
export HELIX_BASE_URL=http://localhost:8001
python examples/labos_client.py health
python examples/labos_client.py plan    # stop at the plan approval gate
python examples/labos_client.py run     # full path + bundle
```
