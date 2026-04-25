# Environment Setup for Helix.AI

## Required Environment Variables

To run Helix.AI, you need to set up the following environment variables:

### 1. Create a `.env` file

Create a `.env` file in the project root directory with the following content:

```bash
# At least one LLM API key is required.
OPENAI_API_KEY=your_openai_api_key_here

# Optional: use DeepSeek as the primary LLM instead of OpenAI.
DEEPSEEK_API_KEY=your_deepseek_api_key_here
```

### 2. Get API Keys

#### OpenAI API Key
1. Go to [OpenAI API](https://platform.openai.com/api-keys)
2. Sign in or create an account
3. Create a new API key
4. Copy the key and add it to your `.env` file

#### DeepSeek API Key (Optional)
1. Go to [DeepSeek Platform](https://platform.deepseek.com/)
2. Sign in or create an account
3. Navigate to API keys section
4. Create a new API key
5. Copy the key and add it to your `.env` file

### 3. LLM selection

- `DEEPSEEK_API_KEY` set → uses DeepSeek as the primary LLM
- Only `OPENAI_API_KEY` set → uses OpenAI
- Neither set → application raises a clear error at startup

**Note:** Helix requires at least one LLM to be available for all classification
decisions (intent, routing, approval, staging). There are no keyword-matching
fallbacks — if the LLM is unavailable the classifier raises explicitly.

### 4. Security Notes

- Never commit your `.env` file to version control
- The `.env` file is already in `.gitignore`
- Keep your API keys secure and don't share them

### 5. Alternative: Set Environment Variables Directly

You can also set environment variables directly in your shell:

```bash
export OPENAI_API_KEY="your_openai_api_key_here"
export DEEPSEEK_API_KEY="your_deepseek_api_key_here"
```

## Development / test variables

| Variable | Purpose |
|----------|---------|
| `HELIX_MOCK_MODE` | Set `1` to disable all LLM calls. Used by the unit-test suite (`pytest.ini` sets this automatically). Classifiers raise `*ClassificationError` rather than silently returning defaults. |
| `HELIX_LLM_ROUTER_FIRST` | Default `1`. Set `0` to enable keyword-based routing branches for local debugging of specific routing paths. Not for production use. |
| `HELIX_AGENT_DISABLED` | Set `1` to skip the BioAgent path entirely (forces fallback router). |
| `HELIX_SANDBOX_HOST_FALLBACK` | Set `1` to skip Docker availability checks and run sandboxed scripts on the host. For local dev only. |

## Session storage (local vs S3)

| Variable | Purpose |
|----------|---------|
| `HELIX_LOCAL_SESSIONS_ONLY` | Set to `1` or `true` to **disable S3** for session initialization (EC2 beta with **only local disk**). Sessions are stored under `HELIX_SESSIONS_DIR` or `./sessions`. |
| `HELIX_SESSIONS_DIR` | Absolute path for session files (e.g. `/var/lib/helix/sessions`). Also accepts `HELIX_STORAGE_DIR` as an alias. |
| `S3_BUCKET_NAME` | Used for S3 session markers / cloud paths when **not** in local-only mode. |

## Troubleshooting

### Error: "No LLM API key configured"
Make sure you have at least one of the following set in your `.env` file:
- `OPENAI_API_KEY` (recommended)
- `DEEPSEEK_API_KEY`

### Error: classifier raises in tests
The unit-test suite runs with `HELIX_MOCK_MODE=1` (set in `pytest.ini`). If a test
triggers a live LLM call, it will raise rather than silently return a default. Use the
autouse fixtures in `tests/unit/backend/conftest.py` or add a `@patch` decorator to
mock the specific classifier your test exercises.

