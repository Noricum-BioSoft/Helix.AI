# Environment Setup for Helix.AI

## Required Environment Variables

To run Helix.AI, you need to set up the following environment variables:

### 1. Create a `.env` file

Create a `.env` file in the project root directory with the following content:

```bash
# OpenAI API Key (required for fallback and some features)
OPENAI_API_KEY=your_openai_api_key_here

# DeepSeek API Key (optional - if not provided, will fallback to OpenAI)
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

### 3. Environment Variable Priority

- If `DEEPSEEK_API_KEY` is set: Uses DeepSeek as the primary LLM
- If `DEEPSEEK_API_KEY` is not set: Falls back to OpenAI
- If neither is set: The application will show an error

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

## Troubleshooting

### Error: "DEEPSEEK_API_KEY must be set"
This means the DeepSeek API key is not configured. Either:
1. Add `DEEPSEEK_API_KEY` to your `.env` file, or
2. The application will automatically fallback to OpenAI if you have `OPENAI_API_KEY` set

### Error: "No API keys found"
Make sure you have at least one of the following set:
- `OPENAI_API_KEY` (recommended)
- `DEEPSEEK_API_KEY` (optional)

