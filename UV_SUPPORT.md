# 🚀 uv Support in Helix.AI

Helix.AI now supports [uv](https://github.com/astral-sh/uv), a fast Python package installer and resolver written in Rust, as an alternative to pip.

## What is uv?

uv is an extremely fast Python package installer and resolver, written in Rust. It is designed as a drop-in replacement for pip and pip-tools.

**Benefits:**
- ⚡ **10-100x faster** than pip
- 🔒 **Deterministic** dependency resolution
- 📦 **Better dependency management**
- 🔄 **Drop-in replacement** for pip

## Installation

### Install uv

```bash
# macOS/Linux
curl -LsSf https://astral.sh/uv/install.sh | sh

# Windows (PowerShell)
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.sh | iex"

# Or using pip
pip install uv
```

## Usage

### Automatic Detection

The `start.sh` script automatically detects if `uv` is installed and uses it if available, falling back to pip if not found.

```bash
./start.sh
```

The script will:
1. Check if `uv` is available
2. Use `uv` if found (faster installs)
3. Fall back to `pip` if `uv` is not available

### Manual Installation

#### Using uv (Recommended)

```bash
cd backend
uv pip install -r requirements.txt
```

#### Using pip (Traditional)

```bash
cd backend
pip install -r requirements.txt
```

## Benefits for Helix.AI

1. **Faster Development**: Much faster dependency installation during development
2. **Better Resolution**: More reliable dependency resolution
3. **Backward Compatible**: Works as a drop-in replacement for pip
4. **Optional**: Script falls back to pip if uv is not installed

## Virtual Environments

uv works with your existing virtual environments:

```bash
# Create virtual environment (if desired)
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Use uv within the virtual environment
uv pip install -r backend/requirements.txt
```

Or let uv manage the environment:

```bash
# uv can also manage virtual environments
cd backend
uv venv  # Creates .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate
uv pip install -r requirements.txt
```

## Migration Notes

- ✅ **No changes required** to `requirements.txt`
- ✅ **Fully backward compatible** with pip
- ✅ **Automatic detection** in startup script
- ✅ **Optional** - works without uv installed

## Troubleshooting

### uv not detected

If the script doesn't detect uv:
1. Verify installation: `uv --version`
2. Ensure uv is in your PATH
3. Restart your terminal

### Installation issues

If you encounter issues with uv:
- The script will automatically fall back to pip
- You can manually use pip if needed
- Check uv documentation: https://github.com/astral-sh/uv

## Future Enhancements

Potential future improvements:
- Lock file generation with `uv pip compile`
- Project management with `uv init`
- Virtual environment management with `uv venv`

## References

- [uv Documentation](https://github.com/astral-sh/uv)
- [uv GitHub Repository](https://github.com/astral-sh/uv)
- [pip vs uv Comparison](https://github.com/astral-sh/uv#vs-pip)




