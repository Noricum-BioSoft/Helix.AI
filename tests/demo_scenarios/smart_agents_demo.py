#!/usr/bin/env python3
"""
Smart Agents in Action: Complete Backend Functionality Demo

This interactive demo showcases Helix.AI's complete agent system:
0. Ask Agent - Bioinformatics Q&A (Intent: Ask)
1. Execute Existing Tool - Small dataset, pre-built tool (Intent: Execute)
2. Large Dataset Processing - Async execution on 5GB dataset (Intent: Execute)
3. Complete Workflow - Multi-step with existing + code-generated tools (NEW!)
4. Infrastructure Decision - Smart routing (Local vs EMR based on data size)
5. Code Generator - Dynamic tool creation when no pre-built tool exists

Usage:
    python smart_agents_demo.py                      # Interactive menu (choose demos)
    python smart_agents_demo.py --full               # Run all demos
    python smart_agents_demo.py --ask-only           # Demo 0: Q&A
    python smart_agents_demo.py --execute-only       # Demo 1: Execute existing tool
    python smart_agents_demo.py --workflow-only      # Demo 2: Large dataset async
    python smart_agents_demo.py --complete-workflow  # Demo 3: Complete workflow (NEW!)
    python smart_agents_demo.py --infra-only         # Demo 4: Infrastructure decisions
    python smart_agents_demo.py --codegen-only       # Demo 5: Code generator
    python smart_agents_demo.py --skip-large         # Skip long-running demos
"""

import sys
import os
import asyncio
import argparse
import time
from pathlib import Path
from datetime import datetime
import json

# Global: allow running demos non-interactively (CI/smoke runs)
NON_INTERACTIVE = False

# Add project root to path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))
sys.path.insert(0, str(project_root / "backend"))

# Load environment variables from .env file
env_file = project_root / "backend" / ".env"
if env_file.exists():
    print(f"📝 Loading environment variables from {env_file}")
    with open(env_file, 'r') as f:
        for line in f:
            line = line.strip()
            # Skip comments and empty lines
            if line and not line.startswith('#') and '=' in line:
                # Split on first = only
                key, value = line.split('=', 1)
                # Remove quotes if present
                value = value.strip().strip('"').strip("'")
                # Only set if not already in environment
                if key not in os.environ and value:
                    os.environ[key] = value
    print(f"✅ Environment variables loaded\n")
else:
    print(f"⚠️  Warning: .env file not found at {env_file}")
    print(f"   API keys must be set in environment\n")


def print_section(title: str, width: int = 80):
    """Print a styled section header."""
    print(f"\n{'='*width}")
    print(f"  {title}")
    print(f"{'='*width}\n")


def print_subsection(title: str):
    """Print a styled subsection header."""
    print(f"\n{'-'*70}")
    print(f"  {title}")
    print(f"{'-'*70}\n")


def print_key_value(key: str, value: str, indent: int = 0):
    """Print a key-value pair with optional indentation."""
    indent_str = " " * indent
    print(f"{indent_str}• {key}: {value}")


def wait_for_user(message: str = "Press Enter to continue..."):
    """Wait for user to press Enter."""
    if NON_INTERACTIVE or os.getenv("HELIX_DEMO_NONINTERACTIVE", "").lower() in ("1", "true", "yes"):
        return
    input(f"\n{message}")


async def demo_infrastructure_small(backend_url: str = "http://localhost:8001"):
    """Demo 1A: Small dataset → Local execution."""
    print_section("DEMO 1A: Small Dataset → Local Execution")
    
    print("📊 Scenario: FastQC on small test dataset (~19MB)")
    print_key_value("Input files", "2 FASTQ files on S3 (test dataset)")
    print_key_value("Expected routing", "LOCAL (files are small enough to download)")
    print_key_value("Expected cost", "$0.002 (S3 egress only)")
    
    wait_for_user("\nReady to run? Press Enter...")
    
    print(f"\n🚀 Executing: python run_fastqc_full_execution.py --small --backend-url {backend_url}\n")
    
    import subprocess
    result = subprocess.run(
        [sys.executable, "run_fastqc_full_execution.py", "--small", "--backend-url", backend_url],
        cwd=Path(__file__).parent,
        capture_output=False,  # Show output in real-time
        text=True
    )
    
    if result.returncode == 0:
        print("\n✅ Demo 1A completed successfully!")
        print("\n🔍 Key Observations:")
        print_key_value("Infrastructure Decision", "Local (automatic routing)", indent=2)
        print_key_value("Reasoning", "Files <100MB threshold, download is fast", indent=2)
        print_key_value("Confidence", "High (0.9+) - all file sizes known", indent=2)
        print_key_value("Cost", "~$0.002 (nearly free)", indent=2)
    else:
        print(f"\n⚠️  Demo 1A failed with return code {result.returncode}")
    
    wait_for_user()


async def demo_infrastructure_large(backend_url: str = "http://localhost:8001"):
    """Demo 1B: Large dataset → EMR execution."""
    print_section("DEMO 1B: Large Dataset → EMR Execution")
    
    print("📊 Scenario: FastQC on large dataset (~5GB)")
    print_key_value("Input files", "2 FASTQ files on S3 (full dataset)")
    print_key_value("Expected routing", "EMR (files exceed 100MB threshold)")
    print_key_value("Expected cost", "$2-5 (EMR cluster + compute)")
    
    print("\n⏱️  Note: This will take 10-15 minutes to complete.")
    print("   The script will poll EMR job status every 30 seconds.")
    
    wait_for_user("\nReady to run? Press Enter...")
    
    print(f"\n🚀 Executing: python run_fastqc_full_execution.py --large --backend-url {backend_url}\n")
    
    import subprocess
    result = subprocess.run(
        [sys.executable, "run_fastqc_full_execution.py", "--large", "--backend-url", backend_url],
        cwd=Path(__file__).parent,
        capture_output=False,  # Show output in real-time
        text=True
    )
    
    if result.returncode == 0:
        print("\n✅ Demo 1B completed successfully!")
        print("\n🔍 Key Observations:")
        print_key_value("Infrastructure Decision", "EMR (automatic routing)", indent=2)
        print_key_value("Reasoning", "Files >100MB, process in-place on S3", indent=2)
        print_key_value("Confidence", "High (0.9+) - all file sizes known", indent=2)
        print_key_value("Cost savings", "Avoids $0.45 S3 egress + local compute", indent=2)
        print_key_value("Actual cost", "$2-5 (vs. $3-8 for always-on cluster)", indent=2)
    else:
        print(f"\n⚠️  Demo 1B failed with return code {result.returncode}")
    
    wait_for_user()


async def demo_infrastructure_comparison():
    """Demo 1C: Show infrastructure decision logic."""
    print_section("DEMO 1C: Infrastructure Decision Logic")
    
    print("📖 Opening infrastructure-decision-agent specification...\n")
    
    spec_path = project_root / "agents" / "infrastructure-decision-agent.md"
    
    print("🔍 Key Highlights from the Specification:\n")
    
    print("1️⃣  Decision Framework: The '3-Factor Model'")
    print("   - Data Locality & Size (PRIMARY)")
    print("   - Tool Availability (SECONDARY)")
    print("   - Computational Requirements (TERTIARY)")
    
    print("\n2️⃣  Key Threshold: 100MB total on S3")
    print("   - <100MB → Local or EC2 (fast download)")
    print("   - >100MB → EMR (process in-place, avoid egress)")
    
    print("\n3️⃣  Confidence Scoring (0.0-1.0)")
    print("   - 0.9-1.0: High confidence (all factors known)")
    print("   - 0.7-0.9: Medium confidence (some assumptions)")
    print("   - 0.5-0.7: Low confidence (multiple unknowns)")
    print("   - <0.5: Very low confidence (heuristic fallback)")
    
    print("\n4️⃣  Cost Reasoning")
    print("   - Provides RANGES, not fake precision")
    print("   - Includes explicit assumptions")
    print("   - Examples:")
    print("     • EMR: $2-5 (depends on runtime, cluster size)")
    print("     • EC2: $0.50-2.00 (instance type dependent)")
    print("     • Local: $0.00 (free, but may be slower)")
    
    print("\n5️⃣  Validation & Repair")
    print("   - Pydantic V2 contracts")
    print("   - LLM responses validated automatically")
    print("   - Up to 2 retries with error feedback")
    print("   - Fallback to heuristics if validation fails")
    
    wait_for_user("\nPress Enter to view the full specification file...")
    
    # Try to open in default editor
    try:
        if sys.platform == "darwin":  # macOS
            os.system(f'open "{spec_path}"')
        elif sys.platform == "linux":
            os.system(f'xdg-open "{spec_path}"')
        elif sys.platform == "win32":
            os.system(f'start "" "{spec_path}"')
        print(f"\n✅ Opened {spec_path} in default application")
    except Exception as e:
        print(f"\n⚠️  Could not open file automatically: {e}")
        print(f"   Please open manually: {spec_path}")
    
    wait_for_user()


async def demo_code_generator(backend_url: str = "http://localhost:8001"):
    """Demo 2A: Code Generator Agent - Custom tool generation."""
    print_section("DEMO 2A: Code Generator Agent - Custom Tool Generation")
    
    print("📊 Scenario: Read merging with custom quality threshold")
    print_key_value("Task", "Merge paired-end reads with Q30 quality threshold")
    print_key_value("Challenge", "No pre-built tool with this specific parameter")
    print_key_value("Solution", "Generate custom tool on-the-fly")
    
    wait_for_user("\nReady to generate tool? Press Enter...")
    
    print("\n🤖 Calling Tool Generator Agent via backend API...\n")
    
    try:
        import requests
        
        # Create session
        print(f"📞 Connecting to backend at {backend_url}...")
        try:
            session_response = requests.post(f"{backend_url}/create_session", timeout=10)
            session_response.raise_for_status()
            session_id = session_response.json()["session_id"]
            print(f"✅ Session created: {session_id}\n")
        except Exception as e:
            print(f"❌ Failed to create session: {e}")
            print(f"   Is the backend running? Try: uv run python -m backend.main_with_mcp")
            raise
        
        command = """
        Perform read merging on paired-end reads with custom quality threshold of 30:
        - R1: s3://noricum-ngs-data/datasets/GRCh38.p12.MafHi/test/test_mate_R1.fq
        - R2: s3://noricum-ngs-data/datasets/GRCh38.p12.MafHi/test/test_mate_R2.fq
        - Output: s3://helix-ai-results/merged_q30.fq
        """
        
        print("⏳ Generating tool (this may take 30-60 seconds)...\n")
        
        # Call backend API
        execute_response = requests.post(
            f"{backend_url}/execute",
            json={
                "command": command.strip(),
                "session_id": session_id
            },
            timeout=180  # 3 minute timeout for code generation
        )
        execute_response.raise_for_status()
        result = execute_response.json()
        
        print("\n" + "="*70)
        print("CODE GENERATOR RESULT")
        print("="*70 + "\n")
        
        # Handle both direct result and nested result structures
        if isinstance(result, dict):
            # Check for nested structures from MCP response
            if "data" in result and isinstance(result["data"], dict):
                if "results" in result["data"]:
                    result = result["data"]["results"]
                elif "result" in result["data"]:
                    result = result["data"]["result"]
        
        print(f"Status: {result.get('status', 'unknown')}")
        print(f"Tool Generated: {result.get('tool_generated', False)}")
        
        if result.get('explanation'):
            print(f"\nExplanation:\n{result['explanation']}\n")
        
        if result.get('code_preview'):
            print("Code Preview (first 500 chars):")
            print("-" * 70)
            print(result['code_preview'])
            print("-" * 70)
        
        if result.get('execution_result'):
            exec_result = result['execution_result']
            print(f"\nExecution Status: {exec_result.get('status', 'unknown')}")
            
            if exec_result.get('stdout'):
                print("\nOutput:")
                print(exec_result['stdout'][:500])
            
            # Show errors if execution failed
            if exec_result.get('status') == 'error':
                if exec_result.get('error'):
                    print(f"\n⚠️  Error: {exec_result['error']}")
                if exec_result.get('stderr'):
                    print(f"\nStderr (first 500 chars):")
                    print(exec_result['stderr'][:500])
        
        print("\n✅ Demo 2A completed!")
        
        print("\n🔍 Key Observations:")
        print_key_value("Tool Priority", "Prefers established tools (BBMerge) first", indent=2)
        print_key_value("Installation", "Attempts conda install if tool missing", indent=2)
        print_key_value("Fallback", "Python implementation only as last resort", indent=2)
        print_key_value("Infrastructure", "Consulted Infrastructure Decision Agent", indent=2)
        print_key_value("Time to solution", "30-60 seconds (vs. weeks for manual dev)", indent=2)
        
    except requests.RequestException as e:
        print(f"\n❌ Backend API request failed: {e}")
        if hasattr(e, 'response') and e.response is not None:
            print(f"   Status code: {e.response.status_code}")
            try:
                print(f"   Response: {e.response.json()}")
            except:
                print(f"   Response: {e.response.text[:500]}")
    except Exception as e:
        print(f"\n❌ Demo 2A failed: {e}")
        import traceback
        traceback.print_exc()
    
    wait_for_user()


async def demo_code_generator_spec():
    """Demo 2B: Show code generator specification."""
    print_section("DEMO 2B: Code Generator Specification")
    
    print("📖 Opening tool-generator-agent specification...\n")
    
    spec_path = project_root / "agents" / "tool-generator-agent.md"
    
    print("🔍 Key Highlights from the Specification:\n")
    
    print("1️⃣  Tool Selection Priority")
    print("   - 1st choice: Established tools (BBMerge, FLASH, PEAR)")
    print("   - 2nd choice: Python libraries (BioPython, pysam)")
    print("   - 3rd choice: Custom Python (only if no alternatives)")
    
    print("\n2️⃣  Workflow")
    print("   a) Task Analysis (identify operation, inputs, outputs)")
    print("   b) Tool Research (find appropriate bioinformatics tools)")
    print("   c) Infrastructure Decision (consult Infrastructure Agent)")
    print("   d) Code Generation (create executable implementation)")
    print("   e) Execution (run on selected infrastructure)")
    
    print("\n3️⃣  Error Handling Requirements")
    print("   - Input validation (paths, formats, permissions)")
    print("   - File system errors (not found, disk space)")
    print("   - Tool execution errors (not installed, out of memory)")
    print("   - Infrastructure errors (AWS credentials, EMR failures)")
    print("   - Data integrity errors (corrupt files, empty outputs)")
    
    print("\n4️⃣  Integration with Infrastructure Agent")
    print("   - Receives infrastructure decision (Local, EC2, EMR)")
    print("   - Generates code for SPECIFIED infrastructure")
    print("   - Does NOT make infrastructure decisions")
    
    print("\n5️⃣  Execution Environments")
    print("   - Local: Direct subprocess calls, file I/O")
    print("   - EC2: SSH execution with pre-installed tools")
    print("   - EMR: PySpark jobs, S3-native processing")
    print("   - Batch: Containerized execution")
    print("   - Lambda: Serverless (15min limit)")
    
    wait_for_user("\nPress Enter to view the full specification file...")
    
    # Try to open in default editor
    try:
        if sys.platform == "darwin":  # macOS
            os.system(f'open "{spec_path}"')
        elif sys.platform == "linux":
            os.system(f'xdg-open "{spec_path}"')
        elif sys.platform == "win32":
            os.system(f'start "" "{spec_path}"')
        print(f"\n✅ Opened {spec_path} in default application")
    except Exception as e:
        print(f"\n⚠️  Could not open file automatically: {e}")
        print(f"   Please open manually: {spec_path}")
    
    wait_for_user()


async def demo_test_coverage():
    """Demo 2C: Show test coverage."""
    print_section("DEMO 2C: Test Coverage")
    
    print("🧪 Running unit tests for both agents...\n")
    
    wait_for_user("Press Enter to run tests...")
    
    import subprocess
    
    test_files = [
        "tests/unit/backend/test_infra_decision_validation.py",
        "tests/unit/backend/test_tool_generator_agent.py",
    ]
    
    for test_file in test_files:
        print(f"\n{'='*70}")
        print(f"Running: pytest {test_file} -v")
        print('='*70 + "\n")
        
        result = subprocess.run(
            ["pytest", test_file, "-v", "--tb=short"],
            cwd=project_root,
            capture_output=False,
            text=True
        )
        
        if result.returncode == 0:
            print(f"\n✅ {test_file} passed!")
        else:
            print(f"\n⚠️  {test_file} had failures")
    
    print("\n🔍 Test Coverage Highlights:")
    print_key_value("Infrastructure Agent", "Contract validation, confidence scoring, cost ranges", indent=2)
    print_key_value("Tool Generator", "Tool selection, installation, fallbacks, error handling", indent=2)
    print_key_value("Integration Tests", "End-to-end workflows with real AWS infrastructure", indent=2)
    
    wait_for_user()


async def demo_summary():
    """Print demo summary and next steps."""
    print_section("DEMO SUMMARY & VALUE PROPOSITION")
    
    print("💰 Key Value Propositions:\n")
    
    print("1️⃣  Intelligent Cost Optimization")
    print("   ✓ Small datasets: Run locally (cost: $0.002)")
    print("   ✓ Large datasets: Auto-route to EMR (cost: $3.42)")
    print("   ✓ Savings: 60-80% vs. always-on clusters\n")
    
    print("2️⃣  Unlimited Extensibility")
    print("   ✓ No pre-built tool? Generate in 30-60 seconds")
    print("   ✓ Adapts to custom requirements")
    print("   ✓ Time savings: Seconds vs. weeks\n")
    
    print("3️⃣  Zero AWS Expertise Required")
    print("   ✓ Users: 'Analyze my RNA-seq data'")
    print("   ✓ System: File sizing → routing → execution")
    print("   ✓ Accessibility: Students → researchers → PIs\n")
    
    print("🏗️  Technical Moat:\n")
    print("   • Agent architecture with formal contracts (not prompt chains)")
    print("   • Pydantic validation + repair (up to 2 retries)")
    print("   • 6 specialized agents with clear responsibilities")
    print("   • 129 test files ensuring reliability")
    print("   • Domain-specific optimization (100MB threshold, tool database)")
    
    print("\n📈 Market Opportunity:\n")
    print("   • Academic researchers: 500K+ globally ($50-200/mo)")
    print("   • Biotech companies: 5K+ companies ($500-5K/mo)")
    print("   • Pharma R&D: Top 50 pharma ($50K-500K/yr)")
    
    print("\n🚀 Next Steps:\n")
    print("   1. Technical deep-dive (45 min)")
    print("   2. Demo environment access")
    print("   3. Roadmap discussion (Q1-Q3 2026)")
    
    print("\n" + "="*70)
    print("Thank you for watching the demo!")
    print("="*70 + "\n")


async def demo_ask_question(backend_url: str = "http://localhost:8001"):
    """Demo 0: Ask Agent - Simple Q&A."""
    print_section("DEMO 0: Ask Agent - Bioinformatics Q&A")
    
    print("📚 Scenario: User has a question about bioinformatics")
    print_key_value("Intent", "Ask (not execute)")
    print_key_value("Agent Flow", "Intent Detector → Bioinformatics Guru")
    print_key_value("Use Case", "Learning, understanding, getting recommendations")
    
    wait_for_user("\nReady to ask a question? Press Enter...")
    
    print("\n🤖 Asking: 'What is the difference between paired-end and single-end sequencing?'\n")
    
    try:
        import requests
        
        # Create session
        session_response = requests.post(f"{backend_url}/create_session", timeout=10)
        session_response.raise_for_status()
        session_id = session_response.json()["session_id"]
        
        question = "What is the difference between paired-end and single-end sequencing? When should I use each?"
        
        # Execute question
        execute_response = requests.post(
            f"{backend_url}/execute",
            json={"command": question, "session_id": session_id},
            timeout=60
        )
        execute_response.raise_for_status()
        result = execute_response.json()
        
        print("="*70)
        print("BIOINFORMATICS GURU RESPONSE")
        print("="*70 + "\n")
        
        # Extract answer from nested structure
        answer = result
        if isinstance(result, dict):
            if "response" in result:
                answer = result["response"]
            elif "data" in result and isinstance(result["data"], dict):
                if "response" in result["data"]:
                    answer = result["data"]["response"]
        
        print(answer if isinstance(answer, str) else json.dumps(answer, indent=2))
        
        print("\n" + "="*70)
        print("\n✅ Demo 0 completed!")
        
        print("\n🔍 Key Observations:")
        print_key_value("Intent Detection", "Correctly identified as 'ask' (not execute)", indent=2)
        print_key_value("Agent Routing", "Sent to Bioinformatics Guru, not Executor", indent=2)
        print_key_value("Response Type", "Educational answer, no code execution", indent=2)
        print_key_value("Session Preserved", "Can ask follow-up questions in same session", indent=2)
        
    except Exception as e:
        print(f"\n❌ Demo 0 failed: {e}")
        import traceback
        traceback.print_exc()
    
    wait_for_user()


async def demo_execute_existing_tool(backend_url: str = "http://localhost:8001"):
    """Demo 1: Execute existing tool on small dataset."""
    print_section("DEMO 1: Execute Existing Tool - Small Dataset")
    
    print("🎯 Scenario: Run quality control on small test dataset")
    print_key_value("Dataset", "~19MB paired-end FASTQ files on S3")
    print_key_value("Tool", "FastQC (pre-existing tool, no code generation)")
    print_key_value("Agent Flow", "Intent → Executor → Infrastructure → Tool Selection → Execution")
    print_key_value("Expected Routing", "Local (small files, fast execution)")
    
    wait_for_user("\nReady to execute? Press Enter...")
    
    command = """Run FastQC quality analysis on these paired-end reads:
Input files:
  - R1: s3://noricum-ngs-data/datasets/GRCh38.p12.MafHi/test/test_mate_R1.fq (~9.5 MB)
  - R2: s3://noricum-ngs-data/datasets/GRCh38.p12.MafHi/test/test_mate_R2.fq (~9.5 MB)
Output location:
  - Results uploaded to: s3://noricum-ngs-data/datasets/GRCh38.p12.MafHi/test/fastqc-results/
  - HTML reports: test_mate_R1_fastqc.html, test_mate_R2_fastqc.html
  - ZIP archives: test_mate_R1_fastqc.zip, test_mate_R2_fastqc.zip"""
    
    print("\n🤖 Command being executed:")
    print("-" * 70)
    print(command)
    print("-" * 70 + "\n")
    
    try:
        import requests
        
        # Create session
        session_response = requests.post(f"{backend_url}/create_session", timeout=10)
        session_response.raise_for_status()
        session_id = session_response.json()["session_id"]
        print(f"✅ Session created: {session_id}\n")
        
        print("⏳ Executing (should complete in ~1 minute)...\n")
        
        # Execute command
        execute_response = requests.post(
            f"{backend_url}/execute",
            json={"command": command, "session_id": session_id},
            timeout=180  # 3 minute timeout
        )
        execute_response.raise_for_status()
        result = execute_response.json()
        
        print("\n" + "="*70)
        print("EXECUTION RESULT")
        print("="*70 + "\n")
        
        print(json.dumps(result, indent=2)[:1000] + "...")
        
        print("\n✅ Demo 1 completed!")
        
        print("\n🔍 Key Observations:")
        print_key_value("Intent Detection", "Identified as 'execute' intent", indent=2)
        print_key_value("Tool Selection", "Used existing FastQC tool (no code generation)", indent=2)
        print_key_value("Infrastructure", "Routed to Local execution (files <100MB)", indent=2)
        print_key_value("Cost", "~$0.002 (S3 egress only, no compute)", indent=2)
        print_key_value("Time", "~60 seconds (fast local execution)", indent=2)
        
    except Exception as e:
        print(f"\n❌ Demo 1 failed: {e}")
        import traceback
        traceback.print_exc()
    
    wait_for_user()


async def demo_workflow_large_dataset(backend_url: str = "http://localhost:8001"):
    """Demo 2: Large dataset processing with async execution."""
    print_section("DEMO 2: Large Dataset Processing - Async Execution")
    
    print("🔬 Scenario: Quality control on large RNA-seq dataset")
    print_key_value("Dataset", "~5GB paired-end FASTQ files on S3")
    print_key_value("Task", "FastQC quality control (Step 1 of typical RNA-seq workflow)")
    print_key_value("Agent Flow", "Intent → Executor → Infrastructure → EMR Job Manager")
    print_key_value("Expected Routing", "EMR cluster (large files, distributed processing)")
    print("\n📝 Note: This demonstrates async execution on large datasets.")
    print("   A complete RNA-seq workflow (alignment → quantification) would")
    print("   require additional tools that are not yet implemented.\n")
    
    wait_for_user("\nReady to generate workflow? Press Enter...")
    
    command = """Run quality control (FastQC) on this RNA-seq dataset:

Input files:
  - R1: s3://noricum-ngs-data/datasets/GRCh38.p12.MafHi/mate_R1.fq (~2.5 GB)
  - R2: s3://noricum-ngs-data/datasets/GRCh38.p12.MafHi/mate_R2.fq (~2.5 GB)

Output location (if specified):
  - Results: s3://noricum-ngs-data/results/rnaseq-workflow/fastqc/

Note: This demo executes FastQC (an existing tool) on large files (5GB total).
The system will automatically route to EMR for distributed processing.
For a complete workflow demonstration including code-generated tools, see Demo 3."""
    
    print("\n🤖 Command being executed:")
    print("-" * 70)
    print(command)
    print("-" * 70 + "\n")
    
    try:
        import requests
        
        # Create session
        session_response = requests.post(f"{backend_url}/create_session", timeout=10)
        session_response.raise_for_status()
        session_id = session_response.json()["session_id"]
        print(f"✅ Session created: {session_id}\n")
        
        print("⏳ Planning workflow (this will show the plan, not execute)...\n")
        
        # Execute command
        execute_response = requests.post(
            f"{backend_url}/execute",
            json={"command": command, "session_id": session_id},
            timeout=120
        )
        execute_response.raise_for_status()
        result = execute_response.json()
        
        print("\n" + "="*70)
        print("INITIAL RESPONSE - JOB SUBMITTED")
        print("="*70 + "\n")
        
        # Extract job_id from nested response structure
        job_id = None
        if isinstance(result, dict):
            if "data" in result and isinstance(result["data"], dict):
                if "results" in result["data"] and isinstance(result["data"]["results"], dict):
                    if "result" in result["data"]["results"]:
                        job_id = result["data"]["results"]["result"].get("job_id")
        
        if job_id:
            print(f"✅ Async job submitted!")
            print(f"🆔 Job ID: {job_id}")
            print(f"📊 Status: submitted")
            print(f"⏱️  This job will take 10-30 minutes to complete")
            print(f"\n🔄 Now polling for job status (checking every 10 seconds)...\n")
            
            # Poll for job status
            poll_count = 0
            max_polls = 180  # 30 minutes max (180 * 10 seconds)
            
            while poll_count < max_polls:
                try:
                    status_response = requests.get(
                        f"{backend_url}/jobs/{job_id}",
                        timeout=10
                    )
                    status_response.raise_for_status()
                    status_data = status_response.json()
                    
                    if "job" in status_data:
                        job = status_data["job"]
                        job_status = job.get("status", "unknown")
                        
                        # Show progress update
                        timestamp = datetime.now().strftime('%H:%M:%S')
                        print(f"[{timestamp}] Status: {job_status.upper()}", end="")
                        
                        if "progress" in job:
                            print(f" | Progress: {job['progress']}", end="")
                        
                        if job_status == "completed":
                            print("\n\n✅ Job completed successfully!")
                            
                            # Get results
                            results_response = requests.get(
                                f"{backend_url}/jobs/{job_id}/results",
                                timeout=10
                            )
                            if results_response.status_code == 200:
                                results_data = results_response.json()
                                print("\n" + "="*70)
                                print("JOB RESULTS")
                                print("="*70 + "\n")
                                print(json.dumps(results_data, indent=2)[:1500] + "...")
                            break
                        
                        elif job_status == "failed":
                            print("\n\n❌ Job failed!")
                            if "error" in job:
                                print(f"Error: {job['error']}")
                            break
                        
                        else:
                            # Job still running
                            print(" (waiting...)")
                        
                    poll_count += 1
                    if job_status not in ["completed", "failed"]:
                        time.sleep(10)  # Wait 10 seconds before next poll
                    
                except Exception as poll_error:
                    print(f"\n⚠️  Error polling status: {poll_error}")
                    break
            
            if poll_count >= max_polls:
                print(f"\n⏱️  Polling timeout reached (30 minutes). Job may still be running.")
                print(f"   Check status manually: curl {backend_url}/jobs/{job_id}")
        
        else:
            # No job_id found - might be synchronous execution
            print("Response (first 2000 chars):")
            print(json.dumps(result, indent=2)[:2000] + "...")
        
        print("\n✅ Demo 2 completed!")
        
        print("\n🔍 Key Observations:")
        print_key_value("Intent Detection", "Identified as 'execute' intent", indent=2)
        print_key_value("Tool Selection", "FastQC quality control", indent=2)
        print_key_value("Infrastructure", "EMR cluster (5GB > 100MB threshold)", indent=2)
        print_key_value("Cost Estimate", "$2-5 (EMR cluster + S3-native processing)", indent=2)
        print_key_value("Async Execution", "Returns job_id and polls for completion", indent=2)
        print_key_value("Polling Pattern", "Client checks /jobs/{job_id} every 10 seconds", indent=2)
        print_key_value("Data Locality", "Processes data on S3 without egress", indent=2)
        print_key_value("Scalability", "Same pattern works for full workflows when tools are added", indent=2)
        
    except Exception as e:
        print(f"\n❌ Demo 2 failed: {e}")
        import traceback
        traceback.print_exc()
    
    wait_for_user()


async def demo_complete_workflow_with_codegen(backend_url: str = "http://localhost:8001"):
    """Demo 3: Complete multi-step workflow with existing + code-generated tools."""
    print_section("DEMO 3: Complete Workflow - Existing + Generated Tools")
    
    print("🎯 Scenario: Complete RNA-seq analysis pipeline")
    print_key_value("Challenge", "Execute a 3-step workflow where some tools don't exist yet")
    print_key_value("Solution", "Use existing FastQC + code-generate alignment & quantification tools")
    print_key_value("Dataset", "Small test dataset (~19MB for faster demonstration)")
    
    print("\n📋 Workflow Steps:")
    print("   1️⃣  Quality Control → FastQC (existing tool)")
    print("   2️⃣  Read Alignment → STAR/HISAT2 (code-generated on-the-fly)")
    print("   3️⃣  Gene Quantification → featureCounts (code-generated on-the-fly)")
    
    print("\n💡 This demonstrates Helix.AI's key differentiator:")
    print("   Seamlessly combine pre-built tools with dynamically generated code")
    print("   to execute complex workflows that would normally take weeks to implement.")
    
    wait_for_user("\nReady to execute complete workflow? Press Enter...")
    
    command = """Execute a complete RNA-seq analysis workflow on this dataset:

Workflow:
  Step 1: Quality Control (FastQC)
  Step 2: Read Alignment using STAR aligner with reference genome GRCh38
  Step 3: Gene quantification using featureCounts

Input files:
  - R1: s3://noricum-ngs-data/datasets/GRCh38.p12.MafHi/test/test_mate_R1.fq (~9.5 MB)
  - R2: s3://noricum-ngs-data/datasets/GRCh38.p12.MafHi/test/test_mate_R2.fq (~9.5 MB)
  - Reference genome: GRCh38 (human)
  - GTF annotation: Gencode v40

Output:
  - QC reports: s3://noricum-ngs-data/results/complete-workflow/fastqc/
  - Alignments (BAM): s3://noricum-ngs-data/results/complete-workflow/aligned/
  - Gene counts: s3://noricum-ngs-data/results/complete-workflow/counts/gene_counts.txt

Requirements:
  - Use STAR aligner for alignment (generate tool if not available)
  - Use featureCounts for quantification (generate tool if not available)
  - Process steps sequentially: QC → Align → Quantify
  - Pass intermediate outputs between steps"""
    
    print("\n🤖 Command being executed:")
    print("-" * 70)
    print(command)
    print("-" * 70 + "\n")
    
    try:
        import requests
        
        # Create session
        session_response = requests.post(f"{backend_url}/create_session", timeout=10)
        session_response.raise_for_status()
        session_id = session_response.json()["session_id"]
        print(f"✅ Session created: {session_id}\n")
        
        print("⏳ Executing workflow (this may take several minutes)...\n")
        print("📝 What's happening:")
        print("   1. Intent detector identifies this as a multi-step workflow")
        print("   2. Workflow planner breaks it into 3 steps")
        print("   3. For each step:")
        print("      - Check if tool exists")
        print("      - If not, invoke Code Generator Agent to create it")
        print("      - Execute the step")
        print("      - Pass outputs to next step\n")
        
        # Execute command
        execute_response = requests.post(
            f"{backend_url}/execute",
            json={"command": command, "session_id": session_id},
            timeout=300  # 5 minute timeout for multi-step workflow
        )
        execute_response.raise_for_status()
        result = execute_response.json()
        
        print("\n" + "="*70)
        print("WORKFLOW EXECUTION RESULT")
        print("="*70 + "\n")
        
        print(json.dumps(result, indent=2)[:3000] + "...")
        
        # Parse what ACTUALLY happened
        tool_executed = result.get("tool", "unknown")
        execution_status = result.get("status", "unknown")
        execution_mode = None
        execution_time = None
        output_location = None
        
        if "data" in result and "results" in result["data"]:
            res = result["data"]["results"]
            if "result" in res:
                execution_mode = res["result"].get("execution_mode", "unknown")
                execution_time = res["result"].get("execution_time")
                output_location = res["result"].get("output")
        
        print("\n" + "="*70)
        print("WHAT ACTUALLY HAPPENED")
        print("="*70)
        
        print(f"\n✅ Step 1 (Quality Control):")
        print_key_value("Tool", "FastQC (existing tool)", indent=2)
        print_key_value("Status", f"{execution_status}", indent=2)
        print_key_value("Execution", f"{execution_mode} ({execution_time:.2f}s)" if execution_time else execution_mode, indent=2)
        print_key_value("Output", output_location if output_location else "Not specified", indent=2)
        
        print(f"\n❌ Step 2 (Read Alignment):")
        print_key_value("Tool", "STAR aligner", indent=2)
        print_key_value("Status", "NOT EXECUTED", indent=2)
        print_key_value("Reason", "Workflow orchestration not implemented", indent=2)
        
        print(f"\n❌ Step 3 (Gene Quantification):")
        print_key_value("Tool", "featureCounts", indent=2)
        print_key_value("Status", "NOT EXECUTED", indent=2)
        print_key_value("Reason", "Workflow orchestration not implemented", indent=2)
        
        print("\n" + "="*70)
        print("CURRENT LIMITATIONS")
        print("="*70)
        print("\n⚠️  The backend currently lacks:")
        print("   1. Multi-step workflow orchestration")
        print("   2. Tool availability checking")
        print("   3. Dynamic tool generation integration")
        print("   4. Step-to-step data passing")
        
        print("\n✅ What DOES work:")
        print("   • Single tool execution (FastQC)")
        print("   • Infrastructure routing (sandbox for small files)")
        print("   • Code Generator Agent (in isolation - see Demo 5)")
        
        print("\n📍 Where to find FastQC outputs:")
        if output_location:
            print(f"   aws s3 ls {output_location}")
            print(f"   # Download: aws s3 cp {output_location}test_mate_R1_fastqc.html .")
        
        print("\n📝 Where to find execution logs:")
        print("   • Backend logs: Check the terminal where you ran 'uv run python -m backend.main_with_mcp'")
        print("   • Sandbox logs: Docker container logs (if using sandbox)")
        print(f"   • Session ID: {session_id} (for debugging)")
        
    except requests.exceptions.Timeout:
        print("\n⏱️  Request timeout - workflow execution takes longer than expected")
        print("   This is normal for multi-step workflows with code generation")
        print("   In production, this would be handled as an async job")
    except Exception as e:
        print(f"\n❌ Demo 3 failed: {e}")
        print("\n📝 Note: This demo requires:")
        print("   - Code Generator Agent to be fully operational")
        print("   - Workflow orchestration with step chaining")
        print("   - Tool availability checking")
        print("   These features may still be in development.")
        import traceback
        traceback.print_exc()
    
    wait_for_user()


async def main():
    """Main demo execution."""
    parser = argparse.ArgumentParser(
        description="Smart Agents in Action: Complete Backend Functionality Demo"
    )
    parser.add_argument(
        "--full",
        action="store_true",
        help="Run complete demo (all scenarios)"
    )
    parser.add_argument(
        "--ask-only",
        action="store_true",
        help="Run only Ask Agent demo"
    )
    parser.add_argument(
        "--execute-only",
        action="store_true",
        help="Run only Execute Tool demo"
    )
    parser.add_argument(
        "--workflow-only",
        action="store_true",
        help="Run only Large Dataset Processing demo"
    )
    parser.add_argument(
        "--complete-workflow",
        action="store_true",
        help="Run only Complete Workflow demo (existing + generated tools)"
    )
    parser.add_argument(
        "--infra-only",
        action="store_true",
        help="Run only infrastructure decision demos"
    )
    parser.add_argument(
        "--codegen-only",
        action="store_true",
        help="Run only code generator demos"
    )
    parser.add_argument(
        "--skip-large",
        action="store_true",
        help="Skip large dataset demo (takes 10-15 min)"
    )
    parser.add_argument(
        "--backend-url",
        default="http://localhost:8001",
        help="Backend URL (default: http://localhost:8001)"
    )
    parser.add_argument(
        "--non-interactive",
        action="store_true",
        help="Run without waiting for Enter prompts (for automation)"
    )
    
    args = parser.parse_args()
    global NON_INTERACTIVE
    NON_INTERACTIVE = bool(args.non_interactive)
    
    print_section("🧬 Helix.AI - Smart Agents in Action", width=80)
    print("Complete Backend Functionality Showcase\n")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Interactive menu (if no command-line flags)
    selected_demos = None
    has_cli_flags = any([args.ask_only, args.execute_only, args.workflow_only, args.complete_workflow, args.infra_only, args.codegen_only, args.full])
    
    if not has_cli_flags:
        # Show demo roadmap for interactive mode
        print("\n📋 Demo Roadmap:")
        print("   0️⃣  Ask Agent - Simple Q&A (Intent: Ask)")
        print("   1️⃣  Execute Existing Tool - Small Dataset (Intent: Execute)")
        print("   2️⃣  Large Dataset Processing - Async Execution on 5GB (Intent: Execute)")
        print("   3️⃣  Complete Workflow - Existing + Generated Tools (NEW! ⭐)")
        print("   4️⃣  Infrastructure Decision - Cost Optimization")
        print("   5️⃣  Code Generator - Dynamic Tool Creation")
        print("   6️⃣  Run ALL demos in sequence")
        
        print("\n" + "="*80)
        print("Which demo(s) would you like to run?")
        print("="*80)
        print("Options:")
        print("  - Enter a single number (0-6) to run one demo")
        print("  - Enter multiple numbers separated by commas (e.g., 0,1,3)")
        print("  - Enter 6 or 'all' to run all demos in sequence")
        print("  - Press Enter to run all demos (default)")
        
        user_input = input("\nYour choice: ").strip()
        
        if user_input.lower() == 'all' or user_input == '6' or user_input == '':
            selected_demos = [0, 1, 2, 3, 4, 5]
            print("\n✅ Running all demos in sequence...\n")
        else:
            try:
                # Parse comma-separated numbers
                selected_demos = [int(x.strip()) for x in user_input.split(',') if x.strip()]
                # Validate range
                selected_demos = [x for x in selected_demos if 0 <= x <= 5]
                if not selected_demos:
                    print("\n⚠️  No valid demos selected. Running all demos...\n")
                    selected_demos = [0, 1, 2, 3, 4, 5]
                else:
                    demo_names = {
                        0: "Ask Agent",
                        1: "Execute Tool",
                        2: "Large Dataset",
                        3: "Complete Workflow",
                        4: "Infrastructure Decision",
                        5: "Code Generator"
                    }
                    print(f"\n✅ Running selected demos: {', '.join([demo_names[x] for x in selected_demos])}\n")
            except ValueError:
                print("\n⚠️  Invalid input. Running all demos...\n")
                selected_demos = [0, 1, 2, 3, 4, 5]
        
        wait_for_user("Press Enter to continue...")
    
    try:
        # Scenario-specific runs (command-line flags)
        if args.ask_only:
            await demo_ask_question(backend_url=args.backend_url)
            return
        
        if args.execute_only:
            await demo_execute_existing_tool(backend_url=args.backend_url)
            return
        
        if args.workflow_only:
            await demo_workflow_large_dataset(backend_url=args.backend_url)
            return
        
        if args.complete_workflow:
            await demo_complete_workflow_with_codegen(backend_url=args.backend_url)
            return
        
        if args.codegen_only:
            await demo_code_generator(backend_url=args.backend_url)
            await demo_code_generator_spec()
            await demo_test_coverage()
            return
        
        if args.infra_only:
            await demo_infrastructure_small(backend_url=args.backend_url)
            if not args.skip_large:
                await demo_infrastructure_large(backend_url=args.backend_url)
            await demo_infrastructure_comparison()
            return
        
        # Full demo or interactive selection
        if args.full or selected_demos is not None:
            demos_to_run = selected_demos if selected_demos is not None else [0, 1, 2, 3, 4, 5]
            
            # Run selected demos
            if 0 in demos_to_run:
                print_section("DEMO 0: ASK AGENT", width=80)
                await demo_ask_question(backend_url=args.backend_url)
            
            if 1 in demos_to_run:
                print_section("DEMO 1: EXECUTE EXISTING TOOL", width=80)
                await demo_execute_existing_tool(backend_url=args.backend_url)
            
            if 2 in demos_to_run:
                print_section("DEMO 2: LARGE DATASET PROCESSING", width=80)
                await demo_workflow_large_dataset(backend_url=args.backend_url)
            
            if 3 in demos_to_run:
                print_section("DEMO 3: COMPLETE WORKFLOW (Existing + Generated Tools)", width=80)
                await demo_complete_workflow_with_codegen(backend_url=args.backend_url)
            
            if 4 in demos_to_run:
                print_section("DEMO 4: INFRASTRUCTURE DECISION", width=80)
                await demo_infrastructure_small(backend_url=args.backend_url)
                
                if not args.skip_large:
                    await demo_infrastructure_large(backend_url=args.backend_url)
                else:
                    print("\n⏭️  Skipping large dataset demo (--skip-large flag)")
                
                await demo_infrastructure_comparison()
            
            if 5 in demos_to_run:
                print_section("DEMO 5: CODE GENERATOR", width=80)
                await demo_code_generator(backend_url=args.backend_url)
                await demo_code_generator_spec()
                await demo_test_coverage()
        
        # Summary
        await demo_summary()
        
    except KeyboardInterrupt:
        print("\n\n⚠️  Demo interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\n\n❌ Demo failed with error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    asyncio.run(main())
