#!/usr/bin/env python3
"""
Proof that read_merging refactoring eliminated redundancy.

This script validates that:
1. No circular routing (read_merging doesn't call broker.execute_tool)
2. Direct execution path exists (calls merge_reads_from_s3 or run_read_merging_raw)
3. Tool-generator-agent is NOT in the execution path for read_merging
"""

import ast
import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(project_root))


def analyze_dispatch_tool():
    """Analyze dispatch_tool to verify refactoring."""
    dispatch_tool_path = project_root / "backend" / "main.py"
    
    with open(dispatch_tool_path, 'r') as f:
        source = f.read()
    
    tree = ast.parse(source)
    
    class ReadMergingAnalyzer(ast.NodeVisitor):
        def __init__(self):
            self.in_read_merging_branch = False
            self.has_broker_execute_tool = False
            self.has_tool_generator = False
            self.has_direct_execution = False
            self.lines = []
            
        def visit_If(self, node):
            # Check if we're in the read_merging elif branch
            if isinstance(node.test, ast.Compare):
                for op in node.test.ops:
                    if isinstance(op, ast.Eq):
                        # Check for tool_name == "read_merging"
                        if (isinstance(node.test.left, ast.Name) and 
                            node.test.left.id == "tool_name"):
                            for comparator in node.test.comparators:
                                if (isinstance(comparator, ast.Constant) and 
                                    comparator.value == "read_merging"):
                                    self.in_read_merging_branch = True
                                    self.lines.append(f"Line {node.lineno}: Found read_merging branch")
                                    self.generic_visit(node)
                                    self.in_read_merging_branch = False
                                    return
                                    
            self.generic_visit(node)
        
        def visit_Call(self, node):
            if self.in_read_merging_branch:
                # Check for broker.execute_tool calls
                if isinstance(node.func, ast.Attribute):
                    if node.func.attr == "execute_tool":
                        self.has_broker_execute_tool = True
                        self.lines.append(f"Line {node.lineno}: Found broker.execute_tool call (CIRCULAR ROUTING)")
                    elif node.func.attr == "generate_and_execute_tool":
                        self.has_tool_generator = True
                        self.lines.append(f"Line {node.lineno}: Found tool-generator-agent call")
                    elif node.func.attr in ("merge_reads_from_s3", "run_read_merging_raw"):
                        self.has_direct_execution = True
                        self.lines.append(f"Line {node.lineno}: Found direct execution: {node.func.attr}")
            self.generic_visit(node)
        
        def visit_Assign(self, node):
            if self.in_read_merging_branch:
                # Check for result assignments from direct calls
                if isinstance(node.value, ast.Call):
                    if isinstance(node.value.func, ast.Attribute):
                        if node.value.func.attr in ("merge_reads_from_s3", "run_read_merging_raw"):
                            self.has_direct_execution = True
                            self.lines.append(f"Line {node.lineno}: Found direct execution assignment: {node.value.func.attr}")
            self.generic_visit(node)
    
    analyzer = ReadMergingAnalyzer()
    analyzer.visit(tree)
    
    return analyzer


def main():
    """Run validation checks."""
    print("=" * 70)
    print("READ_MERGING REFACTORING VALIDATION")
    print("=" * 70)
    print()
    
    analyzer = analyze_dispatch_tool()
    
    print("Analysis Results:")
    print("-" * 70)
    for line in analyzer.lines:
        print(f"  {line}")
    
    if not analyzer.lines:
        print("  ⚠️  No read_merging branch found - check code structure")
    print()
    
    print("Validation Checks:")
    print("-" * 70)
    
    # Check 1: No circular routing
    if analyzer.has_broker_execute_tool:
        print("  ❌ FAIL: Found circular routing (broker.execute_tool in read_merging branch)")
        print("     This should have been removed in the refactor!")
        return False
    else:
        print("  ✅ PASS: No circular routing detected")
    
    # Check 2: No tool-generator-agent
    if analyzer.has_tool_generator:
        print("  ❌ FAIL: Found tool-generator-agent calls")
        print("     Tool-generator-agent should NOT be used for read_merging!")
        return False
    else:
        print("  ✅ PASS: No tool-generator-agent calls detected")
    
    # Check 3: Direct execution exists
    if analyzer.has_direct_execution:
        print("  ✅ PASS: Direct execution path exists (merge_reads_from_s3/run_read_merging_raw)")
    else:
        print("  ❌ FAIL: No direct execution path found")
        print("     Should call merge_reads_from_s3 or run_read_merging_raw directly!")
        return False
    
    print()
    print("=" * 70)
    print("✅ ALL VALIDATION CHECKS PASSED!")
    print("=" * 70)
    print()
    print("Summary:")
    print("  • read_merging now directly calls merge_reads_from_s3/run_read_merging_raw")
    print("  • No redundant tool-generator-agent routing")
    print("  • No circular broker.execute_tool calls")
    print("  • Clean, efficient execution path")
    print()
    
    return True


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)



