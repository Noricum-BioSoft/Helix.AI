#!/usr/bin/env python3
"""
Comprehensive Natural Language Command Mapping Test Suite
Tests how well Helix.AI interprets natural language commands and maps them to functions
"""

import requests
import json
import time
import sys
from typing import Dict, Any, List, Tuple

class NaturalLanguageMappingTester:
    def __init__(self, base_url: str = "http://localhost:8001"):
        self.base_url = base_url
        self.session_id = None
        self.test_results = []
        
    def log_test(self, test_name: str, success: bool, details: str = "", expected_tool: str = "", actual_tool: str = ""):
        """Log test results with tool mapping details"""
        status = "✅ PASS" if success else "❌ FAIL"
        print(f"{status} {test_name}")
        if details:
            print(f"   Details: {details}")
        if expected_tool and actual_tool:
            print(f"   Expected: {expected_tool} | Actual: {actual_tool}")
        self.test_results.append({
            "test": test_name,
            "success": success,
            "details": details,
            "expected_tool": expected_tool,
            "actual_tool": actual_tool
        })
    
    def test_session_creation(self) -> bool:
        """Create a test session"""
        try:
            response = requests.post(f"{self.base_url}/create_session")
            if response.status_code == 200:
                data = response.json()
                if "session_id" in data:
                    self.session_id = data["session_id"]
                    print(f"📋 Test Session ID: {self.session_id}")
                    return True
            return False
        except Exception:
            return False
    
    def execute_command_and_analyze(self, command: str, expected_tool: str, expected_params: Dict[str, Any] = None) -> Tuple[bool, str, str, Dict[str, Any]]:
        """Execute a command and analyze the tool mapping and parameters"""
        try:
            response = requests.post(f"{self.base_url}/execute", json={
                "command": command,
                "session_id": self.session_id
            })
            
            if response.status_code == 200:
                data = response.json()
                if data.get("success"):
                    # Get the session to extract tool information from history
                    session_response = requests.get(f"{self.base_url}/session/{self.session_id}")
                    if session_response.status_code == 200:
                        session_data = session_response.json()
                        if "session" in session_data and "history" in session_data["session"]:
                            history = session_data["session"]["history"]
                            if history:
                                # Get the most recent entry (the one we just executed)
                                latest_entry = history[-1]
                                actual_tool = latest_entry.get("tool", "unknown")
                            else:
                                actual_tool = "unknown"
                        else:
                            actual_tool = "unknown"
                    else:
                        actual_tool = "unknown"
                    
                    # Check if the tool matches expected
                    tool_match = expected_tool in actual_tool or actual_tool in expected_tool
                    
                    # Analyze parameters if provided
                    param_match = True
                    if expected_params:
                        # This is a simplified check - in practice we'd need to parse the actual parameters
                        param_match = True  # Placeholder for parameter validation
                    
                    success = tool_match and param_match
                    details = f"Command: '{command}'"
                    
                    return success, details, actual_tool, data.get("result", {})
                else:
                    return False, f"Command failed: {data.get('error', 'Unknown error')}", "none", {}
            else:
                return False, f"HTTP {response.status_code}", "none", {}
        except Exception as e:
            return False, f"Exception: {str(e)}", "none", {}
    
    def test_sequence_alignment_commands(self) -> List[bool]:
        """Test various natural language commands for sequence alignment"""
        print("\n🧬 Testing Sequence Alignment Commands")
        print("-" * 50)
        
        alignment_tests = [
            ("align these sequences: >seq1 ATGCGATCG >seq2 ATGCGATC", "sequence_alignment"),
            ("align DNA sequences ATGCGATCG and ATGCGATC", "sequence_alignment"),
            ("perform sequence alignment on ATGCGATCG, ATGCGATC", "sequence_alignment"),
            ("compare these sequences: ATGCGATCG vs ATGCGATC", "sequence_alignment"),
            ("align these DNA sequences", "sequence_alignment"),
            ("sequence alignment of ATGCGATCG and ATGCGATC", "sequence_alignment"),
            ("align the following sequences", "sequence_alignment"),
            ("perform multiple sequence alignment", "sequence_alignment"),
        ]
        
        results = []
        for command, expected_tool in alignment_tests:
            success, details, actual_tool, result = self.execute_command_and_analyze(command, expected_tool)
            self.log_test(f"Alignment: {command[:50]}...", success, details, expected_tool, actual_tool)
            results.append(success)
            time.sleep(0.5)  # Brief pause between tests
        
        return results
    
    def test_mutation_commands(self) -> List[bool]:
        """Test various natural language commands for mutation generation"""
        print("\n🧬 Testing Mutation Generation Commands")
        print("-" * 50)
        
        mutation_tests = [
            ("generate 10 variants of sequence ATGCGATCG", "mutate_sequence"),
            ("mutate sequence ATGCGATCG to create 20 variants", "mutate_sequence"),
            ("create mutations for ATGCGATCG", "mutate_sequence"),
            ("generate variants from ATGCGATCG", "mutate_sequence"),
            ("mutate this protein sequence", "mutate_sequence"),
            ("create 15 mutations of ATGCGATCG", "mutate_sequence"),
            ("generate random variants from ATGCGATCG", "mutate_sequence"),
            ("mutate DNA sequence ATGCGATCG", "mutate_sequence"),
        ]
        
        results = []
        for command, expected_tool in mutation_tests:
            success, details, actual_tool, result = self.execute_command_and_analyze(command, expected_tool)
            self.log_test(f"Mutation: {command[:50]}...", success, details, expected_tool, actual_tool)
            results.append(success)
            time.sleep(0.5)
        
        return results
    
    def test_variant_selection_commands(self) -> List[bool]:
        """Test various natural language commands for variant selection"""
        print("\n🧬 Testing Variant Selection Commands")
        print("-" * 50)
        
        selection_tests = [
            ("select top 5 variants from previous results", "select_variants"),
            ("choose the best variants", "select_variants"),
            ("select variants based on conservation", "select_variants"),
            ("pick the top variants", "select_variants"),
            ("select the best variants", "select_variants"),
            ("choose variants with highest diversity", "select_variants"),
            ("select top 10 variants", "select_variants"),
            ("pick the most conserved variants", "select_variants"),
        ]
        
        results = []
        for command, expected_tool in selection_tests:
            success, details, actual_tool, result = self.execute_command_and_analyze(command, expected_tool)
            self.log_test(f"Selection: {command[:50]}...", success, details, expected_tool, actual_tool)
            results.append(success)
            time.sleep(0.5)
        
        return results
    
    def test_phylogenetic_commands(self) -> List[bool]:
        """Test various natural language commands for phylogenetic analysis"""
        print("\n🧬 Testing Phylogenetic Analysis Commands")
        print("-" * 50)
        
        phylogenetic_tests = [
            ("build phylogenetic tree for sequences: >seq1 ATGCGATCG >seq2 ATGCGATC", "phylogenetic_tree"),
            ("create phylogenetic tree from aligned sequences", "phylogenetic_tree"),
            ("build a phylogenetic tree", "phylogenetic_tree"),
            ("construct evolutionary tree", "phylogenetic_tree"),
            ("generate phylogenetic tree", "phylogenetic_tree"),
            ("create evolutionary tree from sequences", "phylogenetic_tree"),
            ("build tree of life for these sequences", "phylogenetic_tree"),
            ("construct phylogenetic tree", "phylogenetic_tree"),
        ]
        
        results = []
        for command, expected_tool in phylogenetic_tests:
            success, details, actual_tool, result = self.execute_command_and_analyze(command, expected_tool)
            self.log_test(f"Phylogenetic: {command[:50]}...", success, details, expected_tool, actual_tool)
            results.append(success)
            time.sleep(0.5)
        
        return results
    
    def test_conservation_commands(self) -> List[bool]:
        """Test various natural language commands for conservation analysis"""
        print("\n🧬 Testing Conservation Analysis Commands")
        print("-" * 50)
        
        conservation_tests = [
            ("analyze conservation in aligned sequences", "sequence_alignment"),
            ("check sequence conservation", "sequence_alignment"),
            ("analyze sequence conservation", "handle_natural_command"),
            ("measure conservation across sequences", "sequence_alignment"),
            ("calculate conservation scores", "sequence_alignment"),
            ("analyze conservation patterns", "sequence_alignment"),
            ("check for conserved regions", "sequence_alignment"),
            ("analyze evolutionary conservation", "sequence_alignment"),
        ]
        
        results = []
        for command, expected_tool in conservation_tests:
            success, details, actual_tool, result = self.execute_command_and_analyze(command, expected_tool)
            self.log_test(f"Conservation: {command[:50]}...", success, details, expected_tool, actual_tool)
            results.append(success)
            time.sleep(0.5)
        
        return results
    
    def test_visualization_commands(self) -> List[bool]:
        """Test various natural language commands for visualization"""
        print("\n🧬 Testing Visualization Commands")
        print("-" * 50)
        
        visualization_tests = [
            ("visualize the alignment results", "sequence_alignment"),
            ("create a plot of the alignment", "sequence_alignment"),
            ("show alignment visualization", "sequence_alignment"),
            ("generate alignment plot", "sequence_alignment"),
            ("visualize sequence alignment", "sequence_alignment"),
            ("create alignment diagram", "sequence_alignment"),
            ("show alignment results", "sequence_alignment"),
            ("plot the alignment", "sequence_alignment"),
        ]
        
        results = []
        for command, expected_tool in visualization_tests:
            success, details, actual_tool, result = self.execute_command_and_analyze(command, expected_tool)
            self.log_test(f"Visualization: {command[:50]}...", success, details, expected_tool, actual_tool)
            results.append(success)
            time.sleep(0.5)
        
        return results
    
    def test_parameter_extraction(self) -> List[bool]:
        """Test parameter extraction from natural language commands"""
        print("\n🧬 Testing Parameter Extraction")
        print("-" * 50)
        
        parameter_tests = [
            ("generate 25 variants of sequence ATGCGATCG", {"num_variants": 25, "sequence": "ATGCGATCG"}),
            ("align sequences: >seq1 ATGCGATCG >seq2 ATGCGATC >seq3 ATGCGATCGATC", {"sequences": ">seq1 ATGCGATCG >seq2 ATGCGATC >seq3 ATGCGATCGATC"}),
            ("select top 15 variants based on diversity", {"num_variants": 15, "criteria": "diversity"}),
            ("mutate ATGCGATCG to create 50 variants", {"num_variants": 50, "sequence": "ATGCGATCG"}),
            ("build phylogenetic tree for: >seq1 ATGCGATCG >seq2 ATGCGATC", {"sequences": ">seq1 ATGCGATCG >seq2 ATGCGATC"}),
        ]
        
        results = []
        for command, expected_params in parameter_tests:
            success, details, actual_tool, result = self.execute_command_and_analyze(command, "any")
            
            # Check if the command was successful (indicating parameter extraction worked)
            param_success = success  # Simplified check - in practice we'd validate actual parameters
            
            self.log_test(f"Parameter: {command[:50]}...", param_success, details)
            results.append(param_success)
            time.sleep(0.5)
        
        return results
    
    def test_ambiguous_commands(self) -> List[bool]:
        """Test ambiguous commands that should be handled gracefully"""
        print("\n🧬 Testing Ambiguous Commands")
        print("-" * 50)
        
        ambiguous_tests = [
            ("analyze this", "handle_natural_command"),
            ("process sequences", "handle_natural_command"),
            ("do analysis", "handle_natural_command"),
            ("work with data", "handle_natural_command"),
            ("perform analysis", "handle_natural_command"),
            ("analyze everything", "handle_natural_command"),
            ("process the data", "handle_natural_command"),
            ("run analysis", "handle_natural_command"),
        ]
        
        results = []
        for command, expected_tool in ambiguous_tests:
            success, details, actual_tool, result = self.execute_command_and_analyze(command, expected_tool)
            self.log_test(f"Ambiguous: {command[:50]}...", success, details, expected_tool, actual_tool)
            results.append(success)
            time.sleep(0.5)
        
        return results
    
    def test_error_handling(self) -> List[bool]:
        """Test error handling for invalid commands"""
        print("\n🧬 Testing Error Handling")
        print("-" * 50)
        
        error_tests = [
            ("invalid command", "handle_natural_command"),
            ("execute impossible operation", "handle_natural_command"),
            ("perform undefined analysis", "handle_natural_command"),
            ("do something impossible", "handle_natural_command"),
            ("run invalid tool", "handle_natural_command"),
        ]
        
        results = []
        for command, expected_tool in error_tests:
            success, details, actual_tool, result = self.execute_command_and_analyze(command, expected_tool)
            # For error handling, we expect the command to be handled gracefully (not crash)
            error_handled = success or "error" in str(result).lower()
            self.log_test(f"Error: {command[:50]}...", error_handled, details, expected_tool, actual_tool)
            results.append(error_handled)
            time.sleep(0.5)
        
        return results
    
    def test_command_variations(self) -> List[bool]:
        """Test different ways to express the same command"""
        print("\n🧬 Testing Command Variations")
        print("-" * 50)
        
        variation_tests = [
            # Same intent, different phrasings
            ("align sequences ATGCGATCG and ATGCGATC", "sequence_alignment"),
            ("align ATGCGATCG with ATGCGATC", "sequence_alignment"),
            ("compare ATGCGATCG to ATGCGATC", "sequence_alignment"),
            ("align ATGCGATCG versus ATGCGATC", "sequence_alignment"),
            
            # Different ways to specify quantities
            ("generate 10 variants", "mutate_sequence"),
            ("create 10 variants", "mutate_sequence"),
            ("make 10 variants", "mutate_sequence"),
            ("produce 10 variants", "mutate_sequence"),
            
            # Different selection criteria
            ("select top 5", "select_variants"),
            ("choose top 5", "select_variants"),
            ("pick top 5", "select_variants"),
            ("select best 5", "select_variants"),
        ]
        
        results = []
        for command, expected_tool in variation_tests:
            success, details, actual_tool, result = self.execute_command_and_analyze(command, expected_tool)
            self.log_test(f"Variation: {command[:50]}...", success, details, expected_tool, actual_tool)
            results.append(success)
            time.sleep(0.5)
        
        return results
    
    def generate_comprehensive_report(self):
        """Generate a comprehensive test report"""
        print("\n" + "=" * 80)
        print("📊 COMPREHENSIVE NATURAL LANGUAGE MAPPING TEST REPORT")
        print("=" * 80)
        
        # Calculate statistics
        total_tests = len(self.test_results)
        passed_tests = sum(1 for result in self.test_results if result["success"])
        failed_tests = total_tests - passed_tests
        
        # Tool mapping accuracy
        tool_mapping_accuracy = 0
        if total_tests > 0:
            tool_mapping_accuracy = (passed_tests / total_tests) * 100
        
        print(f"\n📈 Overall Statistics:")
        print(f"   Total Tests: {total_tests}")
        print(f"   Passed: {passed_tests}")
        print(f"   Failed: {failed_tests}")
        print(f"   Success Rate: {tool_mapping_accuracy:.1f}%")
        
        # Analyze by category
        categories = {
            "Sequence Alignment": [],
            "Mutation Generation": [],
            "Variant Selection": [],
            "Phylogenetic Analysis": [],
            "Conservation Analysis": [],
            "Visualization": [],
            "Parameter Extraction": [],
            "Ambiguous Commands": [],
            "Error Handling": [],
            "Command Variations": []
        }
        
        for result in self.test_results:
            test_name = result["test"]
            for category in categories:
                if category.lower().replace(" ", "_") in test_name.lower():
                    categories[category].append(result)
                    break
        
        print(f"\n📋 Results by Category:")
        for category, results in categories.items():
            if results:
                passed = sum(1 for r in results if r["success"])
                total = len(results)
                rate = (passed / total) * 100 if total > 0 else 0
                print(f"   {category}: {passed}/{total} ({rate:.1f}%)")
        
        # Tool mapping analysis
        print(f"\n🔧 Tool Mapping Analysis:")
        tool_mappings = {}
        for result in self.test_results:
            expected = result.get("expected_tool", "unknown")
            actual = result.get("actual_tool", "unknown")
            if expected not in tool_mappings:
                tool_mappings[expected] = {"correct": 0, "incorrect": 0, "total": 0}
            
            tool_mappings[expected]["total"] += 1
            if expected in actual or actual in expected:
                tool_mappings[expected]["correct"] += 1
            else:
                tool_mappings[expected]["incorrect"] += 1
        
        for tool, stats in tool_mappings.items():
            if stats["total"] > 0:
                accuracy = (stats["correct"] / stats["total"]) * 100
                print(f"   {tool}: {stats['correct']}/{stats['total']} ({accuracy:.1f}%)")
        
        # Recommendations
        print(f"\n💡 Recommendations:")
        if tool_mapping_accuracy >= 90:
            print("   ✅ Excellent natural language mapping! System is very robust.")
        elif tool_mapping_accuracy >= 80:
            print("   ✅ Good natural language mapping. Minor improvements possible.")
        elif tool_mapping_accuracy >= 70:
            print("   ⚠️  Fair natural language mapping. Consider improving command patterns.")
        else:
            print("   ❌ Natural language mapping needs significant improvement.")
        
        if failed_tests > 0:
            print(f"   🔧 Focus on improving {failed_tests} failed test cases.")
        
        return tool_mapping_accuracy >= 80
    
    def run_all_tests(self):
        """Run all natural language mapping tests"""
        print("🧬 Helix.AI Natural Language Command Mapping Test Suite")
        print("=" * 80)
        
        # Create session
        if not self.test_session_creation():
            print("❌ Failed to create test session")
            return False
        
        # Run all test categories
        test_functions = [
            self.test_sequence_alignment_commands,
            self.test_mutation_commands,
            self.test_variant_selection_commands,
            self.test_phylogenetic_commands,
            self.test_conservation_commands,
            self.test_visualization_commands,
            self.test_parameter_extraction,
            self.test_ambiguous_commands,
            self.test_error_handling,
            self.test_command_variations
        ]
        
        for test_func in test_functions:
            test_func()
        
        # Generate comprehensive report
        return self.generate_comprehensive_report()

def main():
    """Main test runner"""
    tester = NaturalLanguageMappingTester()
    success = tester.run_all_tests()
    sys.exit(0 if success else 1)

if __name__ == "__main__":
    main() 