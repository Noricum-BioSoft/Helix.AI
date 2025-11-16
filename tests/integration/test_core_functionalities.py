#!/usr/bin/env python3
"""
Comprehensive Test Suite for Helix.AI Core Functionalities
Tests all common operations identified across workflows
"""

import requests
import json
import time
import sys
from typing import Dict, Any, List

class HelixAITester:
    def __init__(self, base_url: str = "http://localhost:8001"):
        self.base_url = base_url
        self.session_id = None
        self.test_results = []
        
    def log_test(self, test_name: str, success: bool, details: str = ""):
        """Log test results"""
        status = "✅ PASS" if success else "❌ FAIL"
        print(f"{status} {test_name}")
        if details:
            print(f"   Details: {details}")
        self.test_results.append({
            "test": test_name,
            "success": success,
            "details": details
        })
        
    def test_health_check(self) -> bool:
        """Test 1: Health Check"""
        try:
            response = requests.get(f"{self.base_url}/health")
            if response.status_code == 200:
                data = response.json()
                if data.get("status") == "healthy":
                    self.log_test("Health Check", True)
                    return True
                else:
                    self.log_test("Health Check", False, f"Unexpected response: {data}")
                    return False
            else:
                self.log_test("Health Check", False, f"Status code: {response.status_code}")
                return False
        except Exception as e:
            self.log_test("Health Check", False, f"Exception: {str(e)}")
            return False
    
    def test_session_creation(self) -> bool:
        """Test 2: Session Creation"""
        try:
            response = requests.post(f"{self.base_url}/create_session")
            if response.status_code == 200:
                data = response.json()
                if "session_id" in data:
                    self.session_id = data["session_id"]
                    self.log_test("Session Creation", True, f"Session ID: {self.session_id}")
                    return True
                else:
                    self.log_test("Session Creation", False, f"No session_id in response: {data}")
                    return False
            else:
                self.log_test("Session Creation", False, f"Status code: {response.status_code}")
                return False
        except Exception as e:
            self.log_test("Session Creation", False, f"Exception: {str(e)}")
            return False
    
    def test_sequence_alignment(self) -> bool:
        """Test 3: Sequence Alignment (Most Common Operation)"""
        command = "align these sequences: >seq1 ATGCGATCG >seq2 ATGCGATC >seq3 ATGCGATCGATC"
        try:
            response = requests.post(f"{self.base_url}/execute", json={
                "command": command,
                "session_id": self.session_id
            })
            if response.status_code == 200:
                data = response.json()
                if data.get("success") and "alignment" in data.get("result", {}):
                    alignment = data["result"]["alignment"]
                    if len(alignment) >= 2:
                        self.log_test("Sequence Alignment", True, f"Aligned {len(alignment)} sequences")
                        return True
                    else:
                        self.log_test("Sequence Alignment", False, f"Expected 2+ sequences, got {len(alignment)}")
                        return False
                else:
                    self.log_test("Sequence Alignment", False, f"Unexpected response: {data}")
                    return False
            else:
                self.log_test("Sequence Alignment", False, f"Status code: {response.status_code}")
                return False
        except Exception as e:
            self.log_test("Sequence Alignment", False, f"Exception: {str(e)}")
            return False
    
    def test_mutation_generation(self) -> bool:
        """Test 4: Mutation Generation"""
        command = "generate 10 variants of sequence ATGCGATCG"
        try:
            response = requests.post(f"{self.base_url}/execute", json={
                "command": command,
                "session_id": self.session_id
            })
            if response.status_code == 200:
                data = response.json()
                if data.get("success"):
                    result = data.get("result", {})
                    if "mutated_sequences" in result or "variants" in result or "variants" in result.get("statistics", {}):
                        self.log_test("Mutation Generation", True, "Variants generated successfully")
                        return True
                    else:
                        self.log_test("Mutation Generation", False, f"No variants in response: {result}")
                        return False
                else:
                    self.log_test("Mutation Generation", False, f"Unexpected response: {data}")
                    return False
            else:
                self.log_test("Mutation Generation", False, f"Status code: {response.status_code}")
                return False
        except Exception as e:
            self.log_test("Mutation Generation", False, f"Exception: {str(e)}")
            return False
    
    def test_phylogenetic_analysis(self) -> bool:
        """Test 5: Phylogenetic Analysis"""
        command = "build phylogenetic tree for sequences: >seq1 ATGCGATCG >seq2 ATGCGATC >seq3 ATGCGATCGATC"
        try:
            response = requests.post(f"{self.base_url}/execute", json={
                "command": command,
                "session_id": self.session_id
            })
            if response.status_code == 200:
                data = response.json()
                if data.get("success"):
                    self.log_test("Phylogenetic Analysis", True, "Tree analysis completed")
                    return True
                else:
                    self.log_test("Phylogenetic Analysis", False, f"Unexpected response: {data}")
                    return False
            else:
                self.log_test("Phylogenetic Analysis", False, f"Status code: {response.status_code}")
                return False
        except Exception as e:
            self.log_test("Phylogenetic Analysis", False, f"Exception: {str(e)}")
            return False
    
    def test_variant_selection(self) -> bool:
        """Test 6: Variant Selection"""
        command = "select top 5 variants from previous results"
        try:
            response = requests.post(f"{self.base_url}/execute", json={
                "command": command,
                "session_id": self.session_id
            })
            if response.status_code == 200:
                data = response.json()
                if data.get("success"):
                    self.log_test("Variant Selection", True, "Variants selected successfully")
                    return True
                else:
                    self.log_test("Variant Selection", False, f"Unexpected response: {data}")
                    return False
            else:
                self.log_test("Variant Selection", False, f"Status code: {response.status_code}")
                return False
        except Exception as e:
            self.log_test("Variant Selection", False, f"Exception: {str(e)}")
            return False
    
    def test_conservation_analysis(self) -> bool:
        """Test 7: Conservation Analysis"""
        command = "analyze conservation in aligned sequences"
        try:
            response = requests.post(f"{self.base_url}/execute", json={
                "command": command,
                "session_id": self.session_id
            })
            if response.status_code == 200:
                data = response.json()
                if data.get("success"):
                    self.log_test("Conservation Analysis", True, "Conservation analysis completed")
                    return True
                else:
                    self.log_test("Conservation Analysis", False, f"Unexpected response: {data}")
                    return False
            else:
                self.log_test("Conservation Analysis", False, f"Status code: {response.status_code}")
                return False
        except Exception as e:
            self.log_test("Conservation Analysis", False, f"Exception: {str(e)}")
            return False
    
    def test_visualization(self) -> bool:
        """Test 8: Data Visualization"""
        command = "visualize the alignment results"
        try:
            response = requests.post(f"{self.base_url}/execute", json={
                "command": command,
                "session_id": self.session_id
            })
            if response.status_code == 200:
                data = response.json()
                if data.get("success"):
                    self.log_test("Data Visualization", True, "Visualization generated")
                    return True
                else:
                    self.log_test("Data Visualization", False, f"Unexpected response: {data}")
                    return False
            else:
                self.log_test("Data Visualization", False, f"Status code: {response.status_code}")
                return False
        except Exception as e:
            self.log_test("Data Visualization", False, f"Exception: {str(e)}")
            return False
    
    def test_multi_step_workflow(self) -> bool:
        """Test 9: Multi-Step Workflow (Generate → Analyze → Select)"""
        steps = [
            "generate 20 variants of sequence ATGCGATCG",
            "align all the generated variants",
            "select top 5 variants based on conservation"
        ]
        
        try:
            for i, step in enumerate(steps, 1):
                print(f"   Step {i}: {step}")
                response = requests.post(f"{self.base_url}/execute", json={
                    "command": step,
                    "session_id": self.session_id
                })
                if response.status_code != 200:
                    self.log_test("Multi-Step Workflow", False, f"Step {i} failed: {response.status_code}")
                    return False
                time.sleep(1)  # Brief pause between steps
            
            self.log_test("Multi-Step Workflow", True, "All 3 steps completed successfully")
            return True
        except Exception as e:
            self.log_test("Multi-Step Workflow", False, f"Exception: {str(e)}")
            return False
    
    def test_natural_language_commands(self) -> bool:
        """Test 10: Natural Language Command Processing"""
        commands = [
            "align these DNA sequences",
            "mutate this protein sequence",
            "build a phylogenetic tree",
            "analyze sequence conservation",
            "select the best variants"
        ]
        
        success_count = 0
        for command in commands:
            try:
                response = requests.post(f"{self.base_url}/execute", json={
                    "command": command,
                    "session_id": self.session_id
                })
                if response.status_code == 200:
                    data = response.json()
                    if data.get("success"):
                        success_count += 1
                    else:
                        print(f"   Command failed: {command}")
                else:
                    print(f"   Command failed: {command} (Status: {response.status_code})")
            except Exception as e:
                print(f"   Command failed: {command} (Exception: {str(e)})")
        
        success_rate = success_count / len(commands)
        if success_rate >= 0.8:  # 80% success rate threshold
            self.log_test("Natural Language Commands", True, f"{success_count}/{len(commands)} commands successful")
            return True
        else:
            self.log_test("Natural Language Commands", False, f"Only {success_count}/{len(commands)} commands successful")
            return False
    
    def test_session_persistence(self) -> bool:
        """Test 11: Session Persistence"""
        try:
            # Get session info
            response = requests.get(f"{self.base_url}/session/{self.session_id}")
            if response.status_code == 200:
                data = response.json()
                if "session" in data and "history" in data["session"]:
                    self.log_test("Session Persistence", True, "Session data persisted correctly")
                    return True
                else:
                    self.log_test("Session Persistence", False, f"Missing session data: {data}")
                    return False
            else:
                self.log_test("Session Persistence", False, f"Status code: {response.status_code}")
                return False
        except Exception as e:
            self.log_test("Session Persistence", False, f"Exception: {str(e)}")
            return False
    
    def test_error_handling(self) -> bool:
        """Test 12: Error Handling"""
        invalid_commands = [
            "invalid command",
            "execute impossible operation",
            "perform undefined analysis"
        ]
        
        error_count = 0
        for command in invalid_commands:
            try:
                response = requests.post(f"{self.base_url}/execute", json={
                    "command": command,
                    "session_id": self.session_id
                })
                if response.status_code == 200:
                    data = response.json()
                    # Check if it's a successful response with an error inside
                    if not data.get("success") or (data.get("success") and "error" in str(data.get("result", {}))):
                        error_count += 1
                        print(f"   ✅ Error handled gracefully: {command}")
                    else:
                        print(f"   ❌ Unexpected success: {command}")
                else:
                    error_count += 1  # HTTP error is also acceptable
                    print(f"   ✅ HTTP error handled: {command}")
            except Exception as e:
                error_count += 1  # Exception is acceptable for invalid commands
                print(f"   ✅ Exception handled: {command} - {str(e)}")
        
        if error_count >= len(invalid_commands) * 0.8:  # 80% should handle errors gracefully
            self.log_test("Error Handling", True, f"{error_count}/{len(invalid_commands)} errors handled gracefully")
            return True
        else:
            self.log_test("Error Handling", False, f"Only {error_count}/{len(invalid_commands)} errors handled gracefully")
            return False
    
    def run_all_tests(self):
        """Run all tests and generate summary"""
        print("🧬 Helix.AI Core Functionality Test Suite")
        print("=" * 50)
        
        tests = [
            self.test_health_check,
            self.test_session_creation,
            self.test_sequence_alignment,
            self.test_mutation_generation,
            self.test_phylogenetic_analysis,
            self.test_variant_selection,
            self.test_conservation_analysis,
            self.test_visualization,
            self.test_multi_step_workflow,
            self.test_natural_language_commands,
            self.test_session_persistence,
            self.test_error_handling
        ]
        
        passed = 0
        for test in tests:
            if test():
                passed += 1
        
        print("\n" + "=" * 50)
        print(f"📊 Test Summary: {passed}/{len(tests)} tests passed")
        
        if passed == len(tests):
            print("🎉 All core functionalities are working correctly!")
        elif passed >= len(tests) * 0.8:
            print("✅ Most core functionalities are working correctly")
        else:
            print("⚠️  Several core functionalities need attention")
        
        return passed == len(tests)

def main():
    """Main test runner"""
    tester = HelixAITester()
    success = tester.run_all_tests()
    sys.exit(0 if success else 1)

if __name__ == "__main__":
    main()


