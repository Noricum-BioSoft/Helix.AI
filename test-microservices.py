#!/usr/bin/env python3
"""
Test script for DataBloom.AI Microservices Architecture
"""

import asyncio
import httpx
import json
import time
from typing import Dict, Any

class MicroservicesTester:
    def __init__(self):
        self.base_urls = {
            "api_gateway": "http://localhost:8000",
            "workflow_engine": "http://localhost:8001",
            "alignment_service": "http://localhost:8002"
        }
        self.session_id = None
        self.workflow_id = None
    
    async def test_health_checks(self) -> bool:
        """Test health endpoints of all services"""
        print("🏥 Testing health checks...")
        
        async with httpx.AsyncClient() as client:
            for service_name, base_url in self.base_urls.items():
                try:
                    response = await client.get(f"{base_url}/health")
                    if response.status_code == 200:
                        print(f"✅ {service_name}: Healthy")
                    else:
                        print(f"❌ {service_name}: Unhealthy (Status: {response.status_code})")
                        return False
                except Exception as e:
                    print(f"❌ {service_name}: Connection failed - {e}")
                    return False
        
        print("✅ All health checks passed!")
        return True
    
    async def test_session_creation(self) -> bool:
        """Test session creation through API Gateway"""
        print("\n📝 Testing session creation...")
        
        async with httpx.AsyncClient() as client:
            try:
                response = await client.post(f"{self.base_urls['api_gateway']}/api/v1/sessions")
                if response.status_code == 200:
                    data = response.json()
                    self.session_id = data.get("session_id")
                    self.workflow_id = data.get("workflow_id")
                    print(f"✅ Session created: {self.session_id}")
                    print(f"✅ Workflow created: {self.workflow_id}")
                    return True
                else:
                    print(f"❌ Session creation failed: {response.status_code}")
                    return False
            except Exception as e:
                print(f"❌ Session creation error: {e}")
                return False
    
    async def test_workflow_execution(self) -> bool:
        """Test workflow execution"""
        if not self.workflow_id:
            print("❌ No workflow ID available")
            return False
        
        print("\n🔄 Testing workflow execution...")
        
        async with httpx.AsyncClient() as client:
            try:
                # Start workflow
                response = await client.post(
                    f"{self.base_urls['api_gateway']}/api/v1/workflows/{self.workflow_id}/start"
                )
                if response.status_code == 200:
                    print("✅ Workflow started successfully")
                    
                    # Check workflow status
                    await asyncio.sleep(2)  # Wait for workflow to progress
                    status_response = await client.get(
                        f"{self.base_urls['api_gateway']}/api/v1/workflows/{self.workflow_id}"
                    )
                    if status_response.status_code == 200:
                        workflow_data = status_response.json()
                        print(f"✅ Workflow status: {workflow_data.get('workflow', {}).get('status')}")
                        return True
                    else:
                        print(f"❌ Failed to get workflow status: {status_response.status_code}")
                        return False
                else:
                    print(f"❌ Workflow start failed: {response.status_code}")
                    return False
            except Exception as e:
                print(f"❌ Workflow execution error: {e}")
                return False
    
    async def test_tool_execution(self) -> bool:
        """Test direct tool execution"""
        print("\n🔧 Testing tool execution...")
        
        if not self.session_id:
            print("❌ No session ID available")
            return False
        
        async with httpx.AsyncClient() as client:
            try:
                # Test alignment service
                test_sequences = ">seq1\nATGCGATCG\n>seq2\nATGCGATC"
                
                request_data = {
                    "session_id": self.session_id,
                    "workflow_id": self.workflow_id,
                    "step_id": "test_step",
                    "input_data": {
                        "sequences": test_sequences,
                        "algorithm": "clustalw"
                    }
                }
                
                response = await client.post(
                    f"{self.base_urls['api_gateway']}/api/v1/tools/sequence_alignment/execute",
                    json=request_data
                )
                
                if response.status_code == 200:
                    result = response.json()
                    if result.get("success"):
                        print("✅ Tool execution successful")
                        return True
                    else:
                        print(f"❌ Tool execution failed: {result.get('error')}")
                        return False
                else:
                    print(f"❌ Tool execution request failed: {response.status_code}")
                    return False
            except Exception as e:
                print(f"❌ Tool execution error: {e}")
                return False
    
    async def test_service_discovery(self) -> bool:
        """Test service discovery"""
        print("\n🔍 Testing service discovery...")
        
        async with httpx.AsyncClient() as client:
            try:
                response = await client.get(f"{self.base_urls['api_gateway']}/api/v1/services")
                if response.status_code == 200:
                    services = response.json().get("services", {})
                    print(f"✅ Found {len(services)} services:")
                    for service_name, service_url in services.items():
                        print(f"   - {service_name}: {service_url}")
                    return True
                else:
                    print(f"❌ Service discovery failed: {response.status_code}")
                    return False
            except Exception as e:
                print(f"❌ Service discovery error: {e}")
                return False
    
    async def run_all_tests(self):
        """Run all tests"""
        print("🚀 Starting DataBloom.AI Microservices Tests")
        print("=" * 50)
        
        tests = [
            ("Health Checks", self.test_health_checks),
            ("Session Creation", self.test_session_creation),
            ("Service Discovery", self.test_service_discovery),
            ("Tool Execution", self.test_tool_execution),
            ("Workflow Execution", self.test_workflow_execution),
        ]
        
        results = []
        for test_name, test_func in tests:
            try:
                result = await test_func()
                results.append((test_name, result))
            except Exception as e:
                print(f"❌ {test_name} failed with exception: {e}")
                results.append((test_name, False))
        
        # Print summary
        print("\n" + "=" * 50)
        print("📊 Test Results Summary:")
        print("=" * 50)
        
        passed = 0
        total = len(results)
        
        for test_name, result in results:
            status = "✅ PASS" if result else "❌ FAIL"
            print(f"{test_name}: {status}")
            if result:
                passed += 1
        
        print(f"\nOverall: {passed}/{total} tests passed")
        
        if passed == total:
            print("🎉 All tests passed! Microservices architecture is working correctly.")
        else:
            print("⚠️  Some tests failed. Check service logs for details.")
        
        return passed == total

async def main():
    """Main test function"""
    tester = MicroservicesTester()
    success = await tester.run_all_tests()
    
    if success:
        print("\n🎯 Architecture Verification Complete!")
        print("The microservices architecture is ready for use.")
    else:
        print("\n🔧 Please check the service logs and ensure all services are running.")
        print("Use 'docker-compose logs' to view service logs.")

if __name__ == "__main__":
    asyncio.run(main()) 