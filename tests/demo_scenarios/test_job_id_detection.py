#!/usr/bin/env python3
"""
Test that job_id detection works with actual backend response structure.
"""

import json

# Simulated response from your terminal output
response = {
    "version": "1.0",
    "success": True,
    "session_id": "7a141a26-26dc-41d3-9a30-31fd01a13234",
    "data": {
        "results": {
            "status": "submitted",
            "result": {
                "type": "job",
                "status": "submitted",
                "job_id": "7f0222f3-c47d-4269-8c67-52df897a8906",
                "message": "FastQC job submitted."
            }
        }
    },
    "raw_result": {
        "status": "submitted",
        "result": {
            "type": "job",
            "status": "submitted",
            "job_id": "7f0222f3-c47d-4269-8c67-52df897a8906"
        }
    }
}

def detect_job_id(result):
    """Same logic as in the script"""
    job_id = None
    if isinstance(result, dict):
        # Check multiple locations for job_id in order of likelihood
        # Direct top-level
        if "job_id" in result:
            job_id = result["job_id"]
            print(f"✅ Found job_id at: result['job_id']")
        # In result field
        elif "result" in result and isinstance(result["result"], dict):
            if "job_id" in result["result"]:
                job_id = result["result"]["job_id"]
                print(f"✅ Found job_id at: result['result']['job_id']")
        # In data field (direct)
        elif "data" in result and isinstance(result["data"], dict):
            if "job_id" in result["data"]:
                job_id = result["data"]["job_id"]
                print(f"✅ Found job_id at: result['data']['job_id']")
            # In data.results.result (nested structure from MCP response)
            elif "results" in result["data"] and isinstance(result["data"]["results"], dict):
                results = result["data"]["results"]
                if "job_id" in results:
                    job_id = results["job_id"]
                    print(f"✅ Found job_id at: result['data']['results']['job_id']")
                elif "result" in results and isinstance(results["result"], dict):
                    if "job_id" in results["result"]:
                        job_id = results["result"]["job_id"]
                        print(f"✅ Found job_id at: result['data']['results']['result']['job_id']")
        # In raw_result field
        elif "raw_result" in result and isinstance(result["raw_result"], dict):
            raw = result["raw_result"]
            if "job_id" in raw:
                job_id = raw["job_id"]
                print(f"✅ Found job_id at: result['raw_result']['job_id']")
            elif "result" in raw and isinstance(raw["result"], dict):
                if "job_id" in raw["result"]:
                    job_id = raw["result"]["job_id"]
                    print(f"✅ Found job_id at: result['raw_result']['result']['job_id']")
    
    return job_id

if __name__ == "__main__":
    print("Testing job_id detection with actual backend response...\n")
    
    job_id = detect_job_id(response)
    
    if job_id:
        print(f"\n✅ SUCCESS: Detected job_id = {job_id}")
    else:
        print(f"\n❌ FAILED: Could not detect job_id")
        print(f"Response structure: {list(response.keys())}")
        print(f"\nFull response:\n{json.dumps(response, indent=2)}")
