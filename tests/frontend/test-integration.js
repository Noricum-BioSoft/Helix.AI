#!/usr/bin/env node
/**
 * Test script to verify frontend-backend MCP integration
 */

const API_BASE_URL = 'http://localhost:8001';

async function testEndpoint(endpoint, method = 'GET', data = null) {
  try {
    const response = await fetch(`${API_BASE_URL}${endpoint}`, {
      method,
      headers: {
        'Content-Type': 'application/json',
      },
      body: data ? JSON.stringify(data) : undefined,
    });

    if (!response.ok) {
      throw new Error(`HTTP ${response.status}: ${response.statusText}`);
    }

    const result = await response.json();
    return { success: true, data: result };
  } catch (error) {
    return { success: false, error: error.message };
  }
}

async function runTests() {
  console.log('🧪 Testing Helix.AI Frontend-Backend Integration');
  console.log('=' .repeat(60));

  // Test 1: Health Check
  console.log('\n1. Testing server health...');
  const healthResult = await testEndpoint('/health');
  if (healthResult.success) {
    console.log('✅ Health check passed:', healthResult.data);
  } else {
    console.log('❌ Health check failed:', healthResult.error);
  }

  // Test 2: List Tools
  console.log('\n2. Testing tools listing...');
  const toolsResult = await testEndpoint('/tools/list');
  if (toolsResult.success) {
    console.log('✅ Tools listing successful');
    console.log(`   Found ${toolsResult.data.tools.length} tools`);
    toolsResult.data.tools.forEach(tool => {
      console.log(`   - ${tool.name}: ${tool.description}`);
    });
  } else {
    console.log('❌ Tools listing failed:', toolsResult.error);
  }

  // Test 3: Sequence Alignment
  console.log('\n3. Testing sequence alignment...');
  const alignmentResult = await testEndpoint('/tools/sequence-alignment', 'POST', {
    sequences: '>seq1\nACTGTTGAC\n>seq2\nACTGCATCC',
    algorithm: 'clustal'
  });
  if (alignmentResult.success) {
    console.log('✅ Sequence alignment successful');
    console.log('   Result:', alignmentResult.data.result?.text || 'No text output');
  } else {
    console.log('❌ Sequence alignment failed:', alignmentResult.error);
  }

  // Test 4: Sequence Mutation
  console.log('\n4. Testing sequence mutation...');
  const mutationResult = await testEndpoint('/tools/mutate-sequence', 'POST', {
    sequence: 'ACTGTTGAC',
    num_variants: 5,
    mutation_rate: 0.1
  });
  if (mutationResult.success) {
    console.log('✅ Sequence mutation successful');
    console.log('   Result:', mutationResult.data.result?.text || 'No text output');
  } else {
    console.log('❌ Sequence mutation failed:', mutationResult.error);
  }

  // Test 5: Legacy Command Execution
  console.log('\n5. Testing legacy command execution...');
  const legacyResult = await testEndpoint('/execute', 'POST', {
    command: 'align the given sequences'
  });
  if (legacyResult.success) {
    console.log('✅ Legacy command execution successful');
  } else {
    console.log('❌ Legacy command execution failed:', legacyResult.error);
  }

  console.log('\n' + '='.repeat(60));
  console.log('🎉 Integration test completed!');
  
  // Summary
  const tests = [
    { name: 'Health Check', result: healthResult.success },
    { name: 'Tools Listing', result: toolsResult.success },
    { name: 'Sequence Alignment', result: alignmentResult.success },
    { name: 'Sequence Mutation', result: mutationResult.success },
    { name: 'Legacy Commands', result: legacyResult.success }
  ];

  const passed = tests.filter(t => t.result).length;
  const total = tests.length;

  console.log(`\n📊 Test Results: ${passed}/${total} tests passed`);
  
  if (passed === total) {
    console.log('🎯 All tests passed! The frontend-backend integration is working correctly.');
  } else {
    console.log('⚠️  Some tests failed. Please check the backend server and try again.');
  }
}

// Check if fetch is available (Node.js 18+)
if (typeof fetch === 'undefined') {
  console.error('❌ This script requires Node.js 18+ or a fetch polyfill');
  process.exit(1);
}

// Run the tests
runTests().catch(error => {
  console.error('❌ Test execution failed:', error);
  process.exit(1);
}); 