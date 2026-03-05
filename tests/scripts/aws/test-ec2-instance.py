#!/usr/bin/env python3
"""
Test script to verify EC2 instance configuration.
Tests S3 access, Python packages, and bioinformatics tools.
"""

import os
import sys
import subprocess
import shutil
import tempfile
from pathlib import Path

def test_python_version():
    """Test Python version."""
    print("\n" + "=" * 60)
    print("🐍 Testing Python Version")
    print("=" * 60)
    
    import sys
    version = sys.version_info
    print(f"Python version: {version.major}.{version.minor}.{version.micro}")
    
    if version.major == 3 and version.minor >= 8:
        print("✅ Python version is 3.8+ (compatible with boto3)")
        return True
    else:
        print(f"⚠️  Python version is {version.major}.{version.minor} (boto3 recommends 3.8+)")
        return False

def test_python_packages():
    """Test required Python packages."""
    print("\n" + "=" * 60)
    print("📦 Testing Python Packages")
    print("=" * 60)
    
    required_packages = {
        'boto3': 'AWS SDK',
        'botocore': 'AWS SDK core',
        'Bio': 'BioPython',  # biopython is imported as Bio
        'pandas': 'Pandas',
        'numpy': 'NumPy'
    }
    
    all_ok = True
    for package, description in required_packages.items():
        try:
            module = __import__(package)
            version = getattr(module, '__version__', 'unknown')
            print(f"✅ {package} ({description}): {version}")
        except ImportError:
            print(f"❌ {package} ({description}): NOT INSTALLED")
            all_ok = False
    
    return all_ok

def test_aws_credentials():
    """Test AWS credentials."""
    print("\n" + "=" * 60)
    print("🔐 Testing AWS Credentials")
    print("=" * 60)
    
    try:
        import boto3
        from botocore.exceptions import NoCredentialsError, ClientError
    except ImportError:
        print("❌ boto3 not installed - cannot test AWS credentials")
        return False
    
    try:
        # Try to get caller identity
        sts = boto3.client('sts')
        identity = sts.get_caller_identity()
        print(f"✅ AWS credentials found")
        print(f"   Account: {identity.get('Account', 'N/A')}")
        print(f"   User/Role: {identity.get('Arn', 'N/A')}")
        
        # Try S3 access
        s3 = boto3.client('s3')
        # List buckets to test S3 access
        try:
            buckets = s3.list_buckets()
            print(f"✅ S3 access working (can list buckets)")
            if buckets.get('Buckets'):
                print(f"   Found {len(buckets['Buckets'])} bucket(s)")
            return True
        except ClientError as e:
            error_code = e.response.get('Error', {}).get('Code', '')
            if error_code == 'AccessDenied':
                print(f"⚠️  S3 access denied (may need IAM permissions)")
            else:
                print(f"⚠️  S3 access error: {error_code}")
            return False
            
    except NoCredentialsError:
        print("❌ AWS credentials not found")
        print("   The instance needs an IAM role with S3 permissions")
        return False
    except Exception as e:
        print(f"❌ Error testing AWS credentials: {e}")
        return False

def test_bioinformatics_tools():
    """Test bioinformatics tools availability."""
    print("\n" + "=" * 60)
    print("🧬 Testing Bioinformatics Tools")
    print("=" * 60)
    
    tools = {
        'bbmerge.sh': 'BBTools - Read merging',
        'samtools': 'SAMtools - Sequence alignment',
        'bcftools': 'BCFtools - Variant calling',
        'bwa': 'BWA - Sequence alignment',
        'bowtie2': 'Bowtie2 - Sequence alignment',
        'fastqc': 'FastQC - Quality control'
    }
    
    all_ok = True
    for tool, description in tools.items():
        path = shutil.which(tool)
        if path:
            # Try to get version
            try:
                result = subprocess.run(
                    [tool, '--version'] if tool != 'bbmerge.sh' else [tool],
                    capture_output=True,
                    text=True,
                    timeout=5
                )
                version_info = result.stdout.strip()[:50] if result.stdout else "available"
                print(f"✅ {tool} ({description}): {path}")
                if version_info:
                    print(f"   Version info: {version_info}")
            except Exception:
                print(f"✅ {tool} ({description}): {path}")
        else:
            print(f"❌ {tool} ({description}): NOT FOUND")
            all_ok = False
    
    return all_ok

def test_s3_operations():
    """Test S3 upload and download."""
    print("\n" + "=" * 60)
    print("☁️  Testing S3 Operations")
    print("=" * 60)
    
    try:
        import boto3
        from botocore.exceptions import ClientError
        
        s3 = boto3.client('s3')
        region = os.getenv('AWS_REGION', 'us-west-1')
        
        # Use a test bucket (you may need to adjust this)
        test_bucket = os.getenv('HELIX_TEST_BUCKET', 'noricum-ngs-data')
        test_key = 'helix-test/test_file.txt'
        
        # Create test file
        test_content = "Helix.AI EC2 test file\nThis file was created to test S3 access."
        
        print(f"Testing with bucket: {test_bucket}")
        print(f"Test key: {test_key}")
        
        # Test 1: List bucket (read permission)
        print("\n1. Testing S3 list bucket...")
        try:
            s3.head_bucket(Bucket=test_bucket)
            print(f"   ✅ Can access bucket: {test_bucket}")
        except ClientError as e:
            error_code = e.response.get('Error', {}).get('Code', '')
            if error_code == '404':
                print(f"   ⚠️  Bucket {test_bucket} does not exist")
                print(f"   Set HELIX_TEST_BUCKET environment variable to a valid bucket")
                return False
            elif error_code == '403':
                print(f"   ❌ Access denied to bucket {test_bucket}")
                print(f"   Check IAM role permissions")
                return False
            else:
                print(f"   ⚠️  Cannot access bucket: {error_code}")
                return False
        
        # Test 2: Upload
        try:
            with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
                f.write(test_content)
                temp_file = f.name
            
            print("\n2. Testing S3 upload...")
            s3.upload_file(temp_file, test_bucket, test_key)
            print(f"   ✅ Successfully uploaded to s3://{test_bucket}/{test_key}")
            
            # Test 3: Download
            print("\n3. Testing S3 download...")
            download_path = temp_file + '.downloaded'
            s3.download_file(test_bucket, test_key, download_path)
            
            # Verify content
            with open(download_path, 'r') as f:
                downloaded_content = f.read()
            
            if downloaded_content == test_content:
                print(f"   ✅ Successfully downloaded and verified content")
            else:
                print(f"   ⚠️  Downloaded but content mismatch")
            
            # Cleanup
            os.unlink(temp_file)
            os.unlink(download_path)
            
            # Test 4: Delete
            try:
                s3.delete_object(Bucket=test_bucket, Key=test_key)
                print(f"\n4. Testing S3 delete...")
                print(f"   ✅ Successfully deleted test file from S3")
            except:
                print(f"\n4. Testing S3 delete...")
                print(f"   ⚠️  Could not delete test file (may need delete permission)")
            
            return True
            
        except ClientError as e:
            error_code = e.response.get('Error', {}).get('Code', '')
            if error_code == 'AccessDenied':
                print(f"   ❌ Access denied - check IAM role permissions")
            else:
                print(f"   ❌ S3 operation failed: {error_code}")
            return False
            
    except Exception as e:
        print(f"❌ Error testing S3 operations: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_tool_execution():
    """Test executing bioinformatics tools."""
    print("\n" + "=" * 60)
    print("🔧 Testing Tool Execution")
    print("=" * 60)
    
    tools_to_test = [
        ('samtools', ['--version']),
        ('bcftools', ['--version']),
        ('bwa', []),
        ('bowtie2', ['--version']),
        ('fastqc', ['--version']),
    ]
    
    all_ok = True
    for tool, args in tools_to_test:
        tool_path = shutil.which(tool)
        if tool_path:
            try:
                result = subprocess.run(
                    [tool] + args,
                    capture_output=True,
                    text=True,
                    timeout=10
                )
                if result.returncode == 0 or tool == 'bwa':  # bwa may return non-zero for --version
                    version_info = result.stdout.strip()[:80] if result.stdout else "executable"
                    print(f"✅ {tool}: {version_info}")
                else:
                    print(f"⚠️  {tool}: returned exit code {result.returncode}")
                    if result.stderr:
                        print(f"   Error: {result.stderr[:100]}")
            except subprocess.TimeoutExpired:
                print(f"⚠️  {tool}: execution timed out")
            except Exception as e:
                print(f"❌ {tool}: error - {e}")
                all_ok = False
        else:
            print(f"⚠️  {tool}: not found")
    
    # Test bbmerge.sh (special case - may need input files)
    if shutil.which('bbmerge.sh'):
        print(f"✅ bbmerge.sh: found at {shutil.which('bbmerge.sh')}")
    else:
        print(f"⚠️  bbmerge.sh: not found")
        all_ok = False
    
    return all_ok

def main():
    print("=" * 60)
    print("Helix.AI EC2 Instance Test Suite")
    print("=" * 60)
    print("\nThis script tests:")
    print("  - Python version and packages")
    print("  - AWS credentials and S3 access")
    print("  - Bioinformatics tools availability")
    print("  - Tool execution")
    
    results = {}
    
    # Run tests
    results['python_version'] = test_python_version()
    results['python_packages'] = test_python_packages()
    results['aws_credentials'] = test_aws_credentials()
    results['bioinformatics_tools'] = test_bioinformatics_tools()
    results['s3_operations'] = test_s3_operations()
    results['tool_execution'] = test_tool_execution()
    
    # Summary
    print("\n" + "=" * 60)
    print("📊 Test Summary")
    print("=" * 60)
    
    all_passed = True
    for test_name, passed in results.items():
        status = "✅ PASS" if passed else "❌ FAIL"
        print(f"{status}: {test_name.replace('_', ' ').title()}")
        if not passed:
            all_passed = False
    
    print("\n" + "=" * 60)
    if all_passed:
        print("✅ All tests passed! EC2 instance is ready.")
    else:
        print("❌ Some tests failed. See details above.")
        print("\n💡 To fix issues:")
        if not results.get('python_packages'):
            print("   - Run: ./scripts/aws/install-conda-on-instance.sh")
        if not results.get('aws_credentials'):
            print("   - Run: ./scripts/aws/attach-iam-role-to-instance.sh")
        if not results.get('bioinformatics_tools'):
            print("   - Run: ./scripts/aws/install-conda-on-instance.sh")
    print("=" * 60)
    
    return 0 if all_passed else 1

if __name__ == '__main__':
    sys.exit(main())

