#!/usr/bin/env python3
"""
Test script to verify EC2 setup and optionally create an instance.
"""

import os
import sys
import boto3
from pathlib import Path

# Add project root to path
project_root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(project_root / "backend"))

from ec2_executor import EC2Executor, get_ec2_executor

def check_requirements():
    """Check if required environment variables are set."""
    print("🔍 Checking EC2 setup requirements...\n")
    
    required_vars = {
        'HELIX_EC2_KEY_NAME': 'EC2 Key Pair Name',
        'HELIX_EC2_KEY_FILE': 'Path to SSH private key file',
    }
    
    optional_vars = {
        'HELIX_EC2_INSTANCE_ID': 'Existing EC2 Instance ID',
        'HELIX_EC2_INSTANCE_TYPE': 'Instance Type (default: t3.medium)',
        'HELIX_EC2_AMI_ID': 'Custom AMI ID',
        'HELIX_EC2_AUTO_CREATE': 'Auto-create if instance not found (true/false)',
        'AWS_REGION': 'AWS Region (default: us-east-1)',
    }
    
    missing_required = []
    for var, desc in required_vars.items():
        value = os.getenv(var)
        if value:
            print(f"✅ {var}: {value}")
        else:
            print(f"❌ {var}: NOT SET - {desc}")
            missing_required.append(var)
    
    print("\n📋 Optional variables:")
    for var, desc in optional_vars.items():
        value = os.getenv(var)
        if value:
            print(f"   ✅ {var}: {value}")
        else:
            print(f"   ⚪ {var}: not set ({desc})")
    
    # Check if key file exists
    key_file = os.getenv('HELIX_EC2_KEY_FILE')
    if key_file:
        key_path = Path(key_file).expanduser()
        if key_path.exists():
            # Check permissions
            stat = key_path.stat()
            mode = oct(stat.st_mode)[-3:]
            if mode != '600':
                print(f"\n⚠️  Warning: Key file permissions are {mode}, should be 600")
                print(f"   Run: chmod 600 {key_file}")
            else:
                print(f"✅ Key file exists and has correct permissions (600)")
        else:
            print(f"\n❌ Key file not found: {key_file}")
            missing_required.append('HELIX_EC2_KEY_FILE_EXISTS')
    
    # Check AWS credentials
    try:
        session = boto3.Session()
        credentials = session.get_credentials()
        if credentials:
            print(f"\n✅ AWS credentials found")
        else:
            print(f"\n❌ AWS credentials not found")
            missing_required.append('AWS_CREDENTIALS')
    except Exception as e:
        print(f"\n❌ Error checking AWS credentials: {e}")
        missing_required.append('AWS_CREDENTIALS')
    
    return len(missing_required) == 0, missing_required

def test_instance_creation():
    """Test creating an EC2 instance."""
    print("\n🚀 Testing EC2 instance creation...\n")
    
    executor = get_ec2_executor()
    
    # Check for existing instances
    print("🔍 Checking for existing instances...")
    existing = executor._find_existing_instance()
    if existing:
        print(f"✅ Found existing instance: {existing}")
        return existing
    
    # Check configured instance ID
    configured_id = os.getenv('HELIX_EC2_INSTANCE_ID')
    if configured_id:
        print(f"🔍 Checking configured instance: {configured_id}")
        status = executor._check_instance_status(configured_id)
        if status == "running":
            print(f"✅ Configured instance is running: {configured_id}")
            return configured_id
        elif status == "not_found":
            print(f"⚠️  Configured instance does not exist: {configured_id}")
            # Check if auto-create is enabled (accepts 'true', '1', 'yes', 'on', or any non-empty value)
            auto_create_val = os.getenv('HELIX_EC2_AUTO_CREATE', 'false').lower().strip()
            auto_create = auto_create_val in ('true', '1', 'yes', 'on', 'auto-create') or (auto_create_val and auto_create_val != 'false')
            if auto_create:
                print("✅ Auto-create enabled, will create new instance")
            else:
                print("❌ Auto-create disabled. Set HELIX_EC2_AUTO_CREATE=true (or '1', 'yes', 'on') to auto-create")
                return None
        else:
            print(f"⚠️  Configured instance status: {status}")
    
    # Try to get or create instance
    print("\n🔧 Attempting to get or create instance...")
    try:
        instance_id = executor.get_or_create_instance()
        if instance_id:
            print(f"✅ Success! Instance ID: {instance_id}")
            print(f"\n📝 Note: Tools installation takes 10-15 minutes.")
            print(f"   Check status: ssh -i {os.getenv('HELIX_EC2_KEY_FILE')} ec2-user@<instance-ip>")
            print(f"   Then run: cat /opt/helix-tools/install.log")
            return instance_id
        else:
            print("❌ Failed to get or create instance")
            return None
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()
        return None

def main():
    print("=" * 60)
    print("Helix.AI EC2 Setup Test")
    print("=" * 60)
    
    # Check requirements
    all_good, missing = check_requirements()
    
    if not all_good:
        print("\n❌ Missing required configuration:")
        for item in missing:
            print(f"   - {item}")
        print("\n💡 Set the required environment variables and try again.")
        return 1
    
    # Test instance creation
    instance_id = test_instance_creation()
    
    if instance_id:
        print("\n" + "=" * 60)
        print("✅ Setup complete!")
        print("=" * 60)
        print(f"\nInstance ID: {instance_id}")
        print("\nNext steps:")
        print("1. Wait 10-15 minutes for tools to install")
        print("2. Set HELIX_EC2_INSTANCE_ID if not already set:")
        print(f"   export HELIX_EC2_INSTANCE_ID={instance_id}")
        print("3. Enable EC2 execution:")
        print("   export HELIX_USE_EC2=true")
        print("4. Restart your backend server")
        return 0
    else:
        print("\n" + "=" * 60)
        print("❌ Setup incomplete")
        print("=" * 60)
        return 1

if __name__ == '__main__':
    sys.exit(main())

