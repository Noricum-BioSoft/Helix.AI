#!/usr/bin/env python3
"""
Display information about the configured EC2 instance.
"""

import os
import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(project_root / "backend"))

from ec2_executor import get_ec2_executor

def main():
    print("=" * 60)
    print("Helix.AI EC2 Instance Information")
    print("=" * 60)
    
    # Check required environment variables
    key_name = os.getenv('HELIX_EC2_KEY_NAME')
    key_file = os.getenv('HELIX_EC2_KEY_FILE')
    instance_id = os.getenv('HELIX_EC2_INSTANCE_ID')
    
    if not key_name or not key_file:
        print("\n❌ Missing required environment variables:")
        if not key_name:
            print("   - HELIX_EC2_KEY_NAME")
        if not key_file:
            print("   - HELIX_EC2_KEY_FILE")
        print("\n💡 Set these variables and try again.")
        return 1
    
    executor = get_ec2_executor()
    
    # If instance ID is set, show that instance
    if instance_id:
        print(f"\n🔍 Showing configured instance: {instance_id}\n")
        details = executor.get_instance_details(instance_id)
        if details:
            display_instance_info(executor, instance_id)
        else:
            print(f"❌ Could not retrieve instance: {instance_id}")
            print("   The instance may not exist or you may not have permissions.")
            return 1
    else:
        # Try to find existing instances
        print("\n🔍 Searching for Helix.AI instances...\n")
        existing = executor._find_existing_instance()
        if existing:
            print(f"✅ Found instance: {existing}\n")
            display_instance_info(executor, existing)
        else:
            print("❌ No EC2 instance found.")
            print("\n💡 Options:")
            print("   1. Set HELIX_EC2_INSTANCE_ID to a specific instance ID")
            print("   2. Run: python scripts/aws/test-ec2-setup.py")
            print("   3. Create an instance with: python scripts/aws/test-ec2-setup.py")
            return 1
    
    return 0

def display_instance_info(executor, instance_id: str):
    """Display detailed information about an EC2 instance."""
    print("=" * 60)
    print("📊 EC2 Instance Information")
    print("=" * 60)
    
    details = executor.get_instance_details(instance_id)
    if not details:
        print(f"❌ Could not retrieve details for instance: {instance_id}")
        return
    
    print(f"\n🆔 Instance ID: {details.get('InstanceId', 'N/A')}")
    print(f"📦 Instance Type: {details.get('InstanceType', 'N/A')}")
    print(f"🟢 State: {details.get('State', 'N/A')}")
    
    public_ip = details.get('PublicIpAddress')
    private_ip = details.get('PrivateIpAddress')
    
    if public_ip:
        print(f"🌐 Public IP: {public_ip}")
        key_file = os.getenv('HELIX_EC2_KEY_FILE')
        if key_file:
            print(f"\n🔗 SSH Command:")
            print(f"   ssh -i {key_file} ec2-user@{public_ip}")
    else:
        print(f"🌐 Public IP: Not available (instance may be in private subnet)")
    
    if private_ip:
        print(f"🔒 Private IP: {private_ip}")
    
    # Get additional details
    try:
        response = executor.ec2_client.describe_instances(InstanceIds=[instance_id])
        if response['Reservations']:
            instance = response['Reservations'][0]['Instances'][0]
            
            # Launch time
            launch_time = instance.get('LaunchTime')
            if launch_time:
                print(f"⏰ Launched: {launch_time}")
            
            # Tags
            tags = {tag['Key']: tag['Value'] for tag in instance.get('Tags', [])}
            if tags:
                print(f"\n🏷️  Tags:")
                for key, value in tags.items():
                    print(f"   {key}: {value}")
            
            # Security groups
            security_groups = instance.get('SecurityGroups', [])
            if security_groups:
                print(f"\n🔐 Security Groups:")
                for sg in security_groups:
                    print(f"   - {sg.get('GroupName', 'N/A')} ({sg.get('GroupId', 'N/A')})")
            
            # VPC/Subnet
            vpc_id = instance.get('VpcId')
            subnet_id = instance.get('SubnetId')
            if vpc_id:
                print(f"\n🌐 Network:")
                print(f"   VPC: {vpc_id}")
                if subnet_id:
                    print(f"   Subnet: {subnet_id}")
            
            # Check if tools are installed
            if public_ip and details.get('State') == 'running':
                print(f"\n🔧 Tools Installation Status:")
                print(f"   Note: Tools installation takes 10-15 minutes after instance launch")
                key_file = os.getenv('HELIX_EC2_KEY_FILE')
                if key_file:
                    print(f"   Check status: ssh -i {key_file} ec2-user@{public_ip} 'cat /opt/helix-tools/install.log'")
    
    except Exception as e:
        print(f"\n⚠️  Could not retrieve additional details: {e}")
    
    print("=" * 60)

if __name__ == '__main__':
    sys.exit(main())

