#!/usr/bin/env python3
"""
Check SSH access to EC2 instance and security group configuration.
"""

import os
import sys
import boto3
from pathlib import Path

# Add project root to path
project_root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(project_root / "backend"))

def get_my_ip():
    """Get current public IP address."""
    try:
        import urllib.request
        ip = urllib.request.urlopen('https://checkip.amazonaws.com').read().decode('utf-8').strip()
        return ip
    except Exception as e:
        print(f"⚠️  Could not determine your IP: {e}")
        return None

def check_security_group(ec2_client, sg_id, my_ip):
    """Check security group SSH rules."""
    try:
        response = ec2_client.describe_security_groups(GroupIds=[sg_id])
        if not response['SecurityGroups']:
            return False, "Security group not found"
        
        sg = response['SecurityGroups'][0]
        print(f"\n🔐 Security Group: {sg.get('GroupName', 'N/A')} ({sg_id})")
        print(f"   Description: {sg.get('Description', 'N/A')}")
        
        # Check ingress rules
        ingress_rules = sg.get('IpPermissions', [])
        ssh_allowed = False
        ssh_rules = []
        
        for rule in ingress_rules:
            if rule.get('IpProtocol') == 'tcp' and rule.get('FromPort') == 22:
                ssh_rules.append(rule)
                for ip_range in rule.get('IpRanges', []):
                    cidr = ip_range.get('CidrIp', '')
                    if cidr == '0.0.0.0/0':
                        ssh_allowed = True
                        print(f"   ✅ SSH (port 22) is open to: {cidr} (all IPs)")
                    elif my_ip:
                        # Check if my IP is in the CIDR range
                        import ipaddress
                        try:
                            if ipaddress.ip_address(my_ip) in ipaddress.ip_network(cidr, strict=False):
                                ssh_allowed = True
                                print(f"   ✅ SSH (port 22) is open to: {cidr} (includes your IP)")
                            else:
                                print(f"   ⚠️  SSH (port 22) is open to: {cidr} (does NOT include your IP)")
                        except:
                            print(f"   ⚠️  SSH (port 22) is open to: {cidr}")
        
        if not ssh_rules:
            print(f"   ❌ No SSH (port 22) rules found in security group")
            return False, "No SSH rules"
        
        if not ssh_allowed and my_ip:
            return False, f"SSH not allowed from your IP ({my_ip})"
        
        return True, "SSH access configured"
        
    except Exception as e:
        return False, f"Error checking security group: {e}"

def main():
    print("=" * 60)
    print("EC2 SSH Access Checker")
    print("=" * 60)
    
    # Get instance ID
    instance_id = os.getenv('HELIX_EC2_INSTANCE_ID')
    if not instance_id:
        print("\n❌ HELIX_EC2_INSTANCE_ID not set")
        print("   Set it with: export HELIX_EC2_INSTANCE_ID=i-xxxxx")
        return 1
    
    # Get my IP
    my_ip = get_my_ip()
    if my_ip:
        print(f"\n🌐 Your public IP: {my_ip}")
    
    # Connect to EC2
    region = os.getenv('AWS_REGION', 'us-east-1')
    ec2_client = boto3.client('ec2', region_name=region)
    
    # Get instance details
    try:
        response = ec2_client.describe_instances(InstanceIds=[instance_id])
        if not response['Reservations'] or not response['Reservations'][0]['Instances']:
            print(f"\n❌ Instance {instance_id} not found")
            return 1
        
        instance = response['Reservations'][0]['Instances'][0]
        state = instance['State']['Name']
        public_ip = instance.get('PublicIpAddress')
        private_ip = instance.get('PrivateIpAddress')
        security_groups = instance.get('SecurityGroups', [])
        
        print(f"\n🆔 Instance ID: {instance_id}")
        print(f"🟢 State: {state}")
        
        if state != 'running':
            print(f"\n❌ Instance is not running (state: {state})")
            print(f"   Start it with: aws ec2 start-instances --instance-ids {instance_id}")
            return 1
        
        if not public_ip:
            print(f"\n❌ Instance has no public IP address")
            print(f"   This instance may be in a private subnet without a public IP")
            if private_ip:
                print(f"   Private IP: {private_ip}")
            return 1
        
        print(f"🌐 Public IP: {public_ip}")
        if private_ip:
            print(f"🔒 Private IP: {private_ip}")
        
        # Check security groups
        if not security_groups:
            print(f"\n❌ Instance has no security groups")
            return 1
        
        print(f"\n🔍 Checking security group configuration...")
        all_ok = True
        for sg in security_groups:
            sg_id = sg.get('GroupId')
            ok, msg = check_security_group(ec2_client, sg_id, my_ip)
            if not ok:
                all_ok = False
                print(f"\n❌ {msg}")
        
        if all_ok:
            print(f"\n✅ Security group configuration looks good")
            print(f"\n🔗 Try SSH:")
            key_file = os.getenv('HELIX_EC2_KEY_FILE')
            if key_file:
                print(f"   ssh -i {key_file} ec2-user@{public_ip}")
            else:
                print(f"   ssh -i <your-key-file> ec2-user@{public_ip}")
        else:
            print(f"\n💡 To fix SSH access:")
            if my_ip:
                print(f"   1. Add your IP to the security group:")
                print(f"      aws ec2 authorize-security-group-ingress \\")
                print(f"        --group-id {security_groups[0]['GroupId']} \\")
                print(f"        --protocol tcp \\")
                print(f"        --port 22 \\")
                print(f"        --cidr {my_ip}/32")
            else:
                print(f"   1. Add your IP to the security group (replace YOUR_IP):")
                print(f"      aws ec2 authorize-security-group-ingress \\")
                print(f"        --group-id {security_groups[0]['GroupId']} \\")
                print(f"        --protocol tcp \\")
                print(f"        --port 22 \\")
                print(f"        --cidr YOUR_IP/32")
            print(f"\n   2. Or allow SSH from anywhere (less secure, for testing only):")
            print(f"      aws ec2 authorize-security-group-ingress \\")
            print(f"        --group-id {security_groups[0]['GroupId']} \\")
            print(f"        --protocol tcp \\")
            print(f"        --port 22 \\")
            print(f"        --cidr 0.0.0.0/0")
        
        return 0 if all_ok else 1
        
    except Exception as e:
        print(f"\n❌ Error: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == '__main__':
    sys.exit(main())

