#!/usr/bin/env python3
"""
Diagnose CDK deployment errors by checking CloudFormation events
"""

import boto3
import sys
from datetime import datetime

def diagnose_stack_error(stack_name, region):
    cf = boto3.client('cloudformation', region_name=region)
    
    try:
        # Check if stack exists
        response = cf.describe_stacks(StackName=stack_name)
        stack = response['Stacks'][0]
        print(f"✅ Stack exists: {stack_name}")
        print(f"   Status: {stack['StackStatus']}")
        if 'StackStatusReason' in stack:
            print(f"   Reason: {stack['StackStatusReason']}")
        print()
        
    except cf.exceptions.ClientError as e:
        if 'does not exist' in str(e):
            print(f"ℹ️  Stack does not exist: {stack_name}")
            print("   This is normal for first deployment")
            return
        else:
            print(f"❌ Error checking stack: {e}")
            return
    
    # Get recent stack events
    print("📋 Recent Stack Events (most recent first):")
    print("-" * 100)
    
    try:
        events_response = cf.describe_stack_events(StackName=stack_name)
        events = events_response['StackEvents'][:20]  # Last 20 events
        
        for event in events:
            timestamp = event['Timestamp'].strftime('%Y-%m-%d %H:%M:%S')
            resource_type = event.get('ResourceType', 'N/A')
            logical_id = event.get('LogicalResourceId', 'N/A')
            status = event.get('ResourceStatus', 'N/A')
            reason = event.get('ResourceStatusReason', '')
            
            # Highlight failures
            if 'FAILED' in status or 'failed' in reason.lower():
                print(f"❌ {timestamp} | {logical_id}")
                print(f"   Type: {resource_type}")
                print(f"   Status: {status}")
                if reason:
                    print(f"   Reason: {reason}")
                print()
            elif 'COMPLETE' in status:
                print(f"✅ {timestamp} | {logical_id} | {status}")
            else:
                print(f"⏳ {timestamp} | {logical_id} | {status}")
                
    except Exception as e:
        print(f"❌ Error fetching events: {e}")
        
    # Check for existing resources that might conflict
    print("\n🔍 Checking for potential conflicts:")
    print("-" * 100)
    
    ecr = boto3.client('ecr', region_name=region)
    s3 = boto3.client('s3', region_name=region)
    
    # Check ECR repository
    try:
        repos = ecr.describe_repositories(repositoryNames=['helix-ai-backend'])
        if repos['repositories']:
            print(f"⚠️  ECR repository 'helix-ai-backend' already exists")
            print(f"   URI: {repos['repositories'][0]['repositoryUri']}")
            print(f"   This might cause a conflict if not managed by the stack")
    except ecr.exceptions.RepositoryNotFoundException:
        print(f"✅ ECR repository 'helix-ai-backend' does not exist (good for new deployment)")
    except Exception as e:
        print(f"ℹ️  Could not check ECR: {e}")
    
    # Check S3 bucket
    account_id = boto3.client('sts').get_caller_identity()['Account']
    bucket_name = f"helix-ai-frontend-{account_id}-{region}"
    try:
        s3.head_bucket(Bucket=bucket_name)
        print(f"⚠️  S3 bucket '{bucket_name}' already exists")
        print(f"   This might cause a conflict if not managed by the stack")
    except:
        print(f"✅ S3 bucket '{bucket_name}' does not exist (good for new deployment)")

if __name__ == '__main__':
    stack_name = 'HelixAIStack'
    region = 'us-west-1'
    
    print(f"🔍 Diagnosing deployment error for {stack_name} in {region}")
    print("=" * 100)
    print()
    
    diagnose_stack_error(stack_name, region)




