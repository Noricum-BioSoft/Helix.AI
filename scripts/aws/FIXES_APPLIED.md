# Fixes Applied to Resolve SSM Issue

## Problems Found

1. **Missing SSM Permissions**: EC2 instance IAM role didn't have SSM permissions, so Systems Manager couldn't connect
2. **Wrong Package Manager**: User data script used `apt-get` but Amazon Linux 2023 uses `dnf`
3. **Instance Already Created**: Existing instance has the old IAM role attached (roles are immutable once attached)

## Fixes Applied

1. ✅ **Fixed user data script**: Changed `apt-get` to `dnf` for Amazon Linux 2023
2. ✅ **Added SSM permissions**: Added `AmazonSSMManagedInstanceCore` managed policy to EC2 role
3. ✅ **S3 fallback**: User data script now downloads from S3 if ECR pull fails

## Solution

Since the existing instance has the old IAM role attached, we need to **terminate and recreate** the instance to get the new role and updated user data.

**Run this to fix everything:**
```bash
cd /Users/eoberortner/git/Helix.AI/scripts/aws
./terminate-and-recreate-instance.sh deploy.config
```

This will:
1. Update the CDK stack with all fixes
2. Terminate the current instance
3. CloudFormation will auto-create a new instance with:
   - Correct IAM role (with SSM permissions)
   - Fixed user data (dnf, S3 fallback)
   - Will automatically download image from S3 on startup

## Alternative (if you don't want to terminate)

If you prefer not to terminate, you can:
1. Update the stack: `cd infrastructure && cdk deploy`
2. Manually create a new instance with the updated configuration
3. Or wait for SSM to eventually work (but it won't because of missing permissions)

The terminate-and-recreate approach is the cleanest solution.
