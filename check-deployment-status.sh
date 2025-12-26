#!/bin/bash
# Quick script to check CDK deployment status

STACK_NAME="HelixAIStack"
REGION="us-west-1"

echo "📊 Checking Deployment Status"
echo "=============================="
echo ""

# Check CloudFormation Stack Status
echo "🔹 CloudFormation Stack Status:"
aws cloudformation describe-stacks \
  --stack-name $STACK_NAME \
  --region $REGION \
  --query 'Stacks[0].{Status:StackStatus,LastUpdate:LastUpdatedTime}' \
  --output table

echo ""

# Check ECS Service Status
echo "🔹 ECS Service Status:"
SERVICE_NAME=$(aws cloudformation describe-stack-resources \
  --stack-name $STACK_NAME \
  --region $REGION \
  --query 'StackResources[?LogicalResourceId==`ServiceD69D759B`].PhysicalResourceId' \
  --output text 2>/dev/null | cut -d/ -f3)

if [ -n "$SERVICE_NAME" ]; then
  CLUSTER_NAME=$(aws cloudformation describe-stack-resources \
    --stack-name $STACK_NAME \
    --region $REGION \
    --query 'StackResources[?LogicalResourceId==`ClusterEB0386A7`].PhysicalResourceId' \
    --output text 2>/dev/null | cut -d/ -f2)
  
  if [ -n "$CLUSTER_NAME" ]; then
    aws ecs describe-services \
      --cluster $CLUSTER_NAME \
      --services $SERVICE_NAME \
      --region $REGION \
      --query 'services[0].{Status:status,DesiredCount:desiredCount,RunningCount:runningCount,PendingCount:pendingCount,FailedTasks:deployments[0].failedTasks,RolloutState:deployments[0].rolloutState}' \
      --output table
  fi
else
  echo "  Service not found yet"
fi

echo ""

# Check Recent Stack Events (last 5)
echo "🔹 Recent Stack Events:"
aws cloudformation describe-stack-events \
  --stack-name $STACK_NAME \
  --region $REGION \
  --max-items 5 \
  --query 'StackEvents[*].{Time:Timestamp,LogicalId:LogicalResourceId,Status:ResourceStatus}' \
  --output table

echo ""
echo "💡 Tip: Run this script repeatedly to monitor progress"
