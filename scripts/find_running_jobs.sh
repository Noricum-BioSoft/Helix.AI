#!/bin/bash
# Find running EMR jobs (when backend was restarted and job IDs lost)

REGION="${AWS_REGION:-us-west-1}"

echo "=========================================="
echo "Finding Running EMR Jobs"
echo "=========================================="
echo ""

# Find active EMR clusters
echo "Looking for active EMR clusters..."
CLUSTERS=$(aws emr list-clusters \
    --region "$REGION" \
    --active \
    --query 'Clusters[*].[Id,Name,Status.State]' \
    --output text 2>/dev/null)

if [ -z "$CLUSTERS" ]; then
    echo "❌ No active EMR clusters found"
    exit 0
fi

echo "Active clusters:"
echo "$CLUSTERS"
echo ""

# For each cluster, list running/pending steps
while IFS=$'\t' read -r CLUSTER_ID CLUSTER_NAME CLUSTER_STATE; do
    echo "=========================================="
    echo "Cluster: $CLUSTER_NAME ($CLUSTER_ID)"
    echo "State: $CLUSTER_STATE"
    echo "=========================================="
    
    # List steps
    STEPS=$(aws emr list-steps \
        --cluster-id "$CLUSTER_ID" \
        --region "$REGION" \
        --step-states PENDING RUNNING \
        --query 'Steps[*].[Id,Name,Status.State,CreationDateTime]' \
        --output text 2>/dev/null)
    
    if [ -z "$STEPS" ]; then
        echo "No running/pending steps"
    else
        echo ""
        echo "Running/Pending Steps:"
        echo "----------------------------------------------------------------------"
        printf "%-20s %-40s %-12s %s\n" "STEP ID" "NAME" "STATE" "CREATED"
        echo "----------------------------------------------------------------------"
        while IFS=$'\t' read -r STEP_ID STEP_NAME STEP_STATE CREATED_AT; do
            printf "%-20s %-40s %-12s %s\n" "$STEP_ID" "$STEP_NAME" "$STEP_STATE" "$CREATED_AT"
        done <<< "$STEPS"
        
        echo ""
        echo "📋 To check step details:"
        while IFS=$'\t' read -r STEP_ID _; do
            echo "   aws emr describe-step --cluster-id $CLUSTER_ID --step-id $STEP_ID"
        done <<< "$STEPS"
        
        echo ""
        echo "📋 To view step logs:"
        while IFS=$'\t' read -r STEP_ID _; do
            echo "   aws emr get-step-logs --cluster-id $CLUSTER_ID --step-id $STEP_ID"
        done <<< "$STEPS"
    fi
    
    echo ""
done <<< "$CLUSTERS"

# Also list recently completed steps (last 5)
echo "=========================================="
echo "Recently Completed Steps (Last 5)"
echo "=========================================="
echo ""

while IFS=$'\t' read -r CLUSTER_ID CLUSTER_NAME CLUSTER_STATE; do
    COMPLETED=$(aws emr list-steps \
        --cluster-id "$CLUSTER_ID" \
        --region "$REGION" \
        --step-states COMPLETED FAILED CANCELLED \
        --query 'Steps[:5].[Id,Name,Status.State,Status.Timeline.EndDateTime]' \
        --output text 2>/dev/null)
    
    if [ -n "$COMPLETED" ]; then
        echo "Cluster: $CLUSTER_NAME ($CLUSTER_ID)"
        echo "----------------------------------------------------------------------"
        printf "%-20s %-40s %-12s %s\n" "STEP ID" "NAME" "STATE" "ENDED"
        echo "----------------------------------------------------------------------"
        while IFS=$'\t' read -r STEP_ID STEP_NAME STEP_STATE ENDED_AT; do
            printf "%-20s %-40s %-12s %s\n" "$STEP_ID" "$STEP_NAME" "$STEP_STATE" "$ENDED_AT"
        done <<< "$COMPLETED"
        echo ""
    fi
done <<< "$CLUSTERS"
