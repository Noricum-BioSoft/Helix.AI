#!/bin/bash
# Monitor for new EC2 instance creation
source deploy.config
STACK_NAME="${STACK_NAME:-HelixAIStack}"

echo "Waiting for new EC2 instance to be created..."
while true; do
    INSTANCE_ID=$(aws cloudformation describe-stacks \
        --stack-name "${STACK_NAME}" \
        --region "${AWS_REGION}" \
        --query "Stacks[0].Outputs[?OutputKey=='EC2InstanceId'].OutputValue" \
        --output text 2>/dev/null || echo "")
    
    if [ -n "${INSTANCE_ID}" ] && [ "${INSTANCE_ID}" != "None" ]; then
        STATE=$(aws ec2 describe-instances --instance-ids "${INSTANCE_ID}" --region "${AWS_REGION}" --query "Reservations[0].Instances[0].State.Name" --output text 2>/dev/null || echo "not-found")
        echo "Found instance: ${INSTANCE_ID} (State: ${STATE})"
        if [ "${STATE}" = "running" ]; then
            echo "✅ New instance is running!"
            echo ""
            echo "Instance ID: ${INSTANCE_ID}"
            echo "The instance should automatically:"
            echo "  1. Run user data script (with dnf, not apt-get)"
            echo "  2. Download Docker image from S3"
            echo "  3. Start the container"
            echo ""
            echo "Wait 2-3 minutes for user data to complete, then check:"
            echo "  aws ssm start-session --target ${INSTANCE_ID} --region ${AWS_REGION}"
            break
        fi
    fi
    echo "Waiting... (checking every 10 seconds)"
    sleep 10
done
