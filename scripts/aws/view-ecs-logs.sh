#!/bin/bash
# View ECS service logs - convenient script for monitoring backend logs

set -e

LOG_GROUP="/ecs/helix-ai"
REGION="us-west-1"
FOLLOW="${1:-follow}"  # Default to follow mode

if [ "$FOLLOW" = "follow" ] || [ "$FOLLOW" = "-f" ] || [ "$FOLLOW" = "--follow" ]; then
    echo "📊 Following logs from $LOG_GROUP (press Ctrl+C to stop)..."
    echo ""
    aws logs tail "$LOG_GROUP" --follow --region "$REGION" --format short
elif [ "$FOLLOW" = "recent" ] || [ "$FOLLOW" = "-r" ] || [ "$FOLLOW" = "--recent" ]; then
    echo "📊 Showing recent logs from $LOG_GROUP (last 30 minutes)..."
    echo ""
    aws logs tail "$LOG_GROUP" --region "$REGION" --since 30m --format short
elif [ "$FOLLOW" = "filter" ] || [ "$FOLLOW" = "-g" ] || [ "$FOLLOW" = "--grep" ]; then
    echo "📊 Following logs filtered for commands and errors..."
    echo ""
    aws logs tail "$LOG_GROUP" --follow --region "$REGION" --format short | grep -E "\[handle_command\]|command|ERROR|error|tool|Tool|WARNING|warning"
else
    echo "Usage: $0 [follow|recent|filter]"
    echo ""
    echo "  follow  - Follow logs in real-time (default)"
    echo "  recent  - Show recent logs from last 30 minutes"
    echo "  filter  - Follow logs filtered for commands/errors"
    echo ""
    echo "Examples:"
    echo "  $0           # Follow logs (default)"
    echo "  $0 follow    # Follow logs"
    echo "  $0 recent    # Show recent logs"
    echo "  $0 filter    # Follow filtered logs"
fi


