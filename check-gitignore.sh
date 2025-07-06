#!/bin/bash

# Check for unwanted files that should be ignored by git
echo "🔍 Checking for unwanted files in git repository..."
echo "=================================================="

# Check if node_modules exists and is tracked
if [ -d "frontend/node_modules" ]; then
    echo "⚠️  node_modules directory found in frontend/"
    if git ls-files frontend/node_modules >/dev/null 2>&1; then
        echo "❌ node_modules is being tracked by git!"
        echo "   Run: git rm -r --cached frontend/node_modules"
    else
        echo "✅ node_modules is properly ignored"
    fi
else
    echo "✅ No node_modules directory found"
fi

# Check for other common unwanted files
unwanted_files=(
    "frontend/dist"
    "frontend/.vite"
    "frontend/build"
    "backend/__pycache__"
    "backend/*.pyc"
    "backend/venv"
    "backend/.env"
    "backend/logs"
    "*.log"
    ".DS_Store"
)

echo ""
echo "Checking for other unwanted files:"
for file in "${unwanted_files[@]}"; do
    if [ -e "$file" ] || [ -d "$file" ]; then
        if git ls-files "$file" >/dev/null 2>&1; then
            echo "❌ $file is being tracked by git!"
        else
            echo "✅ $file is properly ignored"
        fi
    else
        echo "✅ $file does not exist"
    fi
done

# Check for large files that might be accidentally committed
echo ""
echo "Checking for large files (>10MB):"
find . -type f -size +10M -not -path "./.git/*" -not -path "./node_modules/*" -not -path "./frontend/node_modules/*" 2>/dev/null | while read file; do
    if git ls-files "$file" >/dev/null 2>&1; then
        size=$(du -h "$file" | cut -f1)
        echo "⚠️  Large file tracked: $file ($size)"
    fi
done

# Check gitignore files
echo ""
echo "Checking .gitignore files:"
if [ -f ".gitignore" ]; then
    echo "✅ Root .gitignore exists"
else
    echo "❌ Root .gitignore missing"
fi

if [ -f "frontend/.gitignore" ]; then
    echo "✅ Frontend .gitignore exists"
else
    echo "❌ Frontend .gitignore missing"
fi

if [ -f "backend/.gitignore" ]; then
    echo "✅ Backend .gitignore exists"
else
    echo "❌ Backend .gitignore missing"
fi

echo ""
echo "🎯 Recommendations:"
echo "1. If any files are being tracked that shouldn't be, run:"
echo "   git rm --cached <file>"
echo "2. Commit the .gitignore changes:"
echo "   git add .gitignore frontend/.gitignore backend/.gitignore"
echo "3. Commit the removal of unwanted files:"
echo "   git commit -m 'Remove unwanted files and update .gitignore'"
echo ""
echo "✅ Git ignore check completed!" 