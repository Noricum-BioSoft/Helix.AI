#!/bin/bash
# Setup script for EC2 instance with pre-installed bioinformatics tools
# This script can be run on an EC2 instance to install all necessary tools

set -e

echo "=========================================="
echo "Helix.AI Bioinformatics Tools Setup"
echo "=========================================="
echo ""

# Update system
echo "Step 1: Updating system packages..."
yum update -y

# Install basic dependencies
echo "Step 2: Installing basic dependencies..."
yum install -y wget curl git gcc gcc-c++ make

# Install Miniconda
echo "Step 3: Installing Miniconda..."
cd /tmp
if [ ! -f "/opt/conda/bin/conda" ]; then
    wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    chmod +x Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda
    rm -f Miniconda3-latest-Linux-x86_64.sh
fi

# Add conda to PATH
export PATH="/opt/conda/bin:$PATH"

# Initialize conda
echo "Step 4: Configuring conda..."
/opt/conda/bin/conda init bash
source ~/.bashrc

# Add bioconda channels
echo "Step 5: Adding bioconda channels..."
/opt/conda/bin/conda config --add channels bioconda
/opt/conda/bin/conda config --add channels conda-forge
/opt/conda/bin/conda config --set channel_priority strict

# Install bioinformatics tools
echo "Step 6: Installing bioinformatics tools..."
echo "This may take 10-15 minutes..."

# Core tools for read processing
/opt/conda/bin/conda install -y bbtools  # Includes bbmerge.sh
/opt/conda/bin/conda install -y flash
/opt/conda/bin/conda install -y pear

# Sequence analysis tools
/opt/conda/bin/conda install -y samtools
/opt/conda/bin/conda install -y bcftools
/opt/conda/bin/conda install -y bwa
/opt/conda/bin/conda install -y bowtie2

# Quality control
/opt/conda/bin/conda install -y fastqc

# Python bioinformatics libraries
/opt/conda/bin/conda install -y biopython
/opt/conda/bin/conda install -y pandas
/opt/conda/bin/conda install -y numpy

# Create symlinks for easy access
echo "Step 7: Creating system-wide symlinks..."
ln -sf /opt/conda/bin/bbmerge.sh /usr/local/bin/bbmerge.sh || true
ln -sf /opt/conda/bin/samtools /usr/local/bin/samtools || true
ln -sf /opt/conda/bin/bcftools /usr/local/bin/bcftools || true
ln -sf /opt/conda/bin/bwa /usr/local/bin/bwa || true
ln -sf /opt/conda/bin/bowtie2 /usr/local/bin/bowtie2 || true
ln -sf /opt/conda/bin/fastqc /usr/local/bin/fastqc || true

# Create working directory
echo "Step 8: Setting up working directory..."
mkdir -p /opt/helix-tools
chmod 777 /opt/helix-tools

# Verify installations
echo ""
echo "Step 9: Verifying installations..."
echo "=========================================="

tools=("bbmerge.sh" "samtools" "bcftools" "bwa" "bowtie2" "fastqc")
for tool in "${tools[@]}"; do
    if command -v "$tool" &> /dev/null; then
        echo "✅ $tool: $(which $tool)"
    else
        echo "❌ $tool: NOT FOUND"
    fi
done

# Check Python libraries
echo ""
echo "Python libraries:"
/opt/conda/bin/python -c "import Bio; print('✅ BioPython:', Bio.__version__)" || echo "❌ BioPython: NOT FOUND"
/opt/conda/bin/python -c "import pandas; print('✅ pandas:', pandas.__version__)" || echo "❌ pandas: NOT FOUND"
/opt/conda/bin/python -c "import numpy; print('✅ numpy:', numpy.__version__)" || echo "❌ numpy: NOT FOUND"

echo ""
echo "=========================================="
echo "✅ Setup completed!"
echo "=========================================="
echo ""
echo "Tools are installed in /opt/conda/bin"
echo "Working directory: /opt/helix-tools"
echo ""
echo "To use this instance with Helix.AI:"
echo "1. Note the instance ID"
echo "2. Set environment variable: export HELIX_EC2_INSTANCE_ID=<instance-id>"
echo "3. Set SSH key: export HELIX_EC2_KEY_FILE=~/.ssh/your-key.pem"
echo "4. Enable EC2 execution: export HELIX_USE_EC2=true"

