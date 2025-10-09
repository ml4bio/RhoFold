#!/bin/bash

# Quick test script for RhoFold diffusion module
# This script runs all validation tests

set -e

echo "RhoFold Diffusion Module Quick Test"
echo "===================================="

# Check if conda environment exists
if ! conda info --envs | grep -q "rnafold"; then
    echo "Creating conda environment..."
    ./setup_env.sh
else
    echo "Conda environment 'rnafold' already exists"
fi

# Activate environment
echo "Activating conda environment..."
source $(conda info --base)/etc/profile.d/conda.sh
conda activate rnafold

# Run validation tests
echo "Running validation tests..."
python train/validate_setup.py

echo "Quick test completed!"
