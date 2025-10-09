#!/bin/bash

# Setup script for RhoFold diffusion training environment
# This script creates a conda environment with all necessary dependencies

set -e

echo "Setting up RhoFold diffusion training environment..."

# Create conda environment
conda create -n rnafold python=3.10 -y

# Activate environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate rnafold

# Install PyTorch with CUDA support (for 2080Ti)
echo "Installing PyTorch with CUDA support..."
conda install pytorch torchvision torchaudio pytorch-cuda=11.8 -c pytorch -c nvidia -y

# Install core dependencies
echo "Installing core dependencies..."
conda install -c conda-forge -y \
    numpy \
    scipy \
    pandas \
    matplotlib \
    seaborn \
    tqdm \
    h5py \
    biopython \
    pyyaml \
    ml-collections \
    einops \
    wandb \
    tensorboard

# Install additional dependencies
pip install \
    openmm \
    mdtraj \
    prody \
    biotite \
    esm \
    fair-esm \
    transformers \
    accelerate \
    deepspeed

# Install AMBER tools for structure relaxation
echo "Installing AMBER tools..."
conda install -c conda-forge -y \
    ambermd \
    ambertools

# Verify installation
echo "Verifying installation..."
python -c "import torch; print(f'PyTorch version: {torch.__version__}'); print(f'CUDA available: {torch.cuda.is_available()}'); print(f'CUDA version: {torch.version.cuda}')"
python -c "import openmm; print(f'OpenMM version: {openmm.__version__}')"

echo "Environment setup complete!"
echo "To activate the environment, run: conda activate rnafold"
