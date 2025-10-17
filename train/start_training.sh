#!/bin/bash

# Training script for RhoFold with diffusion module
# Usage: bash start_training.sh [index_file] [batch_size] [epochs]

# Default parameters
INDEX_FILE=${1:-"fold-0_train_ids"}
BATCH_SIZE=${2:-4}
EPOCHS=${3:-100}
DATA_DIR="/home/rcy/code/RhoFold/RNA3D_DATA"
SAVE_DIR="./checkpoints_fold-0"

echo "Starting RhoFold Diffusion Training"
echo "=================================="
echo "Data directory: $DATA_DIR"
echo "Index file: $INDEX_FILE"
echo "Batch size: $BATCH_SIZE"
echo "Epochs: $EPOCHS"
echo "Save directory: $SAVE_DIR"
echo ""

# Check if data directory exists
if [ ! -d "$DATA_DIR" ]; then
    echo "❌ Data directory not found: $DATA_DIR"
    exit 1
fi

# Check if index file exists
INDEX_PATH="$DATA_DIR/list/$INDEX_FILE"
if [ ! -f "$INDEX_PATH" ]; then
    echo "❌ Index file not found: $INDEX_PATH"
    echo "Available index files:"
    ls -la "$DATA_DIR/list/"
    exit 1
fi

# Check if subdirectories exist
for subdir in "seq" "rMSA" "pdb"; do
    if [ ! -d "$DATA_DIR/$subdir" ]; then
        echo "❌ Subdirectory not found: $DATA_DIR/$subdir"
        exit 1
    fi
done

echo "✓ All data directories found"
echo ""

# Activate conda environment
echo "Activating conda environment..."
source ~/miniconda3/etc/profile.d/conda.sh
conda activate rnafold

if [ $? -ne 0 ]; then
    echo "❌ Failed to activate conda environment 'rnafold'"
    echo "Please create the environment first:"
    echo "  conda create -n rnafold python=3.9"
    echo "  conda activate rnafold"
    echo "  conda install pytorch torchvision torchaudio pytorch-cuda=11.8 -c pytorch -c nvidia"
    echo "  pip install ml-collections biotite mdtraj tqdm"
    exit 1
fi

echo "✓ Conda environment activated"
echo ""

# Run training
echo "Starting training..."
python train/train_diffusion_final.py \
    --data_dir "$DATA_DIR" \
    --index_file "$INDEX_FILE" \
    --batch_size "$BATCH_SIZE" \
    --epochs "$EPOCHS" \
    --learning_rate 1e-4 \
    --weight_decay 1e-5 \
    --max_seq_len 512 \
    --msa_depth 128 \
    --save_every 10 \
    --save_dir "$SAVE_DIR"

echo ""
echo "Training completed!"

