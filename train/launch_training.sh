#!/bin/bash

# Launch script for RhoFold diffusion training
# Supports both single GPU and multi-GPU training

set -e

# Default values
DATA_DIR=""
INDEX_FILE=""
EXPERIMENT_NAME="diffusion_test"
BATCH_SIZE=4
EPOCHS=100
LEARNING_RATE=1e-4
MAX_SEQ_LEN=512
MSA_DEPTH=128
WORLD_SIZE=4
USE_WANDB=false
CHECKPOINT_DIR="./checkpoints"

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --data_dir)
            DATA_DIR="$2"
            shift 2
            ;;
        --index_file)
            INDEX_FILE="$2"
            shift 2
            ;;
        --experiment_name)
            EXPERIMENT_NAME="$2"
            shift 2
            ;;
        --batch_size)
            BATCH_SIZE="$2"
            shift 2
            ;;
        --epochs)
            EPOCHS="$2"
            shift 2
            ;;
        --learning_rate)
            LEARNING_RATE="$2"
            shift 2
            ;;
        --max_seq_len)
            MAX_SEQ_LEN="$2"
            shift 2
            ;;
        --msa_depth)
            MSA_DEPTH="$2"
            shift 2
            ;;
        --world_size)
            WORLD_SIZE="$2"
            shift 2
            ;;
        --use_wandb)
            USE_WANDB=true
            shift
            ;;
        --checkpoint_dir)
            CHECKPOINT_DIR="$2"
            shift 2
            ;;
        --help)
            echo "Usage: $0 --data_dir DATA_DIR --index_file INDEX_FILE [OPTIONS]"
            echo ""
            echo "Required arguments:"
            echo "  --data_dir DATA_DIR        Directory containing .seq, .a3m, .pdb files"
            echo "  --index_file INDEX_FILE    Text file with list of sample names"
            echo ""
            echo "Optional arguments:"
            echo "  --experiment_name NAME     Experiment name for logging (default: diffusion_test)"
            echo "  --batch_size SIZE          Batch size per GPU (default: 4)"
            echo "  --epochs EPOCHS            Number of epochs (default: 100)"
            echo "  --learning_rate LR         Learning rate (default: 1e-4)"
            echo "  --max_seq_len LEN          Maximum sequence length (default: 512)"
            echo "  --msa_depth DEPTH          Maximum MSA depth (default: 128)"
            echo "  --world_size SIZE          Number of GPUs to use (default: 4)"
            echo "  --use_wandb                Use Weights & Biases for logging"
            echo "  --checkpoint_dir DIR       Directory to save checkpoints (default: ./checkpoints)"
            echo "  --help                     Show this help message"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

# Check required arguments
if [[ -z "$DATA_DIR" || -z "$INDEX_FILE" ]]; then
    echo "Error: --data_dir and --index_file are required"
    echo "Use --help for usage information"
    exit 1
fi

# Check if data directory exists
if [[ ! -d "$DATA_DIR" ]]; then
    echo "Error: Data directory $DATA_DIR does not exist"
    exit 1
fi

# Check if index file exists
if [[ ! -f "$INDEX_FILE" ]]; then
    echo "Error: Index file $INDEX_FILE does not exist"
    exit 1
fi

# Check if CUDA is available
if ! command -v nvidia-smi &> /dev/null; then
    echo "Error: nvidia-smi not found. Please ensure CUDA is installed"
    exit 1
fi

# Check number of available GPUs
AVAILABLE_GPUS=$(nvidia-smi --list-gpus | wc -l)
if [[ $AVAILABLE_GPUS -lt $WORLD_SIZE ]]; then
    echo "Warning: Only $AVAILABLE_GPUS GPUs available, but $WORLD_SIZE requested"
    WORLD_SIZE=$AVAILABLE_GPUS
fi

echo "Starting RhoFold diffusion training..."
echo "======================================"
echo "Data directory: $DATA_DIR"
echo "Index file: $INDEX_FILE"
echo "Experiment name: $EXPERIMENT_NAME"
echo "Batch size per GPU: $BATCH_SIZE"
echo "Total effective batch size: $((BATCH_SIZE * WORLD_SIZE))"
echo "Epochs: $EPOCHS"
echo "Learning rate: $LEARNING_RATE"
echo "Max sequence length: $MAX_SEQ_LEN"
echo "MSA depth: $MSA_DEPTH"
echo "Number of GPUs: $WORLD_SIZE"
echo "Use wandb: $USE_WANDB"
echo "Checkpoint directory: $CHECKPOINT_DIR"
echo "======================================"

# Create checkpoint directory
mkdir -p "$CHECKPOINT_DIR"

# Build command
CMD="python train_distributed.py"
CMD="$CMD --data_dir $DATA_DIR"
CMD="$CMD --index_file $INDEX_FILE"
CMD="$CMD --experiment_name $EXPERIMENT_NAME"
CMD="$CMD --batch_size $BATCH_SIZE"
CMD="$CMD --epochs $EPOCHS"
CMD="$CMD --learning_rate $LEARNING_RATE"
CMD="$CMD --max_seq_len $MAX_SEQ_LEN"
CMD="$CMD --msa_depth $MSA_DEPTH"
CMD="$CMD --world_size $WORLD_SIZE"
CMD="$CMD --checkpoint_dir $CHECKPOINT_DIR"

if [[ "$USE_WANDB" == "true" ]]; then
    CMD="$CMD --use_wandb"
fi

echo "Running command: $CMD"
echo ""

# Run training
eval $CMD

echo ""
echo "Training completed!"
