# RhoFold Diffusion Training

This directory contains the training pipeline for RhoFold with diffusion module integration.

## Overview

The training system includes:
- **Diffusion Module**: Replaces the E2Eformer with a diffusion-based architecture
- **Data Loading**: Handles .seq, .a3m, and .pdb files
- **Loss Functions**: Comprehensive loss computation including FAPE, distance, secondary structure, etc.
- **Distributed Training**: Multi-GPU support for 4x2080Ti setup
- **Validation**: Comprehensive testing and validation scripts

## Quick Start

### 1. Environment Setup

```bash
# Create and activate conda environment
./setup_env.sh
conda activate rnafold
```

### 2. Run Validation Tests

```bash
# Run all validation tests
./quick_test.sh

# Or run individual tests
python train/validate_setup.py
python train/test_forward.py
```

### 3. Prepare Training Data

Your training data should be organized as follows:
```
data_dir/
├── sample1.seq      # RNA sequence
├── sample1.a3m      # MSA alignment
├── sample1.pdb      # 3D structure
├── sample2.seq
├── sample2.a3m
├── sample2.pdb
└── ...
```

Create an index file listing all sample names:
```
sample1
sample2
sample3
...
```

### 4. Start Training

#### Single GPU Training
```bash
python train/train_diffusion.py \
    --data_dir /path/to/data \
    --index_file /path/to/index.txt \
    --batch_size 4 \
    --epochs 100 \
    --learning_rate 1e-4
```

#### Multi-GPU Training (4x2080Ti)
```bash
./train/launch_training.sh \
    --data_dir /path/to/data \
    --index_file /path/to/index.txt \
    --batch_size 4 \
    --epochs 100 \
    --world_size 4 \
    --use_wandb
```

## File Structure

```
train/
├── data_loader.py           # Data loading pipeline
├── losses.py               # Loss function definitions
├── train_diffusion.py      # Single GPU training
├── train_distributed.py    # Multi-GPU training
├── test_forward.py         # Forward pass validation
├── validate_setup.py       # Comprehensive validation
├── launch_training.sh      # Training launch script
└── README.md              # This file
```

## Key Components

### Diffusion Module

The diffusion module (`rhofold/model/diffusion.py`) includes:
- **SinusoidalPositionalEmbedding**: Time step encoding
- **TimeEmbedding**: Time embedding MLP
- **ResidualBlock**: UNet residual blocks with time conditioning
- **UNet2D**: 2D UNet for pairwise feature processing
- **DiffusionModule**: DDPM-based diffusion process
- **DiffusionStack**: Main stack replacing E2EformerStack

### Data Loading

The data loader (`data_loader.py`) handles:
- **Sequence Loading**: From .seq files
- **MSA Loading**: From .a3m files
- **Structure Loading**: From .pdb files
- **Tokenization**: Using RNA alphabet
- **Batching**: With proper padding and masking

### Loss Functions

The loss system (`losses.py`) includes:
- **MLM Loss**: Masked language modeling
- **Distance Loss**: Pairwise distance prediction
- **Secondary Structure Loss**: Base pairing prediction
- **FAPE Loss**: Frame-aligned point error
- **Clash Loss**: Steric clash penalty
- **pLDDT Loss**: Confidence prediction

## Training Configuration

### Model Parameters

```python
diffusion_config = {
    'time_embed_dim': 256,      # Time embedding dimension
    'num_timesteps': 1000,      # Number of diffusion steps
    'noise_schedule': 'cosine', # Noise schedule type
    'base_channels': 64,        # UNet base channels
}
```

### Training Parameters

```python
training_config = {
    'batch_size': 4,            # Per GPU batch size
    'learning_rate': 1e-4,      # Learning rate
    'weight_decay': 1e-5,       # Weight decay
    'epochs': 100,              # Number of epochs
    'max_seq_len': 512,         # Maximum sequence length
    'msa_depth': 128,           # Maximum MSA depth
}
```

## Validation Tests

### Forward Pass Test
- Tests model forward pass with dummy data
- Validates output shapes and types
- Checks for errors and exceptions

### Loss Computation Test
- Tests all loss function components
- Validates loss values and gradients
- Checks for NaN values

### Gradient Flow Test
- Tests gradient flow through diffusion module
- Validates parameter updates
- Ensures gradients reach diffusion layers

### Memory Usage Test
- Tests memory consumption
- Validates batch size limits
- Checks for memory leaks

## Distributed Training

### Multi-GPU Setup

The distributed training supports:
- **Data Parallelism**: Across multiple GPUs
- **Gradient Synchronization**: Using NCCL backend
- **Checkpoint Saving**: On rank 0 process
- **Logging**: Centralized on rank 0

### Launch Script

```bash
./train/launch_training.sh \
    --data_dir /path/to/data \
    --index_file /path/to/index.txt \
    --world_size 4 \
    --batch_size 4 \
    --epochs 100 \
    --use_wandb
```

## Monitoring

### Weights & Biases Integration

```bash
# Enable wandb logging
--use_wandb

# Set experiment name
--experiment_name "my_experiment"
```

### Checkpoint Saving

Checkpoints are saved in the specified directory:
- `best_model.pth`: Best model based on validation loss
- `checkpoint_epoch_N.pth`: Regular checkpoints

## Troubleshooting

### Common Issues

1. **CUDA Out of Memory**
   - Reduce batch size
   - Reduce sequence length
   - Use gradient checkpointing

2. **Data Loading Errors**
   - Check file formats (.seq, .a3m, .pdb)
   - Verify index file format
   - Check file permissions

3. **Gradient Issues**
   - Check learning rate
   - Verify loss computation
   - Check for NaN values

### Debug Mode

Run with debug information:
```bash
python train/validate_setup.py --debug
```

## Performance Tips

1. **Memory Optimization**
   - Use appropriate batch sizes
   - Enable gradient checkpointing
   - Use mixed precision training

2. **Speed Optimization**
   - Use multiple data loader workers
   - Enable pin memory
   - Use distributed training

3. **Convergence**
   - Monitor loss curves
   - Adjust learning rate schedule
   - Use early stopping

## Example Usage

```bash
# 1. Setup environment
./setup_env.sh
conda activate rnafold

# 2. Run validation
./quick_test.sh

# 3. Prepare data (create your data directory and index file)

# 4. Start training
./train/launch_training.sh \
    --data_dir ./data \
    --index_file ./data/index.txt \
    --experiment_name "diffusion_v1" \
    --batch_size 4 \
    --epochs 100 \
    --use_wandb
```

## Support

For issues and questions:
1. Check the validation tests first
2. Review the error messages
3. Check the troubleshooting section
4. Verify your data format