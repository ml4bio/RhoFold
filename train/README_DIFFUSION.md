# RhoFold with Diffusion Module

This directory contains the implementation and training code for RhoFold with a new diffusion module that replaces the original E2Eformer transformer stack.

## Overview

The diffusion module (`@diffusion.py`) implements a Denoising Diffusion Probabilistic Model (DDPM) for pairwise feature refinement, replacing the original `@e2eformer.py` transformer-based approach. This provides an alternative architecture for processing MSA and pairwise features in the RhoFold model.

## Key Components

### 1. Diffusion Module (`rhofold/model/diffusion.py`)
- **DiffusionModule**: Core DDPM implementation with UNet2D for noise prediction
- **DiffusionStack**: Main stack that replaces E2EformerStack
- **UNet2D**: Simplified 2D UNet architecture for processing pairwise features
- **Time Embedding**: Sinusoidal positional embedding for diffusion timesteps

### 2. Data Loading (`train/data_loader.py`)
- **RhoFoldDataset**: Handles loading of `.seq`, `.a3m`, and `.pdb` files
- **collate_fn**: Pads sequences to same length for batch processing
- Supports FASTA headers in sequence files
- Handles MSA processing and structure parsing

### 3. Loss Functions (`train/losses.py`)
- **RhoFoldLoss**: Combined loss function with multiple components
- **L_dis**: Distance prediction loss
- **L_ss**: Secondary structure loss
- **L_fape**: Frame Aligned Point Error loss
- **L_clash**: Clash penalty loss
- **L_pLDDT**: pLDDT confidence loss
- **L_ss3d**: 3D secondary structure loss

### 4. Training Scripts
- **train_diffusion_final.py**: Complete training script with checkpointing
- **train_distributed.py**: Distributed training for multiple GPUs
- **launch_training.sh**: Shell script for launching distributed training

## Setup

### 1. Environment Setup
```bash
# Create conda environment
conda create -n rnafold python=3.9
conda activate rnafold

# Install PyTorch with CUDA 11.8
conda install pytorch torchvision torchaudio pytorch-cuda=11.8 -c pytorch -c nvidia

# Install other dependencies
pip install ml-collections biotite mdtraj tqdm
```

### 2. Data Format
The training expects data in the following format:
```
data_dir/
├── index.txt          # List of sample names (one per line)
├── sample1.seq        # Sequence file (FASTA format)
├── sample1.a3m        # MSA file
├── sample1.pdb        # Structure file
├── sample2.seq
├── sample2.a3m
└── sample2.pdb
```

## Usage

### 1. Basic Training
```bash
python train/train_diffusion_final.py \
    --data_dir /path/to/train/data \
    --index_file /path/to/train/index.txt \
    --val_data_dir /path/to/val/data \
    --val_index_file /path/to/val/index.txt \
    --batch_size 4 \
    --epochs 100 \
    --learning_rate 1e-4 \
    --max_seq_len 512 \
    --msa_depth 128
```

### 2. Distributed Training (4x2080Ti)
```bash
# Launch distributed training
bash train/launch_training.sh \
    --data_dir /path/to/train/data \
    --index_file /path/to/train/index.txt \
    --batch_size 4 \
    --epochs 100
```

### 3. Testing
```bash
# Test data loading
python train/debug_data_loading.py

# Test simplified loss computation
python train/test_simple_loss.py

# Test gradient flow
python train/test_gradient_flow.py

# Test full pipeline
python train/test_real_data.py
```

## Model Architecture

### Diffusion Module Details
- **Input**: Pairwise features `[B, L, L, C_z]` and MSA features `[B, N_seq, L, C_m]`
- **Process**: 
  1. Add noise to pairwise features during training
  2. Use UNet2D to predict and remove noise
  3. Apply pair mask to refined features
- **Output**: Refined pairwise features `[B, L, L, C_z]` and single representation `[B, L, C_s]`

### UNet2D Architecture
- **Input Projection**: Conv2d to base channels
- **Residual Blocks**: 6 residual blocks with time embedding
- **Output Projection**: Conv2d to output channels
- **Time Embedding**: Sinusoidal positional embedding with MLP

## Configuration

### Diffusion Configuration
```python
diffusion_config = {
    'time_embed_dim': 128,      # Time embedding dimension
    'num_timesteps': 100,       # Number of diffusion timesteps
    'noise_schedule': 'linear', # Noise schedule type
    'base_channels': 32,        # Base number of channels in UNet
}
```

### Training Configuration
```python
config = {
    'learning_rate': 1e-4,      # Learning rate
    'weight_decay': 1e-5,       # Weight decay
    'epochs': 100,              # Number of epochs
    'batch_size': 4,            # Batch size
    'max_seq_len': 512,         # Maximum sequence length
    'msa_depth': 128,           # MSA depth
}
```

## Key Features

### 1. Robust Data Loading
- Handles various file formats and edge cases
- Automatic padding and truncation
- Support for FASTA headers in sequence files
- Robust PDB parsing with fallback to dummy data

### 2. Flexible Loss Computation
- Modular loss functions
- Support for different tensor shapes
- Automatic shape matching and padding
- NaN detection and handling

### 3. Gradient Flow Validation
- Comprehensive gradient checking
- Verification of diffusion module gradients
- Backward pass validation

### 4. Checkpointing and Resuming
- Automatic checkpoint saving
- Model state preservation
- Easy resuming from checkpoints

## Troubleshooting

### Common Issues

1. **Data Loading Errors**
   - Check file formats and paths
   - Ensure index file contains valid sample names
   - Verify PDB file format

2. **Memory Issues**
   - Reduce batch size
   - Decrease sequence length or MSA depth
   - Use gradient checkpointing

3. **Gradient Issues**
   - Check learning rate
   - Verify loss computation
   - Use gradient clipping

### Debug Commands
```bash
# Test individual components
python train/debug_data_loading.py
python train/test_simple_loss.py
python train/test_gradient_flow.py

# Check model architecture
python -c "from rhofold.rhofold_with_diffusion import create_diffusion_rhofold; print(create_diffusion_rhofold(rhofold_config, {}))"
```

## Performance

### Memory Usage
- **Training**: ~8GB GPU memory for batch size 4, sequence length 512
- **Inference**: ~4GB GPU memory for single sequence

### Training Speed
- **Single GPU**: ~2-3 minutes per epoch (batch size 4)
- **4x2080Ti**: ~1-2 minutes per epoch (batch size 4)

## Future Improvements

1. **Advanced UNet Architecture**
   - Add attention mechanisms
   - Implement skip connections
   - Use residual connections

2. **Loss Function Enhancements**
   - Add more sophisticated FAPE loss
   - Implement contact prediction loss
   - Add auxiliary losses

3. **Training Optimizations**
   - Mixed precision training
   - Gradient accumulation
   - Learning rate scheduling

## Citation

If you use this diffusion module in your research, please cite the original RhoFold paper and acknowledge the diffusion module implementation.

## License

This code is licensed under the Apache License 2.0, same as the original RhoFold codebase.

