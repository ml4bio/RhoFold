# RhoFold Diffusion Module

This document describes the diffusion-based module designed to replace the E2Eformer transformer in RhoFold for RNA structure prediction.

## Overview

The diffusion module implements a DDPM (Denoising Diffusion Probabilistic Models) architecture to refine pairwise features in RNA structure prediction. It replaces the traditional transformer-based E2Eformer with a diffusion process that iteratively denoises pairwise representations.

## Architecture

### Core Components

1. **UNet2D**: A 2D U-Net architecture for processing pairwise feature maps
2. **Time Embedding**: Sinusoidal positional embeddings for diffusion timesteps
3. **Residual Blocks**: Time-conditioned residual blocks with group normalization
4. **Diffusion Process**: Forward and reverse diffusion processes for feature refinement

### Key Features

- **Multi-scale Processing**: Downsampling and upsampling paths with skip connections
- **Time Conditioning**: Each residual block is conditioned on diffusion timesteps
- **Noise Schedule**: Configurable linear or cosine noise schedules
- **Memory Efficient**: Optimized for processing large pairwise feature maps

## Module Structure

```
DiffusionModule
├── TimeEmbedding
│   ├── SinusoidalPositionalEmbedding
│   └── MLP layers
├── UNet2D
│   ├── Input projection
│   ├── DownBlocks (3 levels)
│   ├── Bottleneck
│   ├── UpBlocks (3 levels)
│   └── Output projection
└── Diffusion Process
    ├── Forward process (q_sample)
    ├── Reverse process (p_sample)
    └── Noise prediction
```

## Usage

### Basic Usage

```python
from rhofold.model.diffusion import DiffusionStack
from rhofold.config import rhofold_config

# Create diffusion configuration
diffusion_config = {
    'time_embed_dim': 256,
    'num_timesteps': 1000,
    'noise_schedule': 'cosine',
    'base_channels': 64,
}

# Create diffusion stack
diffusion_stack = DiffusionStack(
    c_m=256,  # MSA channel dimension
    c_z=128,  # Pair channel dimension
    c_s=384,  # Single representation dimension
    no_blocks=12,  # Number of diffusion blocks
    diffusion_config=diffusion_config,
)

# Forward pass
msa_fea, pair_fea, single_fea = diffusion_stack(
    m=msa_embedding,
    z=pair_embedding,
    msa_mask=msa_mask,
    pair_mask=pair_mask,
    chunk_size=None,
)
```

### Integration with RhoFold

```python
from rhofold.rhofold_with_diffusion import create_diffusion_rhofold

# Create RhoFold with diffusion module
model = create_diffusion_rhofold(
    config=rhofold_config,
    diffusion_config={
        'time_embed_dim': 256,
        'num_timesteps': 1000,
        'noise_schedule': 'cosine',
        'base_channels': 64,
    }
)
```

## Configuration Parameters

### Diffusion Configuration

- `time_embed_dim`: Dimension of time embeddings (default: 256)
- `num_timesteps`: Number of diffusion timesteps (default: 1000)
- `noise_schedule`: Noise schedule type ('linear' or 'cosine')
- `base_channels`: Base number of channels in UNet (default: 64)

### UNet Architecture

- `in_channels`: Input channels (c_z + 1 for noise channel)
- `out_channels`: Output channels (c_z)
- `num_layers_per_block`: Number of residual layers per block (default: 2)

## Training Process

### Forward Process (Training)

1. Sample random timesteps t ~ Uniform(0, T)
2. Add noise to pairwise features: x_t = sqrt(α_t) * x_0 + sqrt(1-α_t) * ε
3. Predict noise: ε_θ = UNet(x_t, t)
4. Compute loss: L = MSE(ε, ε_θ)

### Reverse Process (Inference)

1. Start with pure noise: x_T ~ N(0, I)
2. Iteratively denoise: x_{t-1} = p_θ(x_{t-1} | x_t)
3. Apply mask and return refined features

## Advantages over Transformer

1. **Better Global Context**: Diffusion process naturally captures long-range dependencies
2. **Iterative Refinement**: Multiple denoising steps allow for gradual feature improvement
3. **Noise Robustness**: Training with noise makes the model more robust
4. **Scalability**: UNet architecture scales better to large sequences
5. **Interpretability**: Each timestep represents a different level of detail

## Implementation Details

### Memory Optimization

- Chunked processing for large sequences
- Gradient checkpointing for memory efficiency
- In-place operations where possible

### Training Stability

- Group normalization instead of batch normalization
- SiLU activation functions
- Proper weight initialization

### Inference Speed

- Configurable number of timesteps
- Optional fast sampling methods
- Batch processing support

## Comparison with AlphaFold3

This implementation draws inspiration from AlphaFold3's diffusion module but is specifically designed for RNA structure prediction:

- **2D Processing**: Focuses on pairwise features rather than 3D coordinates
- **RNA-specific**: Optimized for RNA sequence and structure patterns
- **MSA Integration**: Maintains compatibility with MSA embeddings
- **Efficient Architecture**: Streamlined for RNA-specific requirements

## Future Improvements

1. **Conditional Diffusion**: Incorporate MSA information as conditioning
2. **Adaptive Timesteps**: Dynamic timestep selection based on sequence length
3. **Multi-scale Features**: Process features at multiple resolutions
4. **Attention Mechanisms**: Integrate attention within UNet blocks
5. **Progressive Training**: Start with fewer timesteps and increase during training

## References

- Ho, J., Jain, A., & Abbeel, P. (2020). Denoising diffusion probabilistic models.
- Rombach, R., et al. (2022). High-resolution image synthesis with latent diffusion models.
- AlphaFold3: Protein structure prediction with diffusion models.
