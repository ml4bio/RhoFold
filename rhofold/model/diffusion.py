# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import torch
import torch.nn as nn
import torch.nn.functional as F
import math
from typing import Tuple, Optional, Sequence
from functools import partial

from rhofold.model.primitives import Linear, LayerNorm
from rhofold.utils.tensor_utils import add


class SinusoidalPositionalEmbedding(nn.Module):
    """Sinusoidal positional embedding for time steps."""
    
    def __init__(self, dim):
        super().__init__()
        self.dim = dim
        
    def forward(self, timesteps):
        device = timesteps.device
        half_dim = self.dim // 2
        embeddings = math.log(10000) / (half_dim - 1)
        embeddings = torch.exp(torch.arange(half_dim, device=device) * -embeddings)
        embeddings = timesteps[:, None] * embeddings[None, :]
        embeddings = torch.cat((embeddings.sin(), embeddings.cos()), dim=-1)
        return embeddings


class TimeEmbedding(nn.Module):
    """Time embedding module for diffusion timesteps."""
    
    def __init__(self, time_embed_dim):
        super().__init__()
        self.time_embed_dim = time_embed_dim
        self.time_mlp = nn.Sequential(
            SinusoidalPositionalEmbedding(time_embed_dim),
            nn.Linear(time_embed_dim, time_embed_dim * 4),
            nn.SiLU(),
            nn.Linear(time_embed_dim * 4, time_embed_dim),
        )
        
    def forward(self, timesteps):
        return self.time_mlp(timesteps)


class ResidualBlock(nn.Module):
    """Residual block with time embedding for UNet."""
    
    def __init__(self, in_channels, out_channels, time_embed_dim, dropout=0.1):
        super().__init__()
        self.in_channels = in_channels
        self.out_channels = out_channels
        
        self.norm1 = nn.GroupNorm(8, in_channels)
        self.conv1 = nn.Conv2d(in_channels, out_channels, 3, padding=1)
        
        self.time_proj = nn.Linear(time_embed_dim, out_channels)
        
        self.norm2 = nn.GroupNorm(8, out_channels)
        self.conv2 = nn.Conv2d(out_channels, out_channels, 3, padding=1)
        
        self.dropout = nn.Dropout(dropout)
        
        if in_channels != out_channels:
            self.residual_conv = nn.Conv2d(in_channels, out_channels, 1)
        else:
            self.residual_conv = nn.Identity()
            
    def forward(self, x, time_emb):
        residual = self.residual_conv(x)
        
        # First conv block
        h = self.norm1(x)
        h = F.silu(h)
        h = self.conv1(h)
        
        # Add time embedding
        time_emb = self.time_proj(time_emb)[:, :, None, None]
        h = h + time_emb
        
        # Second conv block
        h = self.norm2(h)
        h = F.silu(h)
        h = self.dropout(h)
        h = self.conv2(h)
        
        return h + residual


class DownBlock(nn.Module):
    """Downsampling block for UNet."""
    
    def __init__(self, in_channels, out_channels, time_embed_dim, num_layers=2):
        super().__init__()
        self.layers = nn.ModuleList([
            ResidualBlock(in_channels if i == 0 else out_channels, 
                         out_channels, time_embed_dim)
            for i in range(num_layers)
        ])
        self.downsample = nn.Conv2d(out_channels, out_channels, 3, stride=2, padding=1)
        
    def forward(self, x, time_emb):
        skip_connections = []
        
        for layer in self.layers:
            x = layer(x, time_emb)
            skip_connections.append(x)
            
        x = self.downsample(x)
        return x, skip_connections


class UpBlock(nn.Module):
    """Upsampling block for UNet."""
    
    def __init__(self, in_channels, out_channels, time_embed_dim, num_layers=2, use_skip=True):
        super().__init__()
        self.use_skip = use_skip
        self.upsample = nn.ConvTranspose2d(in_channels, in_channels, 3, stride=2, padding=1, output_padding=1)
        
        self.layers = nn.ModuleList([
            ResidualBlock(in_channels if i == 0 else out_channels, 
                         out_channels, time_embed_dim)
            for i in range(num_layers)
        ])
        
    def forward(self, x, skip_connections, time_emb):
        x = self.upsample(x)
        
        # Concatenate with skip connections if enabled
        for i, layer in enumerate(self.layers):
            if self.use_skip and i < len(skip_connections):
                skip = skip_connections[-(i+1)]
                x = torch.cat([x, skip], dim=1)
            x = layer(x, time_emb)
            
        return x


class UNet2D(nn.Module):
    """2D UNet for processing pairwise features."""
    
    def __init__(self, in_channels, out_channels, time_embed_dim=256, 
                 base_channels=64, num_layers_per_block=2):
        super().__init__()
        self.in_channels = in_channels
        self.out_channels = out_channels
        self.time_embed_dim = time_embed_dim
        
        # Time embedding
        self.time_embed = TimeEmbedding(time_embed_dim)
        
        # Initial projection
        self.input_proj = nn.Conv2d(in_channels, base_channels, 3, padding=1)
        
        # Simplified architecture - just residual blocks without downsampling/upsampling
        self.blocks = nn.ModuleList([
            ResidualBlock(base_channels, base_channels, time_embed_dim)
            for _ in range(6)  # 6 residual blocks
        ])
        
        # Output projection
        self.output_proj = nn.Conv2d(base_channels, out_channels, 3, padding=1)
        
    def forward(self, x, timesteps):
        # Time embedding
        time_emb = self.time_embed(timesteps)
        
        # Input projection
        x = self.input_proj(x)
        
        # Process through residual blocks
        for block in self.blocks:
            x = block(x, time_emb)
            
        # Output projection
        x = self.output_proj(x)
        
        return x


class DiffusionModule(nn.Module):
    """DDPM-based diffusion module for pairwise feature refinement."""
    
    def __init__(self, c_z, c_s, time_embed_dim=256, num_timesteps=1000, 
                 noise_schedule="linear", base_channels=64):
        super().__init__()
        self.c_z = c_z
        self.c_s = c_s
        self.num_timesteps = num_timesteps
        self.time_embed_dim = time_embed_dim
        
        # UNet for denoising
        self.unet = UNet2D(
            in_channels=c_z + 1,  # +1 for noise channel
            out_channels=c_z,
            time_embed_dim=time_embed_dim,
            base_channels=base_channels
        )
        
        # Noise schedule
        if noise_schedule == "linear":
            self.register_buffer('betas', self.linear_beta_schedule(num_timesteps))
        elif noise_schedule == "cosine":
            self.register_buffer('betas', self.cosine_beta_schedule(num_timesteps))
        else:
            raise ValueError(f"Unknown noise schedule: {noise_schedule}")
            
        self.register_buffer('alphas', 1.0 - self.betas)
        self.register_buffer('alphas_cumprod', torch.cumprod(self.alphas, dim=0))
        self.register_buffer('alphas_cumprod_prev', F.pad(self.alphas_cumprod[:-1], (1, 0), value=1.0))
        
        # Calculations for diffusion q(x_t | x_{t-1}) and others
        self.register_buffer('sqrt_alphas_cumprod', torch.sqrt(self.alphas_cumprod))
        self.register_buffer('sqrt_one_minus_alphas_cumprod', torch.sqrt(1.0 - self.alphas_cumprod))
        
    def linear_beta_schedule(self, timesteps, beta_start=0.0001, beta_end=0.02):
        return torch.linspace(beta_start, beta_end, timesteps)
        
    def cosine_beta_schedule(self, timesteps, s=0.008):
        steps = timesteps + 1
        x = torch.linspace(0, timesteps, steps)
        alphas_cumprod = torch.cos(((x / timesteps) + s) / (1 + s) * torch.pi * 0.5) ** 2
        alphas_cumprod = alphas_cumprod / alphas_cumprod[0]
        betas = 1 - (alphas_cumprod[1:] / alphas_cumprod[:-1])
        return torch.clip(betas, 0.0001, 0.9999)
    
    def extract(self, a, t, x_shape):
        """Extract coefficients for given timesteps."""
        batch_size = t.shape[0]
        # Ensure both tensors are on the same device
        a = a.to(t.device)
        out = a.gather(-1, t)
        return out.reshape(batch_size, *((1,) * (len(x_shape) - 1)))
    
    def q_sample(self, x_start, t, noise=None):
        """Sample from q(x_t | x_0)."""
        if noise is None:
            noise = torch.randn_like(x_start)
            
        sqrt_alphas_cumprod_t = self.extract(self.sqrt_alphas_cumprod, t, x_start.shape)
        sqrt_one_minus_alphas_cumprod_t = self.extract(self.sqrt_one_minus_alphas_cumprod, t, x_start.shape)
        
        return sqrt_alphas_cumprod_t * x_start + sqrt_one_minus_alphas_cumprod_t * noise
    
    def predict_noise(self, x_noisy, t):
        """Predict noise from noisy input."""
        # Convert from [B, H, W, C] to [B, C, H, W] for UNet
        x_noisy = x_noisy.permute(0, 3, 1, 2)
        
        # Add noise channel
        noise_channel = torch.randn_like(x_noisy[:, :1])
        x_input = torch.cat([x_noisy, noise_channel], dim=1)
        
        # Predict noise
        predicted_noise = self.unet(x_input, t)
        
        # Convert back to [B, H, W, C]
        predicted_noise = predicted_noise.permute(0, 2, 3, 1)
        
        return predicted_noise
    
    def p_sample(self, x, t, t_index):
        """Sample from p(x_{t-1} | x_t)."""
        betas_t = self.extract(self.betas, t, x.shape)
        sqrt_one_minus_alphas_cumprod_t = self.extract(self.sqrt_one_minus_alphas_cumprod, t, x.shape)
        sqrt_recip_alphas_t = self.extract(torch.sqrt(1.0 / self.alphas), t, x.shape)
        
        # Predict noise
        predicted_noise = self.predict_noise(x, t)
        
        # Compute mean
        model_mean = sqrt_recip_alphas_t * (x - betas_t * predicted_noise / sqrt_one_minus_alphas_cumprod_t)
        
        if t_index == 0:
            return model_mean
        else:
            posterior_variance_t = self.extract(self.betas, t, x.shape)
            noise = torch.randn_like(x)
            return model_mean + torch.sqrt(posterior_variance_t) * noise
    
    def p_sample_loop(self, shape, device):
        """Sample from the diffusion process."""
        b = shape[0]
        # Start from pure noise
        img = torch.randn(shape, device=device)
        
        for i in reversed(range(self.num_timesteps)):
            t = torch.full((b,), i, device=device, dtype=torch.long)
            img = self.p_sample(img, t, i)
            
        return img
    
    def forward(self, z, msa_fea, msa_mask, pair_mask, chunk_size=None, 
                inplace_safe=False, _mask_trans=True, _attn_chunk_size=None, 
                _offload_inference=False, **kwargs):
        """
        Forward pass for diffusion module.
        
        Args:
            z: [*, N_res, N_res, C_z] pair embedding
            msa_fea: [*, N_seq, N_res, C_m] MSA embedding (used for conditioning)
            msa_mask: [*, N_seq, N_res] MSA mask
            pair_mask: [*, N_res, N_res] pair mask
            chunk_size: Chunk size for processing
            inplace_safe: Whether to use in-place operations
            _mask_trans: Whether to apply mask in transitions
            _attn_chunk_size: Chunk size for attention
            _offload_inference: Whether to offload for inference
            **kwargs: Additional arguments
            
        Returns:
            m: [*, N_seq, N_res, C_m] MSA embedding (unchanged)
            z: [*, N_res, N_res, C_z] refined pair embedding
            s: [*, N_res, C_s] single embedding (from first row of MSA)
        """
        device = z.device
        batch_size = z.shape[0]
        
        # During training, add noise to pairwise features
        if self.training:
            # Sample random timesteps (ensure they're within valid range)
            t = torch.randint(0, self.num_timesteps, (batch_size,), device=device)
            # Clamp to ensure valid range
            t = torch.clamp(t, 0, self.num_timesteps - 1)
            
            # Add noise to pairwise features
            noise = torch.randn_like(z)
            z_noisy = self.q_sample(z, t, noise)
            
            # Predict and remove noise
            predicted_noise = self.predict_noise(z_noisy, t)
            z_refined = z_noisy - predicted_noise
        else:
            # During inference, use full denoising process
            z_refined = self.p_sample_loop(z.shape, device)
        
        # Apply pair mask
        if pair_mask is not None:
            z_refined = z_refined * pair_mask.unsqueeze(-1)
        
        # Generate single representation from MSA features
        s = msa_fea[..., 0, :, :]  # Take first row (query sequence)
        
        # Project to single representation dimension
        if s.shape[-1] != self.c_s:
            s = F.linear(s, torch.randn(self.c_s, s.shape[-1], device=device))
        
        return msa_fea, z_refined, s


class DiffusionStack(nn.Module):
    """
    Main diffusion stack that replaces E2EformerStack.
    
    This module processes MSA and pairwise features using diffusion-based refinement.
    """
    
    def __init__(
        self,
        c_m: int,
        c_z: int,
        c_hidden_msa_att: int,
        c_hidden_opm: int,
        c_hidden_mul: int,
        c_hidden_pair_att: int,
        c_s: int,
        no_heads_msa: int,
        no_heads_pair: int,
        no_blocks: int,
        transition_n: int,
        blocks_per_ckpt: int,
        tune_chunk_size: bool = False,
        diffusion_config: Optional[dict] = None,
        **kwargs,
    ):
        """
        Args:
            c_m: MSA channel dimension
            c_z: Pair channel dimension
            c_hidden_msa_att: Hidden dimension in MSA attention (unused in diffusion)
            c_hidden_opm: Hidden dimension in outer product mean (unused in diffusion)
            c_hidden_mul: Hidden dimension in multiplicative updates (unused in diffusion)
            c_hidden_pair_att: Hidden dimension in triangular attention (unused in diffusion)
            c_s: Channel dimension of the output "single" embedding
            no_heads_msa: Number of heads used for MSA attention (unused in diffusion)
            no_heads_pair: Number of heads used for pair attention (unused in diffusion)
            no_blocks: Number of diffusion blocks in the stack
            transition_n: Factor for transition layers (unused in diffusion)
            blocks_per_ckpt: Number of blocks per checkpoint (unused in diffusion)
            tune_chunk_size: Whether to tune chunk size (unused in diffusion)
            diffusion_config: Configuration for diffusion module
            **kwargs: Additional arguments
        """
        super(DiffusionStack, self).__init__()
        
        self.c_m = c_m
        self.c_z = c_z
        self.c_s = c_s
        self.no_blocks = no_blocks
        
        # Default diffusion configuration
        if diffusion_config is None:
            diffusion_config = {
                'time_embed_dim': 256,
                'num_timesteps': 1000,
                'noise_schedule': 'cosine',
                'base_channels': 64,
            }
        
        # Create diffusion modules
        self.diffusion_blocks = nn.ModuleList()
        for _ in range(no_blocks):
            diffusion_module = DiffusionModule(
                c_z=c_z,
                c_s=c_s,
                **diffusion_config
            )
            self.diffusion_blocks.append(diffusion_module)
        
        # Final projection for single representation
        self.single_proj = Linear(c_m, c_s)
        
    def forward(
        self,
        m: torch.Tensor,
        z: torch.Tensor,
        msa_mask: torch.Tensor,
        pair_mask: torch.Tensor,
        chunk_size: int,
        inplace_safe: bool = False,
        _mask_trans: bool = True,
    ) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        """
        Forward pass through diffusion stack.
        
        Args:
            m: [*, N_seq, N_res, C_m] MSA embedding
            z: [*, N_res, N_res, C_z] pair embedding
            msa_mask: [*, N_seq, N_res] MSA mask
            pair_mask: [*, N_res, N_res] pair mask
            chunk_size: Inference-time subbatch size
            inplace_safe: Whether to use in-place operations
            _mask_trans: Whether to apply mask in transitions
            
        Returns:
            m: [*, N_seq, N_res, C_m] MSA embedding
            z: [*, N_res, N_res, C_z] refined pair embedding
            s: [*, N_res, C_s] single embedding
        """
        
        # Process through diffusion blocks
        for diffusion_block in self.diffusion_blocks:
            m, z, s = diffusion_block(
                z=z,
                msa_fea=m,
                msa_mask=msa_mask,
                pair_mask=pair_mask,
                chunk_size=chunk_size,
                inplace_safe=inplace_safe,
                _mask_trans=_mask_trans,
            )
        
        # Final single representation projection
        s = self.single_proj(m[..., 0, :, :])
        
        return m, z, s
