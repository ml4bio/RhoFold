"""
Debug test to understand the data flow and shapes.
"""

import os
import sys
import torch
import numpy as np
from pathlib import Path

# Add parent directory to path
sys.path.append(str(Path(__file__).parent.parent))

from rhofold.config import rhofold_config
from rhofold.rhofold_with_diffusion import create_diffusion_rhofold


def test_data_shapes():
    """Test with correct data shapes."""
    print("Testing data shapes...")
    
    try:
        # Create model with smaller config
        diffusion_config = {
            'time_embed_dim': 128,
            'num_timesteps': 100,
            'noise_schedule': 'linear',
            'base_channels': 32,
        }
        
        model = create_diffusion_rhofold(rhofold_config, diffusion_config)
        model.eval()
        
        device = torch.device('cpu')
        model.to(device)
        
        # Create dummy data with correct shapes
        batch_size = 1
        seq_len = 50
        msa_depth = 3
        c_m = 256  # MSA feature dimension
        c_z = 128  # Pairwise feature dimension
        
        # MSA tokens: [batch, msa_depth, seq_len]
        tokens = torch.randint(0, 4, (batch_size, msa_depth, seq_len), dtype=torch.long, device=device)
        
        # RNA-FM tokens: [batch, seq_len]
        rna_fm_tokens = torch.randint(0, 4, (batch_size, seq_len), dtype=torch.long, device=device)
        
        # Sequence
        sequence = ["AUCG" * (seq_len // 4)] * batch_size
        
        print(f"Input shapes:")
        print(f"  tokens: {tokens.shape}")
        print(f"  rna_fm_tokens: {rna_fm_tokens.shape}")
        print(f"  sequence: {len(sequence)} sequences of length {len(sequence[0])}")
        
        # Test MSA embedder separately
        print("\nTesting MSA embedder...")
        msa_fea, pair_fea = model.msa_embedder.forward(
            tokens=tokens,
            rna_fm_tokens=rna_fm_tokens,
            is_BKL=True
        )
        
        print(f"MSA embedder outputs:")
        print(f"  msa_fea: {msa_fea.shape}")
        print(f"  pair_fea: {pair_fea.shape}")
        
        # Test diffusion module separately
        print("\nTesting diffusion module...")
        if hasattr(model.e2eformer, 'diffusion_blocks'):
            diffusion_block = model.e2eformer.diffusion_blocks[0]
            
            # Create proper masks
            msa_mask = torch.ones(msa_fea.shape[:3], device=device)
            pair_mask = torch.ones(pair_fea.shape[:3], device=device)
            
            print(f"Diffusion inputs:")
            print(f"  z (pair_fea): {pair_fea.shape}")
            print(f"  msa_fea: {msa_fea.shape}")
            print(f"  msa_mask: {msa_mask.shape}")
            print(f"  pair_mask: {pair_mask.shape}")
            
            # Test diffusion forward
            with torch.no_grad():
                m_out, z_out, s_out = diffusion_block(
                    z=pair_fea,
                    msa_fea=msa_fea,
                    msa_mask=msa_mask,
                    pair_mask=pair_mask,
                    chunk_size=None,
                    inplace_safe=False,
                    _mask_trans=True,
                )
                
                print(f"Diffusion outputs:")
                print(f"  m_out: {m_out.shape}")
                print(f"  z_out: {z_out.shape}")
                print(f"  s_out: {s_out.shape}")
        
        return True
        
    except Exception as e:
        print(f"✗ Test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    """Run debug test."""
    print("RhoFold Diffusion Debug Test")
    print("=" * 40)
    
    success = test_data_shapes()
    
    if success:
        print("\n✓ Debug test passed!")
    else:
        print("\n✗ Debug test failed!")
    
    return success


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
