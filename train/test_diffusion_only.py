"""
Test only the diffusion module functionality without structure module.
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


def test_diffusion_forward():
    """Test diffusion module forward pass."""
    print("Testing diffusion module forward pass...")
    
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
        
        # MSA tokens: [batch, msa_depth, seq_len]
        tokens = torch.randint(0, 4, (batch_size, msa_depth, seq_len), dtype=torch.long, device=device)
        
        # RNA-FM tokens: [batch, seq_len]
        rna_fm_tokens = torch.randint(0, 4, (batch_size, seq_len), dtype=torch.long, device=device)
        
        # Sequence
        sequence = ["AUCG" * (seq_len // 4)] * batch_size
        
        print(f"Input shapes:")
        print(f"  tokens: {tokens.shape}")
        print(f"  rna_fm_tokens: {rna_fm_tokens.shape}")
        
        # Test MSA embedder
        with torch.no_grad():
            msa_fea, pair_fea = model.msa_embedder.forward(
                tokens=tokens,
                rna_fm_tokens=rna_fm_tokens,
                is_BKL=True
            )
            
            print(f"MSA embedder outputs:")
            print(f"  msa_fea: {msa_fea.shape}")
            print(f"  pair_fea: {pair_fea.shape}")
            
            # Test diffusion stack
            msa_fea_out, pair_fea_out, single_fea_out = model.e2eformer(
                m=msa_fea,
                z=pair_fea,
                msa_mask=torch.ones(msa_fea.shape[:3]).to(device),
                pair_mask=torch.ones(pair_fea.shape[:3]).to(device),
                chunk_size=None,
            )
            
            print(f"Diffusion stack outputs:")
            print(f"  msa_fea_out: {msa_fea_out.shape}")
            print(f"  pair_fea_out: {pair_fea_out.shape}")
            print(f"  single_fea_out: {single_fea_out.shape}")
            
            # Check if diffusion output has correct shape
            if pair_fea_out.shape == pair_fea.shape:
                print("✓ Diffusion output shape is correct")
            else:
                print(f"✗ Diffusion output shape is incorrect: {pair_fea_out.shape} vs {pair_fea.shape}")
                return False
            
            return True
            
    except Exception as e:
        print(f"✗ Forward pass failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_diffusion_gradient():
    """Test diffusion module gradient flow."""
    print("Testing diffusion module gradient flow...")
    
    try:
        # Create model with smaller config
        diffusion_config = {
            'time_embed_dim': 128,
            'num_timesteps': 100,
            'noise_schedule': 'linear',
            'base_channels': 32,
        }
        
        model = create_diffusion_rhofold(rhofold_config, diffusion_config)
        model.train()
        
        device = torch.device('cpu')
        model.to(device)
        
        # Create dummy data
        batch_size = 1
        seq_len = 50
        msa_depth = 3
        
        # MSA tokens: [batch, msa_depth, seq_len]
        tokens = torch.randint(0, 4, (batch_size, msa_depth, seq_len), dtype=torch.long, device=device)
        
        # RNA-FM tokens: [batch, seq_len]
        rna_fm_tokens = torch.randint(0, 4, (batch_size, seq_len), dtype=torch.long, device=device)
        
        # Sequence
        sequence = ["AUCG" * (seq_len // 4)] * batch_size
        
        # Forward pass
        msa_fea, pair_fea = model.msa_embedder.forward(
            tokens=tokens,
            rna_fm_tokens=rna_fm_tokens,
            is_BKL=True
        )
        
        msa_fea_out, pair_fea_out, single_fea_out = model.e2eformer(
            m=msa_fea,
            z=pair_fea,
            msa_mask=torch.ones(msa_fea.shape[:3]).to(device),
            pair_mask=torch.ones(pair_fea.shape[:3]).to(device),
            chunk_size=None,
        )
        
        # Create a simple loss (sum of all outputs)
        total_loss = pair_fea_out.sum() + single_fea_out.sum()
        
        # Backward pass
        total_loss.backward()
        print("✓ Backward pass successful")
        
        # Check gradients in diffusion module
        if hasattr(model.e2eformer, 'diffusion_blocks'):
            diffusion_grads = []
            for block in model.e2eformer.diffusion_blocks:
                for param in block.parameters():
                    if param.grad is not None:
                        diffusion_grads.append(param.grad.abs().mean().item())
            
            if diffusion_grads:
                avg_grad = np.mean(diffusion_grads)
                print(f"  Average gradient magnitude in diffusion: {avg_grad:.6f}")
                
                if avg_grad > 1e-8:
                    print("✓ Gradients are flowing to diffusion module")
                    return True
                else:
                    print("✗ Gradients are too small in diffusion module")
                    return False
            else:
                print("✗ No gradients found in diffusion module")
                return False
        else:
            print("✗ No diffusion module found")
            return False
            
    except Exception as e:
        print(f"✗ Gradient flow test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    """Run diffusion-only tests."""
    print("RhoFold Diffusion Module Test")
    print("=" * 40)
    
    tests = [
        ("Diffusion Forward Pass", test_diffusion_forward),
        ("Diffusion Gradient Flow", test_diffusion_gradient),
    ]
    
    results = []
    
    for test_name, test_func in tests:
        print(f"\n{test_name}:")
        print("-" * 20)
        try:
            success = test_func()
            results.append((test_name, success))
        except Exception as e:
            print(f"✗ {test_name} failed with exception: {e}")
            results.append((test_name, False))
    
    # Summary
    print("\n" + "=" * 40)
    print("TEST SUMMARY:")
    print("=" * 40)
    
    passed = 0
    total = len(results)
    
    for test_name, success in results:
        status = "✓ PASS" if success else "✗ FAIL"
        print(f"{test_name:20} {status}")
        if success:
            passed += 1
    
    print("-" * 40)
    print(f"Passed: {passed}/{total}")
    
    if passed == total:
        print("\n🎉 All diffusion tests passed!")
        print("\nThe diffusion module is working correctly!")
        print("Key achievements:")
        print("✓ Model creation successful")
        print("✓ Forward pass works with correct shapes")
        print("✓ Gradients flow properly through diffusion module")
        print("✓ Output shapes are correct (B, L, L, C)")
        print("\nNext steps:")
        print("1. Fix data loading issues (PDB parsing)")
        print("2. Test with real data")
        print("3. Run full training")
    else:
        print(f"\n❌ {total - passed} tests failed.")
    
    return passed == total


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
