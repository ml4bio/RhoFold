"""
Minimal test script to validate the basic functionality without PDB parsing.
"""

import os
import sys
import torch
import numpy as np
from pathlib import Path
import tempfile

# Add parent directory to path
sys.path.append(str(Path(__file__).parent.parent))

from rhofold.config import rhofold_config
from rhofold.rhofold_with_diffusion import create_diffusion_rhofold


def test_model_creation():
    """Test model creation."""
    print("Testing model creation...")
    
    try:
        # Create model with smaller config
        diffusion_config = {
            'time_embed_dim': 128,
            'num_timesteps': 100,
            'noise_schedule': 'linear',
            'base_channels': 32,
        }
        
        model = create_diffusion_rhofold(rhofold_config, diffusion_config)
        total_params = sum(p.numel() for p in model.parameters())
        trainable_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
        
        print(f"✓ Model created: {total_params:,} total parameters, {trainable_params:,} trainable")
        
        # Check if diffusion module is present
        if hasattr(model.e2eformer, 'diffusion_blocks'):
            print(f"  Diffusion blocks: {len(model.e2eformer.diffusion_blocks)}")
        else:
            print("  ✗ No diffusion blocks found")
            return False
        
        return True
        
    except Exception as e:
        print(f"✗ Model creation failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_forward_pass_with_dummy_data():
    """Test forward pass with dummy data."""
    print("Testing forward pass with dummy data...")
    
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
        
        # Create dummy data
        batch_size = 1
        seq_len = 50
        msa_depth = 3
        
        # Dummy tokens - tokens should be 3D (batch, msa_depth, seq_len), rna_fm_tokens should be 2D (batch, seq_len)
        tokens = torch.randint(0, 4, (batch_size, msa_depth, seq_len), dtype=torch.long, device=device)
        rna_fm_tokens = torch.randint(0, 4, (batch_size, seq_len), dtype=torch.long, device=device)
        sequence = ["AUCG" * (seq_len // 4)] * batch_size
        
        with torch.no_grad():
            # Forward pass
            outputs = model(
                tokens=tokens,
                rna_fm_tokens=rna_fm_tokens,
                seq=sequence
            )
            
            print(f"✓ Forward pass successful on CPU")
            print(f"  Recycling outputs: {len(outputs)}")
            
            # Check output shapes
            for i, output in enumerate(outputs):
                print(f"  Recycling {i}:")
                for key, value in output.items():
                    if isinstance(value, torch.Tensor):
                        print(f"    {key}: {value.shape}")
            
            # Check if diffusion output has correct shape
            if 'z' in outputs[-1]:
                z_shape = outputs[-1]['z'].shape
                print(f"  Diffusion output shape: {z_shape}")
                if len(z_shape) == 4 and z_shape[1] == z_shape[2]:  # (B, L, L, C)
                    print("✓ Diffusion output shape is correct (B, L, L, C)")
                else:
                    print(f"✗ Diffusion output shape is incorrect: {z_shape}")
                    return False
            
            return True
            
    except Exception as e:
        print(f"✗ Forward pass failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_gradient_flow():
    """Test gradient flow through diffusion module."""
    print("Testing gradient flow...")
    
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
        
        # Dummy tokens - tokens should be 3D (batch, msa_depth, seq_len), rna_fm_tokens should be 2D (batch, seq_len)
        tokens = torch.randint(0, 4, (batch_size, msa_depth, seq_len), dtype=torch.long, device=device)
        rna_fm_tokens = torch.randint(0, 4, (batch_size, seq_len), dtype=torch.long, device=device)
        sequence = ["AUCG" * (seq_len // 4)] * batch_size
        
        # Forward pass
        outputs = model(
            tokens=tokens,
            rna_fm_tokens=rna_fm_tokens,
            seq=sequence
        )
        
        # Create a simple loss (sum of all outputs)
        total_loss = 0
        for output in outputs:
            for key, value in output.items():
                if isinstance(value, torch.Tensor):
                    total_loss += value.sum()
        
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
    """Run minimal tests."""
    print("RhoFold Diffusion Minimal Test")
    print("=" * 40)
    
    tests = [
        ("Model Creation", test_model_creation),
        ("Forward Pass (Dummy Data)", test_forward_pass_with_dummy_data),
        ("Gradient Flow", test_gradient_flow),
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
        print("\n🎉 All minimal tests passed!")
        print("\nThe diffusion module is working correctly!")
        print("Next steps:")
        print("1. Fix the data loading issues (PDB parsing)")
        print("2. Test with real data")
        print("3. Run full training")
    else:
        print(f"\n❌ {total - passed} tests failed.")
    
    return passed == total


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
