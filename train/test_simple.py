"""
Simple test script to validate the basic functionality without complex PDB parsing.
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
from train.data_loader import create_dataloader
from train.losses import RhoFoldLoss


def create_simple_test_data(data_dir, num_samples=2, seq_len=50):
    """Create simple test data without complex PDB parsing."""
    print(f"Creating simple test data in {data_dir}...")
    
    os.makedirs(data_dir, exist_ok=True)
    
    # Create test sequences
    sequences = []
    for i in range(num_samples):
        seq = "AUCG" * (seq_len // 4)
        sequences.append(seq)
    
    # Create index file
    with open(os.path.join(data_dir, "index.txt"), "w") as f:
        for i in range(num_samples):
            f.write(f"test_sample_{i}\n")
    
    # Create files for each sample
    for i in range(num_samples):
        sample_name = f"test_sample_{i}"
        seq = sequences[i]
        
        # Create .seq file
        with open(os.path.join(data_dir, f"{sample_name}.seq"), "w") as f:
            f.write(seq)
        
        # Create .a3m file with multiple sequences
        with open(os.path.join(data_dir, f"{sample_name}.a3m"), "w") as f:
            f.write(f">{sample_name}_1\n{seq}\n")
            f.write(f">{sample_name}_2\n{seq}\n")
            f.write(f">{sample_name}_3\n{seq}\n")
        
        # Create simple .pdb file
        with open(os.path.join(data_dir, f"{sample_name}.pdb"), "w") as f:
            f.write("HEADER    TEST RNA STRUCTURE\n")
            # Create simple coordinates for each residue
            atom_id = 1
            for res_id in range(1, len(seq) + 1):
                x = res_id * 3.4
                y = 0.0
                z = 0.0
                f.write(f"ATOM  {atom_id:5d}  P   A{res_id:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           P\n")
                atom_id += 1
            f.write("END\n")
    
    print(f"✓ Created {num_samples} test samples")


def test_basic_functionality():
    """Test basic functionality without complex data loading."""
    print("Testing basic functionality...")
    
    try:
        with tempfile.TemporaryDirectory() as temp_dir:
            create_simple_test_data(temp_dir, num_samples=2, seq_len=50)
            
            # Test data loader
            dataloader = create_dataloader(
                data_dir=temp_dir,
                index_file=os.path.join(temp_dir, "index.txt"),
                batch_size=1,
                num_workers=0,
                max_seq_len=50,
                msa_depth=3,
                shuffle=False,
            )
            
            print(f"✓ DataLoader created with {len(dataloader)} batches")
            
            # Test loading a batch
            batch = next(iter(dataloader))
            
            # Check batch structure
            required_keys = ['tokens', 'msa_tokens', 'msa_mask', 'pair_mask', 'sequence', 'atom_coords', 'atom_types', 'atom_mask', 'target']
            for key in required_keys:
                if key not in batch:
                    print(f"✗ Missing key in batch: {key}")
                    return False
            
            print("✓ Batch structure correct")
            
            # Check tensor shapes
            print(f"  Tokens shape: {batch['tokens'].shape}")
            print(f"  MSA tokens shape: {batch['msa_tokens'].shape}")
            print(f"  Atom coords shape: {batch['atom_coords'].shape}")
            
            return True
            
    except Exception as e:
        print(f"✗ Data loading failed: {e}")
        import traceback
        traceback.print_exc()
        return False


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


def test_forward_pass_cpu():
    """Test forward pass on CPU only."""
    print("Testing forward pass on CPU...")
    
    try:
        with tempfile.TemporaryDirectory() as temp_dir:
            create_simple_test_data(temp_dir, num_samples=2, seq_len=50)
            
            # Create data loader
            dataloader = create_dataloader(
                data_dir=temp_dir,
                index_file=os.path.join(temp_dir, "index.txt"),
                batch_size=1,
                num_workers=0,
                max_seq_len=50,
                msa_depth=3,
                shuffle=False,
            )
            
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
            
            with torch.no_grad():
                batch = next(iter(dataloader))
                
                # Move to device
                for key in batch:
                    if isinstance(batch[key], torch.Tensor):
                        batch[key] = batch[key].to(device)
                    elif isinstance(batch[key], dict):
                        for subkey in batch[key]:
                            if isinstance(batch[key][subkey], torch.Tensor):
                                batch[key][subkey] = batch[key][subkey].to(device)
                
                # Forward pass
                outputs = model(
                    tokens=batch['tokens'],
                    rna_fm_tokens=batch['msa_tokens'],
                    seq=batch['sequence']
                )
                
                print(f"✓ Forward pass successful on CPU")
                print(f"  Recycling outputs: {len(outputs)}")
                
                # Check output shapes
                for i, output in enumerate(outputs):
                    print(f"  Recycling {i}:")
                    for key, value in output.items():
                        if isinstance(value, torch.Tensor):
                            print(f"    {key}: {value.shape}")
            
            return True
            
    except Exception as e:
        print(f"✗ Forward pass failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    """Run basic tests."""
    print("RhoFold Diffusion Basic Test")
    print("=" * 40)
    
    tests = [
        ("Model Creation", test_model_creation),
        ("Data Loading", test_basic_functionality),
        ("Forward Pass (CPU)", test_forward_pass_cpu),
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
        print("\n🎉 All basic tests passed!")
    else:
        print(f"\n❌ {total - passed} tests failed.")
    
    return passed == total


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
