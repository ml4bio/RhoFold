"""
Comprehensive validation script for RhoFold diffusion training setup.
This script validates all components before starting training.
"""

import os
import sys
import torch
import torch.nn as nn
import numpy as np
from pathlib import Path
import tempfile
import subprocess
import time

# Add parent directory to path
sys.path.append(str(Path(__file__).parent.parent))

from rhofold.config import rhofold_config
from rhofold.rhofold_with_diffusion import create_diffusion_rhofold
from train.data_loader import create_dataloader
from train.losses import RhoFoldLoss


def check_environment():
    """Check if the environment is properly set up."""
    print("Checking environment...")
    
    # Check Python version
    python_version = sys.version_info
    if python_version.major != 3 or python_version.minor < 8:
        print(f"✗ Python version {python_version.major}.{python_version.minor} is not supported. Need Python 3.8+")
        return False
    else:
        print(f"✓ Python version {python_version.major}.{python_version.minor}")
    
    # Check PyTorch
    try:
        import torch
        print(f"✓ PyTorch version: {torch.__version__}")
        
        # Check CUDA
        if torch.cuda.is_available():
            print(f"✓ CUDA available: {torch.cuda.device_count()} GPUs")
            for i in range(torch.cuda.device_count()):
                print(f"  GPU {i}: {torch.cuda.get_device_name(i)}")
        else:
            print("✗ CUDA not available")
            return False
    except ImportError:
        print("✗ PyTorch not installed")
        return False
    
    # Check other dependencies
    required_packages = [
        'numpy', 'scipy', 'pandas', 'matplotlib', 'seaborn', 'tqdm',
        'h5py', 'Bio', 'yaml', 'ml_collections', 'einops',
        'openmm', 'mdtraj', 'prody', 'biotite'
    ]
    
    missing_packages = []
    for package in required_packages:
        try:
            __import__(package)
            print(f"✓ {package}")
        except ImportError:
            missing_packages.append(package)
            print(f"✗ {package}")
    
    if missing_packages:
        print(f"Missing packages: {missing_packages}")
        return False
    
    return True


def create_test_data(data_dir, num_samples=4, seq_len=100):
    """Create test data for validation."""
    print(f"Creating test data in {data_dir}...")
    
    os.makedirs(data_dir, exist_ok=True)
    
    # Create test sequences
    sequences = []
    for i in range(num_samples):
        # Create varied sequences
        if i % 2 == 0:
            seq = "AUCG" * (seq_len // 4)
        else:
            seq = "GCAU" * (seq_len // 4)
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
        
        # Create .pdb file with realistic RNA structure
        with open(os.path.join(data_dir, f"{sample_name}.pdb"), "w") as f:
            f.write("HEADER    TEST RNA STRUCTURE\n")
            f.write("ATOM      1  P   A   1      20.154  16.967  23.862  1.00 20.00           P\n")
            f.write("ATOM      2  O1P A   1      19.123  17.123  24.123  1.00 20.00           O\n")
            f.write("ATOM      3  O2P A   1      21.123  17.123  24.123  1.00 20.00           O\n")
            f.write("ATOM      4  O5' A   1      20.123  15.123  22.123  1.00 20.00           O\n")
            f.write("ATOM      5  C5' A   1      19.123  14.123  21.123  1.00 20.00           C\n")
            f.write("ATOM      6  C4' A   1      18.123  13.123  20.123  1.00 20.00           C\n")
            f.write("ATOM      7  O4' A   1      17.123  12.123  19.123  1.00 20.00           O\n")
            f.write("ATOM      8  C1' A   1      16.123  11.123  18.123  1.00 20.00           C\n")
            f.write("ATOM      9  N1  A   1      15.123  10.123  17.123  1.00 20.00           N\n")
            f.write("ATOM     10  C2  A   1      14.123   9.123  16.123  1.00 20.00           C\n")
            f.write("ATOM     11  N3  A   1      13.123   8.123  15.123  1.00 20.00           N\n")
            f.write("ATOM     12  C4  A   1      12.123   7.123  14.123  1.00 20.00           C\n")
            f.write("ATOM     13  C5  A   1      11.123   6.123  13.123  1.00 20.00           C\n")
            f.write("ATOM     14  C6  A   1      10.123   5.123  12.123  1.00 20.00           C\n")
            f.write("ATOM     15  N6  A   1       9.123   4.123  11.123  1.00 20.00           N\n")
            f.write("ATOM     16  N7  A   1       8.123   3.123  10.123  1.00 20.00           N\n")
            f.write("ATOM     17  C8  A   1       7.123   2.123   9.123  1.00 20.00           C\n")
            f.write("ATOM     18  N9  A   1       6.123   1.123   8.123  1.00 20.00           N\n")
            f.write("ATOM     19  C2' A   1       5.123   0.123   7.123  1.00 20.00           C\n")
            f.write("ATOM     20  C3' A   1       4.123  -0.877   6.123  1.00 20.00           C\n")
            f.write("ATOM     21  O3' A   1       3.123  -1.877   5.123  1.00 20.00           O\n")
            f.write("ATOM     22  H1  A   1       2.123  -2.877   4.123  1.00 20.00           H\n")
            f.write("ATOM     23  H2  A   1       1.123  -3.877   3.123  1.00 20.00           H\n")
            f.write("END\n")
    
    print(f"✓ Created {num_samples} test samples")


def test_model_creation():
    """Test model creation and parameter counting."""
    print("Testing model creation...")
    
    try:
        # Test different diffusion configurations
        configs = [
            {
                'time_embed_dim': 128,
                'num_timesteps': 500,
                'noise_schedule': 'linear',
                'base_channels': 32,
            },
            {
                'time_embed_dim': 256,
                'num_timesteps': 1000,
                'noise_schedule': 'cosine',
                'base_channels': 64,
            }
        ]
        
        for i, config in enumerate(configs):
            model = create_diffusion_rhofold(rhofold_config, diffusion_config=config)
            total_params = sum(p.numel() for p in model.parameters())
            trainable_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
            
            print(f"✓ Config {i+1}: {total_params:,} total parameters, {trainable_params:,} trainable")
            
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


def test_data_loading():
    """Test data loading pipeline."""
    print("Testing data loading...")
    
    try:
        with tempfile.TemporaryDirectory() as temp_dir:
            create_test_data(temp_dir, num_samples=4, seq_len=100)
            
            # Test data loader
            dataloader = create_dataloader(
                data_dir=temp_dir,
                index_file=os.path.join(temp_dir, "index.txt"),
                batch_size=2,
                num_workers=0,
                max_seq_len=100,
                msa_depth=3,
                shuffle=True,
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


def test_forward_pass():
    """Test forward pass with different batch sizes."""
    print("Testing forward pass...")
    
    try:
        with tempfile.TemporaryDirectory() as temp_dir:
            create_test_data(temp_dir, num_samples=4, seq_len=100)
            
            # Create data loader
            dataloader = create_dataloader(
                data_dir=temp_dir,
                index_file=os.path.join(temp_dir, "index.txt"),
                batch_size=2,
                num_workers=0,
                max_seq_len=100,
                msa_depth=3,
                shuffle=False,
            )
            
            # Create model
            diffusion_config = {
                'time_embed_dim': 256,
                'num_timesteps': 1000,
                'noise_schedule': 'cosine',
                'base_channels': 64,
            }
            
            model = create_diffusion_rhofold(rhofold_config, diffusion_config)
            model.eval()
            
            # Test on different devices
            devices = ['cpu']
            if torch.cuda.is_available():
                devices.append('cuda')
            
            for device_name in devices:
                device = torch.device(device_name)
                model.to(device)
                
                print(f"  Testing on {device_name}...")
                
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
                    
                    print(f"    ✓ Forward pass successful on {device_name}")
                    print(f"    Recycling outputs: {len(outputs)}")
                    
                    # Check output shapes
                    for i, output in enumerate(outputs):
                        print(f"    Recycling {i}:")
                        for key, value in output.items():
                            if isinstance(value, torch.Tensor):
                                print(f"      {key}: {value.shape}")
            
            return True
            
    except Exception as e:
        print(f"✗ Forward pass failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_loss_computation():
    """Test loss computation."""
    print("Testing loss computation...")
    
    try:
        with tempfile.TemporaryDirectory() as temp_dir:
            create_test_data(temp_dir, num_samples=4, seq_len=100)
            
            # Create data loader
            dataloader = create_dataloader(
                data_dir=temp_dir,
                index_file=os.path.join(temp_dir, "index.txt"),
                batch_size=2,
                num_workers=0,
                max_seq_len=100,
                msa_depth=3,
                shuffle=False,
            )
            
            # Create model
            diffusion_config = {
                'time_embed_dim': 256,
                'num_timesteps': 1000,
                'noise_schedule': 'cosine',
                'base_channels': 64,
            }
            
            model = create_diffusion_rhofold(rhofold_config, diffusion_config)
            model.eval()
            
            # Create loss function
            criterion = RhoFoldLoss()
            
            device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
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
                
                # Use the last recycling output for loss computation
                output = outputs[-1]
                
                # Prepare predictions for loss computation
                pred = prepare_predictions(output, batch, device)
                
                # Compute loss
                total_loss, loss_dict = criterion(pred, batch['target'])
                
                print(f"✓ Loss computation successful")
                print(f"  Total loss: {total_loss.item():.6f}")
                print("  Individual losses:")
                for key, value in loss_dict.items():
                    print(f"    {key}: {value.item():.6f}")
                
                # Check for NaN values
                if torch.isnan(total_loss):
                    print("✗ Loss contains NaN values")
                    return False
                
                return True
            
    except Exception as e:
        print(f"✗ Loss computation failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_gradient_flow():
    """Test gradient flow through diffusion module."""
    print("Testing gradient flow...")
    
    try:
        with tempfile.TemporaryDirectory() as temp_dir:
            create_test_data(temp_dir, num_samples=4, seq_len=100)
            
            # Create data loader
            dataloader = create_dataloader(
                data_dir=temp_dir,
                index_file=os.path.join(temp_dir, "index.txt"),
                batch_size=2,
                num_workers=0,
                max_seq_len=100,
                msa_depth=3,
                shuffle=False,
            )
            
            # Create model
            diffusion_config = {
                'time_embed_dim': 256,
                'num_timesteps': 1000,
                'noise_schedule': 'cosine',
                'base_channels': 64,
            }
            
            model = create_diffusion_rhofold(rhofold_config, diffusion_config)
            model.train()
            
            # Create loss function
            criterion = RhoFoldLoss()
            
            device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
            model.to(device)
            
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
            
            # Use the last recycling output for loss computation
            output = outputs[-1]
            
            # Prepare predictions for loss computation
            pred = prepare_predictions(output, batch, device)
            
            # Compute loss
            total_loss, loss_dict = criterion(pred, batch['target'])
            
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


def test_memory_usage():
    """Test memory usage with different batch sizes."""
    print("Testing memory usage...")
    
    if not torch.cuda.is_available():
        print("Skipping memory test (no CUDA)")
        return True
    
    try:
        with tempfile.TemporaryDirectory() as temp_dir:
            create_test_data(temp_dir, num_samples=8, seq_len=200)
            
            # Create data loader
            dataloader = create_dataloader(
                data_dir=temp_dir,
                index_file=os.path.join(temp_dir, "index.txt"),
                batch_size=4,
                num_workers=0,
                max_seq_len=200,
                msa_depth=4,
                shuffle=False,
            )
            
            # Create model
            diffusion_config = {
                'time_embed_dim': 256,
                'num_timesteps': 1000,
                'noise_schedule': 'cosine',
                'base_channels': 64,
            }
            
            model = create_diffusion_rhofold(rhofold_config, diffusion_config)
            model.eval()
            
            device = torch.device('cuda')
            model.to(device)
            
            # Clear cache
            torch.cuda.empty_cache()
            initial_memory = torch.cuda.memory_allocated()
            
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
                
                peak_memory = torch.cuda.max_memory_allocated()
                memory_used = (peak_memory - initial_memory) / 1024**3  # GB
                
                print(f"✓ Memory usage: {memory_used:.2f} GB")
                print(f"  Peak memory: {peak_memory / 1024**3:.2f} GB")
                
                if memory_used > 8.0:  # 8GB threshold
                    print("⚠ Warning: High memory usage detected")
                
                return True
            
    except Exception as e:
        print(f"✗ Memory test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def prepare_predictions(output, batch, device):
    """Prepare model output for loss computation."""
    pred = {}
    
    # Distance predictions
    if 'p' in output and 'c4_' in output and 'n' in output:
        # Combine distance predictions
        distance_logits = torch.stack([
            output['p'], output['c4_'], output['n']
        ], dim=1)  # (B, 3, L, L, n_bins)
        pred['distance_logits'] = distance_logits
    
    # Secondary structure predictions
    if 'ss' in output:
        pred['ss_logits'] = output['ss']
    
    # pLDDT predictions
    if 'plddt' in output:
        pred['plddt_logits'] = output['plddt']
    
    # Atom coordinates
    if 'cords_c1\'' in output:
        pred['atom_coords'] = output['cords_c1\''][-1]  # Use last recycling
        # Create dummy atom types (in practice, these would be computed)
        B, L = pred['atom_coords'].shape[:2]
        pred['atom_types'] = torch.zeros(B, L, 23, dtype=torch.long, device=device)
    
    # Create dummy pseudo coordinates for ss3d loss
    if 'atom_coords' in pred:
        B, L, N, _ = pred['atom_coords'].shape
        pred['pseudo_coords'] = pred['atom_coords'][:, :, :4, :]  # Use first 4 atoms
    
    # Create dummy frames for FAPE loss
    if 'atom_coords' in pred:
        coords = pred['atom_coords'][:, :, :3, :]  # Use first 3 atoms
        pred['pred_frames'] = build_local_frame(coords)
    
    return pred


def build_local_frame(coords):
    """Build local residue frames from coordinates."""
    v1 = coords[:, :, 1] - coords[:, :, 0]  # C1' - C4'
    v2 = coords[:, :, 2] - coords[:, :, 0]  # N1/N9 - C4'
    
    # Normalize
    x_axis = v1 / (torch.norm(v1, dim=-1, keepdim=True) + 1e-6)
    
    # Orthogonalize
    proj = (v2 * x_axis).sum(dim=-1, keepdim=True) * x_axis
    y_axis = v2 - proj
    y_axis = y_axis / (torch.norm(y_axis, dim=-1, keepdim=True) + 1e-6)
    
    # Cross product
    z_axis = torch.cross(x_axis, y_axis, dim=-1)
    z_axis = z_axis / (torch.norm(z_axis, dim=-1, keepdim=True) + 1e-6)
    
    # Stack to rotation matrix
    R = torch.stack([x_axis, y_axis, z_axis], dim=-1)
    t = coords[:, :, 0]  # Origin
    
    return {'R': R, 't': t}


def main():
    """Run all validation tests."""
    print("RhoFold Diffusion Training Setup Validation")
    print("=" * 50)
    
    tests = [
        ("Environment Check", check_environment),
        ("Model Creation", test_model_creation),
        ("Data Loading", test_data_loading),
        ("Forward Pass", test_forward_pass),
        ("Loss Computation", test_loss_computation),
        ("Gradient Flow", test_gradient_flow),
        ("Memory Usage", test_memory_usage),
    ]
    
    results = []
    
    for test_name, test_func in tests:
        print(f"\n{test_name}:")
        print("-" * 30)
        try:
            success = test_func()
            results.append((test_name, success))
        except Exception as e:
            print(f"✗ {test_name} failed with exception: {e}")
            results.append((test_name, False))
    
    # Summary
    print("\n" + "=" * 50)
    print("VALIDATION SUMMARY:")
    print("=" * 50)
    
    passed = 0
    total = len(results)
    
    for test_name, success in results:
        status = "✓ PASS" if success else "✗ FAIL"
        print(f"{test_name:20} {status}")
        if success:
            passed += 1
    
    print("-" * 50)
    print(f"Passed: {passed}/{total}")
    
    if passed == total:
        print("\n🎉 All tests passed! The setup is ready for training.")
        print("\nNext steps:")
        print("1. Prepare your training data in the required format")
        print("2. Create an index file listing your samples")
        print("3. Run: python train/test_forward.py (for single GPU test)")
        print("4. Run: ./train/launch_training.sh --data_dir YOUR_DATA --index_file YOUR_INDEX")
    else:
        print(f"\n❌ {total - passed} tests failed. Please fix the issues before training.")
    
    return passed == total


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
