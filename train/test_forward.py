"""
Test script to validate forward pass and gradient flow for RhoFold diffusion module.
This script runs a small batch test to ensure everything works correctly.
"""

import os
import sys
import torch
import torch.nn as nn
import numpy as np
from pathlib import Path
import tempfile

# Add parent directory to path
sys.path.append(str(Path(__file__).parent.parent))

from rhofold.config import rhofold_config
from rhofold.rhofold_with_diffusion import create_diffusion_rhofold
from train.data_loader import create_dataloader
from train.losses import RhoFoldLoss


def create_dummy_data(data_dir, num_samples=2, seq_len=50):
    """Create dummy data for testing."""
    os.makedirs(data_dir, exist_ok=True)
    
    # Create dummy sequences
    sequences = ["AUCG" * (seq_len // 4) for _ in range(num_samples)]
    
    # Create index file
    with open(os.path.join(data_dir, "index.txt"), "w") as f:
        for i in range(num_samples):
            f.write(f"sample_{i}\n")
    
    # Create dummy files for each sample
    for i in range(num_samples):
        sample_name = f"sample_{i}"
        seq = sequences[i]
        
        # Create .seq file
        with open(os.path.join(data_dir, f"{sample_name}.seq"), "w") as f:
            f.write(seq)
        
        # Create .a3m file
        with open(os.path.join(data_dir, f"{sample_name}.a3m"), "w") as f:
            f.write(f">{sample_name}\n{seq}\n")
            f.write(f">{sample_name}_2\n{seq}\n")
        
        # Create .pdb file
        with open(os.path.join(data_dir, f"{sample_name}.pdb"), "w") as f:
            f.write("HEADER    TEST STRUCTURE\n")
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


def test_forward_pass():
    """Test forward pass with dummy data."""
    print("Testing forward pass...")
    
    # Create dummy data
    with tempfile.TemporaryDirectory() as temp_dir:
        create_dummy_data(temp_dir, num_samples=2, seq_len=50)
        
        # Create data loader
        dataloader = create_dataloader(
            data_dir=temp_dir,
            index_file=os.path.join(temp_dir, "index.txt"),
            batch_size=2,
            num_workers=0,
            max_seq_len=50,
            msa_depth=2,
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
        
        # Get a batch
        batch = next(iter(dataloader))
        
        # Move to device
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        model.to(device)
        
        for key in batch:
            if isinstance(batch[key], torch.Tensor):
                batch[key] = batch[key].to(device)
            elif isinstance(batch[key], dict):
                for subkey in batch[key]:
                    if isinstance(batch[key][subkey], torch.Tensor):
                        batch[key][subkey] = batch[key][subkey].to(device)
        
        # Forward pass
        try:
            with torch.no_grad():
                outputs = model(
                    tokens=batch['tokens'],
                    rna_fm_tokens=batch['msa_tokens'],
                    seq=batch['sequence']
                )
            
            print("✓ Forward pass successful")
            print(f"  Number of recycling outputs: {len(outputs)}")
            
            # Check output shapes
            for i, output in enumerate(outputs):
                print(f"  Recycling {i}:")
                for key, value in output.items():
                    if isinstance(value, torch.Tensor):
                        print(f"    {key}: {value.shape}")
                    elif isinstance(value, dict):
                        for subkey, subvalue in value.items():
                            if isinstance(subvalue, torch.Tensor):
                                print(f"    {key}.{subkey}: {subvalue.shape}")
            
            # Check diffusion output shape
            if hasattr(model.e2eformer, 'diffusion_blocks'):
                print("  Diffusion module detected")
                print(f"  Number of diffusion blocks: {len(model.e2eformer.diffusion_blocks)}")
            
            return True
            
        except Exception as e:
            print(f"✗ Forward pass failed: {e}")
            import traceback
            traceback.print_exc()
            return False


def test_loss_computation():
    """Test loss computation."""
    print("Testing loss computation...")
    
    # Create dummy data
    with tempfile.TemporaryDirectory() as temp_dir:
        create_dummy_data(temp_dir, num_samples=2, seq_len=50)
        
        # Create data loader
        dataloader = create_dataloader(
            data_dir=temp_dir,
            index_file=os.path.join(temp_dir, "index.txt"),
            batch_size=2,
            num_workers=0,
            max_seq_len=50,
            msa_depth=2,
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
        
        # Get a batch
        batch = next(iter(dataloader))
        
        # Move to device
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        model.to(device)
        
        for key in batch:
            if isinstance(batch[key], torch.Tensor):
                batch[key] = batch[key].to(device)
            elif isinstance(batch[key], dict):
                for subkey in batch[key]:
                    if isinstance(batch[key][subkey], torch.Tensor):
                        batch[key][subkey] = batch[key][subkey].to(device)
        
        # Forward pass
        with torch.no_grad():
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
        try:
            total_loss, loss_dict = criterion(pred, batch['target'])
            
            print("✓ Loss computation successful")
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
    
    # Create dummy data
    with tempfile.TemporaryDirectory() as temp_dir:
        create_dummy_data(temp_dir, num_samples=2, seq_len=50)
        
        # Create data loader
        dataloader = create_dataloader(
            data_dir=temp_dir,
            index_file=os.path.join(temp_dir, "index.txt"),
            batch_size=2,
            num_workers=0,
            max_seq_len=50,
            msa_depth=2,
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
        
        # Get a batch
        batch = next(iter(dataloader))
        
        # Move to device
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        model.to(device)
        
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
        try:
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
            print(f"✗ Backward pass failed: {e}")
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
    """Run all tests."""
    print("Running RhoFold diffusion module tests...")
    print("=" * 50)
    
    # Test forward pass
    success1 = test_forward_pass()
    print()
    
    # Test loss computation
    success2 = test_loss_computation()
    print()
    
    # Test gradient flow
    success3 = test_gradient_flow()
    print()
    
    # Summary
    print("=" * 50)
    print("Test Summary:")
    print(f"  Forward pass: {'✓ PASS' if success1 else '✗ FAIL'}")
    print(f"  Loss computation: {'✓ PASS' if success2 else '✗ FAIL'}")
    print(f"  Gradient flow: {'✓ PASS' if success3 else '✗ FAIL'}")
    
    if success1 and success2 and success3:
        print("\n🎉 All tests passed! The diffusion module is ready for training.")
    else:
        print("\n❌ Some tests failed. Please check the errors above.")


if __name__ == "__main__":
    main()
