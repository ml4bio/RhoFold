"""
Test script using real data from diffusion_test folder.
Tests the complete pipeline with real RNA structure data.
"""

import os
import sys
import torch
import torch.nn as nn
import numpy as np
from pathlib import Path

# Add parent directory to path
sys.path.append(str(Path(__file__).parent.parent))

from rhofold.config import rhofold_config
from rhofold.rhofold_with_diffusion import create_diffusion_rhofold
from train.data_loader import create_dataloader
from train.losses import RhoFoldLoss


def test_real_data_pipeline():
    """Test complete pipeline with real data."""
    print("Testing with real data from diffusion_test folder...")
    
    try:
        # Use real data directory
        data_dir = "/home/rcy/code/RhoFold/example/input/diffusion_test"
        index_file = os.path.join(data_dir, "index.txt")
        
        # Create index file for the real data
        with open(index_file, "w") as f:
            f.write("1a1t_B\n")
        
        # Create data loader
        dataloader = create_dataloader(
            data_dir=data_dir,
            index_file=index_file,
            batch_size=1,
            num_workers=0,
            max_seq_len=50,
            msa_depth=3,
            shuffle=False,
        )
        
        print(f"✓ DataLoader created with {len(dataloader)} batches")
        
        # Create model with diffusion module
        diffusion_config = {
            'time_embed_dim': 128,
            'num_timesteps': 100,
            'noise_schedule': 'linear',
            'base_channels': 32,
        }
        
        model = create_diffusion_rhofold(rhofold_config, diffusion_config)
        model.eval()
        
        device = torch.device('cpu')  # Use CPU to avoid memory issues
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
            
            # Convert msa_tokens from 3D to 2D for RNA-FM
            rna_fm_tokens = batch['msa_tokens'][:, 0, :]  # Take first row [batch, seq_len]
            
            # Convert sequence from list to string
            seq_str = batch['sequence'][0] if isinstance(batch['sequence'], list) else batch['sequence']
            
            print(f"✓ Processing sequence: {seq_str}")
            print(f"  Sequence length: {len(seq_str)}")
            
            # Full forward pass through the model
            outputs = model(
                tokens=batch['tokens'],
                rna_fm_tokens=rna_fm_tokens,
                seq=seq_str
            )
            
            print(f"✓ Full forward pass successful")
            print(f"  Number of recycling outputs: {len(outputs)}")
            
            # Check the last recycling output (final prediction)
            final_output = outputs[-1]
            print(f"  Final output keys: {list(final_output.keys())}")
            
            # Check key outputs from structure module
            if 'plddt' in final_output:
                plddt = final_output['plddt']
                if isinstance(plddt, tuple):
                    print(f"  pLDDT (tuple): {len(plddt)} elements, shape: {plddt[0].shape}")
                else:
                    print(f"  pLDDT shape: {plddt.shape}")
            
            coords_key = 'cords_c1\''
            if coords_key in final_output:
                coords = final_output[coords_key]
                if isinstance(coords, list):
                    print(f"  Coordinates (list): {len(coords)} elements, shape: {coords[0].shape}")
                elif isinstance(coords, tuple):
                    print(f"  Coordinates (tuple): {len(coords)} elements, shape: {coords[0].shape}")
                else:
                    print(f"  Coordinates shape: {coords.shape}")
            
            # Check distance predictions
            distance_keys = ['p', 'c4_', 'n']
            for key in distance_keys:
                if key in final_output:
                    print(f"  {key} shape: {final_output[key].shape}")
            
            # Check secondary structure predictions
            if 'ss' in final_output:
                print(f"  SS shape: {final_output['ss'].shape}")
            
            return True
            
    except Exception as e:
        print(f"✗ Real data test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_loss_computation_real_data():
    """Test loss computation with real data."""
    print("Testing loss computation with real data...")
    
    try:
        # Use real data directory
        data_dir = "/home/rcy/code/RhoFold/example/input/diffusion_test"
        index_file = os.path.join(data_dir, "index.txt")
        
        # Create index file for the real data
        with open(index_file, "w") as f:
            f.write("1a1t_B\n")
        
        # Create data loader
        dataloader = create_dataloader(
            data_dir=data_dir,
            index_file=index_file,
            batch_size=1,
            num_workers=0,
            max_seq_len=50,
            msa_depth=3,
            shuffle=False,
        )
        
        # Create model with diffusion module
        diffusion_config = {
            'time_embed_dim': 128,
            'num_timesteps': 100,
            'noise_schedule': 'linear',
            'base_channels': 32,
        }
        
        model = create_diffusion_rhofold(rhofold_config, diffusion_config)
        model.eval()
        
        # Create loss function
        criterion = RhoFoldLoss()
        
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
            
            # Convert msa_tokens from 3D to 2D for RNA-FM
            rna_fm_tokens = batch['msa_tokens'][:, 0, :]
            
            # Convert sequence from list to string
            seq_str = batch['sequence'][0] if isinstance(batch['sequence'], list) else batch['sequence']
            
            # Forward pass
            outputs = model(
                tokens=batch['tokens'],
                rna_fm_tokens=rna_fm_tokens,
                seq=seq_str
            )
            
            # Use the last recycling output for loss computation
            output = outputs[-1]
            
            # Prepare predictions for loss computation
            pred = prepare_predictions(output, batch, device)
            
            print(f"✓ Predictions prepared:")
            for key, value in pred.items():
                if isinstance(value, torch.Tensor):
                    print(f"  {key}: {value.shape}")
                else:
                    print(f"  {key}: {type(value)}")
            
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


def test_gradient_flow_real_data():
    """Test gradient flow with real data."""
    print("Testing gradient flow with real data...")
    
    try:
        # Use real data directory
        data_dir = "/home/rcy/code/RhoFold/example/input/diffusion_test"
        index_file = os.path.join(data_dir, "index.txt")
        
        # Create index file for the real data
        with open(index_file, "w") as f:
            f.write("1a1t_B\n")
        
        # Create data loader
        dataloader = create_dataloader(
            data_dir=data_dir,
            index_file=index_file,
            batch_size=1,
            num_workers=0,
            max_seq_len=50,
            msa_depth=3,
            shuffle=False,
        )
        
        # Create model with diffusion module
        diffusion_config = {
            'time_embed_dim': 128,
            'num_timesteps': 100,
            'noise_schedule': 'linear',
            'base_channels': 32,
        }
        
        model = create_diffusion_rhofold(rhofold_config, diffusion_config)
        model.train()  # Set to training mode for gradient computation
        
        # Create loss function
        criterion = RhoFoldLoss()
        
        device = torch.device('cpu')
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
        
        # Convert msa_tokens from 3D to 2D for RNA-FM
        rna_fm_tokens = batch['msa_tokens'][:, 0, :]
        
        # Convert sequence from list to string
        seq_str = batch['sequence'][0] if isinstance(batch['sequence'], list) else batch['sequence']
        
        # Forward pass
        outputs = model(
            tokens=batch['tokens'],
            rna_fm_tokens=rna_fm_tokens,
            seq=seq_str
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


def prepare_predictions(output, batch, device):
    """Prepare model output for loss computation."""
    pred = {}
    
    # Distance predictions
    if 'p' in output and 'c4_' in output and 'n' in output:
        # The outputs are (B, n_bins, L, L), we need (B, L, L, n_bins)
        p_logits = output['p'].permute(0, 2, 3, 1)  # (B, L, L, n_bins)
        c4_logits = output['c4_'].permute(0, 2, 3, 1)  # (B, L, L, n_bins)
        n_logits = output['n'].permute(0, 2, 3, 1)  # (B, L, L, n_bins)
        
        # Stack to (B, 3, L, L, n_bins) - keep as 5D tensor
        distance_logits = torch.stack([p_logits, c4_logits, n_logits], dim=1)
        
        pred['distance_logits'] = distance_logits
        print(f"DEBUG: distance_logits shape: {distance_logits.shape}")
    
    # Secondary structure predictions
    if 'ss' in output:
        ss_logits = output['ss']
        # Remove extra dimension if present: (B, 1, L, L) -> (B, L, L)
        if ss_logits.dim() == 4 and ss_logits.shape[1] == 1:
            ss_logits = ss_logits.squeeze(1)
        pred['ss_logits'] = ss_logits
    
    # pLDDT predictions
    if 'plddt' in output:
        plddt = output['plddt']
        if isinstance(plddt, tuple):
            pred['plddt_logits'] = plddt[0]  # Take first element
        else:
            pred['plddt_logits'] = plddt
    
    # Atom coordinates
    coords_key = 'cords_c1\''
    if coords_key in output:
        coords = output[coords_key]
        if isinstance(coords, list):
            coords = coords[0]  # Take first element
        elif isinstance(coords, tuple):
            coords = coords[0]  # Take first element
        
        # Ensure coords has shape (B, L, N, 3) for loss functions
        if coords.dim() == 3:  # (B, L, 3) -> (B, L, 1, 3)
            coords = coords.unsqueeze(2)
        elif coords.dim() == 2:  # (L, 3) -> (1, L, 1, 3)
            coords = coords.unsqueeze(0).unsqueeze(2)
        
        pred['atom_coords'] = coords
        # Create dummy atom types and mask (in practice, these would be computed)
        B, L, N = pred['atom_coords'].shape[:3]
        pred['atom_types'] = torch.zeros(B, L, N, dtype=torch.long, device=device)
        pred['atom_mask'] = torch.ones(B, L, N, dtype=torch.bool, device=device)
    
    # Create dummy pseudo coordinates for ss3d loss
    if 'atom_coords' in pred:
        if pred['atom_coords'].dim() == 4:
            B, L, N, _ = pred['atom_coords'].shape
            pred['pseudo_coords'] = pred['atom_coords'][:, :, :4, :]  # Use first 4 atoms
        else:
            # Handle 3D case
            B, L, _ = pred['atom_coords'].shape
            pred['pseudo_coords'] = pred['atom_coords'][:, :, :3]  # Use first 3 coordinates
    
    # Create dummy frames for FAPE loss
    if 'atom_coords' in pred:
        if pred['atom_coords'].dim() == 4:
            B, L, N, _ = pred['atom_coords'].shape
            if N >= 3:
                coords = pred['atom_coords'][:, :, :3, :]  # Use first 3 atoms
            else:
                # Pad with dummy atoms if we don't have enough
                coords = pred['atom_coords']
                if N == 1:
                    # Duplicate the single atom to create 3 atoms
                    coords = coords.repeat(1, 1, 3, 1)
                elif N == 2:
                    # Add a third atom
                    coords = torch.cat([coords, coords[:, :, :1, :]], dim=2)
            pred['pred_frames'] = build_local_frame(coords)
        else:
            # Handle 3D case - create dummy 4D coordinates
            B, L, _ = pred['atom_coords'].shape
            coords = pred['atom_coords'].unsqueeze(2).expand(B, L, 3, 3)  # Expand to 4D
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
    """Run tests with real data."""
    print("RhoFold Real Data Pipeline Test")
    print("=" * 40)
    
    tests = [
        ("Real Data Pipeline", test_real_data_pipeline),
        ("Loss Computation", test_loss_computation_real_data),
        ("Gradient Flow", test_gradient_flow_real_data),
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
    print("\n" + "=" * 40)
    print("REAL DATA TEST SUMMARY:")
    print("=" * 40)
    
    passed = 0
    total = len(results)
    
    for test_name, success in results:
        status = "✓ PASS" if success else "✗ FAIL"
        print(f"{test_name:25} {status}")
        if success:
            passed += 1
    
    print("-" * 40)
    print(f"Passed: {passed}/{total}")
    
    if passed == total:
        print("\n🎉 All real data tests passed!")
        print("\nDiffusion module successfully integrated with RhoFold using real data!")
        print("Key achievements:")
        print("✓ Real data loading works")
        print("✓ Complete forward pass works")
        print("✓ Loss computation works")
        print("✓ Gradient flow works")
        print("\nThe model is ready for training with real data!")
    else:
        print(f"\n❌ {total - passed} tests failed.")
        print("Please fix the issues before training.")
    
    return passed == total


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
