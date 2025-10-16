"""
Complete pipeline test for RhoFold with diffusion module.
Tests the integration of diffusion module with structure module and loss computation.
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


def create_test_data_with_structure(data_dir, num_samples=2, seq_len=50):
    """Create test data with proper structure for full pipeline testing."""
    print(f"Creating test data with structure in {data_dir}...")
    
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
        
        # Create proper .pdb file with realistic RNA structure
        with open(os.path.join(data_dir, f"{sample_name}.pdb"), "w") as f:
            f.write("HEADER    TEST RNA STRUCTURE\n")
            atom_id = 1
            for res_id in range(1, len(seq) + 1):
                # Create coordinates for each residue with proper spacing
                x = res_id * 3.4
                y = 0.0
                z = 0.0
                # Use proper PDB format with correct spacing
                f.write(f"ATOM  {atom_id:5d}  P   A{res_id:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           P\n")
                atom_id += 1
            f.write("END\n")
    
    print(f"✓ Created {num_samples} test samples")


def test_structure_module_integration():
    """Test structure module integration with diffusion module."""
    print("Testing structure module integration...")
    
    try:
        with tempfile.TemporaryDirectory() as temp_dir:
            create_test_data_with_structure(temp_dir, num_samples=1, seq_len=50)
            
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
                
                # Full forward pass through the model
                # Convert sequence from list to string
                seq_str = batch['sequence'][0] if isinstance(batch['sequence'], list) else batch['sequence']
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
        print(f"✗ Structure module integration failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_loss_computation():
    """Test complete loss computation with diffusion module."""
    print("Testing loss computation...")
    
    try:
        with tempfile.TemporaryDirectory() as temp_dir:
            create_test_data_with_structure(temp_dir, num_samples=1, seq_len=50)
            
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
                
                # Forward pass
                # Convert sequence from list to string
                seq_str = batch['sequence'][0] if isinstance(batch['sequence'], list) else batch['sequence']
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
    """Test gradient flow through the complete model."""
    print("Testing gradient flow...")
    
    try:
        with tempfile.TemporaryDirectory() as temp_dir:
            create_test_data_with_structure(temp_dir, num_samples=1, seq_len=50)
            
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
            
            # Forward pass
            # Convert sequence from list to string
            seq_str = batch['sequence'][0] if isinstance(batch['sequence'], list) else batch['sequence']
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
        # Combine distance predictions
        distance_logits = torch.stack([
            output['p'], output['c4_'], output['n']
        ], dim=1)  # (B, 3, L, L, n_bins)
        pred['distance_logits'] = distance_logits
        print(f"DEBUG: distance_logits shape: {distance_logits.shape}")
    
    # Secondary structure predictions
    if 'ss' in output:
        pred['ss_logits'] = output['ss']
    
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
            pred['atom_coords'] = coords[0]  # Take first element
        elif isinstance(coords, tuple):
            pred['atom_coords'] = coords[0]  # Take first element
        else:
            pred['atom_coords'] = coords
        # Create dummy atom types (in practice, these would be computed)
        B, L = pred['atom_coords'].shape[:2]
        pred['atom_types'] = torch.zeros(B, L, 23, dtype=torch.long, device=device)
    
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
            coords = pred['atom_coords'][:, :, :3, :]  # Use first 3 atoms
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
    """Run complete pipeline tests."""
    print("RhoFold Complete Pipeline Test")
    print("=" * 40)
    
    tests = [
        ("Structure Module Integration", test_structure_module_integration),
        ("Loss Computation", test_loss_computation),
        ("Gradient Flow", test_gradient_flow),
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
    print("PIPELINE TEST SUMMARY:")
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
        print("\n🎉 All pipeline tests passed!")
        print("\nDiffusion module successfully integrated with RhoFold!")
        print("Key achievements:")
        print("✓ Diffusion module works with structure module")
        print("✓ Complete forward pass works")
        print("✓ Loss computation works")
        print("✓ Gradient flow works")
        print("\nThe model is ready for training!")
    else:
        print(f"\n❌ {total - passed} tests failed.")
        print("Please fix the issues before training.")
    
    return passed == total


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
