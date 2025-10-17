"""
Debug all tensors in collate function to find the problematic one.
"""

import os
import sys
import torch
from pathlib import Path

# Add parent directory to path
sys.path.append(str(Path(__file__).parent.parent))

from train.data_loader import RhoFoldDataset


def debug_all_tensors():
    """Debug all tensors in collate function."""
    print("Debugging all tensors...")
    
    # Create dataset directly
    dataset = RhoFoldDataset(
        data_dir="/home/rcy/code/RhoFold/RNA3D_DATA",
        index_file="fold-0_train_ids",
        max_seq_len=50,
        msa_depth=5,
    )
    
    # Load a few samples manually
    samples = []
    for i in range(2):
        sample = dataset[i]
        samples.append(sample)
        print(f"Sample {i}: {sample['sample_name']}")
    
    # Get max sequence length
    max_len = max(item['tokens'].shape[1] for item in samples)
    print(f"Max sequence length: {max_len}")
    
    # Process each key in the batch
    collated = {}
    
    for key in samples[0].keys():
        print(f"\nProcessing key: {key}")
        
        if key == 'target':
            # Handle target dict separately
            collated[key] = {}
            for target_key in samples[0]['target'].keys():
                print(f"  Processing target key: {target_key}")
                
                if target_key in ['true_frames']:
                    # Handle frames dict
                    print(f"    Handling true_frames dict")
                    collated[key][target_key] = {
                        'R': torch.stack([item['target'][target_key]['R'] for item in samples]),
                        't': torch.stack([item['target'][target_key]['t'] for item in samples])
                    }
                else:
                    # Pad other tensors
                    tensors = [item['target'][target_key] for item in samples]
                    print(f"    Tensors: {[t.shape for t in tensors]}")
                    
                    # Convert numpy arrays to torch tensors if needed
                    tensors = [torch.tensor(t) if not isinstance(t, torch.Tensor) else t for t in tensors]
                    
                    if tensors[0].dim() == 1:
                        print(f"    1D tensor processing")
                        if tensors[0].shape[0] == max_len:
                            # Sequence length tensor
                            print(f"      Sequence length tensor")
                            padded = torch.zeros(len(samples), max_len, dtype=tensors[0].dtype)
                            for i, tensor in enumerate(tensors):
                                if len(tensor) > 0:
                                    padded[i, :len(tensor)] = tensor
                            collated[key][target_key] = padded
                        else:
                            # Flattened tensor (like atom_coords) - need to pad to same length
                            print(f"      Flattened tensor processing")
                            if target_key == 'atom_coords':
                                # atom_coords is (L*N, 3) - need to pad to max length
                                max_atoms = max(tensor.shape[0] for tensor in tensors)
                                padded = torch.zeros(len(samples), max_atoms, 3, dtype=tensors[0].dtype)
                                for i, tensor in enumerate(tensors):
                                    if tensor.shape[0] > 0:
                                        padded[i, :tensor.shape[0]] = tensor
                                collated[key][target_key] = padded
                            elif target_key == 'atom_mask':
                                # atom_mask is (L*N,) - need to pad to max length
                                max_atoms = max(tensor.shape[0] for tensor in tensors)
                                padded = torch.zeros(len(samples), max_atoms, dtype=tensors[0].dtype)
                                for i, tensor in enumerate(tensors):
                                    if tensor.shape[0] > 0:
                                        padded[i, :tensor.shape[0]] = tensor
                                collated[key][target_key] = padded
                            elif target_key in ['masked_tokens', 'plddt_true', 'res_mask']:
                                # These are 1D tensors that need padding to max sequence length
                                print(f"      Processing {target_key} as 1D sequence tensor")
                                padded = torch.zeros(len(samples), max_len, dtype=tensors[0].dtype)
                                for i, tensor in enumerate(tensors):
                                    print(f"        Tensor {i}: shape={tensor.shape}, len={len(tensor)}")
                                    if len(tensor) > 0:
                                        print(f"        Attempting: padded[{i}, :{len(tensor)}] = tensor")
                                        try:
                                            padded[i, :len(tensor)] = tensor
                                            print(f"        ✓ Success")
                                        except Exception as e:
                                            print(f"        ❌ Error: {e}")
                                            print(f"        padded[{i}, :{len(tensor)}] shape: {padded[i, :len(tensor)].shape}")
                                            print(f"        tensor shape: {tensor.shape}")
                                            raise
                                collated[key][target_key] = padded
                            elif target_key == 'basepair_list':
                                # Special handling for basepair_list which can be empty
                                print(f"      Processing basepair_list")
                                max_basepairs = max(tensor.shape[0] for tensor in tensors)
                                if max_basepairs > 0:
                                    padded = torch.zeros(len(samples), max_basepairs, 2, dtype=tensors[0].dtype)
                                    for i, tensor in enumerate(tensors):
                                        if tensor.shape[0] > 0:
                                            padded[i, :tensor.shape[0]] = tensor
                                    collated[key][target_key] = padded
                                else:
                                    # All basepair_list are empty, create empty tensor
                                    collated[key][target_key] = torch.zeros(len(samples), 0, 2, dtype=tensors[0].dtype)
                            else:
                                # Other flattened tensors - check if they can be stacked
                                print(f"      Other flattened tensor: {target_key}")
                                try:
                                    collated[key][target_key] = torch.stack(tensors)
                                    print(f"      ✓ Stacked successfully")
                                except RuntimeError as e:
                                    print(f"      ❌ Stacking failed: {e}")
                                    # If stacking fails, pad to same length
                                    max_len_tensor = max(tensor.shape[0] for tensor in tensors)
                                    if tensors[0].dim() == 1:
                                        padded = torch.zeros(len(samples), max_len_tensor, dtype=tensors[0].dtype)
                                        for i, tensor in enumerate(tensors):
                                            if len(tensor) > 0:
                                                padded[i, :len(tensor)] = tensor
                                        collated[key][target_key] = padded
                                        print(f"      ✓ Padded successfully")
                                    else:
                                        # Fallback: just stack what we can
                                        collated[key][target_key] = torch.stack(tensors)
                    elif tensors[0].dim() == 2:
                        # 2D tensors - check if it's square (pairwise) or rectangular
                        if tensors[0].shape[0] == tensors[0].shape[1]:
                            # Square matrix (pairwise)
                            print(f"    Square matrix (pairwise)")
                            padded = torch.zeros(len(samples), max_len, max_len, dtype=tensors[0].dtype)
                            for i, tensor in enumerate(tensors):
                                padded[i, :tensor.shape[0], :tensor.shape[1]] = tensor
                            collated[key][target_key] = padded
                        else:
                            # Rectangular matrix - need to check if it's atom_coords related
                            print(f"    Rectangular matrix: {tensors[0].shape}")
                            if target_key in ['atom_coords', 'atom_mask']:
                                # These are flattened tensors that need padding
                                print(f"      Processing {target_key} as flattened tensor")
                                max_atoms = max(tensor.shape[0] for tensor in tensors)
                                if target_key == 'atom_coords':
                                    padded = torch.zeros(len(samples), max_atoms, 3, dtype=tensors[0].dtype)
                                else:
                                    padded = torch.zeros(len(samples), max_atoms, dtype=tensors[0].dtype)
                                for i, tensor in enumerate(tensors):
                                    padded[i, :tensor.shape[0]] = tensor
                                collated[key][target_key] = padded
                            else:
                                # Other rectangular matrices - just stack
                                print(f"      Other rectangular matrix: {target_key}")
                                collated[key][target_key] = torch.stack(tensors)
                    else:
                        print(f"    Higher dimensional tensor: {tensors[0].dim()}D")
                        collated[key][target_key] = torch.stack(tensors)
        else:
            print(f"  Non-target key: {key}")
            # Handle other keys
            collated[key] = [item[key] for item in samples]
    
    print("\n✓ Collate function completed successfully!")
    return collated


def main():
    """Run the comprehensive debug."""
    print("Comprehensive Collate Debug Test")
    print("=" * 40)
    
    try:
        result = debug_all_tensors()
        print(f"\nFinal result keys: {list(result.keys())}")
    except Exception as e:
        print(f"❌ Debug failed: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()
