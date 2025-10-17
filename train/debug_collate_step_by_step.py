"""
Step-by-step debug script to identify the exact problematic tensor.
"""

import os
import sys
import torch
from pathlib import Path

# Add parent directory to path
sys.path.append(str(Path(__file__).parent.parent))

from train.data_loader import RhoFoldDataset


def debug_collate_step_by_step():
    """Debug collate function step by step."""
    print("Step-by-step collate debugging...")
    
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
    
    # Manually implement collate logic to debug
    print("\nManually implementing collate logic...")
    
    # Get max sequence length
    max_len = max(item['tokens'].shape[1] for item in samples)
    print(f"Max sequence length: {max_len}")
    
    # Process target dict
    target_key = 'masked_tokens'
    tensors = [item['target'][target_key] for item in samples]
    print(f"\nProcessing {target_key}:")
    for i, tensor in enumerate(tensors):
        print(f"  Tensor {i}: shape={tensor.shape}, len={len(tensor)}")
    
    # Check if this is the problematic case
    if target_key in ['masked_tokens', 'plddt_true', 'res_mask']:
        print(f"  This is a 1D tensor that needs padding to max sequence length")
        padded = torch.zeros(len(samples), max_len, dtype=tensors[0].dtype)
        print(f"  Created padded tensor: {padded.shape}")
        
        for i, tensor in enumerate(tensors):
            print(f"  Processing tensor {i}: shape={tensor.shape}, len={len(tensor)}")
            if len(tensor) > 0:
                print(f"    Attempting: padded[{i}, :{len(tensor)}] = tensor")
                try:
                    padded[i, :len(tensor)] = tensor
                    print(f"    ✓ Success")
                except Exception as e:
                    print(f"    ❌ Error: {e}")
                    print(f"    padded[{i}, :{len(tensor)}] shape: {padded[i, :len(tensor)].shape}")
                    print(f"    tensor shape: {tensor.shape}")
                    break
            else:
                print(f"    Skipping empty tensor")


def main():
    """Run the step-by-step debug."""
    print("Step-by-Step Collate Debug Test")
    print("=" * 40)
    
    try:
        debug_collate_step_by_step()
    except Exception as e:
        print(f"❌ Debug failed: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()
