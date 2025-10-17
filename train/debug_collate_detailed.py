"""
Detailed debug script to identify the exact problematic tensor.
"""

import os
import sys
import torch
from pathlib import Path

# Add parent directory to path
sys.path.append(str(Path(__file__).parent.parent))

from train.data_loader import RhoFoldDataset, collate_fn


def debug_collate_detailed():
    """Debug collate function with detailed logging."""
    print("Detailed collate debugging...")
    
    # Create dataset directly
    dataset = RhoFoldDataset(
        data_dir="/home/rcy/code/RhoFold/RNA3D_DATA",
        index_file="fold-0_train_ids",
        max_seq_len=50,
        msa_depth=5,
    )
    
    print(f"Dataset length: {len(dataset)}")
    
    # Load a few samples manually
    samples = []
    for i in range(2):
        try:
            sample = dataset[i]
            samples.append(sample)
            print(f"Sample {i}: {sample['sample_name']}")
            print(f"  Sequence length: {len(sample['sequence'])}")
            print(f"  Tokens shape: {sample['tokens'].shape}")
            print(f"  MSA tokens shape: {sample['msa_tokens'].shape}")
            print(f"  Valid MSA: {sample.get('valid_msa', 'Not found')}")
            
            # Check target tensors
            if 'target' in sample:
                for key, value in sample['target'].items():
                    if isinstance(value, torch.Tensor):
                        print(f"  Target {key}: {value.shape} ({value.dtype})")
                    elif isinstance(value, dict):
                        print(f"  Target {key}: dict with keys {list(value.keys())}")
                        for subkey, subvalue in value.items():
                            if isinstance(subvalue, torch.Tensor):
                                print(f"    {subkey}: {subvalue.shape} ({subvalue.dtype})")
            
        except Exception as e:
            print(f"Error loading sample {i}: {e}")
            continue
    
    if len(samples) < 2:
        print("Not enough samples loaded for collate test")
        return
    
    # Test collate function
    print("\nTesting collate function...")
    try:
        collated = collate_fn(samples)
        print("✓ Collate function successful!")
        
        # Print collated shapes
        for key, value in collated.items():
            if isinstance(value, torch.Tensor):
                print(f"Collated {key}: {value.shape} ({value.dtype})")
            elif isinstance(value, dict):
                print(f"Collated {key}: dict with keys {list(value.keys())}")
                for subkey, subvalue in value.items():
                    if isinstance(subvalue, torch.Tensor):
                        print(f"  {subkey}: {subvalue.shape} ({subvalue.dtype})")
        
    except Exception as e:
        print(f"❌ Collate function failed: {e}")
        import traceback
        traceback.print_exc()


def main():
    """Run the detailed debug."""
    print("Detailed Collate Debug Test")
    print("=" * 40)
    
    try:
        debug_collate_detailed()
    except Exception as e:
        print(f"❌ Debug failed: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()
