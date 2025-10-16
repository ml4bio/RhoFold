"""
Debug script to understand data loading issues.
"""

import os
import sys
import torch
import numpy as np
from pathlib import Path

# Add parent directory to path
sys.path.append(str(Path(__file__).parent.parent))

from train.data_loader import RhoFoldDataset


def debug_data_loading():
    """Debug data loading step by step."""
    print("Debugging data loading...")
    
    # Use real data directory
    data_dir = "/home/rcy/code/RhoFold/example/input/diffusion_test"
    index_file = os.path.join(data_dir, "index.txt")
    
    # Create index file for the real data
    with open(index_file, "w") as f:
        f.write("1a1t_B\n")
    
    # Create dataset
    dataset = RhoFoldDataset(
        data_dir=data_dir,
        index_file=index_file,
        max_seq_len=50,
        msa_depth=3,
    )
    
    print(f"Dataset length: {len(dataset)}")
    
    # Load a single sample
    try:
        sample = dataset[0]
        print("✓ Single sample loaded successfully")
        print(f"Sample keys: {list(sample.keys())}")
        
        for key, value in sample.items():
            if isinstance(value, torch.Tensor):
                print(f"  {key}: {value.shape} ({value.dtype})")
            elif isinstance(value, dict):
                print(f"  {key}: dict with keys {list(value.keys())}")
                for subkey, subvalue in value.items():
                    if isinstance(subvalue, torch.Tensor):
                        print(f"    {subkey}: {subvalue.shape} ({subvalue.dtype})")
            else:
                print(f"  {key}: {type(value)} = {value}")
        
        # Test collate function
        from train.data_loader import collate_fn
        batch = collate_fn([sample])
        print("\n✓ Collate function successful")
        print(f"Batch keys: {list(batch.keys())}")
        
        for key, value in batch.items():
            if isinstance(value, torch.Tensor):
                print(f"  {key}: {value.shape} ({value.dtype})")
            elif isinstance(value, dict):
                print(f"  {key}: dict with keys {list(value.keys())}")
                for subkey, subvalue in value.items():
                    if isinstance(subvalue, torch.Tensor):
                        print(f"    {subkey}: {subvalue.shape} ({subvalue.dtype})")
            else:
                print(f"  {key}: {type(value)}")
        
        return True
        
    except Exception as e:
        print(f"✗ Error loading sample: {e}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == "__main__":
    success = debug_data_loading()
    if success:
        print("\n🎉 Data loading debug successful!")
    else:
        print("\n❌ Data loading debug failed.")
    sys.exit(0 if success else 1)
