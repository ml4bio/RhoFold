"""
Debug script to identify the problematic tensor in collate_fn.
"""

import os
import sys
import torch
from pathlib import Path

# Add parent directory to path
sys.path.append(str(Path(__file__).parent.parent))

from train.data_loader import create_dataloader


def debug_collate():
    """Debug collate function to find problematic tensors."""
    print("Debugging collate function...")
    
    # Create data loader with very small batch
    dataloader = create_dataloader(
        data_dir="/home/rcy/code/RhoFold/RNA3D_DATA",
        index_file="fold-0_train_ids",
        batch_size=2,
        num_workers=0,
        max_seq_len=50,
        msa_depth=5,
        shuffle=False,
    )
    
    print(f"Dataset length: {len(dataloader.dataset)}")
    
    # Test loading a single batch
    try:
        batch = next(iter(dataloader))
        print("✓ Batch loaded successfully!")
        
        # Print all tensor shapes
        for key, value in batch.items():
            if isinstance(value, torch.Tensor):
                print(f"{key}: {value.shape} ({value.dtype})")
            elif isinstance(value, dict):
                print(f"{key}: dict with keys {list(value.keys())}")
                for subkey, subvalue in value.items():
                    if isinstance(subvalue, torch.Tensor):
                        print(f"  {subkey}: {subvalue.shape} ({subvalue.dtype})")
                    elif isinstance(subvalue, dict):
                        print(f"  {subkey}: dict with keys {list(subvalue.keys())}")
                        for subsubkey, subsubvalue in subvalue.items():
                            if isinstance(subsubvalue, torch.Tensor):
                                print(f"    {subsubkey}: {subsubvalue.shape} ({subsubvalue.dtype})")
            else:
                print(f"{key}: {type(value)}")
        
    except Exception as e:
        print(f"❌ Error loading batch: {e}")
        import traceback
        traceback.print_exc()


def main():
    """Run the debug."""
    print("Collate Debug Test")
    print("=" * 30)
    
    try:
        debug_collate()
    except Exception as e:
        print(f"❌ Debug failed: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()
