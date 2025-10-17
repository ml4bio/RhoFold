"""
Simple test for collate function with minimal data.
"""

import os
import sys
import torch
from pathlib import Path

# Add parent directory to path
sys.path.append(str(Path(__file__).parent.parent))

from train.data_loader import create_dataloader


def test_simple_collate():
    """Test collate function with very simple data."""
    print("Testing simple collate function...")
    
    # Create data loader with very small batch and sequence length
    dataloader = create_dataloader(
        data_dir="/home/rcy/code/RhoFold/RNA3D_DATA",
        index_file="fold-0_train_ids",
        batch_size=1,  # Single sample
        num_workers=0,
        max_seq_len=20,  # Very small sequence length
        msa_depth=3,     # Very small MSA depth
        shuffle=False,
    )
    
    print(f"Dataset length: {len(dataloader.dataset)}")
    
    # Test loading a single batch
    try:
        batch = next(iter(dataloader))
        print("✓ Single batch loaded successfully!")
        
        # Print tensor shapes
        for key, value in batch.items():
            if isinstance(value, torch.Tensor):
                print(f"{key}: {value.shape} ({value.dtype})")
            elif isinstance(value, dict):
                print(f"{key}: dict with keys {list(value.keys())}")
                for subkey, subvalue in value.items():
                    if isinstance(subvalue, torch.Tensor):
                        print(f"  {subkey}: {subvalue.shape} ({subvalue.dtype})")
        
    except Exception as e:
        print(f"❌ Error loading batch: {e}")
        import traceback
        traceback.print_exc()


def main():
    """Run the test."""
    print("Simple Collate Test")
    print("=" * 30)
    
    try:
        test_simple_collate()
        print("\n🎉 Test passed!")
    except Exception as e:
        print(f"\n❌ Test failed: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()
