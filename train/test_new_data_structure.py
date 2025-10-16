"""
Test script for the new data structure with subdirectories.
"""

import os
import sys
import torch
import numpy as np
from pathlib import Path

# Add parent directory to path
sys.path.append(str(Path(__file__).parent.parent))

from train.data_loader import create_dataloader


def test_new_data_structure():
    """Test data loading with the new directory structure."""
    print("Testing new data structure...")
    
    # Check if the data directory exists
    data_dir = "/home/rcy/code/RhoFold/RNA3D_DATA"
    if not os.path.exists(data_dir):
        print(f"❌ Data directory not found: {data_dir}")
        return False
    
    # Check subdirectories
    subdirs = ["list", "pdb", "rMSA", "seq"]
    for subdir in subdirs:
        subdir_path = os.path.join(data_dir, subdir)
        if not os.path.exists(subdir_path):
            print(f"❌ Subdirectory not found: {subdir_path}")
            return False
        else:
            print(f"✓ Found subdirectory: {subdir_path}")
    
    # List files in each subdirectory
    for subdir in subdirs:
        subdir_path = os.path.join(data_dir, subdir)
        files = os.listdir(subdir_path)
        print(f"  {subdir}/ contains {len(files)} files")
        if len(files) > 0:
            print(f"    Sample files: {files[:3]}...")  # Show first 3 files
    
    # Look for index files in list/ directory
    list_dir = os.path.join(data_dir, "list")
    index_files = [f for f in os.listdir(list_dir) if not f.endswith('.txt')]
    print(f"\nIndex files found in list/: {index_files}")
    
    if not index_files:
        print("❌ No index files found in list/ directory")
        return False
    
    # Use the first index file for testing
    index_file = index_files[0]
    print(f"Using index file: {index_file}")
    
    # Test data loading
    try:
        dataloader = create_dataloader(
            data_dir=data_dir,
            index_file=index_file,
            batch_size=1,
            num_workers=0,
            max_seq_len=100,  # Small for testing
            msa_depth=10,     # Small for testing
            shuffle=False,
        )
        
        print(f"✓ DataLoader created successfully")
        print(f"  Dataset length: {len(dataloader.dataset)}")
        
        # Test loading a single batch
        batch = next(iter(dataloader))
        print(f"✓ Successfully loaded batch")
        print(f"  Batch keys: {list(batch.keys())}")
        
        # Print tensor shapes
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
        print(f"❌ Data loading failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    """Run the test."""
    print("New Data Structure Test")
    print("=" * 50)
    
    success = test_new_data_structure()
    
    if success:
        print("\n🎉 New data structure test passed!")
        print("You can now use the following command to train:")
        print(f"python train/train_diffusion_final.py \\")
        print(f"    --data_dir /home/rcy/code/RhoFold/RNA3D_DATA \\")
        print(f"    --index_file <your_index_file> \\")
        print(f"    --batch_size 4 \\")
        print(f"    --epochs 100")
    else:
        print("\n❌ New data structure test failed.")
        print("Please check your data directory structure.")
    
    return success


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
