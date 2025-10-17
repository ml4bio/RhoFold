"""
Test script for MSA handling logic with three scenarios:
1. MSA and seq length different - use alignment
2. Missing MSA data - skip MSA module
3. Normal data - process normally
"""

import os
import sys
import torch
import numpy as np
from pathlib import Path

# Add parent directory to path
sys.path.append(str(Path(__file__).parent.parent))

from train.data_loader import create_dataloader


def test_msa_handling():
    """Test MSA handling with different scenarios."""
    print("Testing MSA handling logic...")
    
    # Test with real data
    data_dir = "/home/rcy/code/RhoFold/RNA3D_DATA"
    index_file = "fold-0_train_ids"
    
    # Create dataloader with small batch size for testing
    dataloader = create_dataloader(
        data_dir=data_dir,
        index_file=index_file,
        batch_size=2,
        num_workers=0,
        max_seq_len=100,
        msa_depth=10,
        shuffle=False,
    )
    
    print(f"Dataset length: {len(dataloader.dataset)}")
    
    # Test loading a few batches
    for batch_idx, batch in enumerate(dataloader):
        print(f"\n--- Batch {batch_idx + 1} ---")
        
        # Check MSA validity
        valid_msa = batch.get('valid_msa', None)
        if valid_msa is not None:
            print(f"MSA validity: {valid_msa.tolist()}")
            print(f"Valid MSA samples: {valid_msa.sum().item()}/{len(valid_msa)}")
        else:
            print("No MSA validity information")
        
        # Check tensor shapes
        print(f"Tokens shape: {batch['tokens'].shape}")
        print(f"MSA tokens shape: {batch['msa_tokens'].shape}")
        print(f"Sequence lengths: {[len(seq) for seq in batch['sequence']]}")
        
        # Check if MSA tokens are zero-filled for invalid samples
        if valid_msa is not None:
            for i, is_valid in enumerate(valid_msa.tolist()):
                msa_tokens_sample = batch['msa_tokens'][i]
                if not is_valid:
                    # Should be mostly zeros for invalid MSA
                    non_zero_ratio = (msa_tokens_sample != 0).float().mean().item()
                    print(f"  Sample {i} (invalid MSA): non-zero ratio = {non_zero_ratio:.3f}")
                else:
                    # Should have some non-zero values for valid MSA
                    non_zero_ratio = (msa_tokens_sample != 0).float().mean().item()
                    print(f"  Sample {i} (valid MSA): non-zero ratio = {non_zero_ratio:.3f}")
        
        # Only test first few batches
        if batch_idx >= 2:
            break
    
    print("\n✓ MSA handling test completed!")


def test_alignment_function():
    """Test the MSA alignment function directly."""
    print("\nTesting MSA alignment function...")
    
    from train.utils_a3m import read_a3m, align_msa_to_seq
    
    # Test with a simple case
    test_seq = "AUCG" * 5  # 20 bases
    test_msa = [
        ("seq1", "AUCG" * 5),  # Same length
        ("seq2", "AUCG" * 4 + "A"),  # Different length
        ("seq3", "AUCG" * 6),  # Longer
    ]
    
    # Test alignment
    try:
        aligned = align_msa_to_seq(test_seq, test_msa)
        if aligned is not None:
            print(f"✓ Alignment successful: {aligned.shape}")
            print(f"  Target length: {len(test_seq)}")
            print(f"  Aligned shape: {aligned.shape}")
        else:
            print("❌ Alignment failed")
    except Exception as e:
        print(f"❌ Alignment error: {e}")


def main():
    """Run all tests."""
    print("MSA Handling Test Suite")
    print("=" * 50)
    
    # Test alignment function
    test_alignment_function()
    
    # Test data loading with MSA handling
    test_msa_handling()
    
    print("\n🎉 All tests completed!")


if __name__ == "__main__":
    main()
