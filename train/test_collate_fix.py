"""
Test script to verify the collate_fn fix for different sequence lengths.
"""

import os
import sys
import torch
from pathlib import Path

# Add parent directory to path
sys.path.append(str(Path(__file__).parent.parent))

from train.data_loader import create_dataloader


def test_collate_fix():
    """Test collate_fn with different sequence lengths."""
    print("Testing collate_fn fix...")
    
    # Test with real data
    data_dir = "/home/rcy/code/RhoFold/RNA3D_DATA"
    index_file = "fold-0_train_ids"
    
    # Create dataloader with small batch size for testing
    dataloader = create_dataloader(
        data_dir=data_dir,
        index_file=index_file,
        batch_size=4,
        num_workers=0,  # Use single worker to avoid multiprocessing issues
        max_seq_len=200,  # Allow longer sequences
        msa_depth=10,
        shuffle=False,
    )
    
    print(f"Dataset length: {len(dataloader.dataset)}")
    
    # Test loading a few batches
    for batch_idx, batch in enumerate(dataloader):
        print(f"\n--- Batch {batch_idx + 1} ---")
        
        # Check tensor shapes
        print(f"Tokens shape: {batch['tokens'].shape}")
        print(f"MSA tokens shape: {batch['msa_tokens'].shape}")
        print(f"Atom coords shape: {batch['atom_coords'].shape}")
        print(f"Atom types shape: {batch['atom_types'].shape}")
        print(f"Atom mask shape: {batch['atom_mask'].shape}")
        
        # Check sequence lengths
        seq_lengths = [len(seq) for seq in batch['sequence']]
        print(f"Sequence lengths: {seq_lengths}")
        
        # Check MSA validity
        if 'valid_msa' in batch:
            valid_msa = batch['valid_msa']
            print(f"MSA validity: {valid_msa.tolist()}")
            print(f"Valid MSA samples: {valid_msa.sum().item()}/{len(valid_msa)}")
        
        # Only test first few batches
        if batch_idx >= 2:
            break
    
    print("\n✓ Collate fix test completed!")


def main():
    """Run the test."""
    print("Collate Function Fix Test")
    print("=" * 40)
    
    try:
        test_collate_fix()
        print("\n🎉 Test passed! Collate function is working correctly.")
    except Exception as e:
        print(f"\n❌ Test failed: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()
