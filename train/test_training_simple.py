"""
Simple training test with minimal configuration.
"""

import os
import sys
import torch
from pathlib import Path

# Add parent directory to path
sys.path.append(str(Path(__file__).parent.parent))

from train.data_loader import create_dataloader
from rhofold.config import rhofold_config
from rhofold.rhofold_with_diffusion import create_diffusion_rhofold


def test_training_simple():
    """Test simple training with minimal configuration."""
    print("Testing simple training...")
    
    # Create data loader with single worker
    dataloader = create_dataloader(
        data_dir="/home/rcy/code/RhoFold/RNA3D_DATA",
        index_file="fold-0_train_ids",
        batch_size=2,
        num_workers=0,  # Single worker to avoid multiprocessing issues
        max_seq_len=50,  # Small sequence length
        msa_depth=5,     # Small MSA depth
        shuffle=False,
    )
    
    print(f"Dataset length: {len(dataloader.dataset)}")
    
    # Create model
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Using device: {device}")
    
    diffusion_config = {
        'time_embed_dim': 64,
        'num_timesteps': 10,
        'noise_schedule': 'linear',
        'base_channels': 16,
    }
    
    model = create_diffusion_rhofold(rhofold_config, diffusion_config)
    model.to(device)
    
    print(f"Model created and moved to {device}")
    
    # Test a few batches
    for batch_idx, batch in enumerate(dataloader):
        print(f"\n--- Batch {batch_idx + 1} ---")
        
        # Move batch to device
        for key in batch:
            if isinstance(batch[key], torch.Tensor):
                batch[key] = batch[key].to(device)
            elif isinstance(batch[key], dict):
                for subkey in batch[key]:
                    if isinstance(batch[key][subkey], torch.Tensor):
                        batch[key][subkey] = batch[key][subkey].to(device)
                    elif isinstance(batch[key][subkey], dict):
                        for subsubkey in batch[key][subkey]:
                            if isinstance(batch[key][subkey][subsubkey], torch.Tensor):
                                batch[key][subkey][subsubkey] = batch[key][subkey][subsubkey].to(device)
        
        # Check tensor shapes and devices
        print(f"Tokens shape: {batch['tokens'].shape}, device: {batch['tokens'].device}")
        print(f"MSA tokens shape: {batch['msa_tokens'].shape}, device: {batch['msa_tokens'].device}")
        print(f"Valid MSA: {batch.get('valid_msa', 'Not found')}")
        
        # Convert inputs for model
        rna_fm_tokens = batch['msa_tokens'][:, 0, :]
        seq_str = batch['sequence'][0] if isinstance(batch['sequence'], list) else batch['sequence']
        valid_msa = batch.get('valid_msa', torch.ones(batch['tokens'].shape[0], dtype=torch.bool))
        
        try:
            # Forward pass
            with torch.no_grad():
                outputs = model(
                    tokens=batch['tokens'],
                    rna_fm_tokens=rna_fm_tokens,
                    seq=seq_str,
                    valid_msa=valid_msa
                )
            
            print(f"✓ Forward pass successful")
            print(f"  Number of outputs: {len(outputs)}")
            if outputs:
                output = outputs[-1]
                print(f"  Output keys: {list(output.keys())}")
                for key, value in output.items():
                    if isinstance(value, torch.Tensor):
                        print(f"    {key}: {value.shape}, device: {value.device}")
                    elif isinstance(value, (list, tuple)):
                        print(f"    {key}: {type(value)} with {len(value)} elements")
                        if value and isinstance(value[0], torch.Tensor):
                            print(f"      First element shape: {value[0].shape}, device: {value[0].device}")
            
        except Exception as e:
            print(f"❌ Forward pass failed: {e}")
            import traceback
            traceback.print_exc()
            break
        
        # Only test first few batches
        if batch_idx >= 2:
            break
    
    print("\n✓ Simple training test completed!")


def main():
    """Run the test."""
    print("Simple Training Test")
    print("=" * 30)
    
    try:
        test_training_simple()
        print("\n🎉 Test passed! Training setup is working correctly.")
    except Exception as e:
        print(f"\n❌ Test failed: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()
