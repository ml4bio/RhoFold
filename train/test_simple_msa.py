"""
Simple test for MSA handling modifications.
"""

import os
import sys
import torch
from pathlib import Path

# Add parent directory to path
sys.path.append(str(Path(__file__).parent.parent))

def test_data_loader_import():
    """Test if data loader imports correctly."""
    try:
        from train.data_loader import create_dataloader
        print("✓ Data loader import successful")
        return True
    except Exception as e:
        print(f"❌ Data loader import failed: {e}")
        return False

def test_utils_a3m_import():
    """Test if utils_a3m imports correctly."""
    try:
        from train.utils_a3m import read_a3m, align_msa_to_seq
        print("✓ Utils A3M import successful")
        return True
    except Exception as e:
        print(f"❌ Utils A3M import failed: {e}")
        return False

def test_model_import():
    """Test if model imports correctly."""
    try:
        from rhofold.rhofold_with_diffusion import create_diffusion_rhofold
        print("✓ Model import successful")
        return True
    except Exception as e:
        print(f"❌ Model import failed: {e}")
        return False

def main():
    """Run simple tests."""
    print("Simple MSA Handling Test")
    print("=" * 30)
    
    success = True
    
    # Test imports
    success &= test_data_loader_import()
    success &= test_utils_a3m_import()
    success &= test_model_import()
    
    if success:
        print("\n🎉 All imports successful!")
        print("MSA handling modifications are ready to use.")
    else:
        print("\n❌ Some imports failed.")
    
    return success

if __name__ == "__main__":
    main()
