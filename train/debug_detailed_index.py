"""
详细检查所有可能的索引问题
"""

import os
import sys
import torch
from pathlib import Path

# Add parent directory to path
sys.path.append(str(Path(__file__).parent.parent))

from train.data_loader import RhoFoldDataset, collate_fn
from torch.utils.data import DataLoader


def check_all_indices():
    """检查所有可能的索引问题"""
    print("🔍 详细检查所有索引...")
    
    # 创建dataset
    dataset = RhoFoldDataset(
        data_dir="/home/rcy/code/RhoFold/RNA3D_DATA",
        index_file="fold-0_train_ids",
        max_seq_len=50,
        msa_depth=5,
    )
    
    # 检查前10个样本
    for i in range(min(10, len(dataset))):
        try:
            sample = dataset[i]
            print(f"\n--- Sample {i} ---")
            
            # 检查每个tensor的索引范围
            for key, value in sample.items():
                if isinstance(value, torch.Tensor) and value.dtype in (torch.int32, torch.int64):
                    min_val = value.min().item()
                    max_val = value.max().item()
                    print(f"{key}: min={min_val}, max={max_val}, shape={value.shape}")
                    
                    # 检查是否有负值
                    if min_val < 0:
                        print(f"  ❌ {key}: contains negative values!")
                    
                    # 检查是否有可能越界的值
                    if key == 'tokens' or key == 'msa_tokens':
                        if max_val >= 10:  # MSA vocab size
                            print(f"  ❌ {key}: max value {max_val} >= MSA vocab size 10")
                    elif key == 'rna_fm_tokens':
                        if max_val >= 30:  # RNA-FM vocab size
                            print(f"  ❌ {key}: max value {max_val} >= RNA-FM vocab size 30")
                    elif key in ['atom_types']:
                        if max_val >= 23:  # Max atom types
                            print(f"  ❌ {key}: max value {max_val} >= max atom types 23")
                    elif key in ['distance_bins', 'ss_true']:
                        if max_val >= 64:  # Distance/SS bins
                            print(f"  ❌ {key}: max value {max_val} >= max bins 64")
                    elif key in ['plddt_true']:
                        if max_val > 1.0:  # pLDDT should be 0-1
                            print(f"  ❌ {key}: max value {max_val} > 1.0")
                    elif key in ['res_mask', 'pair_mask', 'atom_mask']:
                        if max_val > 1:  # Masks should be 0-1
                            print(f"  ❌ {key}: max value {max_val} > 1")
                            
        except Exception as e:
            print(f"❌ Sample {i} failed: {e}")
            import traceback
            traceback.print_exc()


def check_model_inputs():
    """检查模型输入是否有问题"""
    print("\n🔍 检查模型输入...")
    
    # 创建dataloader
    dataset = RhoFoldDataset(
        data_dir="/home/rcy/code/RhoFold/RNA3D_DATA",
        index_file="fold-0_train_ids",
        max_seq_len=50,
        msa_depth=5,
    )
    
    loader = DataLoader(
        dataset, 
        batch_size=2, 
        collate_fn=collate_fn,
        num_workers=0,
        shuffle=False
    )
    
    # 检查第一个batch
    try:
        batch = next(iter(loader))
        print("✓ Batch loaded successfully")
        
        # 检查关键输入
        for key in ['tokens', 'msa_tokens', 'rna_fm_tokens']:
            if key in batch:
                tensor = batch[key]
                print(f"{key}: shape={tensor.shape}, dtype={tensor.dtype}")
                if tensor.dtype in (torch.int32, torch.int64):
                    min_val = tensor.min().item()
                    max_val = tensor.max().item()
                    print(f"  min={min_val}, max={max_val}")
                    
                    # 检查索引范围
                    if key in ['tokens', 'msa_tokens']:
                        if max_val >= 10:
                            print(f"  ❌ {key}: max value {max_val} >= MSA vocab size 10")
                    elif key == 'rna_fm_tokens':
                        if max_val >= 30:
                            print(f"  ❌ {key}: max value {max_val} >= RNA-FM vocab size 30")
                            
    except Exception as e:
        print(f"❌ Batch loading failed: {e}")
        import traceback
        traceback.print_exc()


def main():
    """运行详细检查"""
    print("🔍 详细索引检查")
    print("=" * 50)
    
    try:
        check_all_indices()
        check_model_inputs()
        
        print("\n" + "=" * 50)
        print("✅ 详细检查完成！")
        
    except Exception as e:
        print(f"❌ 检查失败: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()
