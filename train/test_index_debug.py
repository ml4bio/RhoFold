"""
使用debug_index_check.py来定位CUDA错误
"""

import os
import sys
import torch
from pathlib import Path

# Add parent directory to path
sys.path.append(str(Path(__file__).parent.parent))

from train.data_loader import RhoFoldDataset, collate_fn
from torch.utils.data import DataLoader
from train.debug_index_check import scan_dataset, scan_dataloader


def test_dataset_scan():
    """方式一：检查单个Dataset"""
    print("=" * 60)
    print("方式一：检查单个 Dataset")
    print("=" * 60)
    
    # 创建dataset
    dataset = RhoFoldDataset(
        data_dir="/home/rcy/code/RhoFold/RNA3D_DATA",
        index_file="fold-0_train_ids",
        max_seq_len=50,
        msa_depth=5,
    )
    
    print(f"Dataset length: {len(dataset)}")
    
    # 扫描前100个样本
    scan_dataset(dataset, vocab_size_seq=30, vocab_size_msa=17, num_samples=100)


def test_dataloader_scan():
    """方式二：检查DataLoader（含padding后的张量）"""
    print("\n" + "=" * 60)
    print("方式二：检查 DataLoader（含 padding 后的张量）")
    print("=" * 60)
    
    # 创建dataset和dataloader
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
        num_workers=0,  # 使用单进程避免多进程问题
        shuffle=False
    )
    
    print(f"Dataset length: {len(dataset)}")
    
    # 扫描前5个batch
    scan_dataloader(loader, vocab_size_seq=30, vocab_size_msa=17, num_batches=5)


def main():
    """运行两种检查方式"""
    print("🔍 使用debug_index_check.py定位CUDA错误")
    print("=" * 60)
    
    try:
        # 方式一：检查单个Dataset
        test_dataset_scan()
        
        # 方式二：检查DataLoader
        test_dataloader_scan()
        
        print("\n" + "=" * 60)
        print("✅ 所有检查完成！")
        print("请查看上面的错误信息来定位问题所在。")
        
    except Exception as e:
        print(f"❌ 检查过程中出现错误: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()
