"""
使用CPU运行训练来获取详细错误信息
"""

import os
import sys
import torch
from pathlib import Path

# Add parent directory to path
sys.path.append(str(Path(__file__).parent.parent))

from train.data_loader import create_dataloader
from train.train_diffusion_final import DiffusionTrainer


def test_cpu_training():
    """在CPU上运行训练来获取详细错误信息"""
    print("🔍 CPU训练测试...")
    
    # 创建数据加载器
    train_dataloader = create_dataloader(
        data_dir="/home/rcy/code/RhoFold/RNA3D_DATA",
        index_file="fold-0_train_ids",
        batch_size=1,  # 使用小batch
        num_workers=0,
        max_seq_len=30,  # 使用更小的序列长度
        msa_depth=3,
        shuffle=False,
    )
    
    val_dataloader = create_dataloader(
        data_dir="/home/rcy/code/RhoFold/RNA3D_DATA",
        index_file="fold-0_train_ids",
        batch_size=1,
        num_workers=0,
        max_seq_len=30,
        msa_depth=3,
        shuffle=False,
    )
    
    # 创建训练器（使用CPU）
    config = {
        'diffusion': {
            'num_timesteps': 1000,
            'noise_schedule': 'linear',
            'time_embed_dim': 256,
            'base_channels': 64,
        },
        'learning_rate': 1e-4,
        'weight_decay': 1e-5,
        'epochs': 1,
    }
    
    trainer = DiffusionTrainer(config)
    
    print("✓ 训练器创建成功")
    
    # 尝试运行一个训练步骤
    try:
        print("开始训练测试...")
        trainer.train(train_dataloader, val_dataloader, "./checkpoints")
    except Exception as e:
        print(f"❌ 训练失败: {e}")
        import traceback
        traceback.print_exc()


def main():
    """运行CPU训练测试"""
    print("🔍 CPU训练测试")
    print("=" * 50)
    
    try:
        test_cpu_training()
    except Exception as e:
        print(f"❌ 测试失败: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()