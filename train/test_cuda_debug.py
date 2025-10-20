"""
使用CUDA_LAUNCH_BLOCKING=1进行调试
"""

import os
import sys
import torch
from pathlib import Path

# Add parent directory to path
sys.path.append(str(Path(__file__).parent.parent))

from train.data_loader import create_dataloader
from train.train_diffusion_final import DiffusionTrainer


def test_cuda_debug():
    """使用CUDA_LAUNCH_BLOCKING=1进行调试"""
    print("🔍 CUDA调试测试...")
    
    # 创建数据加载器 - 使用更小的参数
    train_dataloader = create_dataloader(
        data_dir="/home/rcy/code/RhoFold/RNA3D_DATA",
        index_file="fold-0_train_ids",
        batch_size=1,  # 最小batch size
        num_workers=0,
        max_seq_len=20,  # 更小的序列长度
        msa_depth=3,     # 更小的MSA深度
        shuffle=False,
    )
    
    val_dataloader = create_dataloader(
        data_dir="/home/rcy/code/RhoFold/RNA3D_DATA",
        index_file="fold-0_train_ids",
        batch_size=1,
        num_workers=0,
        max_seq_len=20,
        msa_depth=3,
        shuffle=False,
    )
    
    # 创建训练器
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
    
    # 尝试运行几个训练步骤
    try:
        print("开始CUDA调试训练...")
        print("注意：使用CUDA_LAUNCH_BLOCKING=1，错误信息会更详细")
        
        # 只运行几个batch来测试
        for epoch in range(1):
            print(f"\nEpoch {epoch}:")
            for batch_idx, batch in enumerate(train_dataloader):
                if batch_idx >= 5:  # 只测试前5个batch
                    break
                    
                print(f"  Batch {batch_idx}:")
                try:
                    # 打印batch中的字段
                    print(f"    Batch keys: {list(batch.keys())}")
                    
                    # 移动数据到设备
                    for key in batch:
                        if isinstance(batch[key], torch.Tensor):
                            batch[key] = batch[key].to(trainer.device)
                        elif isinstance(batch[key], dict):
                            for subkey in batch[key]:
                                if isinstance(batch[key][subkey], torch.Tensor):
                                    batch[key][subkey] = batch[key][subkey].to(trainer.device)
                                elif isinstance(batch[key][subkey], dict):
                                    for subsubkey in batch[key][subkey]:
                                        if isinstance(batch[key][subkey][subsubkey], torch.Tensor):
                                            batch[key][subkey][subsubkey] = batch[key][subkey][subsubkey].to(trainer.device)
                    
                    # 准备输入 - 检查实际可用的字段
                    if 'rna_fm_tokens' in batch:
                        rna_fm_tokens = batch['rna_fm_tokens'].squeeze(1)  # [B, L]
                    else:
                        # 如果没有rna_fm_tokens，使用tokens作为替代
                        print(f"    Warning: rna_fm_tokens not found, using tokens as substitute")
                        rna_fm_tokens = batch['tokens'][:, 0, :]  # 取第一个MSA序列
                    
                    seq_str = batch['seq'][0] if 'seq' in batch else "A" * batch['tokens'].shape[2]  # 取第一个序列或生成默认序列
                    valid_msa = batch.get('valid_msa', torch.ones(batch['tokens'].shape[0], dtype=torch.bool))
                    
                    # 前向传播
                    outputs = trainer.model(
                        tokens=batch['tokens'],
                        rna_fm_tokens=rna_fm_tokens,
                        seq=seq_str,
                        valid_msa=valid_msa
                    )
                    
                    print(f"    ✓ 前向传播成功")
                    
                    # 计算损失
                    loss = trainer.criterion(outputs, batch['target'])
                    print(f"    ✓ 损失计算成功: {loss.item():.4f}")
                    
                    # 反向传播
                    loss.backward()
                    print(f"    ✓ 反向传播成功")
                    
                    # 清除梯度
                    trainer.optimizer.zero_grad()
                    
                except Exception as e:
                    print(f"    ❌ Batch {batch_idx} 失败: {e}")
                    import traceback
                    traceback.print_exc()
                    break
                    
        print("\n✅ CUDA调试测试完成！")
        
    except Exception as e:
        print(f"❌ 训练失败: {e}")
        import traceback
        traceback.print_exc()


def main():
    """运行CUDA调试测试"""
    print("🔍 CUDA调试测试")
    print("=" * 50)
    
    try:
        test_cuda_debug()
    except Exception as e:
        print(f"❌ 测试失败: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()
