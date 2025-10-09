# RhoFold Diffusion Training Setup

## 项目概述

本项目为RhoFold模型集成了扩散模块，用于RNA结构预测的端到端训练。扩散模块替代了原始的E2Eformer，通过DDPM（Denoising Diffusion Probabilistic Models）方法改进配对特征的表示学习。

## 核心组件

### 1. 扩散模块 (`rhofold/model/diffusion.py`)
- **SinusoidalPositionalEmbedding**: 时间步编码
- **TimeEmbedding**: 时间嵌入MLP
- **ResidualBlock**: 带时间条件的UNet残差块
- **UNet2D**: 2D UNet用于配对特征处理
- **DiffusionModule**: 基于DDPM的扩散过程
- **DiffusionStack**: 替代E2EformerStack的主栈

### 2. 数据加载 (`train/data_loader.py`)
- 处理.seq、.a3m、.pdb文件格式
- 支持RNA序列tokenization
- 自动填充和掩码处理
- 分布式训练支持

### 3. 损失函数 (`train/losses.py`)
- **MLM Loss**: 掩码语言建模
- **Distance Loss**: 配对距离预测
- **Secondary Structure Loss**: 碱基配对预测
- **FAPE Loss**: 帧对齐点误差
- **Clash Loss**: 立体冲突惩罚
- **pLDDT Loss**: 置信度预测

### 4. 训练脚本
- **单GPU训练**: `train/train_diffusion.py`
- **多GPU训练**: `train/train_distributed.py`
- **验证测试**: `train/validate_setup.py`
- **前向测试**: `train/test_forward.py`

## 快速开始

### 1. 环境设置

```bash
# 创建conda环境
./setup_env.sh

# 激活环境
conda activate rnafold
```

### 2. 运行验证测试

```bash
# 运行所有验证测试
./quick_test.sh

# 或运行单个测试
python train/validate_setup.py
python train/test_forward.py
```

### 3. 准备训练数据

数据格式要求：
```
data_dir/
├── sample1.seq      # RNA序列
├── sample1.a3m      # MSA比对
├── sample1.pdb      # 3D结构
├── sample2.seq
├── sample2.a3m
├── sample2.pdb
└── ...
```

创建索引文件：
```
sample1
sample2
sample3
...
```

### 4. 生成示例数据（用于测试）

```bash
python train/create_sample_data.py \
    --output_dir ./sample_data \
    --num_samples 20 \
    --seq_length 100 \
    --msa_depth 5
```

### 5. 开始训练

#### 单GPU训练
```bash
python train/train_diffusion.py \
    --data_dir ./sample_data \
    --index_file ./sample_data/index.txt \
    --batch_size 4 \
    --epochs 10 \
    --learning_rate 1e-4
```

#### 多GPU训练（4x2080Ti）
```bash
./train/launch_training.sh \
    --data_dir ./sample_data \
    --index_file ./sample_data/index.txt \
    --batch_size 4 \
    --epochs 10 \
    --world_size 4 \
    --use_wandb
```

## 验证测试

### 前向传播测试
- 测试模型前向传播
- 验证输出形状和类型
- 检查错误和异常

### 损失计算测试
- 测试所有损失函数组件
- 验证损失值和梯度
- 检查NaN值

### 梯度流测试
- 测试扩散模块的梯度流
- 验证参数更新
- 确保梯度到达扩散层

### 内存使用测试
- 测试内存消耗
- 验证批次大小限制
- 检查内存泄漏

## 分布式训练

### 多GPU设置

分布式训练支持：
- **数据并行**: 跨多个GPU
- **梯度同步**: 使用NCCL后端
- **检查点保存**: 在rank 0进程
- **日志记录**: 集中在rank 0

### 启动脚本

```bash
./train/launch_training.sh \
    --data_dir /path/to/data \
    --index_file /path/to/index.txt \
    --world_size 4 \
    --batch_size 4 \
    --epochs 100 \
    --use_wandb
```

## 配置参数

### 模型参数

```python
diffusion_config = {
    'time_embed_dim': 256,      # 时间嵌入维度
    'num_timesteps': 1000,      # 扩散步数
    'noise_schedule': 'cosine', # 噪声调度类型
    'base_channels': 64,        # UNet基础通道数
}
```

### 训练参数

```python
training_config = {
    'batch_size': 4,            # 每GPU批次大小
    'learning_rate': 1e-4,      # 学习率
    'weight_decay': 1e-5,       # 权重衰减
    'epochs': 100,              # 训练轮数
    'max_seq_len': 512,         # 最大序列长度
    'msa_depth': 128,           # 最大MSA深度
}
```

## 监控和日志

### Weights & Biases集成

```bash
# 启用wandb日志
--use_wandb

# 设置实验名称
--experiment_name "my_experiment"
```

### 检查点保存

检查点保存在指定目录：
- `best_model.pth`: 基于验证损失的最佳模型
- `checkpoint_epoch_N.pth`: 常规检查点

## 性能优化

### 内存优化
- 使用适当的批次大小
- 启用梯度检查点
- 使用混合精度训练

### 速度优化
- 使用多个数据加载器工作进程
- 启用pin memory
- 使用分布式训练

### 收敛优化
- 监控损失曲线
- 调整学习率调度
- 使用早停

## 故障排除

### 常见问题

1. **CUDA内存不足**
   - 减少批次大小
   - 减少序列长度
   - 使用梯度检查点

2. **数据加载错误**
   - 检查文件格式（.seq, .a3m, .pdb）
   - 验证索引文件格式
   - 检查文件权限

3. **梯度问题**
   - 检查学习率
   - 验证损失计算
   - 检查NaN值

### 调试模式

运行调试信息：
```bash
python train/validate_setup.py --debug
```

## 文件结构

```
RhoFold/
├── rhofold/
│   ├── model/
│   │   ├── diffusion.py          # 扩散模块
│   │   ├── e2eformer.py          # 原始E2Eformer
│   │   └── ...
│   ├── rhofold_with_diffusion.py # 集成扩散模块的主模型
│   └── ...
├── train/
│   ├── data_loader.py            # 数据加载
│   ├── losses.py                 # 损失函数
│   ├── train_diffusion.py        # 单GPU训练
│   ├── train_distributed.py      # 多GPU训练
│   ├── test_forward.py           # 前向测试
│   ├── validate_setup.py         # 验证测试
│   ├── create_sample_data.py     # 示例数据生成
│   ├── launch_training.sh        # 训练启动脚本
│   └── README.md                 # 详细文档
├── setup_env.sh                  # 环境设置脚本
├── quick_test.sh                 # 快速测试脚本
└── TRAINING_SETUP.md             # 本文档
```

## 使用示例

```bash
# 1. 设置环境
./setup_env.sh
conda activate rnafold

# 2. 运行验证
./quick_test.sh

# 3. 生成示例数据
python train/create_sample_data.py --output_dir ./data --num_samples 50

# 4. 开始训练
./train/launch_training.sh \
    --data_dir ./data \
    --index_file ./data/index.txt \
    --experiment_name "diffusion_v1" \
    --batch_size 4 \
    --epochs 100 \
    --use_wandb
```

## 技术细节

### 扩散过程

1. **前向过程**: 向配对特征添加噪声
2. **反向过程**: 使用UNet预测并去除噪声
3. **时间嵌入**: 将时间步信息编码到网络中
4. **条件生成**: 使用MSA特征作为条件

### 损失函数组合

- **MLM**: 序列级别的掩码语言建模
- **Distance**: 配对距离预测（3个距离类型）
- **Secondary Structure**: 碱基配对预测
- **FAPE**: 3D结构对齐误差
- **Clash**: 立体冲突惩罚
- **pLDDT**: 置信度预测

### 分布式训练

- 使用PyTorch DDP进行数据并行
- NCCL后端用于GPU间通信
- 梯度同步和参数更新
- 检查点和日志管理

## 支持

如有问题和疑问：
1. 首先运行验证测试
2. 查看错误信息
3. 参考故障排除部分
4. 验证数据格式

## 总结

本训练系统提供了完整的RhoFold扩散模块训练管道，包括：

✅ **环境设置**: 自动化的conda环境配置
✅ **数据加载**: 支持.seq/.a3m/.pdb格式
✅ **模型集成**: 扩散模块无缝集成
✅ **损失函数**: 全面的损失计算
✅ **分布式训练**: 4x2080Ti多GPU支持
✅ **验证测试**: 全面的测试和验证
✅ **监控日志**: Weights & Biases集成
✅ **文档说明**: 详细的使用文档

现在您可以开始训练新的RhoFold扩散模型了！
