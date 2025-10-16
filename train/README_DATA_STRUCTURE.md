# RhoFold Diffusion Training - Data Structure Guide

## 数据目录结构

您的训练数据应该按照以下结构组织：

```
RNA3D_DATA/
├── list/           # 索引文件目录
│   ├── fold-0_train_ids
│   ├── fold-1_train_ids
│   ├── valid_fold-0
│   └── ...
├── seq/            # 序列文件目录
│   ├── 1a1t_B.seq
│   ├── 2noq_A.seq
│   └── ...
├── rMSA/           # MSA文件目录
│   ├── 1a1t_B.a3m
│   ├── 2noq_A.a3m
│   └── ...
└── pdb/            # 结构文件目录
    ├── 1a1t_B.pdb
    ├── 2noq_A.pdb
    └── ...
```

## 文件命名规则

- **索引文件**: 存储在 `list/` 目录下，无 `.txt` 后缀，每行一个样本名称
- **序列文件**: `{sample_name}.seq` 存储在 `seq/` 目录
- **MSA文件**: `{sample_name}.a3m` 存储在 `rMSA/` 目录  
- **结构文件**: `{sample_name}.pdb` 存储在 `pdb/` 目录

## 使用方法

### 1. 测试数据加载

```bash
cd /home/rcy/code/RhoFold
conda activate rnafold
python train/test_new_data_structure.py
```

### 2. 开始训练

#### 方法1: 使用启动脚本（推荐）

```bash
# 使用默认参数
bash train/start_training.sh

# 指定参数
bash train/start_training.sh fold-0_train_ids 8 200
```

#### 方法2: 直接运行训练脚本

```bash
python train/train_diffusion_final.py \
    --data_dir /home/rcy/code/RhoFold/RNA3D_DATA \
    --index_file fold-0_train_ids \
    --batch_size 4 \
    --epochs 100 \
    --learning_rate 1e-4 \
    --max_seq_len 512 \
    --msa_depth 128
```

### 3. 分布式训练（4x2080Ti）

```bash
# 启动分布式训练
bash train/launch_training.sh \
    --data_dir /home/rcy/code/RhoFold/RNA3D_DATA \
    --index_file fold-0_train_ids \
    --batch_size 16 \
    --epochs 100
```

## 参数说明

- `--data_dir`: 数据根目录路径
- `--index_file`: 索引文件名（不包含路径，会自动在 `list/` 目录下查找）
- `--batch_size`: 批次大小
- `--epochs`: 训练轮数
- `--max_seq_len`: 最大序列长度（默认512）
- `--msa_depth`: MSA最大深度（默认128）
- `--learning_rate`: 学习率（默认1e-4）
- `--save_dir`: 模型保存目录（默认./checkpoints）

## 注意事项

1. **环境要求**: 确保已激活 `rnafold` conda环境
2. **数据完整性**: 确保每个样本都有对应的 `.seq`、`.a3m`、`.pdb` 三个文件
3. **文件格式**: 
   - `.seq` 文件可以是纯序列或FASTA格式（会自动跳过头部）
   - `.a3m` 文件为标准MSA格式
   - `.pdb` 文件为标准PDB格式
4. **内存使用**: 根据GPU显存调整 `batch_size`，2080Ti建议使用4-8
5. **序列长度**: 超过 `max_seq_len` 的序列会被截断

## 故障排除

### 常见错误

1. **文件未找到**: 检查文件路径和命名是否正确
2. **序列长度不匹配**: MSA序列长度与主序列长度不一致时会自动调整
3. **内存不足**: 减小 `batch_size` 或 `max_seq_len`
4. **环境问题**: 确保conda环境正确激活且依赖包已安装

### 调试命令

```bash
# 检查数据目录结构
ls -la /home/rcy/code/RhoFold/RNA3D_DATA/*/

# 检查索引文件内容
head -5 /home/rcy/code/RhoFold/RNA3D_DATA/list/fold-0_train_ids

# 检查样本文件是否存在
sample_name=$(head -1 /home/rcy/code/RhoFold/RNA3D_DATA/list/fold-0_train_ids)
echo "Sample: $sample_name"
ls -la /home/rcy/code/RhoFold/RNA3D_DATA/seq/${sample_name}.seq
ls -la /home/rcy/code/RhoFold/RNA3D_DATA/rMSA/${sample_name}.a3m
ls -la /home/rcy/code/RhoFold/RNA3D_DATA/pdb/${sample_name}.pdb
```

