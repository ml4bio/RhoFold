# MSA Handling Modifications

## 概述

根据您的要求，我已经修改了数据加载和模型处理逻辑，以支持三种不同的MSA数据情况：

1. **MSA和序列长度不同** - 使用对齐算法将MSA对齐到给定序列
2. **缺失MSA数据** - 跳过MSA模块，使用零向量填充
3. **正常数据** - 保持原有处理流程

## 修改的文件

### 1. `train/data_loader.py`

#### 新增功能：
- **`_load_msa_with_alignment()`**: 使用`utils_a3m.py`中的对齐函数处理MSA
- **`valid_msa`标志**: 在数据样本中标记MSA是否有效
- **三种处理逻辑**:
  - 有效MSA + 对齐成功 → 正常处理
  - 无效MSA或对齐失败 → 零向量填充
  - 缺失MSA文件 → 零向量填充

#### 关键修改：
```python
# 在__getitem__中
if msa_path.exists():
    try:
        msa_sequences = self._load_msa_with_alignment(msa_path, sequence)
        valid_msa = True
    except Exception as e:
        msa_sequences = None
        valid_msa = False
else:
    msa_sequences = None
    valid_msa = False

# 在_process_sample中
if valid_msa and msa_sequences is not None:
    # 正常处理MSA
    msa_tokens = process_msa_sequences(msa_sequences)
else:
    # 零向量填充
    msa_tokens = torch.zeros(self.msa_depth, seq_len, dtype=torch.long)
```

### 2. `train/train_diffusion_final.py`

#### 新增功能：
- **`valid_msa`参数传递**: 从batch中提取MSA有效性标志
- **模型调用更新**: 将`valid_msa`传递给模型

#### 关键修改：
```python
# 获取MSA有效性标志
valid_msa = batch.get('valid_msa', torch.ones(batch['tokens'].shape[0], dtype=torch.bool))

# 传递给模型
outputs = self.model(
    tokens=batch['tokens'],
    rna_fm_tokens=rna_fm_tokens,
    seq=seq_str,
    valid_msa=valid_msa
)
```

### 3. `rhofold/rhofold_with_diffusion.py`

#### 新增功能：
- **`valid_msa`参数支持**: 在forward方法中接受MSA有效性标志
- **MSA处理逻辑**: 对无效MSA样本进行特殊处理

#### 关键修改：
```python
def forward_one_cycle(self, tokens, rna_fm_tokens, recycling_inputs, seq, valid_msa=None):
    # ... 正常MSA处理 ...
    
    # 处理无效MSA样本
    if valid_msa is not None:
        invalid_msa_mask = ~valid_msa
        if invalid_msa_mask.any():
            # 对无效MSA样本，将MSA特征置零
            # 模型将完全依赖RNA-FM序列嵌入
            msa_fea[invalid_msa_mask] = 0.0
```

### 4. `train/utils_a3m.py`

#### 功能：
- **`read_a3m()`**: 读取A3M格式文件
- **`align_msa_to_seq()`**: 将MSA序列对齐到给定RNA序列
- **对齐算法**: 使用简化的Needleman-Wunsch算法处理长度不匹配

## 使用方法

### 1. 数据准备

确保您的数据按以下结构组织：
```
RNA3D_DATA/
├── list/           # 索引文件
├── seq/            # 序列文件
├── rMSA/           # MSA文件（.a3m格式）
└── pdb/            # 结构文件
```

### 2. 训练命令

```bash
# 基本训练
python train/train_diffusion_final.py \
    --data_dir /home/rcy/code/RhoFold/RNA3D_DATA \
    --index_file fold-0_train_ids \
    --batch_size 4 \
    --epochs 100

# 使用启动脚本
bash train/start_training.sh fold-0_train_ids 4 100
```

### 3. 测试命令

```bash
# 测试数据加载
python train/test_new_data_structure.py

# 测试MSA处理
python train/test_msa_handling.py

# 简单导入测试
python train/test_simple_msa.py
```

## 处理逻辑详解

### 情况1: MSA和序列长度不同

```python
# 使用对齐算法
aligned_msa = align_msa_to_seq(sequence, msa_path)
if aligned_msa is not None:
    # 对齐成功，正常处理
    msa_sequences = convert_aligned_to_sequences(aligned_msa)
    valid_msa = True
else:
    # 对齐失败，使用零向量
    msa_sequences = None
    valid_msa = False
```

### 情况2: 缺失MSA数据

```python
if not msa_path.exists():
    # 文件不存在，使用零向量
    msa_sequences = None
    valid_msa = False
    print(f"Warning: MSA file not found for {sample_name}")
```

### 情况3: 正常数据

```python
if valid_msa and msa_sequences is not None:
    # 正常处理MSA序列
    msa_tokens = process_msa_sequences(msa_sequences)
else:
    # 零向量填充
    msa_tokens = torch.zeros(msa_depth, seq_len, dtype=torch.long)
```

## 模型行为

### 有效MSA样本
- 使用完整的MSA嵌入
- 结合RNA-FM序列嵌入
- 正常的前向传播

### 无效MSA样本
- MSA嵌入被置零
- 完全依赖RNA-FM序列嵌入
- 模型仍能正常训练和推理

## 注意事项

1. **对齐算法**: 使用简化的Needleman-Wunsch算法，对于复杂情况可能需要改进
2. **性能影响**: 对齐过程会增加数据加载时间
3. **内存使用**: 零向量填充不会显著增加内存使用
4. **训练稳定性**: 无效MSA样本仍能参与训练，有助于模型鲁棒性

## 故障排除

### 常见问题

1. **对齐失败**: 检查序列格式和MSA文件完整性
2. **内存不足**: 减小batch_size或max_seq_len
3. **导入错误**: 确保所有依赖包已安装

### 调试命令

```bash
# 检查数据目录结构
ls -la /home/rcy/code/RhoFold/RNA3D_DATA/*/

# 检查特定样本的MSA文件
sample_name=$(head -1 /home/rcy/code/RhoFold/RNA3D_DATA/list/fold-0_train_ids)
ls -la /home/rcy/code/RhoFold/RNA3D_DATA/rMSA/${sample_name}.a3m

# 运行简单测试
python train/test_simple_msa.py
```

## 总结

这些修改确保了训练流程能够处理各种MSA数据情况，提高了模型的鲁棒性和训练稳定性。无论MSA数据是否完整，模型都能正常训练和推理。
