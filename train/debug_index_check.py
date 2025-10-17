import torch
from torch.utils.data import DataLoader

def check_tensor(name, tensor, vocab_size=None):
    """检查 tensor 是否存在非法索引"""
    if not torch.is_tensor(tensor):
        return []

    errors = []
    if tensor.dtype in (torch.int32, torch.int64):
        if (tensor < 0).any():
            errors.append(f"{name}: contains negative indices")
        if vocab_size is not None and (tensor >= vocab_size).any():
            maxv = tensor.max().item()
            errors.append(f"{name}: max index {maxv} >= vocab_size {vocab_size}")
    if torch.isnan(tensor).any():
        errors.append(f"{name}: contains NaN")
    if torch.isinf(tensor).any():
        errors.append(f"{name}: contains Inf")

    return errors


def scan_dataset(dataset, vocab_size_seq=30, vocab_size_msa=6, num_samples=100):
    """
    对 RNADataset 的前 num_samples 个样本进行扫描。
    vocab_size_seq/MSA 分别对应 RNA-FM 和 align_msa_to_seq 的词表大小。
    """
    print(f"🔍 Checking first {num_samples} samples...")
    for i in range(min(num_samples, len(dataset))):
        try:
            sample = dataset[i]
        except Exception as e:
            print(f"[!] Sample {i} load failed: {e}")
            continue

        for key, value in sample.items():
            if key in ("seq", "seq_tokens"):
                errs = check_tensor(f"sample[{i}].{key}", value, vocab_size_seq)
            elif key.startswith("msa"):
                errs = check_tensor(f"sample[{i}].{key}", value, vocab_size_msa)
            else:
                errs = check_tensor(f"sample[{i}].{key}", value)
            for e in errs:
                print(f"[ERROR] {e}")
    print("✅ Scan finished.")


def scan_dataloader(dataloader, vocab_size_seq=30, vocab_size_msa=6, num_batches=10):
    """
    对 DataLoader 的前 num_batches 个 batch 扫描。
    """
    print(f"🔍 Checking first {num_batches} batches...")
    for bi, batch in enumerate(dataloader):
        if bi >= num_batches:
            break
        for key, value in batch.items():
            if key in ("seq", "seq_tokens"):
                errs = check_tensor(f"batch[{bi}].{key}", value, vocab_size_seq)
            elif key.startswith("msa"):
                errs = check_tensor(f"batch[{bi}].{key}", value, vocab_size_msa)
            else:
                errs = check_tensor(f"batch[{bi}].{key}", value)
            for e in errs:
                print(f"[ERROR] {e}")
    print("✅ Scan finished.")
