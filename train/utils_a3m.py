# utils_a3m.py
import os
import numpy as np

def read_a3m(a3m_path):
    """读取 A3M 文件，返回序列列表 [(header, seq), ...]"""
    seqs, header, seq = [], None, []
    with open(a3m_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    seqs.append((header, "".join(seq)))
                header = line[1:]
                seq = []
            else:
                seq.append(line)
        if header is not None:
            seqs.append((header, "".join(seq)))
    return seqs


def align_msa_to_seq(seq, a3m_input, alphabet="ACGUN"):
    """
    对齐 A3M 文件中的序列到给定 RNA 序列。
    如果 A3M 中 query 长度与 seq 不一致，则根据 gap 信息重建对齐。

    参数：
        seq: str, 原始 RNA 序列（无 gap）
        a3m_input: str, A3M 文件路径 或 已读入的 [(header, seq), ...]
        alphabet: str, RNA 碱基表

    返回：
        aligned_msa: np.ndarray, shape=(N, L), 每行等长，与 seq 一一对齐。
                     若对齐失败返回 None
    """
    if isinstance(a3m_input, str):
        if not os.path.exists(a3m_input):
            raise FileNotFoundError(f"A3M file not found: {a3m_input}")
        msa_list = read_a3m(a3m_input)
    else:
        msa_list = a3m_input

    if len(msa_list) == 0:
        return None

    query_header, query_seq = msa_list[0]
    seq_clean = "".join([c for c in query_seq if c.isalpha() and c.upper() in alphabet])
    
    # Step 1: 检查 query 与目标 seq 是否完全匹配
    if seq_clean.upper() == seq.upper():
        # query alignment 可直接用
        msa = [s for _, s in msa_list]
    else:
        # 如果不匹配，则尝试用对齐算法修复（简易 Needleman-Wunsch）
        from difflib import SequenceMatcher
        sm = SequenceMatcher(None, seq_clean.upper(), seq.upper())
        blocks = sm.get_opcodes()
        # 重新构造 alignment mapping
        rebuilt = []
        i, j = 0, 0
        for tag, i1, i2, j1, j2 in blocks:
            if tag == "equal":
                rebuilt.append(query_seq[i1:i2])
                i, j = i2, j2
            elif tag == "insert":
                # 在 query 缺失 -> 插入 gaps
                rebuilt.append("-" * (j2 - j1))
                j = j2
            elif tag == "delete":
                # query 多余 -> 删除这些列
                i = i2
            elif tag == "replace":
                # 用 gap 填充差异
                len_diff = (j2 - j1) - (i2 - i1)
                if len_diff > 0:
                    rebuilt.append(query_seq[i1:i2] + "-" * len_diff)
                elif len_diff < 0:
                    rebuilt.append(query_seq[i1:i2 + len_diff])
                else:
                    rebuilt.append(query_seq[i1:i2])
                i, j = i2, j2
        rebuilt_query = "".join(rebuilt)
        # replace query in msa_list
        msa_list[0] = (query_header, rebuilt_query)
        msa = [s for _, s in msa_list]

    # Step 2: 清理非 ACGUN 的字符为 '-'
    msa = [
        "".join([c if c.upper() in alphabet + "-" else "-" for c in s])
        for s in msa
    ]

    # Step 3: 确保所有行长度与目标 seq 一致
    target_len = len(seq)
    fixed = []
    for s in msa:
        s = s.replace(".", "-")  # 替换不常见 gap 符号
        if len(s) < target_len:
            s = s + "-" * (target_len - len(s))
        elif len(s) > target_len:
            s = s[:target_len]
        fixed.append(s)
    msa = fixed

    # Step 4: 转为 numpy 数组 (N,L)
    alphabet_dict = {a: i+1 for i, a in enumerate(alphabet)}
    alphabet_dict["-"] = 0
    arr = np.array([[alphabet_dict.get(c.upper(), 0) for c in row] for row in msa], dtype=np.int64)

    return arr