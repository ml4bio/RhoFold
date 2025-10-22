#!/bin/bash
# ==========================================================
# RhoFold Git 版本管理与实验分支创建脚本
# 适用于科研训练环境 (单GPU / 多GPU)
# ==========================================================

# 项目路径（自动检测当前目录）
PROJECT_DIR=$(pwd)
PROJECT_NAME=$(basename "$PROJECT_DIR")

# 分支名称（新实验版本）
NEW_BRANCH="train_improvement"

# 实验结果目录
RESULTS_DIR="${PROJECT_DIR}/results"
mkdir -p "$RESULTS_DIR"

echo "📁 项目路径: $PROJECT_DIR"
echo "🌿 新实验分支: $NEW_BRANCH"
echo "📂 实验结果目录: $RESULTS_DIR"
echo "----------------------------------------------------------"

# ==========================================================
# 1️⃣ 初始化 Git 仓库（如未初始化）
# ==========================================================
if [ ! -d "$PROJECT_DIR/.git" ]; then
    echo "🧩 初始化 Git 仓库..."
    git init
    git add .
    git commit -m "Initial commit: baseline RhoFold version"
else
    echo "✅ 已检测到现有 Git 仓库"
fi

# ==========================================================
# 2️⃣ 创建 .gitignore 文件
# ==========================================================
cat > "${PROJECT_DIR}/.gitignore" << 'EOF'
# 忽略大型文件与临时结果
__pycache__/
*.pyc
*.pkl
*.pt
*.bin
*.ckpt
*.log
*.tmp
*.swp
.DS_Store

# 数据与模型目录
RNA3D_DATA/
checkpoints/
checkpoints_fold-0/
results/
pretrained/
envs/
database/
EOF
echo "✅ 已创建 .gitignore"

# ==========================================================
# 3️⃣ 保存当前稳定版本为 v1_baseline
# ==========================================================
echo "🧱 保存当前版本快照: v1_baseline"
git add .
git commit -m "Save baseline version before diffusion refactor"
git tag -f v1_baseline

# ==========================================================
# 4️⃣ 创建新分支进行实验改动
# ==========================================================
git checkout -b "$NEW_BRANCH"
echo "✅ 新分支 '$NEW_BRANCH' 创建完成"

# ==========================================================
# 5️⃣ 创建新的实验结果目录
# ==========================================================
EXP_DIR="${RESULTS_DIR}/$(date +'%Y%m%d')_${NEW_BRANCH}"
mkdir -p "$EXP_DIR/logs" "$EXP_DIR/checkpoints"
echo "✅ 已创建实验输出路径: $EXP_DIR"

# ==========================================================
# 6️⃣ 记录当前 Git 版本信息
# ==========================================================
GIT_HASH=$(git rev-parse --short HEAD)
echo "Commit: $GIT_HASH" > "$EXP_DIR/version.txt"
git log -1 --pretty=fuller >> "$EXP_DIR/version.txt"
echo "✅ 当前版本信息已保存到 $EXP_DIR/version.txt"

echo "----------------------------------------------------------"
echo "✨ 初始化完成！"
echo "当前分支: $(git branch --show-current)"
echo "Git 版本: $GIT_HASH"
echo "结果目录: $EXP_DIR"
echo "----------------------------------------------------------"
echo "下一步建议："
echo "1️⃣ 在 $NEW_BRANCH 上修改 train/train_diffusion_final.py"
echo "2️⃣ 训练时将结果输出到: $EXP_DIR"
echo "3️⃣ 完成后执行: git add . && git commit -m 'Finish diffusion version'"
echo "4️⃣ 可随时切回旧版本: git checkout v1_baseline"
echo "----------------------------------------------------------"