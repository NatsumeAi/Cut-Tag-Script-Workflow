# 全自动 Cut&Tag 分析流程 (My-CutTag-Pipeline)

这是一个基于 `bash` 脚本的全自动 Cut&Tag 分析流程。它整合了从 FASTQ 质控、Spike-in 归一化、基因组比对、Peak Calling 到 Motif 分析的完整步骤。

## 亮点

* **全自动**：仅需提供 FASTQ、FASTA 和 GFF 文件，脚本自动处理所有中间步骤。
* **Spike-in 归一化**：自动使用 Spike-in 基因组比对结果计算归一化因子。
* **健壮性 GFF 解析**：自动从 GFF/GTF 文件中提取唯一的 TSS 位置，能兼容多种不规范的 GFF 格式。
* **断点续传**：流程使用检查点文件，如果中途失败，修复后可自动从失败步骤恢复。

---

## 1. 依赖安装 (Conda)

本项目依赖的所有软件均通过 Conda 管理。请使用以下命令一键创建所需环境：

```bash
# 1. 克隆本项目 (或下载 .yml 文件)
# git clone [https://github.com/YourUsername/My-CutTag-Pipeline.git](https://github.com/YourUsername/My-CutTag-Pipeline.git)
# cd My-CutTag-Pipeline

# 2. 使用 environment.yml 文件创建环境
conda env create -f environment.yml

# 3. 激活环境
conda activate CutTag
