# 自动 Cut&Tag 分析流程

这是一个基于 `bash` 脚本的全自动 Cut&Tag 分析流程。它整合了从 FASTQ 质控、Spike-in 归一化、基因组比对、Peak Calling 到 Motif 分析的完整步骤。


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
