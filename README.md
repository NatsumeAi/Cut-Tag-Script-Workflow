# 自动 Cut&Tag 分析流程

这是一个基于 `bash` 脚本的全自动 Cut&Tag 分析流程。它整合了从 FASTQ 质控、Spike-in 归一化、基因组比对、Peak Calling 到 Motif 分析的完整步骤。


## 1. 依赖安装 (Conda)

本项目依赖的所有软件均通过 Conda 管理。请使用以下命令一键创建所需环境：

```bash
# 1. 克隆本项目 (或下载 .yml 文件)
# git clone https://github.com/NatsumeAi/Cut-Tag-Script-Workflow.git
# cd Cut-Tag-Script-Workflow

# 2. 使用 environment.yml 文件创建环境
conda env create -f environment.yml

# 3. 激活环境
conda activate CutTag
```
## 2. 运行
```bash
Example：
./CutTag.sh -a WT_2_1.fq.gz -b WT_2_2.fq.gz -d WT_IgG_2_1.fq.gz -e WT_IgG_2_2.fq.gz -n WT2vsIgG2 -f SWO-M.genome.fa -g SWO-M_251011.gff3 -s SpikeIn.fa
Help:
./CutTag.sh -h
```
