# Cut&Tag Analysis Pipline

This is a fully automated cut & tag analysis process based on the 'bash' script. It integrates the complete steps from FASTQ quality control, Spike-in normalization, genome alignment, Peak Calling to Motif analysis.
这是一个基于 `bash` 脚本的全自动 Cut&Tag 分析流程。它整合了从 FASTQ 质控、Spike-in 归一化、基因组比对、Peak Calling 到 Motif 分析的完整步骤。


## 1. Dependencies (Conda)

```bash
# git clone https://github.com/NatsumeAi/Cut-Tag-Script-Workflow.git
# cd Cut-Tag-Script-Workflow

conda env create -f environment.yml

conda activate CutTag

```
## 2. 运行
```bash
# Example：
./CutTag.sh -a WT_2_1.fq.gz -b WT_2_2.fq.gz -d WT_IgG_2_1.fq.gz -e WT_IgG_2_2.fq.gz -n WT2vsIgG2 -f SWO-M.genome.fa -g SWO-M_251011.gff3 -s SpikeIn.fa
# Help:
./CutTag.sh -h
```
