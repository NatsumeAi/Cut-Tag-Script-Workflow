#!/bin/bash

# ==============================================================================
# 全自动 C-U-T-&-T-a-g 分析流程 (v5 - 完全自动化构建索引)
#
# 功能:
# 1. 零手动索引: 自动根据 FASTA 构建 目标基因组 和 Spike-in 的 Bowtie2 索引。
# 2. 命令行全配置: 通过参数指定所有输入 (FASTQ, FASTA, GFF)。
# 3. 环境自适应: 直接使用 Conda 环境中的软件 (picard, homer, bowtie2, etc.)。
# 4. 自动 Spike-in 归一化: 自动计算归一化因子并抽样。
# 5. 完整下游: 包含 Peak Calling, 注释, Metaplot, Motif 分析。
#
# ==============================================================================

# --- 1. 严格模式 ---
set -euo pipefail

# --- 2. 帮助信息 ---
usage() {
    echo "用法: $0 [数据参数] [基因组参数]"
    echo ""
    echo "数据参数:"
    echo "  -a [必需] Treatment 样本 Read 1 (e.g., WT_1_R1.fq.gz)"
    echo "  -b [必需] Treatment 样本 Read 2 (e.g., WT_1_R2.fq.gz)"
    echo "  -d [必需] Control (IgG) 样本 Read 1 (e.g., IgG_1_R1.fq.gz)"
    echo "  -e [必需] Control (IgG) 样本 Read 2 (e.g., IgG_1_R2.fq.gz)"
    echo "  -n [必需] 输出文件名前缀 (e.g., WT1_vs_IgG)"
    echo ""
    echo "基因组与 Spike-in 参数 (均为必需):"
    echo "  -f  [必需] 目标基因组 FASTA 文件 (e.g., genome.fa)"
    echo "  -g  [必需] 目标基因组 GFF/GTF 文件 (e.g., genes.gff)"
    echo "  -s  [必需] Spike-in 基因组 FASTA 文件 (e.g., E_coli.fa)"
    echo ""
    echo "其他:"
    echo "  -h  显示此帮助信息"
    echo ""
    exit 1
}

# ==============================================================================
# --- 3. 运行参数配置 (仅保留通用参数) ---
# ==============================================================================

THREADS=16
BIN_SIZE=10
TSS_UPSTREAM=2000
TSS_DOWNSTREAM=2000

# ==============================================================================

# --- 4. 解析输入参数 ---
TREAT_R1=""
TREAT_R2=""
CTRL_R1=""
CTRL_R2=""
OUTPUT_PREFIX=""
GENOME_FA=""
GENOME_GFF=""
SPIKEIN_FA=""

while getopts ":a:b:d:e:n:f:g:s:h" opt; do
    case ${opt} in
        a) TREAT_R1=$OPTARG ;;
        b) TREAT_R2=$OPTARG ;;
        d) CTRL_R1=$OPTARG ;;
        e) CTRL_R2=$OPTARG ;;
        n)  OUTPUT_PREFIX=$OPTARG ;;
        f)  GENOME_FA=$(realpath "$OPTARG") ;;
        g)  GENOME_GFF=$(realpath "$OPTARG") ;;
        s)  SPIKEIN_FA=$(realpath "$OPTARG") ;;
        h)  usage ;;
        \?) echo "无效的选项: -$OPTARG" >&2; usage ;;
        :)  echo "选项 -$OPTARG 需要一个参数。" >&2; usage ;;
    esac
done

# 检查必需参数
if [ -z "$TREAT_R1" ] || [ -z "$TREAT_R2" ] || [ -z "$CTRL_R1" ] || [ -z "$CTRL_R2" ] \
   || [ -z "$OUTPUT_PREFIX" ] || [ -z "$GENOME_FA" ] || [ -z "$GENOME_GFF" ] || [ -z "$SPIKEIN_FA" ]; then
    echo "错误: 缺少必需的参数 (请检查是否提供了 -s Spike-in FASTA)。"
    usage
fi

# --- 5. 阶段 0: 检查依赖并自动准备所有索引文件 ---
echo "--- [$(date)] 阶段 0: 环境检查与索引构建 ---"

# 检查关键软件
for tool in samtools bowtie2-build bc fastp seqkit macs2 bedtools computeMatrix picard findMotifsGenome.pl; do
    command -v $tool >/dev/null 2>&1 || { echo >&2 "错误: 未找到 '$tool'。请确保已激活 Conda 环境。"; exit 1; }
done
echo "所有软件依赖检查通过。"

# A. 准备目标基因组 FASTA 索引 (.fai)
if [ ! -f "${GENOME_FA}.fai" ]; then
    echo "构建目标基因组 FASTA 索引..."
    samtools faidx "$GENOME_FA"
fi

# B. 计算目标基因组大小
echo "计算目标基因组有效大小..."
CALCULATED_GENOME_SIZE=$(awk '{s+=$2} END {print s}' "${GENOME_FA}.fai")
echo "目标基因组大小: $CALCULATED_GENOME_SIZE bp"

# C. 准备目标基因组 Bowtie2 索引
# 直接使用 FASTA 文件路径作为索引前缀
GENOME_IDX_PREFIX="$GENOME_FA"
if [ ! -f "${GENOME_IDX_PREFIX}.1.bt2" ]; then
    echo "未找到目标基因组 Bowtie2 索引，正在构建 (这可能需要一些时间)..."
    bowtie2-build --threads $THREADS "$GENOME_FA" "$GENOME_IDX_PREFIX" > /dev/null
    echo "目标基因组索引构建完成。"
else
    echo "已找到目标基因组索引。"
fi

# D. 准备 Spike-in Bowtie2 索引
SPIKEIN_IDX_PREFIX="$SPIKEIN_FA"
if [ ! -f "${SPIKEIN_IDX_PREFIX}.1.bt2" ]; then
    echo "未找到 Spike-in Bowtie2 索引，正在构建..."
    bowtie2-build --threads $THREADS "$SPIKEIN_FA" "$SPIKEIN_IDX_PREFIX" > /dev/null
    echo "Spike-in 索引构建完成。"
else
    echo "已找到 Spike-in 索引。"
fi

# E. 准备 TSS BED 文件
GENOME_TSS_BED="${GENOME_GFF}.tss.bed"
if [ ! -f "$GENOME_TSS_BED" ]; then
    echo "正在从 GFF 生成 TSS BED 文件..."
    awk -v FS="\t" -v OFS="\t" '
        # (A) BEGIN 块: 定义一个包含所有常见"基因"特征词的"集合"
        #     这使得脚本无需修改就能适应不同的 GFF 注释来源
        BEGIN { 
            valid_types["gene"] = 1;
            valid_types["transcript"] = 1;
            valid_types["mRNA"] = 1;
            valid_types["lnc_RNA"] = 1;
            valid_types["miRNA"] = 1;
            valid_types["rRNA"] = 1;
            valid_types["tRNA"] = 1;
        }
        
        # (B) 主体块: 
        #     1. 跳过注释行 (以 # 开头)
        #     2. 确保 GFF 至少有 9 列 (NF >= 9)
        #     3. 检查第 3 列 ($3) 是否在我们定义的 "集合" 中
        !/^#/ && NF >= 9 && ($3 in valid_types) {
            
            # (C) 健壮性核心："清洗"所有关键字段
            
            # 清理染色体名 ($1) 和 链($7): 
            # 删除所有可能的空白字符 (空格, \t, \n, \r)
            gsub(/[[:space:]]+/, "", $1);
            gsub(/[[:space:]]+/, "", $7);
            
            # 清理属性 ($9): 
            # 将所有空白字符 (包括 \n 和 \r) 替换为分号 ";"
            # 这是防止换行符破坏输出格式的关键！
            gsub(/[[:space:]\r\n]+/, ";", $9);

            # (D) 打印 TSS
            #     只有当被"清洗"过的 $7 精确等于 "+" 或 "-" 时才打印
            if ($7 == "+") {
                # 打印: [chr] [start] [end] [name] [score] [strand]
                print $1, $4-1, $4, $9, "0", $7 
            } 
            else if ($7 == "-") {
                # 打印: [chr] [start] [end] [name] [score] [strand]
                print $1, $5-1, $5, $9, "0", $7
            }
        }
    ' "$GENOME_GFF" | sort -k1,1 -k2,2n -k3,3n -k6,6 -u > "$GENOME_TSS_BED"
fi
echo "--- [$(date)] 准备工作完成 ---"


# --- 6. 目录与路径设置 ---
TREAT_ID="${OUTPUT_PREFIX}_Treatment"
CTRL_ID="${OUTPUT_PREFIX}_Control"

# 创建目录结构
DIR_QC="1_QualityControl"
DIR_SPIKEIN="2_SpikeIn_Alignment"
DIR_NORM="3_Normalized_Fastq"
DIR_GENOME="4_Genome_Alignment"
DIR_BIGWIG="5_BigWig"
DIR_PEAKS="6_PeakCalling"
DIR_ANNO="7_Annotation"
DIR_METAPLOT="8_Metaplot"
DIR_MOTIF="9_Motif"
DIR_CHECKPOINT="0_Checkpoints"

mkdir -p $DIR_QC $DIR_SPIKEIN $DIR_NORM $DIR_GENOME $DIR_BIGWIG \
         $DIR_PEAKS $DIR_ANNO $DIR_METAPLOT $DIR_MOTIF $DIR_CHECKPOINT \
         $DIR_PEAKS $DIR_ANNO $DIR_METAPLOT $DIR_MOTIF

# 定义关键文件路径
TREAT_R1_CLEAN="${DIR_QC}/${TREAT_ID}_R1.clean.fq.gz"
TREAT_R2_CLEAN="${DIR_QC}/${TREAT_ID}_R2.clean.fq.gz"
CTRL_R1_CLEAN="${DIR_QC}/${CTRL_ID}_R1.clean.fq.gz"
CTRL_R2_CLEAN="${DIR_QC}/${CTRL_ID}_R2.clean.fq.gz"

TREAT_SPK_BAM="${DIR_SPIKEIN}/${TREAT_ID}-spikein.bam"
TREAT_SPK_SORT_BAM="${DIR_SPIKEIN}/${TREAT_ID}-spikein.sort.bam"
CTRL_SPK_BAM="${DIR_SPIKEIN}/${CTRL_ID}-spikein.bam"
CTRL_SPK_SORT_BAM="${DIR_SPIKEIN}/${CTRL_ID}-spikein.sort.bam"

TREAT_R1_NORM_FQ="${DIR_NORM}/${TREAT_ID}_R1.clean.spk-norm.fq.gz"
TREAT_R2_NORM_FQ="${DIR_NORM}/${TREAT_ID}_R2.clean.spk-norm.fq.gz"
CTRL_R1_NORM_FQ="${DIR_NORM}/${CTRL_ID}_R1.clean.spk-norm.fq.gz"
CTRL_R2_NORM_FQ="${DIR_NORM}/${CTRL_ID}_R2.clean.spk-norm.fq.gz"

TREAT_GENOME_BAM="${DIR_GENOME}/${TREAT_ID}.bam"
TREAT_GENOME_SORT_BAM="${DIR_GENOME}/${TREAT_ID}.sort.bam"
TREAT_GENOME_RMUP_BAM="${DIR_GENOME}/${TREAT_ID}.sort.rmdup.bam"
CTRL_GENOME_BAM="${DIR_GENOME}/${CTRL_ID}.bam"
CTRL_GENOME_SORT_BAM="${DIR_GENOME}/${CTRL_ID}.sort.bam"
CTRL_GENOME_RMUP_BAM="${DIR_GENOME}/${CTRL_ID}.sort.rmdup.bam"

TREAT_BIGWIG_FILE="${DIR_BIGWIG}/${TREAT_ID}.bigWig"
CTRL_BIGWIG_FILE="${DIR_BIGWIG}/${CTRL_ID}.bigWig"
PEAK_FILE="${DIR_PEAKS}/${OUTPUT_PREFIX}_peaks.narrowPeak"
PEAK_ANNO_FILE="${DIR_ANNO}/${OUTPUT_PREFIX}_peaks.annotated.txt"
MATRIX_FILE="${DIR_METAPLOT}/${OUTPUT_PREFIX}.TSS.gz"


# --- 7. 处理函数 ---

# [函数] 质控
run_qc() {
    local R1_IN=$1; local R2_IN=$2; local R1_OUT=$3; local R2_OUT=$4; local SAMPLE_ID=$5
    echo "--- [步骤 1] Fastp 质控: $SAMPLE_ID ---"
    fastp -i "$R1_IN" -o "$R1_OUT" -I "$R2_IN" -O "$R2_OUT" \
          -h "${DIR_QC}/${SAMPLE_ID}.fastp.html" -j "${DIR_QC}/${SAMPLE_ID}.fastp.json" -p $THREADS >/dev/null 2>&1
}

# [函数] Spike-in 比对与计数 (使用 SPIKEIN_IDX_PREFIX)
get_spikein_count() {
    local R1_CLEAN=$1; local R2_CLEAN=$2; local SPK_BAM=$3; local SPK_SORT_BAM=$4; local SAMPLE_ID=$5
    echo "--- [步骤 2] Spike-in 比对: $SAMPLE_ID ---" >&2
    # 使用阶段0生成的索引
    bowtie2 -p $THREADS -x "$SPIKEIN_IDX_PREFIX" -1 "$R1_CLEAN" -2 "$R2_CLEAN" -S "$SPK_BAM" >/dev/null 2>&1
    samtools sort -@ $THREADS -O bam -o "$SPK_SORT_BAM" "$SPK_BAM"
    samtools index "$SPK_SORT_BAM"
    rm "$SPK_BAM"
    
    local COUNT=$(samtools view -c -F 4 "$SPK_SORT_BAM")
    echo "$COUNT" # 仅返回数字用于变量赋值
}

# [函数] 归一化及基因组分析
process_genome() {
    local R1_CLEAN=$1; local R2_CLEAN=$2; local R1_NORM_FQ=$3; local R2_NORM_FQ=$4
    local NORM_FACTOR=$5; local GENOME_BAM=$6; local GENOME_SORT_BAM=$7; local GENOME_RMUP_BAM=$8
    local BIGWIG_FILE=$9; local SAMPLE_ID=${10}

    echo "--- [步骤 4] 归一化抽样: $SAMPLE_ID (Factor: $NORM_FACTOR) ---"
    seqkit sample -p "$NORM_FACTOR" "$R1_CLEAN" -o "$R1_NORM_FQ" >/dev/null 2>&1
    seqkit grep -f <(zcat "$R1_NORM_FQ" | awk 'NR%4==1' | cut -d ' ' -f 1 | sed 's/@//') "$R2_CLEAN" | gzip > "$R2_NORM_FQ"

    echo "--- [步骤 5] 基因组比对: $SAMPLE_ID ---"
    bowtie2 -p $THREADS -x "$GENOME_IDX_PREFIX" -1 "$R1_NORM_FQ" -2 "$R2_NORM_FQ" \
            --rg-id "$SAMPLE_ID" --rg "SM:$SAMPLE_ID" \
            -S "$GENOME_BAM" >/dev/null 2>&1
    samtools sort -@ $THREADS -O bam -o "$GENOME_SORT_BAM" "$GENOME_BAM"
    samtools index "$GENOME_SORT_BAM"
    rm "$GENOME_BAM"

    echo "--- [步骤 6] 去重 (Picard): $SAMPLE_ID ---"
    picard  MarkDuplicates --REMOVE_DUPLICATES true -I "$GENOME_SORT_BAM" -O "$GENOME_RMUP_BAM" -M "${DIR_GENOME}/${SAMPLE_ID}_markdup_metrics.txt" > "${DIR_GENOME}/${SAMPLE_ID}.rmdup.log" 2>&1
    samtools index "$GENOME_RMUP_BAM"
    rm "$GENOME_SORT_BAM" "${GENOME_SORT_BAM}.bai"

    echo "--- [步骤 9-Pre] 生成 BigWig: $SAMPLE_ID ---"
    bamCoverage --bam "$GENOME_RMUP_BAM" --binSize $BIN_SIZE --outFileName "$BIGWIG_FILE" \
                --outFileFormat bigwig --ignoreDuplicates --normalizeUsing RPGC \
                --numberOfProcessors $THREADS --effectiveGenomeSize $CALCULATED_GENOME_SIZE >/dev/null 2>&1
}


# --- 8. 主流程执行 ---

# --- 定义检查点文件 ---
FLAG_QC="${DIR_CHECKPOINT}/${OUTPUT_PREFIX}.1_qc.done"
FLAG_SPIKEIN="${DIR_CHECKPOINT}/${OUTPUT_PREFIX}.2_spikein.done"
FLAG_GENOME="${DIR_CHECKPOINT}/${OUTPUT_PREFIX}.3_genome.done"
FLAG_PEAKCALL="${DIR_CHECKPOINT}/${OUTPUT_PREFIX}.4_peakcall.done"
FLAG_ANNO="${DIR_CHECKPOINT}/${OUTPUT_PREFIX}.5_annotation.done"
FLAG_METAPLOT="${DIR_CHECKPOINT}/${OUTPUT_PREFIX}.6_metaplot.done"
FLAG_MOTIF="${DIR_CHECKPOINT}/${OUTPUT_PREFIX}.7_motif.done"

# 定义 Spike-in 计数缓存文件
COUNT_TREAT_CACHE="${DIR_CHECKPOINT}/${OUTPUT_PREFIX}.treat.count"
COUNT_CTRL_CACHE="${DIR_CHECKPOINT}/${OUTPUT_PREFIX}.ctrl.count"


echo "--- [$(date)] 阶段 A-1: 质控 ---"
if [ ! -f "$FLAG_QC" ]; then
    run_qc "$TREAT_R1" "$TREAT_R2" "$TREAT_R1_CLEAN" "$TREAT_R2_CLEAN" "$TREAT_ID" &
    run_qc "$CTRL_R1" "$CTRL_R2" "$CTRL_R1_CLEAN" "$CTRL_R2_CLEAN" "$CTRL_ID" &
    wait
    touch "$FLAG_QC" # <--- 成功后创建标记
else
    echo ">>> [跳过] 步骤 1: Fastp 质控 (已完成)"
fi


echo "--- [$(date)] 阶段 A-2: Spike-in 计数 ---"
if [ ! -f "$FLAG_SPIKEIN" ]; then
    # [修正] 直接将函数的 stdout (即 COUNT) 重定向到缓存文件
    get_spikein_count "$TREAT_R1_CLEAN" "$TREAT_R2_CLEAN" "$TREAT_SPK_BAM" "$TREAT_SPK_SORT_BAM" "$TREAT_ID" > "$COUNT_TREAT_CACHE" &
    PID_T=$!
    get_spikein_count "$CTRL_R1_CLEAN" "$CTRL_R2_CLEAN" "$CTRL_SPK_BAM" "$CTRL_SPK_SORT_BAM" "$CTRL_ID" > "$COUNT_CTRL_CACHE" &
    PID_C=$!
    
    # 等待两个后台作业完成
    wait $PID_T
    wait $PID_C

    # [修正] 现在从缓存文件中安全地读取计数值
    TREAT_SPK_COUNT=$(cat "$COUNT_TREAT_CACHE")
    CTRL_SPK_COUNT=$(cat "$COUNT_CTRL_CACHE")
    
    # 检查计数值是否为空 (防止函数失败)
    if [ -z "$TREAT_SPK_COUNT" ] || [ -z "$CTRL_SPK_COUNT" ]; then
        echo "错误: Spike-in 计数失败，缓存文件为空。" >&2
        exit 1
    fi
    
    touch "$FLAG_SPIKEIN" # <--- 成功后创建标记
else
    echo ">>> [跳过] 步骤 2: Spike-in 比对 (已完成)"
    # <--- [重要] 从缓存文件读回计数值
    TREAT_SPK_COUNT=$(cat "$COUNT_TREAT_CACHE")
    CTRL_SPK_COUNT=$(cat "$COUNT_CTRL_CACHE")
fi
echo ">>> Spike-in Reads - Treat: $TREAT_SPK_COUNT | Ctrl: $CTRL_SPK_COUNT"


echo "--- [$(date)] 阶段 B: 计算归一化因子 ---"
MIN_COUNT=$( (echo "$TREAT_SPK_COUNT"; echo "$CTRL_SPK_COUNT") | sort -n | head -n 1 )
if [ "$MIN_COUNT" -eq 0 ]; then echo "错误: Spike-in 计数为 0，无法归一化。"; exit 1; fi
TREAT_FACTOR=$(echo "scale=6; $MIN_COUNT / $TREAT_SPK_COUNT" | bc)
CTRL_FACTOR=$(echo "scale=6; $MIN_COUNT / $CTRL_SPK_COUNT" | bc)
echo ">>> 归一化因子 - Treat: $TREAT_FACTOR | Ctrl: $CTRL_FACTOR"


echo "--- [$(date)] 阶段 C: 基因组比对与处理 ---"
if [ ! -f "$FLAG_GENOME" ]; then
    process_genome "$TREAT_R1_CLEAN" "$TREAT_R2_CLEAN" "$TREAT_R1_NORM_FQ" "$TREAT_R2_NORM_FQ" \
                   "$TREAT_FACTOR" "$TREAT_GENOME_BAM" "$TREAT_GENOME_SORT_BAM" "$TREAT_GENOME_RMUP_BAM" \
                   "$TREAT_BIGWIG_FILE" "$TREAT_ID" &
    PID_GT=$!
    process_genome "$CTRL_R1_CLEAN" "$CTRL_R2_CLEAN" "$CTRL_R1_NORM_FQ" "$CTRL_R2_NORM_FQ" \
                   "$CTRL_FACTOR" "$CTRL_GENOME_BAM" "$CTRL_GENOME_SORT_BAM" "$CTRL_GENOME_RMUP_BAM" \
                   "$CTRL_BIGWIG_FILE" "$CTRL_ID" &
    PID_GC=$!
    wait $PID_GT
    wait $PID_GC
    touch "$FLAG_GENOME" # <--- 成功后创建标记
else
    echo ">>> [跳过] 步骤 4-6: 基因组比对、去重、BigWig (已完成)"
fi


echo "--- [$(date)] 阶段 D: 下游分析 ---"

if [ ! -f "$FLAG_PEAKCALL" ]; then
    echo "[步骤 7] Peak Calling (MACS2)"
    macs2 callpeak -t "$TREAT_GENOME_RMUP_BAM" -c "$CTRL_GENOME_RMUP_BAM" \
                   -f BAMPE -g $CALCULATED_GENOME_SIZE --bdg --outdir "$DIR_PEAKS" -n "$OUTPUT_PREFIX" \
                   > "${DIR_PEAKS}/${OUTPUT_PREFIX}-macs2.log" 2>&1
    touch "$FLAG_PEAKCALL" # <--- 成功后创建标记
else
    echo ">>> [跳过] 步骤 7: Peak Calling (MACS2) (已完成)"
fi


if [ ! -f "$FLAG_ANNO" ]; then
    echo "[步骤 8] Peak 注释"
    bedtools intersect -a "$PEAK_FILE" \
                       -b <(awk -v FS="\t" '{if($3=="gene" || $3=="transcript" || $3=="chromosome"){pass}else{print $0}}' "$GENOME_GFF") \
                       -wa -wb -loj | uniq -w 40 > "$PEAK_ANNO_FILE"
    touch "$FLAG_ANNO" # <--- 成功后创建标记
else
    echo ">>> [跳过] 步骤 8: Peak 注释 (已完成)"
fi


if [ ! -f "$FLAG_METAPLOT" ]; then
    echo "[步骤 9] Metaplot 绘制"
    computeMatrix reference-point -S "$TREAT_BIGWIG_FILE" "$CTRL_BIGWIG_FILE" \
        -R "$GENOME_TSS_BED" --referencePoint TSS -a $TSS_UPSTREAM -b $TSS_DOWNSTREAM \
        -out "$MATRIX_FILE" -p $THREADS
    plotProfile -m "$MATRIX_FILE" --perGroup --numPlotsPerRow 2 \
        -out "${DIR_METAPLOT}/${OUTPUT_PREFIX}.TSS.profile.pdf"
    plotHeatmap -m "$MATRIX_FILE" --colorList 'white,red' 'white,blue' \
        -out "${DIR_METAPLOT}/${OUTPUT_PREFIX}.TSS.heatmap.pdf"
    touch "$FLAG_METAPLOT" # <--- 成功后创建标记
else
    echo ">>> [跳过] 步骤 9: Metaplot 绘制 (已完成)"
fi


if [ ! -f "$FLAG_MOTIF" ]; then
    echo "[步骤 10] Motif 鉴定 (HOMER)"
    
    # 1. ...
    RUN_DIR_MOTIF="${DIR_MOTIF}/${OUTPUT_PREFIX}"  # <--- 已删除 "local"
    
    # 2. ...
    mkdir -p "$RUN_DIR_MOTIF"

    # 3. ...
    findMotifsGenome.pl "$PEAK_FILE" "$GENOME_FA" "$RUN_DIR_MOTIF" -p $THREADS \
                        > "${RUN_DIR_MOTIF}/${OUTPUT_PREFIX}-homer.log" 2>&1
                        
    touch "$FLAG_MOTIF" # <--- 成功后创建标记
else
    echo ">>> [跳过] 步骤 10: Motif 鉴定 (HOMER) (已完成)"
fi

echo "--- [SUCCESS] 流程完成 ---"