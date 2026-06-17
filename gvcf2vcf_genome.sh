#!/usr/bin/env bash
set -Eeuo pipefail
shopt -s nullglob

########################################
# gvcf2vcf_genome_gpt_3m6d.sh
#
# 功能: 对 12 个样本（2 亲本 + 10 混池）的 gVCF 做联合分型，输出全基因组 VCF
# 流程: 按 contig 逐个运行 GenomicsDBImport → GenotypeGVCFs，最终 GatherVcfs 合并
# 输入: sample<tab>path 格式的列表文件（每行一个样本）
# 生成列表示例:
#   for f in *.g.vcf.gz; do
#     sample=$(basename "$f" .g.vcf.gz)
#     echo -e "${sample}\t$(realpath "$f")"
#   done > sample.list
########################################

########################################
# Configuration (defaults)
########################################

DEFAULT_LIST="/data2/chenh/bsa/F2_3WORKSPACE/gvcf/F2_3_hunchi/PP64_REF_gvcf/gvcf_sample.list"
LIST_FILE="$DEFAULT_LIST"

qunti_name=F2_3
REF_NAME=PP64

REFERENCE="/data2/chenh/bsa/REF/PP64_Cgla_REF/PP64.genome.fa"
OUTDIR="/data2/chenh/bsa/F2_3WORKSPACE/vcf/260615_PP64_REF_vcf"

INCLUDE_NON_VARIANT_SITES="no"   # 是否包括非变异位点（yes=all-sites，文件极大）
MAX_JOBS=3                        # 同时处理的 contig 数（NFS 场景下 3 并发最稳）
JAVA_XMX_IMPORT="20g"             # GenomicsDBImport 的堆内存
JAVA_XMX_GENOTYPE="28g"           # GenotypeGVCFs 的堆内存（分型比导入更吃内存）
JAVA_XMX_GATHER="12g"             # GatherVcfs 最终合并的堆内存
READER_THREADS=4                   # 读取 gVCF 的线程数（共享文件系统不宜过大）
USE_SHARED_POSIXFS_OPTIMIZATIONS="yes"  # NFS/Lustre 优化开关，建议保持开启

########################################
# CLI 参数覆写（可选）
########################################

[[ $# -ge 1 ]] && LIST_FILE="$1"
[[ $# -ge 2 ]] && REFERENCE="$2"
[[ $# -ge 3 ]] && OUTDIR="$3"
[[ $# -gt 3 ]] && {
    echo "Usage: $0 [gvcf_list] [ref.fa] [outdir]" >&2
    exit 1
}

# Remove trailing slashes
REFERENCE="${REFERENCE%/}"
OUTDIR="${OUTDIR%/}"

########################################
# 派生路径（根据上面的配置自动计算）
########################################

REF_DIR="$(dirname "$REFERENCE")"
REF_BASE="$(basename "$REFERENCE")"
REF_STEM="${REF_DIR}/${REF_BASE%.*}"
REF_DICT="${REF_STEM}.dict"

RUN_PREFIX="$(basename "$LIST_FILE")"
RUN_PREFIX="${RUN_PREFIX%.gvcf.list}"
RUN_PREFIX="${RUN_PREFIX%.list}"

LOG_DIR="${OUTDIR}/logs"
STATUS_DIR="${OUTDIR}/status"
PER_INTERVAL_ROOT="${OUTDIR}/per_interval"
TMP_ROOT="${OUTDIR}/tmp"
DB_ROOT="${OUTDIR}/genomicsdb"

SAMPLE_MAP="${OUTDIR}/${RUN_PREFIX}.sample_map"
INTERVAL_META="${OUTDIR}/${RUN_PREFIX}.intervals.tsv"
SHARD_LIST="${OUTDIR}/${RUN_PREFIX}.shards.list"
FINAL_VCF="${OUTDIR}/${qunti_name}_${REF_NAME}.genome.joint.vcf.gz"

TIMESTAMP="$(date '+%Y%m%d_%H%M%S')"
MASTER_LOG="${LOG_DIR}/${RUN_PREFIX}.run_${TIMESTAMP}.log"

########################################
# Functions
########################################

log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*"; }
die() { log "FATAL: $*"; exit 1; }
need_cmd() { command -v "$1" >/dev/null 2>&1 || die "not found: $1"; }
has_index() { [[ -f "$1.tbi" || -f "$1.idx" || -f "$1.csi" ]]; }

########################################
# GatherVcfs 合并步骤（独立函数，便于单独调用）
########################################

run_merge() {
    log "Building shard list..."
    : > "$SHARD_LIST"
    while IFS=$'\t' read -r idx interval len safe; do
        tag="${idx}_${safe}"
        shard="${PER_INTERVAL_ROOT}/${tag}/${tag}.joint.vcf.gz"
        [[ -s "$shard" ]] || die "Missing shard: $shard"
        has_index "$shard" || gatk IndexFeatureFile -I "$shard"
        echo "$shard" >> "$SHARD_LIST"
    done < "$INTERVAL_META"

    local total
    total=$(wc -l < "$SHARD_LIST")
    log "Shards: $total"

    mapfile -t shards < "$SHARD_LIST"

    if [[ ${#shards[@]} -eq 1 ]]; then
        cp -f "${shards[0]}" "$FINAL_VCF"
    else
        local gather_args=()
        for v in "${shards[@]}"; do gather_args+=(-I "$v"); done
        log "Running GatherVcfs..."
        gatk --java-options "-Xms${JAVA_XMX_GATHER} -Xmx${JAVA_XMX_GATHER} -Djava.io.tmpdir=${SCRATCH_TMP}/gather" \
            GatherVcfs "${gather_args[@]}" -O "$FINAL_VCF"
    fi

    has_index "$FINAL_VCF" || gatk IndexFeatureFile -I "$FINAL_VCF"
    log "Final VCF: $FINAL_VCF"
}

########################################
# 环境准备
########################################

mkdir -p "$LOG_DIR" "$STATUS_DIR" "$PER_INTERVAL_ROOT" "$TMP_ROOT" "$DB_ROOT"

# 优先使用本地内存盘 /dev/shm 做 scratch，避免 NFS/Lustre IO 瓶颈
# 若 /dev/shm 不可用则回退到 OUTDIR/tmp
SCRATCH=""
for d in "/dev/shm/${USER}/${RUN_PREFIX}" "$TMP_ROOT"; do
    mkdir -p "$d" 2>/dev/null && SCRATCH="$d" && break
done
SCRATCH_DB="${SCRATCH}/genomicsdb"
SCRATCH_TMP="${SCRATCH}/tmp"
mkdir -p "$SCRATCH_DB" "$SCRATCH_TMP"

exec > >(tee -a "$MASTER_LOG") 2>&1
trap 'die "Abnormal termination. Log: $MASTER_LOG"' ERR

need_cmd gatk awk sed

[[ -f "$LIST_FILE" ]]  || die "List not found: $LIST_FILE"
[[ -f "$REFERENCE" ]]  || die "Reference not found: $REFERENCE"
[[ -f "${REFERENCE}.fai" ]] || die "Missing .fai: ${REFERENCE}.fai"
[[ -f "$REF_DICT" ]]   || die "Missing .dict: $REF_DICT"

log "=== GVCF->VCF joint calling ==="
log "List: $LIST_FILE | Ref: $REFERENCE | Out: $OUTDIR"
log "Samples: from list | Jobs: $MAX_JOBS | NonVarSites: $INCLUDE_NON_VARIANT_SITES"

########################################
# 解析输入列表（格式: sample<tab>path）
#
# 跳过 # 开头的注释行和空行，验证文件存在且有索引，
# 检查重复样本名。
########################################

: > "$SAMPLE_MAP"
TOTAL_SAMPLES=0
declare -A SEEN_SAMPLES=()

while IFS=$'\t' read -r sample path rest; do
    [[ -n "$sample" && "$sample" != '#'* ]] || continue
    [[ -n "$path" && -z "${rest:-}" ]] && [[ "$sample" != '#' ]] || {
        die "Bad format (need sample<tab>path): line starting with '${sample}'"
    }

    [[ -z "${SEEN_SAMPLES[$sample]:-}" ]] || die "Duplicate sample: $sample"
    SEEN_SAMPLES["$sample"]=1

    [[ -f "$path" ]]    || die "gVCF not found for $sample: $path"
    has_index "$path"   || die "Missing index for $sample: $path"

    printf '%s\t%s\n' "$sample" "$path" >> "$SAMPLE_MAP"
    ((++TOTAL_SAMPLES))
done < "$LIST_FILE"

[[ "$TOTAL_SAMPLES" -gt 0 ]] || die "No valid entries in: $LIST_FILE"
log "Samples loaded: $TOTAL_SAMPLES"

########################################
# 从参考基因组的 .fai 生成 interval 清单
# 格式: 序号(补零)  contig名  长度  安全文件名
# 按 .fai 顺序逐个处理，保证最终 VCF 的 contig 顺序与参考一致
########################################

awk 'BEGIN{OFS="\t"} {
    safe=$1; gsub(/[^A-Za-z0-9._-]/,"_",safe)
    printf "%05d\t%s\t%s\t%s\n", NR, $1, $2, safe
}' "${REFERENCE}.fai" > "$INTERVAL_META"

TOTAL_INTERVALS=$(wc -l < "$INTERVAL_META")
[[ "$TOTAL_INTERVALS" -gt 0 ]] || die "No contigs in .fai"
log "Contigs: $TOTAL_INTERVALS"

########################################
# 处理单个 interval：GenomicsDBImport → GenotypeGVCFs
#
# 输入: 从 SAMPLE_MAP 读取所有 gVCF，按 interval 区域导入
#       到 GenomicsDB workspace，再联合分型输出 VCF
# 断点续跑: 若已有完整输出 VCF 则跳过，支持中断后重跑
########################################

run_interval() {
    local idx="$1" interval="$2" len="$3" safe="$4"
    local tag="${idx}_${safe}"
    local db_dir="${SCRATCH_DB}/${tag}.workspace"
    local out_vcf="${PER_INTERVAL_ROOT}/${tag}/${tag}.joint.vcf.gz"
    local tmp_dir="${SCRATCH_TMP}/${tag}"
    local ok_flag="${STATUS_DIR}/${tag}.ok"
    local log_file="${LOG_DIR}/${tag}.log"

    rm -f "$ok_flag" "${STATUS_DIR}/${tag}.fail"

    {
        log "[$tag] interval=$interval len=$len"
        mkdir -p "$(dirname "$out_vcf")" "$tmp_dir"

        # Resume check
        if [[ -s "$out_vcf" ]] && has_index "$out_vcf"; then
            log "[$tag] Already complete, skipping"
            touch "$ok_flag"
            return 0
        fi

        # Clean partial data
        [[ -d "$db_dir" ]] && rm -rf "$db_dir"
        rm -f "$out_vcf" "${out_vcf}.tbi" "${out_vcf}.idx" "${out_vcf}.csi"

        local extra=()
        [[ "$USE_SHARED_POSIXFS_OPTIMIZATIONS" == "yes" ]] && \
            extra+=(--genomicsdb-shared-posixfs-optimizations true)

        log "[$tag] GenomicsDBImport"
        gatk --java-options "-Xms${JAVA_XMX_IMPORT} -Xmx${JAVA_XMX_IMPORT} -Djava.io.tmpdir=${tmp_dir}" \
            GenomicsDBImport \
            -R "$REFERENCE" \
            --sample-name-map "$SAMPLE_MAP" \
            -L "$interval" \
            --genomicsdb-workspace-path "$db_dir" \
            --reader-threads "$READER_THREADS" \
            "${extra[@]}" \
            --tmp-dir "$tmp_dir" || return 1

        log "[$tag] GenotypeGVCFs"
        local geno_extra=()
        [[ "$INCLUDE_NON_VARIANT_SITES" == "yes" ]] && geno_extra+=(--include-non-variant-sites)
        [[ "$USE_SHARED_POSIXFS_OPTIMIZATIONS" == "yes" ]] && \
            geno_extra+=(--genomicsdb-shared-posixfs-optimizations true)

        gatk --java-options "-Xms${JAVA_XMX_GENOTYPE} -Xmx${JAVA_XMX_GENOTYPE} -Djava.io.tmpdir=${tmp_dir}" \
            GenotypeGVCFs \
            -R "$REFERENCE" \
            -V "gendb://${db_dir}" \
            -L "$interval" \
            "${geno_extra[@]}" \
            -O "$out_vcf" \
            --tmp-dir "$tmp_dir" || return 1

        has_index "$out_vcf" || gatk IndexFeatureFile -I "$out_vcf" || return 1

        # 清理 scratch 空间（workspace 很大，完成后立即释放）
        rm -rf "$db_dir" "$tmp_dir"
        touch "$ok_flag"
        log "[$tag] Done"

    } > >(tee -a "$log_file") 2>&1 || {
        touch "${STATUS_DIR}/${tag}.fail"
        log "[$tag] FAILED — see: $log_file"
        return 1
    }
}

########################################
# 主循环：逐个（或并发）处理所有 interval
#
# MAX_JOBS=1 时串行执行，>1 时用 & 后台并发 + jobs 节流
# 所有 interval 完成后检查 .fail 标记文件
#
# 注意: wait 前忽略 set -e，因为 bash 5.x 中 wait 可能
# 因非关键原因返回非零（如后台子进程的 pipe 退出码），
# 导致脚本在到达 .fail 检查之前就被终止。
# 后续的显式 .fail / .ok 检查才是正确的错误判断。
########################################

log "Processing $TOTAL_INTERVALS intervals (MAX_JOBS=$MAX_JOBS)"

if [[ "$MAX_JOBS" -le 1 ]]; then
    while IFS=$'\t' read -r idx interval len safe; do
        run_interval "$idx" "$interval" "$len" "$safe"
    done < "$INTERVAL_META"
else
    while IFS=$'\t' read -r idx interval len safe; do
        while [[ $(jobs -rp | wc -l) -ge "$MAX_JOBS" ]]; do sleep 5; done
        run_interval "$idx" "$interval" "$len" "$safe" &
    done < "$INTERVAL_META"
    wait || true   # 忽略 wait 退出码，改用下面的 .fail 文件判断
fi

# 检查是否有 interval 失败；有则列出并终止，不进入合并步骤
fail_flags=("$STATUS_DIR"/*.fail)
if [[ ${#fail_flags[@]} -gt 0 ]]; then
    log "Failed intervals:"
    for f in "${fail_flags[@]}"; do log "  $(basename "${f%.fail}")"; done
    die "Aborting — check logs above"
fi

ok_flags=("$STATUS_DIR"/*.ok)
log "OK: ${#ok_flags[@]}/$TOTAL_INTERVALS intervals"
[[ ${#ok_flags[@]} -eq "$TOTAL_INTERVALS" ]] || die "Mismatch in interval count"

########################################
# 合并所有 interval 的分片 VCF → 全基因组 joint VCF
########################################

run_merge

log "=== DONE ==="
log "Final VCF: $FINAL_VCF"
log "Log: $MASTER_LOG"

mv $FINAL_VCF ${OUTDIR}/${qunti_name}_${REF_NAME}.genome.joint.vcf.gz