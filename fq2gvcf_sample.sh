#!/usr/bin/env bash
set -Eeuo pipefail

usage() {
    cat <<'EOF'
Usage:
  bash fq2gvcf_sample.sh SAMPLE_LIST REFERENCE.fa FASTQ_DIR BAM_DIR GVCF_DIR [PARABRICKS_IMAGE]

Required arguments:
  SAMPLE_LIST       One sample ID per line
  REFERENCE.fa      Reference genome FASTA
  FASTQ_DIR         Directory containing {sample}_1_clean.fq.gz and {sample}_2_clean.fq.gz
  BAM_DIR           Output directory for BAM files
  GVCF_DIR          Output directory for GVCF files and gvcf_sample.list

Optional:
  PARABRICKS_IMAGE  Default: nvcr.io/nvidia/clara/clara-parabricks:4.1.1-1

Example:
  bash fq2gvcf_sample.sh sample.list reference.fa ./fastq ./bam ./gvcf
EOF
}

die() {
    echo "[ERROR] $*" >&2
    exit 1
}

need_cmd() {
    command -v "$1" >/dev/null 2>&1 || die "Required command not found: $1"
}

[[ $# -ge 5 && $# -le 6 ]] || { usage >&2; exit 2; }

SAMPLE_LIST="$1"
REFERENCE="$2"
FASTQ_DIR="$3"
BAM_DIR="$4"
GVCF_DIR="$5"
PARABRICKS_IMAGE="${6:-nvcr.io/nvidia/clara/clara-parabricks:4.1.1-1}"

[[ -s "$SAMPLE_LIST" ]] || die "Sample list is missing or empty: $SAMPLE_LIST"
[[ -s "$REFERENCE" ]] || die "Reference FASTA is missing or empty: $REFERENCE"
[[ -d "$FASTQ_DIR" ]] || die "FASTQ directory does not exist: $FASTQ_DIR"

need_cmd docker
need_cmd nvidia-smi
need_cmd bwa
need_cmd samtools
need_cmd gatk
need_cmd realpath

REFERENCE="$(realpath "$REFERENCE")"
FASTQ_DIR="$(realpath "$FASTQ_DIR")"
mkdir -p "$BAM_DIR" "$GVCF_DIR"
BAM_DIR="$(realpath "$BAM_DIR")"
GVCF_DIR="$(realpath "$GVCF_DIR")"
REF_DIR="$(dirname "$REFERENCE")"
REF_BASE="$(basename "$REFERENCE")"
REF_STEM="${REFERENCE%.*}"
REF_DICT="${REF_STEM}.dict"

# Read sample IDs safely: one non-empty ID per line, no whitespace or path separators.
samples=()
declare -A seen=()
while IFS= read -r sample || [[ -n "$sample" ]]; do
    [[ -z "$sample" || "$sample" == #* ]] && continue
    [[ "$sample" =~ ^[A-Za-z0-9._-]+$ ]] || die "Invalid sample ID '$sample'. Use letters, digits, dot, underscore, or hyphen only."
    [[ -z "${seen[$sample]:-}" ]] || die "Duplicate sample ID: $sample"
    seen["$sample"]=1
    samples+=("$sample")
done < "$SAMPLE_LIST"
[[ ${#samples[@]} -gt 0 ]] || die "No valid sample IDs found in $SAMPLE_LIST"

echo "[INFO] Preparing reference indexes..."
[[ -s "${REFERENCE}.fai" ]] || samtools faidx "$REFERENCE"
[[ -s "$REF_DICT" ]] || gatk CreateSequenceDictionary -R "$REFERENCE"
[[ -s "${REF_STEM}.bwt" ]] || bwa index "$REFERENCE"

echo "[INFO] Checking GPU and Parabricks container..."
nvidia-smi >/dev/null
docker image inspect "$PARABRICKS_IMAGE" >/dev/null 2>&1 || \
    die "Parabricks image is not available locally: $PARABRICKS_IMAGE. Authenticate to NVIDIA NGC and pull the image first."

SAMPLE_MAP_TMP="${GVCF_DIR}/gvcf_sample.list.tmp"
SAMPLE_MAP="${GVCF_DIR}/gvcf_sample.list"
: > "$SAMPLE_MAP_TMP"

for sample in "${samples[@]}"; do
    fq1="${FASTQ_DIR}/${sample}_1_clean.fq.gz"
    fq2="${FASTQ_DIR}/${sample}_2_clean.fq.gz"
    bam_file="${BAM_DIR}/${sample}.bam"
    gvcf_file="${GVCF_DIR}/${sample}.g.vcf.gz"

    echo
    echo "[INFO] Sample: $sample"

    [[ -s "$fq1" && -s "$fq2" ]] || die "FASTQ pair is incomplete for $sample: $fq1 / $fq2"

    if [[ -s "$gvcf_file" ]]; then
        echo "[SKIP] Existing GVCF found: $gvcf_file"
        if [[ ! -s "${gvcf_file}.tbi" && ! -s "${gvcf_file}.idx" && ! -s "${gvcf_file}.csi" ]]; then
            echo "[INFO] Creating GVCF index..."
            gatk IndexFeatureFile -I "$gvcf_file"
        fi
    else
        if [[ -s "$bam_file" ]]; then
            samtools quickcheck -v "$bam_file" || die "Existing BAM failed samtools quickcheck: $bam_file"
            echo "[SKIP] Existing BAM passed quickcheck: $bam_file"
        else
            echo "[RUN] Parabricks fq2bam: $sample"
            docker run --gpus all --rm \
                --volume "$FASTQ_DIR:/input:ro" \
                --volume "$REF_DIR:/reference:ro" \
                --volume "$BAM_DIR:/bam" \
                "$PARABRICKS_IMAGE" pbrun fq2bam \
                --ref "/reference/$REF_BASE" \
                --in-fq "/input/${sample}_1_clean.fq.gz" "/input/${sample}_2_clean.fq.gz" "@RG\tID:${sample}\tLB:${sample}\tPL:ILLUMINA\tSM:${sample}\tPU:${sample}" \
                --out-bam "/bam/${sample}.bam"
            [[ -s "$bam_file" ]] || die "fq2bam did not create a non-empty BAM: $bam_file"
            samtools quickcheck -v "$bam_file" || die "BAM failed samtools quickcheck: $bam_file"
        fi

        echo "[RUN] Parabricks HaplotypeCaller: $sample"
        docker run --gpus all --rm \
            --volume "$REF_DIR:/reference:ro" \
            --volume "$BAM_DIR:/bam:ro" \
            --volume "$GVCF_DIR:/gvcf" \
            "$PARABRICKS_IMAGE" pbrun haplotypecaller \
            --ref "/reference/$REF_BASE" \
            --in-bam "/bam/${sample}.bam" \
            --gvcf \
            --out-variants "/gvcf/${sample}.g.vcf.gz"
        [[ -s "$gvcf_file" ]] || die "HaplotypeCaller did not create a non-empty GVCF: $gvcf_file"

        if [[ ! -s "${gvcf_file}.tbi" && ! -s "${gvcf_file}.idx" && ! -s "${gvcf_file}.csi" ]]; then
            echo "[INFO] Creating GVCF index..."
            gatk IndexFeatureFile -I "$gvcf_file"
        fi
    fi

    [[ -s "${gvcf_file}.tbi" || -s "${gvcf_file}.idx" || -s "${gvcf_file}.csi" ]] || die "GVCF index is missing: $gvcf_file"
    printf '%s\\t%s\\n' "$sample" "$gvcf_file" >> "$SAMPLE_MAP_TMP"
    echo "[DONE] $sample"
done

mv "$SAMPLE_MAP_TMP" "$SAMPLE_MAP"
echo
echo "[DONE] FASTQ-to-GVCF stage completed."
echo "[INFO] Sample map: $SAMPLE_MAP"
