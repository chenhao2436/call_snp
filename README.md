# call_snp

[English](README_EN.md) | **简体中文**

**call_snp** 是面向短读长测序数据的 SNP 变异检测流程，目标路线为：

```text
Paired-end FASTQ
      |
      v
Parabricks fq2bam
      |
      v
BAM
      |
      v
Parabricks HaplotypeCaller (--gvcf)
      |
      v
Per-sample GVCF
      |
      v
GATK GenomicsDBImport
      |
      v
GATK GenotypeGVCFs
      |
      v
GatherVcfs
      |
      v
Genome-wide joint-called VCF
```

<badge color="warning">Alpha：需要修复并完成验证</badge>

仓库目前包含两个 Bash 脚本：

1. `fq2gvcf_sample.sh`：使用 NVIDIA Parabricks 对双端 FASTQ 进行比对并生成单样本 GVCF。
2. `gvcf2vcf_genome.sh`：按参考基因组 contig 导入多个 GVCF，执行联合分型，并合并为全基因组 VCF。

两个阶段可以分开运行，例如在 GPU 服务器完成 FASTQ 到 GVCF，再在另一台计算服务器完成联合分型。

**当前状态：** 阶段一已完成一轮脚本修正，包括通用路径参数、参考索引生成、BAM 挂载路径、Parabricks read group 参数、GVCF 输出参数和输入检查。阶段二也修复了最终文件重复移动的问题，并允许自定义群体/参考标签。不过本次仍未在 GPU 服务器上执行完整端到端测试，建议先用少量样本验证后再进行正式批量分析。

## 1. 环境要求

### 阶段一：FASTQ → GVCF

- Linux 服务器。
- Bash。
- Docker Engine。
- NVIDIA GPU、驱动和 NVIDIA Container Toolkit。
- NVIDIA Clara Parabricks 镜像。当前脚本使用：
  `nvcr.io/nvidia/clara/clara-parabricks:4.1.1-1`。
- BWA：建立参考基因组索引。
- SAMtools：建立 FASTA 索引。
- GATK：建立参考序列字典。

### 阶段二：GVCF → 全基因组 VCF

- Linux 服务器，Bash 4 或更高版本。
- GATK 4：`GenomicsDBImport`、`GenotypeGVCFs`、`GatherVcfs`、`IndexFeatureFile`。
- 与 GATK 版本兼容的 Java。
- 常用系统命令：`awk`、`sed`、`tee`、`realpath`、`wc`、`date`、`mkdir`。

仓库的 `environment.yml` 提供宿主机 Conda 环境。Parabricks 通过 Docker 镜像运行，不由 pip 或 Conda 安装。NVIDIA GPU 驱动、Docker 和 Container Toolkit 也需要在宿主机单独配置。

## 2. 安装

克隆仓库：

```bash
git clone https://github.com/chenhao2436/call_snp.git
cd call_snp
```

创建 Conda 环境：

```bash
conda env create -f environment.yml
conda activate call-snp
```

检查软件：

```bash
bash --version
bwa 2>&1 | head
samtools --version
gatk --version
docker --version
```

检查 GPU 容器是否可用：

```bash
nvidia-smi
docker run --rm --gpus all \
  nvcr.io/nvidia/clara/clara-parabricks:4.1.1-1 nvidia-smi
```

Parabricks 镜像可能需要 NVIDIA NGC 权限或事先登录认证，请按 NVIDIA 当前的容器访问要求配置。

## 3. 输入文件

### 3.1 参考基因组

准备同一个版本的参考基因组及其索引：

```text
reference/
├── reference.fa
├── reference.fa.fai       # samtools faidx 生成
├── reference.dict        # gatk CreateSequenceDictionary 生成
└── BWA index files        # bwa index 生成
```

生成索引的正确命令：

```bash
bwa index reference.fa
samtools faidx reference.fa
gatk CreateSequenceDictionary -R reference.fa
```

同一次联合分型中的所有样本必须使用相同参考组装版本和 contig 命名方式。不要把不同版本的参考基因组混入同一个 sample map。

### 3.2 双端 FASTQ

当前阶段一脚本预期文件名格式：

```text
{sample}_1_clean.fq.gz
{sample}_2_clean.fq.gz
```

样本名称不要包含空格或 TAB。每个样本的 R1 和 R2 必须一一对应。

### 3.3 阶段一样本列表

每行一个样本 ID，例如：

```text
Parent_A
Parent_B
Pool_Long
Pool_Round
```

当前脚本没有可靠处理空行和注释行的逻辑，请仅使用干净的样本 ID 列表。

### 3.4 阶段二 GVCF sample map

需要一个两列、以 TAB 分隔的文件：

```text
Parent_A    /absolute/path/to/Parent_A.g.vcf.gz
Parent_B    /absolute/path/to/Parent_B.g.vcf.gz
Pool_Long   /absolute/path/to/Pool_Long.g.vcf.gz
Pool_Round  /absolute/path/to/Pool_Round.g.vcf.gz
```

样本名与文件绝对路径之间必须是真正的 TAB，而不是空格。每个路径必须指向存在、可读且已经建立索引的 GVCF；样本名不能重复。

在 GVCF 所在目录生成 sample map：

```bash
for f in *.g.vcf.gz; do
    sample=$(basename "$f" .g.vcf.gz)
    printf '%s\\t%s\\n' "$sample" "$(realpath "$f")"
done > gvcf_sample.list
```

运行联合分型前，先检查生成的列表和样本顺序。

## 4. 阶段一：FASTQ → GVCF

### 4.1 使用命令行传入路径

当前脚本已改为显式参数，不再要求修改脚本中的服务器绝对路径。参数顺序如下：

| 参数 | 说明 |
|---|---|
| `SAMPLE_LIST` | 每行一个样本 ID，可使用空行和 `#` 注释行 |
| `REFERENCE.fa` | 参考基因组 FASTA |
| `FASTQ_DIR` | FASTQ 文件所在目录 |
| `BAM_DIR` | BAM 输出目录 |
| `GVCF_DIR` | GVCF 和 sample map 输出目录 |
| `PARABRICKS_IMAGE` | 可选，Parabricks Docker 镜像，默认 `nvcr.io/nvidia/clara/clara-parabricks:4.1.1-1` |

样本 ID 仅允许字母、数字、点、下划线和连字符。脚本会检查重复样本、FASTQ 配对文件、必要软件、GPU、Docker 镜像和输出文件。

Parabricks 镜像需要预先拉取到本机；若需要认证，请先登录 NVIDIA NGC：

```bash
docker pull nvcr.io/nvidia/clara/clara-parabricks:4.1.1-1
```

### 4.2 运行命令

**只有在修复本 README 所列问题后，才运行：**

```bash
bash fq2gvcf_sample.sh \
  /path/to/sample.list \
  /path/to/reference.fa \
  /path/to/fastq \
  /path/to/bam \
  /path/to/gvcf
```

预期每个样本执行：

- Parabricks `fq2bam`：双端 FASTQ → 比对后的 BAM。
- Parabricks `haplotypecaller --gvcf`：BAM → 单样本 GVCF。

预期文件结构：

```text
bam/
├── Parent_A.bam
└── Parent_B.bam

gvcf/
├── Parent_A.g.vcf.gz
├── Parent_B.g.vcf.gz
└── gvcf_sample.list
```

具体路径由脚本配置决定。

## 5. 阶段二：GVCF → 全基因组联合分型 VCF

### 5.1 检查 GVCF

每个 GVCF 必须完整、可读、基于同一参考基因组，并且具有索引。当前脚本识别 `.tbi`、`.idx` 或 `.csi`。

缺少索引时可使用：

```bash
gatk IndexFeatureFile -I sample.g.vcf.gz
```

### 5.2 运行联合分型

脚本支持三个位置参数：

```bash
bash gvcf2vcf_genome.sh [gvcf_list] [reference.fa] [output_directory]
```

示例：

```bash
bash gvcf2vcf_genome.sh \
  /data/project/gvcf/gvcf_sample.list \
  /data/reference/reference.fa \
  /data/project/joint_vcf
```

如不提供参数，脚本会使用配置区中的默认路径。为了可移植性，建议显式传入三个参数。

### 5.3 资源参数

`gvcf2vcf_genome.sh` 中当前主要参数如下：

| 参数 | 含义 | 当前默认值 |
|---|---|---:|
| `MAX_JOBS` | 并行处理的 contig 数量 | 3 |
| `JAVA_XMX_IMPORT` | GenomicsDBImport 最大堆内存 | 20g |
| `JAVA_XMX_GENOTYPE` | GenotypeGVCFs 最大堆内存 | 28g |
| `JAVA_XMX_GATHER` | GatherVcfs 最大堆内存 | 12g |
| `READER_THREADS` | GVCF reader 线程数 | 4 |
| `INCLUDE_NON_VARIANT_SITES` | 是否包含非变异位点 | no |

根据样本量、可用内存、CPU 和存储性能调整。上述数值只是当前脚本默认配置，并不是适用于所有服务器的最佳设置。可通过环境变量 `QUANTI_NAME` 和 `REFERENCE_NAME` 修改最终 VCF 文件名前缀，例如 `QUANTI_NAME=F3 REFERENCE_NAME=Zunla bash gvcf2vcf_genome.sh list.txt ref.fa outdir`。

### 5.4 输出文件

当前默认输出文件名由脚本中的 `qunti_name` 和 `REF_NAME` 拼接：

```text
F2_3_PP64.genome.joint.vcf.gz
```

输出目录还包含总日志、每个 contig 的日志、分片 VCF、状态标记、临时目录、sample map 和 interval 元信息。

## 6. 输出目录结构

```text
output_directory/
├── logs/
│   ├── *.run_*.log
│   └── *.log
├── status/
│   ├── *.ok
│   └── *.fail
├── per_interval/
│   └── <interval>/
│       └── <interval>.joint.vcf.gz
├── genomicsdb/
├── tmp/
├── *.sample_map
├── *.intervals.tsv
├── *.shards.list
└── <population>_<reference>.genome.joint.vcf.gz
```

阶段二按参考 FASTA 的 `.fai` 顺序逐个处理 contig，并使用状态文件跳过已完成分片。重跑前仍应检查日志和状态文件，不能只根据某个分片文件存在就认定任务完整成功。

## 7. 资源规划与质量控制

大样本 GVCF 联合分型会消耗大量内存、临时空间和存储 I/O。运行前：

- 估算 GVCF 总体积和临时磁盘空间。
- 确认输出文件系统剩余空间。
- 初次运行将 `MAX_JOBS` 设为较小值，监测资源后再提高并发。
- 尽可能使用本地高速盘存放临时 GenomicsDB workspace。
- 避免在同一共享存储上同时运行多个大型联合分型任务。

脚本优先尝试使用 `/dev/shm` 作为 scratch，不可用时回退到输出目录的 `tmp`。使用前请确认空间足够。

每次分析建议记录：

- 仓库 Git commit 或 release tag。
- 参考基因组名称、版本和校验和。
- Parabricks 镜像 tag 和 digest。
- GATK、BWA、SAMtools、Java、Docker 和 GPU 驱动版本。
- 样本列表、sample map。
- 关键资源参数和完整日志。

下游分析前至少检查：

1. 预期的每个样本都有非空 GVCF。
2. 每个 GVCF 都有有效索引。
3. sample map 中样本名唯一且数量正确。
4. 最终 VCF 非空、可读且有索引。
5. VCF header 中包含预期样本列。
6. contig 名称与顺序和参考基因组一致。
7. 样本数和变异统计符合预期。

当前仓库只负责变异检测和联合分型，尚未实现变异过滤、SNP-only 提取、BSA、GWAS、功能注释或表达分析。

## 8. 已知问题与验证状态

本次产品化过程中已经修复以下源码问题：

| 文件 | 已完成修改 |
|---|---|
| `fq2gvcf_sample.sh` | 删除会被 Bash 当作命令执行的裸文本；改用 `samtools faidx`；修正 BAM 目录挂载；修正 Parabricks read group 和 `--out-variants` 参数；改为显式路径参数；增加输入、GPU、工具和输出检查 |
| `gvcf2vcf_genome.sh` | 删除将最终 VCF 移动到自身路径的冗余 `mv`；群体名和参考标签可通过 `QUANTI_NAME`、`REFERENCE_NAME` 环境变量设置 |

仍需完成的工作：

| 优先级 | 文件 | 未完成事项 |
|---|---|---|
| P1 | `fq2gvcf_sample.sh` | 在实际 Parabricks 4.1.1 GPU 环境进行小样本 dry-run / smoke test 和真实数据验证 |
| P1 | `gvcf2vcf_genome.sh` | 在目标 GATK 版本和目标服务器上验证 GenomicsDBImport、并发、断点重跑和 GatherVcfs |
| P1 | 两阶段 | 增加可公开分发的小型测试数据和自动化回归测试 |
| P2 | 两阶段 | 进一步统一配置文件、版本记录和运行报告 |

**验证声明：** 当前修改基于源码审查及 Parabricks 官方命令参数文档；本轮没有在配备 GPU 的服务器上完成端到端变异检测，也没有进行生物学结果验证。不要将静态修复等同于已通过实际运行测试。

## 9. 故障排查

### Docker 无法访问 GPU

```bash
nvidia-smi
docker run --rm --gpus all <parabricks-image> nvidia-smi
```

### GATK 提示参考基因组索引或字典缺失

```bash
samtools faidx reference.fa
gatk CreateSequenceDictionary -R reference.fa
```

### GATK 提示 GVCF 索引缺失

```bash
gatk IndexFeatureFile -I sample.g.vcf.gz
```

### 某条 contig 联合分型失败

检查 `logs/` 中对应日志和 `status/` 下的 `.fail` 文件；确认 GVCF 可读、索引完整、参考基因组一致，并检查 scratch 空间和 Java 内存设置。

### 找不到最终 VCF

只有所有 contig 成功完成后才会进入合并步骤。检查总日志、每个 interval 日志和状态目录，不要把部分分片当作最终结果。

## 10. 开发路线

- [x] FASTQ-to-GVCF 初始脚本（需要修复和验证）。
- [x] GVCF-to-全基因组 VCF 初始脚本（需要验证）。
- [ ] 将绝对路径替换为统一配置接口。
- [ ] 两阶段增加统一的输入预检和日志。
- [ ] 增加小型可复现测试数据。
- [ ] 增加版本标签与发布记录。
- [ ] 接入 BSA 可视化与交互式绘图模块。
- [ ] 接入转录组基因 ID 转换和绘图脚本。
- [ ] 接入变异-表达整合与功能注释模块。

后续模块属于规划，不代表当前仓库已实现这些功能。

## 引用与反馈

使用本流程时，请引用仓库地址和具体 Git commit，并注明实际使用的软件版本。

- Repository: https://github.com/chenhao2436/call_snp
- NVIDIA Parabricks: https://docs.nvidia.com/clara/parabricks/
- GATK: https://gatk.broadinstitute.org/
- BWA: https://github.com/lh3/bwa
- SAMtools: https://www.htslib.org/

提交 Issue 时请提供运行命令、日志片段、软件版本和输入目录结构。不要在公开 Issue 上传人类基因组数据、访问凭证或敏感文件。
