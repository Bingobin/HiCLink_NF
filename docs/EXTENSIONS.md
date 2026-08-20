# HiCLink_NF 优化与扩展路线

## 当前版本已经解决的问题

1. 移除 reference、chromosome sizes、hg38、10 kb 与 Mustache 路径的硬编码。
2. 用统一 samplesheet 区分 Hi-C/HiChIP，并保留 condition、replicate 和 control。
3. 把 reads、pairs、matrix、feature calling 与报告拆成 DSL2 process。
4. pairtools 显式完成 parse、sort、dedup 和 `UU` 筛选，输出去重及 contact QC。
5. `.mcool` 一次生成多分辨率并平衡；下游参数检查防止调用不存在的分辨率。
6. HiChIP 可以优先使用外部 peak；没有 peak 时才从去重、MAPQ 过滤后的 BAM 调 MACS2。
7. 运行自动保存 validated samplesheet、DAG、trace、timeline、execution report 和 MultiQC。
8. HiChIP 默认不调用全局 compartment/TAD，避免把抗体富集偏倚当作无偏三维结构；需要时必须显式开启。

## 优化建议

### P0：正式数据运行前

- 用一组已有结果的样本做回归验证：valid pairs 数、cis/trans 比、duplicate rate、每个分辨率的 contacts 和稀疏度应与旧流程在可解释范围内一致。
- 为实验协议确认 `walks_policy`、MAPQ、restriction-fragment 处理和 read trimming。当前默认值适合常规 paired-end Hi-C/HiChIP，但不能替代协议核对。
- 根据集群节点和样本深度重新标定 `large/xlarge` 的 CPU、内存、时间；pairtools sort 和 cooler zoomify 是主要 I/O/临时空间热点。
- 固定实际验证过的 Conda lock file 或构建单一版本化容器；当前 YAML 使用最低版本范围，便于首次求解，但不是严格锁版本环境。
- FitHiChIP 官方镜像先用小样本验证挂载路径和输出版本，再提交整批数据。

### P1：性能和稳定性

- 加入按 FASTQ 大小自动估算资源的 dynamic directive，并为 exit 137/143 递增内存或时间。
- 在大数据模式下把 BWA → pairtools 融合为可选 streaming process，减少中间 BAM 写盘；HiChIP 仍单独保留 peak BAM 分支。
- 支持 lane-level 输入和 `pairtools merge`，先按 lane QC，再合并到 biological sample。
- 为 scratch 目录、pairtools sort 临时目录和 publish mode 提供站点 profile，避免共享文件系统成为瓶颈。
- 生成工具版本表和参数快照，并对 reference FASTA/chrom sizes/phasing track 记录 checksum。

## 可扩展分析模块

### 1. 生物学重复与一致性

- replicate-level depth、cis/trans、distance decay、library complexity；
- HiCRep、GenomeDISCO 或 QuASAR-Rep 一致性；
- 在 QC 合格后合并 pairs/matrix，同时保留每个重复的独立结果；
- downsampling 到相同有效 contacts 后再比较 feature recovery。

这是差异分析前最高优先级模块。没有生物学重复时，差异 loop/contact 只能作为描述性结果。

### 2. 差异 contact、TAD、compartment 和 loop

- replicate-aware differential contacts：diffHic、multiHiCcompare 或基于 count matrix 的 edgeR/DESeq2 设计；
- compartment switching：统一 eigenvector 方向后比较 E1，并输出 A→B/B→A；
- TAD boundary gain/loss 与 insulation strength 变化；
- HiChIP differential loops：构建 union anchors/loops、提取每个重复 counts、做条件设计和 FDR 校正；
- 可保留旧流程的 differential Mustache，作为无重复或探索性补充，而不是替代 replicate-aware 统计。

### 3. APA/pileup 与 feature-level QC

- 对 loops、TAD corners、CTCF sites、promoter-enhancer anchors 做 cooltools pileup；
- 报告 APA score、中心富集及 distance-matched background；
- 比较 peak-overlap、loop anchor enrichment 和重复间 feature recovery。

### 4. 基因与调控注释

- loop anchors 交集 promoter、enhancer、CTCF、ATAC/ChIP peaks；
- 输出 promoter-enhancer、enhancer-enhancer、promoter-promoter 类别；
- 基因分配必须保留直接 promoter contact 与最近基因推断的区别；
- 可接 RNA-seq DEG/GSEA，但不把结构关联直接解释成表达调控因果。

### 5. 可视化与交换格式

- `.hic` 导出用于 Juicebox；
- HiGlass/Resgen 浏览配置和 track hub；
- bigWig/bedGraph 导出 coverage、E1、insulation；
- 自动生成每个样本的标准 chromosome snapshot 和比较图。

### 6. 特殊结构模块

- stripe calling（例如 Stripenn）；
- structural variation、copy-number-aware normalization、neo-loop；
- allele-specific/haplotype-resolved contacts；
- trans contact、inter-chromosomal hubs；
- Micro-C、PLAC-seq、Capture-C 通过 assay adapter 扩展，而不是在现有 `hic/hichip` 分支中硬塞协议差异。

## 建议的版本节奏

- `0.1`：当前统一单样本处理、基础结构特征与 HiChIP 分支；
- `0.2`：lane/replicate merge、深度归一、一致性和回归测试；
- `0.3`：replicate-aware differential contacts/loops 与 APA；
- `0.4`：调控注释、浏览器交付和可复现容器锁定；
- `1.0`：真实 Hi-C/HiChIP 基准数据集、CI smoke test、集群 profile 和版本化 schema 全部稳定。
