import re

# 您提供的目标基因列表
genes_of_interest = {
    "BRD4", "USP7", "RUNX1", "PRKDC", "TNPO1", "TACC1", 
    "TRIM35", "UBTF", "MIR155HG", "MLLT1", "SH3GL1", 
    "ATP6V1C1", "HRH2", "WT1"
}

# 文件路径
gtf_file = "/md01/nieyg/ref/genome/hg19/gencode.v32lift37.annotation.gtf"
bed_file = "promoters.bed"
dist = 1000  # TSS 上下游扩展距离 (1000bp = 1kb)

with open(gtf_file, 'r') as f_in, open(bed_file, 'w') as f_out:
    for line in f_in:
        if line.startswith('#'): 
            continue
        parts = line.strip().split('\t')
        
        # 只看 gene 级别的注释来确定 TSS
        if parts[2] != 'gene': 
            continue
        
        info = parts[8]
        # 提取 gene_name
        match = re.search(r'gene_name "(.*?)"', info)
        if match:
            gene_name = match.group(1)
            if gene_name in genes_of_interest:
                chrom = parts[0]
                strand = parts[6]
                start = int(parts[3])
                end = int(parts[4])
                
                # 判断正负链以确定真正的转录起始位点 (TSS)
                if strand == '+':
                    tss = start
                else:
                    tss = end
                    
                # 计算启动子区域 (注意防止坐标出现负数)
                p_start = max(0, tss - dist)
                p_end = tss + dist
                
                # 输出标准的 6 列 BED 格式: chrom, start, end, name, score, strand
                f_out.write(f"{chrom}\t{p_start}\t{p_end}\t{gene_name}\t.\t{strand}\n")

print(f"成功提取 {len(genes_of_interest)} 个基因的启动子坐标至 {bed_file}")


findMotifsGenome.pl promoters.bed  hg19 PRKDC_relatedgene_motifDir_Homer -len 8,10,12  

# 1. 加载所需的包
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(ggpubr) # 用于添加显著性 P 值

# 2. 读取数据 (假定文件路径依然是 normalized_counts.txt)
counts <- read.table("normalized_counts.txt", header = TRUE, row.names = 1, check.names = TRUE)

target_genes <- c("ERG", "ETS1", "ELF1")
genes_present <- intersect(target_genes, rownames(counts))

if(length(genes_present) == 0) {
  stop("在文件中没有找到 ERG, ETS1, 或 ELF1 基因。请检查您的数据。")
}

filtered_counts <- counts[genes_present, , drop = FALSE]

# 3. 数据格式转换
filtered_counts$Gene <- rownames(filtered_counts)
long_data <- pivot_longer(filtered_counts, 
                          cols = -Gene, 
                          names_to = "Sample", 
                          values_to = "Expression")

# 4. 解析样本信息并强制指定顺序 (因子化)
long_data <- long_data %>%
  mutate(
    CellLine = str_split_i(Sample, "\\.", 1),
    Treatment = str_split_i(Sample, "\\.", 2)
  ) %>%
  mutate(
    # 强制指定组别顺序：U 放前面，T 放后面
    Treatment = factor(Treatment, levels = c("U", "T")),
    # 强制指定细胞系顺序：U937 -> MV411 -> MOLM13
    CellLine = factor(CellLine, levels = c("U937", "MV411", "MOLM13"))
  )

# 5. 计算均值和误差棒
summary_data <- long_data %>%
  group_by(Gene, CellLine, Treatment) %>%
  summarise(
    Mean = mean(Expression),
    SD = sd(Expression),
    N = n(),
    SE = SD / sqrt(N),
    .groups = 'drop'
  )

# 6. 绘图
p <- ggplot() +
  # A. 绘制柱状图 (使用均值数据)
  geom_bar(data = summary_data, aes(x = CellLine, y = Mean, fill = Treatment),
           stat = "identity", position = position_dodge(0.8), width = 0.7, color = "black") +
           
  # B. 添加误差棒 (使用均值数据)
  geom_errorbar(data = summary_data, aes(x = CellLine, ymin = Mean - SE, ymax = Mean + SE, group = Treatment), 
                position = position_dodge(0.8), width = 0.25) +
                
  # C. 添加独立样本点 (使用原始长数据 long_data)
  geom_point(data = long_data, aes(x = CellLine, y = Expression, fill = Treatment),
             position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
             shape = 21, color = "black", size = 2, alpha = 0.7) +
             
  # 按基因分面
  facet_wrap(~ Gene, scales = "free_y") +
  
  # D. 添加显著性 P 值 (使用 t.test，针对同一细胞系内的 U 和 T 进行比较)
  stat_compare_means(data = long_data, aes(x = CellLine, y = Expression, group = Treatment),
                     method = "t.test", label = "p.format") +
                     
  # 设置主题与颜色
  theme_classic() +
  scale_fill_manual(values = c("U" = "#56B4E9", "T" = "#E69F00")) + 
  labs(
    title = "",
    x = "",
    y = "Normalized Counts",
    fill = "Group"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    strip.text = element_text(face = "bold", size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold")
  )

# 显示图表
print(p)

ggsave("Gene_Expression_Barplot.pdf", plot = p, width = 8, height = 5)

# -s 表示考虑正负链方向进行提取
# -name 表示将 bed 文件中的 gene_name 作为 fasta 的序列 ID
grep 'gene_name "PRKDC"' /md01/nieyg/ref/genome/hg19/gencode.v32lift37.annotation.gtf | awk -F'\t' 'BEGIN{OFS="\t"} $3=="gene" {tss=($7=="+")?$4:$5; print $1, tss-1000, tss+1000, "PRKDC", ".", $7}' > prkdc.bed

bedtools getfasta -fi /md01/nieyg/ref/genome/hg19/hg19.fa -bed  prkdc.bed -s -name -fo  prkdc_promoters.fasta


fimo --oc ETS1_in_PRKDC_promoter --thresh 1.0E-2 MA0098.1.meme prkdc_promoters.fasta
fimo --oc ELF1_in_PRKDC_promoter --thresh 1.0E-2 MA0473.1.meme prkdc_promoters.fasta
fimo --oc ERG_in_PRKDC_promoter --thresh 1.0E-2 MA0474.2.meme prkdc_promoters.fasta








