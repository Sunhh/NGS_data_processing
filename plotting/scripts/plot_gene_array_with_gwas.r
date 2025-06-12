#!/usr/bin/env Rscript

# 加载所需库
suppressPackageStartupMessages({
  library(optparse)
  library(readxl)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(gggenes)
  library(patchwork)
})

# 命令行参数
option_list <- list(
  make_option(c("-g", "--gene_tsv"), type = "character", help = "Path to gene list file (TSV)"),
  make_option(c("-w", "--gwas_tsv"), type = "character", help = "Path to GWAS/SV file (TSV)"),
  make_option(c("-r", "--region"), type = "character", help = "Genomic region: chr:start-end"),
	make_option(c("-c", "--cutP"), type = "numeric", default = 6, help = "-log10(P) for GWAS cutoff"),
  make_option(c("-o", "--output"), type = "character", default = "gwas_gene_plot.pdf", help = "Output file name [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))


if (0) {
 # Test run
 opt <- list(
 gene_tsv = 'gene_list.tsv',
 gwas_tsv = 'gwas_Pvalues.tsv',
 cutP     = 5.96657624451305,
 region    = 'CLV01_Chr03:31445138-31699695',
 output    = 'test_plot.pdf'
 )
} else {
  if (is.null(opt$gene_tsv) || is.null(opt$gwas_tsv) || is.null(opt$region)) {
    cat("\nMissing required arguments.\n\n")
    print_help(OptionParser(option_list = option_list))
    quit(status = 1)
  }
  if (is.null(opt$cutP)) {opt$cutP <- 6}
  if (is.null(opt$output)) {opt$output <- 'test_plot.pdf'}
}
 

gene_tsv <- opt$gene_tsv
gwas_tsv <- opt$gwas_tsv
cutP  <- opt$cutP
region <- opt$region
output_file <- opt$output

# Example gene_tsv file:
#   geneID              chrID        start      end       strand  highlight   ahrd
#   CLV01C03G016270.1   CLV01_Chr03  31445138   31445434  -                   protein RALF-like 9
#   CLV01C03G016280.1   CLV01_Chr03  31451332   31454875  +       Y           calmodulin-binding protein 60 A-like

# Example gwas_tsv file:
#   chr          pos        SV.length   P.value    gwas.ID
#   CLV01_Chr03  31521699   -48         8.54E-08   3_31521699_v14414099_48
#   CLV01_Chr04  9294168    .           1.51E-07   4_9294168_1



# Parse target region
region_parts <- strsplit(region, "[:-]")[[1]]
tgt_chr <- region_parts[1]
tgt_start <- as.numeric(region_parts[2])
tgt_end <- as.numeric(region_parts[3])

# Load gene list.
### genes <- read_excel(gene_file)
genes <- read.delim(gene_tsv, stringsAsFactors = FALSE)
### Subset genes.
genes <- genes %>%
  filter(chrID == tgt_chr, start <= tgt_end, end >= tgt_start) %>%
  mutate(
    use_strand = ifelse(strand == "+", "forward", "reverse"),
    use_fill   = ifelse(is.na(highlight), "other", ifelse(highlight == "Y", "highlighted", "other"))
  ) %>%
  rename(use_gene = geneID, use_start = start, use_end = end)

# Load in GWAS P values.
gwas <- if (grepl("\\.tsv$", gwas_tsv)) {
  read.delim(gwas_tsv, stringsAsFactors = FALSE)
} else {
  read.csv(gwas_tsv, stringsAsFactors = FALSE)
}

gwas <- gwas %>%
  filter(chr == tgt_chr, pos >= tgt_start, pos <= tgt_end, !is.na(P.value)) %>%
  mutate(
    use_logP = -log10(as.numeric(P.value)),
    use_VARtype = ifelse(SV.length == "." | is.na(SV.length) | SV.length == "" , "SNP", "SV")
  )

# 曼哈顿图
manhattan <- ggplot() +
  geom_point(data = gwas %>% filter(use_VARtype == "SNP"),
             aes(x = pos / 1e6, y = use_logP),
             color = "#999999", shape = 16, size = 5) +
  geom_point(data = gwas %>% filter(use_VARtype == "SV"),
             aes(x = pos / 1e6, y = use_logP),
             color = "#ff0000", shape = 2, size = 5) +
  geom_hline(yintercept = cutP, linetype = "dashed", color = "green") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme_minimal() +
  labs(title = "", y = expression(-log[10](italic(P))), x = NULL) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.line.x = element_line(),
		axis.ticks.x = element_line(),
    axis.line.y = element_line(),
		axis.ticks.y = element_line(),
    legend.position = "none"
  )

# print(manhattan)

# 基因排列图
### Set Yxx to distinguish forward/reverse genes.
genes$use_ymin <- ifelse(genes$use_strand == "forward", 1, 0.7);
genes$use_ymax <- ifelse(genes$use_strand == "forward", 1.3, 1);

gene_plot <- ggplot(genes) +
  geom_rect(aes(
    xmin = use_start / 1e6,
    xmax = use_end / 1e6,
    ymin = use_ymin,
    ymax = use_ymax,
    fill = use_fill), color = "black") +
  annotate("segment",
    x    = min(genes$use_start) / 1e6,
    xend = max(genes$use_end) / 1e6,
    y = 1,
    yend = 1,
    color = "black",
    linewidth = 0.5) +
  scale_fill_manual(values = c("highlighted" = "#ffd1f0", "other" = "#dcebbc")) +
  coord_cartesian(ylim = c(0.6, 1.4)) + # 控制 y 范围
  theme_void() +
  theme(
    panel.grid.major = element_blank(),  # 去除主网格线
    panel.grid.minor = element_blank(),  # 去除次网格线
    panel.border = element_blank(),      # 去除边框（避免干扰视觉）
    legend.position = "none"
  )

# print(gene_plot)
# geom_gene_arrow(arrowhead_height = unit(4, "mm"), arrowhead_width = unit(2, "mm"), arrow_body_height = unit(5, "mm")) +



# 合并图形
final_plot <- manhattan / gene_plot + plot_layout(heights = c(3, 1)) ; # print(final_plot)

# 保存图形
ggsave(output_file, final_plot, width = 180, height = 90, units = "mm", dpi= 300, device = cairo_pdf)
