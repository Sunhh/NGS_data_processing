#!/usr/bin/env Rscript
# RNA-seq volcano plot, Nature 系列期刊排版标准：
#   - 尺寸：单栏 89 mm / 1.5 栏 120 mm / 双栏 183 mm
#   - 字体：Arial（默认，PDF 嵌入 Arial）或 Helvetica，字号 5-7 pt
#   - 线宽 0.25-0.5 pt；PDF 为矢量输出，位图输出 600 dpi
# 用法示例：
#   Rscript volcano_plot.r -i avgExp_fdr.tsv \
#     --control 2 --treatment 3 --fdr 5 -o volcano_CPMMV.pdf

suppressPackageStartupMessages({
  library(argparse);
  library(data.table);
  library(ggplot2);
});

# 系统中某字体文件所在目录（供 ghostscript 嵌入字体）
find_font_dir <- function(family) {
  p <- character(0);
  if (requireNamespace("systemfonts", quietly = TRUE)) {
    sf <- systemfonts::system_fonts();
    p <- sf$path[sf$family == family];
  } else if (nzchar(Sys.which("fc-list"))) {
    p <- sub(":.*$", "", system2("fc-list", c(shQuote(family), "file"), stdout = TRUE));
  }
  p <- p[file.exists(p)];
  if (length(p) == 0) NULL else unique(dirname(p));
};

# 在 R 的 pdf 设备中登记 Arial：字宽与 Helvetica 相同，借用 Helvetica 的 AFM，只把字体名改为 ArialMT 等，
# 输出后再由 ghostscript 从系统字体目录嵌入真正的 Arial
register_arial <- function() {
  if ("Arial" %in% names(pdfFonts())) return(invisible(NULL));
  src <- c("Helvetica", "Helvetica-Bold", "Helvetica-Oblique", "Helvetica-BoldOblique");
  dst <- c("ArialMT", "Arial-BoldMT", "Arial-ItalicMT", "Arial-BoldItalicMT");
  afm <- file.path(tempdir(), paste0(dst, ".afm"));
  for (i in seq_along(src)) {
    x <- readLines(gzfile(file.path(system.file("afm", package = "grDevices"), paste0(src[i], ".afm.gz"))));
    writeLines(sub("^FontName .*", paste("FontName", dst[i]), x), afm[i]);
  }
  pdfFonts(Arial = Type1Font("Arial", c(afm, pdfFonts("Helvetica")[[1]]$metrics[5])));
};

# --font 为 Arial 时检查 gs 与字体，缺少时改用 Helvetica；返回 list(font, dir)
resolve_font <- function(font) {
  dir <- NULL;
  if (font == "Arial") {
    dir <- find_font_dir("Arial");
    if (!nzchar(Sys.which("gs")) || is.null(dir)) {
      message("[警告] 没有 ghostscript (gs) 或系统中没有 Arial 字体，--font 改用 Helvetica");
      font <- "Helvetica";
      dir <- NULL;
    } else {
      register_arial();
    }
  }
  list(font = font, dir = dir);
};

parser <- ArgumentParser(description = "RNA-seq volcano plot (log2FC = log2(treatment/control))");
parser$add_argument("-i", "--input", default = "avgExp_fdr.tsv",
  help = "输入表格（tsv）[default: %(default)s]");
parser$add_argument("-o", "--output", default = "volcano.pdf",
  help = "输出文件，扩展名决定格式：pdf/png/tiff [default: %(default)s]");
parser$add_argument("--skip_header", dest = "skip_header", action = "store_true", default = TRUE,
  help = "第一行是标题行，读入时忽略 [default]");
parser$add_argument("--no_skip_header", dest = "skip_header", action = "store_false",
  help = "文件没有标题行，第一行就是数据");
parser$add_argument("--gene", type = "integer", default = 1,
  help = "基因 ID 所在列号 [default: %(default)s]");
parser$add_argument("--control", type = "integer", default = 2,
  help = "control 表达量所在列号 [default: %(default)s]");
parser$add_argument("--treatment", type = "integer", default = 3,
  help = "treatment 表达量所在列号 [default: %(default)s]");
parser$add_argument("--fc_col", type = "integer", default = NULL,
  help = "treatment/control 值所在列号；指定后直接用该列，不再由 --control / --treatment 计算 [default: 不指定]");
parser$add_argument("--fc_type", default = "ratio", choices = c("ratio", "log2"),
  help = "--fc_col 列的数值类型：ratio = treatment/control 比值，脚本再取 log2；log2 = 已经是 log2(treatment/control) [default: %(default)s]");
parser$add_argument("--fdr", type = "integer", default = 5,
  help = "FDR 所在列号 [default: %(default)s]");
parser$add_argument("--min_value", type = "double", default = 0.01,
  help = "计算 fold change 时 control/treatment 的最小值，小于该值的按该值计；指定 --fc_col 时不用 [default: %(default)s]");
parser$add_argument("--fc_cut", type = "double", default = 1,
  help = "|log2(fold change)| 阈值 [default: %(default)s]");
parser$add_argument("--fdr_cut", type = "double", default = 0.05,
  help = "FDR 阈值 [default: %(default)s]");
parser$add_argument("--include_equal", nargs = "?", const = "both", default = "none",
  choices = c("none", "both", "fc", "fdr"),
  help = "判定 DEG 时阈值是否包含等号：none = FDR < fdr_cut 且 |log2FC| > fc_cut；both（只写 --include_equal 时）= FDR <= fdr_cut 且 |log2FC| >= fc_cut；fc / fdr = 只有该项包含等号 [default: %(default)s]");
parser$add_argument("--xlim", type = "double", nargs = 2, default = NULL,
  metavar = c("XMIN", "XMAX"), help = "X 轴起止坐标 [default: 根据数据自动、对称]");
parser$add_argument("--ylim", type = "double", nargs = 2, default = NULL,
  metavar = c("YMIN", "YMAX"), help = "Y 轴起止坐标 [default: 0 到数据最大值]");
parser$add_argument("--xstep", type = "double", default = NULL,
  help = "X 轴 tick 间隔 [default: 自动]");
parser$add_argument("--ystep", type = "double", default = NULL,
  help = "Y 轴 tick 间隔 [default: 自动]");
parser$add_argument("--col_ns", default = "#6BAED6",
  help = "Not significant 点颜色 [default: %(default)s]");
parser$add_argument("--col_up", default = "#F28E2B",
  help = "Up-regulated 点颜色 [default: %(default)s]");
parser$add_argument("--col_down", default = "#2CA25F",
  help = "Down-regulated 点颜色 [default: %(default)s]");
parser$add_argument("--col_line", default = "#2166AC",
  help = "阈值虚线颜色 [default: %(default)s]");
parser$add_argument("--point_size", type = "double", default = 0.6,
  help = "点大小（mm）[default: %(default)s]");
parser$add_argument("--width", type = "double", default = 89,
  help = "图宽（mm），Nature 单栏 89 / 双栏 183 [default: %(default)s]");
parser$add_argument("--height", type = "double", default = 80,
  help = "图高（mm）[default: %(default)s]");
parser$add_argument("--font_size", type = "double", default = 7,
  help = "基础字号（pt），Nature 要求 5-7 [default: %(default)s]");
parser$add_argument("--font", default = "Arial", choices = c("Arial", "Helvetica"),
  help = "图中文字字体：Arial = pdf 中字体名为 ArialMT，并用 ghostscript 嵌入系统中的 Arial 字体（需装有 gs 和 Arial 字体，缺少时自动改用 Helvetica）；Helvetica = R 自带 pdf 设备，嵌入 Nimbus Sans 作为 Helvetica 的替代字体 [default: %(default)s]");
args <- parser$parse_args();
fnt <- resolve_font(args$font);
args$font <- fnt$font;

# ---------- 读入数据 ----------
# 全部按字符读入，列号定位；是否跳过第一行由 --skip_header/--no_skip_header 决定
dt <- fread(args$input, sep = "\t", header = FALSE, skip = if (args$skip_header) 1 else 0,
  colClasses = "character", data.table = TRUE);
use_fc <- !is.null(args$fc_col);
need_cols <- if (use_fc) c(args$gene, args$fc_col, args$fdr) else c(args$gene, args$control, args$treatment, args$fdr);
bad_cols <- need_cols[need_cols < 1 | need_cols > ncol(dt)];
if (length(bad_cols) > 0) {
  stop("列号超出范围: ", paste(bad_cols, collapse = ", "), "（输入文件共 ", ncol(dt), " 列）");
}
if (args$skip_header) {
  hdr <- strsplit(readLines(args$input, n = 1), "\t")[[1]];
  if (use_fc) {
    message("[INFO] gene: ", hdr[args$gene], "; treatment/control (", args$fc_type, "): ", hdr[args$fc_col],
      "; FDR: ", hdr[args$fdr]);
  } else {
    message("[INFO] gene: ", hdr[args$gene], "; control: ", hdr[args$control],
      "; treatment: ", hdr[args$treatment], "; FDR: ", hdr[args$fdr]);
  }
}

dt <- dt[, ..need_cols];
if (use_fc) {
  setnames(dt, c("gene", "fc", "fdr"));
  dt[, `:=`(fc = as.numeric(fc), fdr = as.numeric(fdr))];
  n_neg <- if (args$fc_type == "ratio") dt[!is.na(fc) & fc < 0, .N] else 0;
  if (n_neg > 0) {
    message("[INFO] treatment/control 比值为负的行 ", n_neg, " 个，按 NA 处理");
    dt[fc < 0, fc := NA];
  }
} else {
  setnames(dt, c("gene", "ctl", "trt", "fdr"));
  dt[, `:=`(ctl = as.numeric(ctl), trt = as.numeric(trt), fdr = as.numeric(fdr))];
}
ok <- complete.cases(dt[, !"gene"]);
if (any(!ok)) message("[INFO] 去除含 NA 的行: ", sum(!ok));
dt <- dt[ok];

# ---------- 计算 log2FC 和 -log10(FDR) ----------
if (use_fc) {
  dt[, log2fc := if (args$fc_type == "ratio") log2(fc) else fc];
} else {
  dt[, log2fc := log2(pmax(trt, args$min_value) / pmax(ctl, args$min_value))];
}
# 比值为 0 或 Inf 时 log2FC 为 -Inf / Inf：照常判定 DEG，作图时压到 X 轴边界上
n_inf <- dt[is.infinite(log2fc), .N];
if (n_inf > 0) message("[INFO] log2FC 为 ±Inf 的基因 ", n_inf, " 个，照常判定 DEG，作图时压到 X 轴边界上");
# FDR = 0 时 -log10 为 Inf，用最小非零 FDR 代替
min_fdr <- dt[fdr > 0, min(fdr)];
n_zero <- dt[fdr <= 0, .N];
if (n_zero > 0) message("[INFO] FDR = 0 的基因 ", n_zero, " 个，按最小非零 FDR (", min_fdr, ") 计");
dt[, neglog10fdr := -log10(pmax(fdr, min_fdr))];

lv <- c("Not significant", "Up-regulated", "Down-regulated");
dt[, group := "Not significant"];
fdr_ok <- if (args$include_equal %in% c("both", "fdr")) dt$fdr <= args$fdr_cut else dt$fdr < args$fdr_cut;
fc_eq <- args$include_equal %in% c("both", "fc");
up_ok <- if (fc_eq) dt$log2fc >= args$fc_cut else dt$log2fc > args$fc_cut;
down_ok <- if (fc_eq) dt$log2fc <= -args$fc_cut else dt$log2fc < -args$fc_cut;
dt[fdr_ok & up_ok, group := "Up-regulated"];
dt[fdr_ok & down_ok, group := "Down-regulated"];
dt[, group := factor(group, levels = lv)];
n_grp <- table(dt$group);
message("[INFO] ", paste(names(n_grp), n_grp, sep = ": ", collapse = "; "));

# ---------- 坐标轴：起止和 tick，tick 总是经过 0 ----------
auto_step <- function(lim) {
  b <- pretty(lim, n = 5);
  return(b[2] - b[1]);
};
make_breaks <- function(lim, step) {
  # 从 step 的整数倍开始，保证 0 在 tick 上
  return(seq(ceiling(lim[1] / step - 1e-8) * step, floor(lim[2] / step + 1e-8) * step, by = step));
};

if (is.null(args$xlim)) {
  xmax <- max(abs(dt$log2fc[is.finite(dt$log2fc)]), args$fc_cut) * 1.05;
  if (is.null(args$xstep)) args$xstep <- auto_step(c(-xmax, xmax));
  xmax <- ceiling(xmax / args$xstep) * args$xstep;
  args$xlim <- c(-xmax, xmax);
}
if (is.null(args$ylim)) {
  ymax <- max(dt$neglog10fdr, -log10(args$fdr_cut)) * 1.05;
  ymin <- -0.05 * ymax;
  if (is.null(args$ystep)) args$ystep <- auto_step(c(ymin, ymax));
  ymax <- ceiling(ymax / args$ystep) * args$ystep;
  args$ylim <- c(ymin, ymax);
}
if (is.null(args$xstep)) args$xstep <- auto_step(args$xlim);
if (is.null(args$ystep)) args$ystep <- auto_step(args$ylim);
x_breaks <- make_breaks(args$xlim, args$xstep);
y_breaks <- make_breaks(args$ylim, args$ystep);

# 超出坐标范围的点压到边界上显示
n_oob <- dt[log2fc < args$xlim[1] | log2fc > args$xlim[2] | neglog10fdr < args$ylim[1] | neglog10fdr > args$ylim[2], .N];
if (n_oob > 0) message("[INFO] 超出坐标范围的点 ", n_oob, " 个，已压到边界上显示");
dt[, `:=`(
  plot_x = pmin(pmax(log2fc, args$xlim[1]), args$xlim[2]),
  plot_y = pmin(pmax(neglog10fdr, args$ylim[1]), args$ylim[2])
)];

# ---------- 作图 ----------
lab <- c(
  "Up-regulated" = sprintf("Up-regulated (n=%d)", n_grp[["Up-regulated"]]),
  "Down-regulated" = sprintf("Down-regulated (n=%d)", n_grp[["Down-regulated"]]),
  "Not significant" = "Not significant"
);
# R 的 pdf 设备把 "-" 画成减号，pdf 输出时改用软连字符 (U+00AD) 才显示为连字符
if (tolower(tools::file_ext(args$output)) == "pdf") lab <- gsub("-", "\u00ad", lab, fixed = TRUE);
cols <- c("Up-regulated" = args$col_up, "Down-regulated" = args$col_down, "Not significant" = args$col_ns);

pt2mm <- 1 / ggplot2::.pt;   # 线宽以 pt 指定
fs <- args$font_size;
# 先画 Not significant，显著点画在上层
dt <- dt[order(group != "Not significant")];

p <- ggplot(dt, aes(x = plot_x, y = plot_y, colour = group)) +
  geom_point(size = args$point_size, stroke = 0, shape = 16) +
  geom_vline(xintercept = c(-args$fc_cut, args$fc_cut), linetype = "dashed",
    colour = args$col_line, linewidth = 0.5 * pt2mm) +
  geom_hline(yintercept = -log10(args$fdr_cut), linetype = "dashed",
    colour = args$col_line, linewidth = 0.5 * pt2mm) +
  scale_colour_manual(values = cols, labels = lab, breaks = lv, name = NULL, drop = FALSE) +
  scale_x_continuous(limits = args$xlim, breaks = x_breaks, expand = expansion(mult = 0.02)) +
  scale_y_continuous(limits = args$ylim, breaks = y_breaks, expand = expansion(mult = c(0, 0.02))) +
  labs(x = expression(log[2] * "(fold change)"), y = expression("-" * log[10] * "(FDR)")) +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) +
  theme_classic(base_size = fs, base_family = args$font, base_line_size = 0.5 * pt2mm) +
  theme(
    axis.text = element_text(size = fs - 1, colour = "black"),
    axis.title = element_text(size = fs, colour = "black"),
    axis.line = element_line(linewidth = 0.5 * pt2mm, colour = "black"),
    axis.ticks = element_line(linewidth = 0.5 * pt2mm, colour = "black"),
    axis.ticks.length = unit(1, "mm"),
    legend.text = element_text(size = fs - 1, colour = "black"),
    legend.position = "inside",
    legend.position.inside = c(1, 1),
    legend.justification = c(1, 1),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.key.size = unit(3, "mm"),
    legend.spacing.y = unit(0, "mm"),
    plot.margin = margin(2, 2, 2, 2, "mm")
  );

# ---------- 输出 ----------
ext <- tolower(tools::file_ext(args$output));
if (ext == "pdf") {
  # useDingbats = FALSE 保证点在 Illustrator 等软件中正常显示
  ggsave(args$output, p, width = args$width, height = args$height, units = "mm",
    device = "pdf", useDingbats = FALSE, family = args$font, encoding = "WinAnsi.enc");
  # 期刊要求 PDF 嵌入字体，需要 ghostscript
  if (nzchar(Sys.which("gs"))) {
    tmp <- paste0(args$output, ".tmp.pdf");
    embedFonts(args$output, outfile = tmp, fontpaths = if (is.null(fnt$dir)) character(0) else fnt$dir);
    invisible(file.rename(tmp, args$output));
  } else {
    message("[WARN] 未找到 ghostscript，PDF 字体未嵌入");
  }
} else if (ext == "png") {
  ggsave(args$output, p, width = args$width, height = args$height, units = "mm", dpi = 600);
} else if (ext %in% c("tiff", "tif")) {
  ggsave(args$output, p, width = args$width, height = args$height, units = "mm", dpi = 600,
    device = "tiff", compression = "lzw");
} else {
  stop("不支持的输出格式: ", ext, "（可用 pdf/png/tiff）");
}
message("[INFO] 输出: ", args$output);
