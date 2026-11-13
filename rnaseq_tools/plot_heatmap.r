#!/usr/bin/env Rscript
# 表达量热图：合并 plot_expr_heatmap.r 与 plot_heatmap_by_geneList.pl 的功能
#   - 4 种数据处理方式 log2_none / log2_row / raw_none / raw_row，每种单独输出文件
#   - 可选 gene list；不指定时对全部基因作图
#   - preset：pub（Nature 系列期刊发表）与 screen（放大浏览）
#   - 聚类树与分组色带用 ggplot2 + gtable 绘制，各部件尺寸以 mm 精确控制

suppressPackageStartupMessages({
  library(argparse);
  library(data.table);
  library(ggplot2);
  library(gtable);
  library(grid);
  library(RColorBrewer);
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

# 用 ghostscript 嵌入字体。fontDir 不为空（Arial 等 TrueType 字体）时完整嵌入：gs 只嵌入子集时把 TrueType
# 写成无编码表的符号字体，Illustrator 打开后换用系统字体、按符号编码取字形而显示乱码；完整嵌入则为标准的 WinAnsi 字体
embed_pdf_fonts <- function(file, fontDir = NULL) {
  tmp <- paste0(file, ".tmp.pdf");
  if (is.null(fontDir)) {
    embedFonts(file, outfile = tmp);
  } else {
    embedFonts(file, outfile = tmp, fontpaths = fontDir, options = "-dSubsetFonts=false");
    clear_cap_flags(tmp);
  }
  invisible(file.rename(tmp, file));
};

# gs 有时给字体描述的 /Flags 误加 AllCap、SmallCap 标志，原位清除（数字后补空格，字节数不变，不影响交叉引用表）
clear_cap_flags <- function(file) {
  x <- readBin(file, "raw", file.size(file));
  for (p in grepRaw("/Flags [0-9]+", x, all = TRUE)) {
    s <- p + 7;
    e <- s;
    while (e < length(x) && x[e + 1] >= as.raw(48) && x[e + 1] <= as.raw(57)) e <- e + 1;
    v <- as.integer(rawToChar(x[s:e]));
    nv <- bitwAnd(v, bitwNot(65536L + 131072L));
    if (nv != v) x[s:e] <- charToRaw(formatC(nv, width = -(e - s + 1)));
  }
  writeBin(x, file);
};

##############################
# 命令行参数
##############################
parser <- ArgumentParser(description = "Expression heatmap with optional gene list (ggplot2).");
parser$add_argument("-i", "--input", default = NULL,
  help = "【必需】表达量表格 (tsv)：第 1 列为基因 ID，其余列为样本");
parser$add_argument("-g", "--gene_list", default = NULL,
  help = "基因列表文件，取每行第 1 列（忽略空行和 # 开头的行）；不指定时对全部基因作图");
parser$add_argument("-o", "--out_prefix", default = NULL,
  help = "输出文件前缀，输出为 <prefix>_<mode>.<ext> [default: 输入文件路径去掉扩展名]");
parser$add_argument("--skip_header", dest = "skip_header", action = "store_true", default = TRUE,
  help = "输入文件第一行为标题行（样本名取自标题行）[default]");
parser$add_argument("--no_skip_header", dest = "skip_header", action = "store_false",
  help = "输入文件没有标题行，样本名记为 Col2, Col3, ...");
parser$add_argument("--modes", default = "log2_none,log2_row,raw_none,raw_row",
  help = "数据处理方式，逗号分隔，可选 log2_none, log2_row, raw_none, raw_row [default: %(default)s]");
parser$add_argument("--pseudo", type = "double", default = 1,
  help = "log2(x + pseudo) 的伪计数 [default: %(default)s]");
parser$add_argument("--fill_log2", action = "store_true", default = FALSE,
  help = "热图涂色统一用 log2(x + pseudo)（不缩放），聚类仍按 --modes 对应的数值；不加时涂色与聚类数值相同");
parser$add_argument("--format", default = "pdf",
  help = "输出图片格式，逗号分隔，可选 pdf, png, tiff [default: %(default)s]");
parser$add_argument("--preset", default = "pub", choices = c("pub", "screen"),
  help = "pub: Nature 系列期刊发表规格；screen: 按格子大小自动放大画布，便于浏览 [default: %(default)s]");
parser$add_argument("--palette", default = "RdBu",
  help = "RColorBrewer 色板名；diverging 类色板自动反转（低值蓝、高值红）[default: %(default)s]");
parser$add_argument("--colors", default = NULL,
  help = "自定义渐变色（低 -> 高），逗号分隔，如 \"#2166AC,#F7F7F7,#B2182B\"；指定后忽略 --palette");
parser$add_argument("--na_col", default = "grey80",
  help = "缺失值颜色 [default: %(default)s]");
parser$add_argument("--no_cluster_rows", action = "store_true", default = FALSE,
  help = "不对基因聚类，按 gene list（无 gene list 时按输入文件）顺序排列");
parser$add_argument("--no_cluster_cols", action = "store_true", default = FALSE,
  help = "不对样本聚类，按输入文件列顺序排列");
parser$add_argument("--dist", default = "euclidean",
  choices = c("euclidean", "maximum", "manhattan", "canberra", "minkowski", "correlation"),
  help = "聚类距离 [default: %(default)s]");
parser$add_argument("--hclust_method", default = "complete",
  choices = c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"),
  help = "hclust 聚类方法 [default: %(default)s]");
parser$add_argument("--k", type = "integer", default = NULL,
  help = "将基因聚类树切成 k 组：输出 <prefix>_<mode>.grp，并在热图左侧画分组色带");
parser$add_argument("--show_rownames", default = "auto", choices = c("auto", "yes", "no"),
  help = "是否显示基因名；auto: pub 基因数 <= 60、screen 基因数 <= 300 时显示 [default: %(default)s]");
parser$add_argument("--col_angle", type = "double", default = NULL,
  help = "样本名旋转角度 [default: 样本数 <= 8 且名字不超过 8 个字符时为 0，否则为 90]");
parser$add_argument("--width", type = "double", default = NULL,
  help = "画布宽度 (mm) [default: pub 89；screen 按格子大小自动计算]");
parser$add_argument("--height", type = "double", default = NULL,
  help = "画布高度 (mm) [default: pub 120；screen 按格子大小自动计算]");
parser$add_argument("--font_size", type = "double", default = NULL,
  help = "基础字号 (pt)，刻度文字比它小 1 pt [default: pub 7；screen 8]");
parser$add_argument("--dpi", type = "integer", default = NULL,
  help = "png/tiff 分辨率 [default: pub 600；screen 150]");
parser$add_argument("--font", default = "Arial", choices = c("Arial", "Helvetica"),
  help = "图中文字字体：Arial = pdf 中字体名为 ArialMT，并用 ghostscript 嵌入系统中的 Arial 字体（需装有 gs 和 Arial 字体，缺少时自动改用 Helvetica）；Helvetica = R 自带 pdf 设备，嵌入 Nimbus Sans 作为 Helvetica 的替代字体 [default: %(default)s]");
args <- parser$parse_args();
fnt <- resolve_font(args$font);
args$font <- fnt$font;
if (is.null(args$input)) {
  parser$print_help();
  message("\n[错误] 必须用 -i/--input 指定表达量文件");
  quit(save = "no", status = 1);
}

##############################
# preset（长度单位均为 mm）
##############################
preset <- switch(args$preset,
  pub = list(width = 89, height = 120, font_size = 7, dpi = 600,
    cell_w = NA, cell_h = NA, row_dendro = 12, col_dendro = 6, strip = 2,
    bar_w = 2.5, bar_h = 15, grp_key = 2.5, title = FALSE, rowname_max = 60),
  screen = list(width = NA, height = NA, font_size = 8, dpi = 150,
    cell_w = 15, cell_h = 3.5, row_dendro = 40, col_dendro = 15, strip = 4,
    bar_w = 4, bar_h = 30, grp_key = 4, title = TRUE, rowname_max = 300));
if (!is.null(args$width)) preset$width <- args$width;
if (!is.null(args$height)) preset$height <- args$height;
if (!is.null(args$font_size)) preset$font_size <- args$font_size;
if (!is.null(args$dpi)) preset$dpi <- args$dpi;

fs <- preset$font_size;
afs <- fs - 1;
lw <- 0.5 / .pt;
gap <- 0.5;

modes <- trimws(strsplit(args$modes, ",")[[1]]);
bad <- setdiff(modes, c("log2_none", "log2_row", "raw_none", "raw_row"));
if (length(bad) > 0) stop("未知的 --modes: ", paste(bad, collapse = ", "));
fmts <- tolower(trimws(strsplit(args$format, ",")[[1]]));
fmts[fmts == "tif"] <- "tiff";
bad <- setdiff(fmts, c("pdf", "png", "tiff"));
if (length(bad) > 0) stop("未知的 --format: ", paste(bad, collapse = ", "));

if (is.null(args$out_prefix)) args$out_prefix <- sub("\\.[^./]*$", "", args$input);

##############################
# 配色
##############################
if (!is.null(args$colors)) {
  pal_base <- trimws(strsplit(args$colors, ",")[[1]]);
  ok <- vapply(pal_base, function(cc) !inherits(try(col2rgb(cc), silent = TRUE), "try-error"), logical(1));
  if (!all(ok)) stop("无法识别的颜色: ", paste(pal_base[!ok], collapse = ", "));
  if (length(pal_base) < 2) stop("--colors 至少需要 2 种颜色");
} else {
  if (!args$palette %in% rownames(brewer.pal.info)) stop("未知的 RColorBrewer 色板: ", args$palette);
  info <- brewer.pal.info[args$palette, ];
  pal_base <- brewer.pal(min(9, info$maxcolors), args$palette);
  if (info$category == "div") pal_base <- rev(pal_base);
}
pal <- colorRampPalette(pal_base)(255);

##############################
# 读入表达量与 gene list
##############################
d <- fread(args$input, sep = "\t", header = args$skip_header,
  colClasses = "character", na.strings = c("", "NA"));
if (ncol(d) < 2) stop("输入文件至少需要 2 列（基因 ID + 样本）");
if (!args$skip_header) setnames(d, c("Gene_ID", paste0("Col", seq(2, ncol(d)))));

raw_chr <- as.matrix(d[, -1, with = FALSE]);
expr <- suppressWarnings(matrix(as.numeric(raw_chr), nrow = nrow(raw_chr),
  dimnames = list(d[[1]], colnames(raw_chr))));
nBad <- sum(is.na(expr) & !is.na(raw_chr));
if (nBad > 0) message(sprintf("[警告] %d 个非数值单元格已按 NA 处理", nBad));

dup <- duplicated(rownames(expr));
if (any(dup)) {
  message(sprintf("[警告] 表达量文件中有 %d 个重复基因 ID，只保留第一次出现的行", sum(dup)));
  expr <- expr[!dup, , drop = FALSE];
}

if (!is.null(args$gene_list)) {
  gl <- trimws(readLines(args$gene_list, warn = FALSE));
  gl <- sub("[[:space:]].*$", "", gl);
  gl <- unique(gl[gl != "" & !startsWith(gl, "#")]);
  miss <- setdiff(gl, rownames(expr));
  if (length(miss) > 0) {
    message(sprintf("[警告] gene list 中有 %d 个基因不在表达量文件中，已跳过: %s%s",
      length(miss), paste(head(miss, 10), collapse = ", "), if (length(miss) > 10) ", ..." else ""));
  }
  expr <- expr[intersect(gl, rownames(expr)), , drop = FALSE];
}

allNA <- rowSums(!is.na(expr)) == 0;
if (any(allNA)) {
  message(sprintf("[警告] %d 个基因全部为缺失值，已去除", sum(allNA)));
  expr <- expr[!allNA, , drop = FALSE];
}
if (nrow(expr) == 0) stop("没有可作图的基因");
message(sprintf("作图基因数: %d；样本数: %d", nrow(expr), ncol(expr)));

if ((any(grepl("^log2", modes)) || args$fill_log2) && any(expr <= -args$pseudo, na.rm = TRUE)) {
  message(sprintf("[警告] 存在 <= -%g 的值，log2(x + %g) 后为 NA", args$pseudo, args$pseudo));
}

##############################
# 函数
##############################
# 数据处理：log2 转换与按行 z-score；按行缩放时去掉标准差为 0 或无法计算的行
transform_mat <- function(x, mode) {
  if (startsWith(mode, "log2")) {
    x <- log2(x + args$pseudo);
    x[is.nan(x)] <- NA;
  }
  if (endsWith(mode, "_row")) {
    s <- apply(x, 1, sd, na.rm = TRUE);
    keep <- !is.na(s) & s > 0;
    if (any(!keep)) message(sprintf("  [%s] 去除 %d 个标准差为 0 的基因", mode, sum(!keep)));
    x <- x[keep, , drop = FALSE];
    x <- (x - rowMeans(x, na.rm = TRUE)) / s[keep];
  }
  x[rowSums(!is.na(x)) > 0, , drop = FALSE];
};

calc_dist <- function(x, method) {
  if (method == "correlation") {
    dd <- as.dist(1 - cor(t(x), use = "pairwise.complete.obs"));
  } else {
    dd <- dist(x, method = method);
  }
  if (anyNA(dd)) stop("距离矩阵含 NA（缺失值过多或相关系数无法计算），请换用其它 --dist 或检查数据");
  dd;
};

# 聚类树线段：pos 为叶节点位置（1..n，按 hc$order），h 为高度
dendro_seg <- function(hc) {
  n <- length(hc$order);
  leafPos <- integer(n);
  leafPos[hc$order] <- seq_len(n);
  nodePos <- numeric(n - 1);
  getPos <- function(j) if (j < 0) leafPos[-j] else nodePos[j];
  getHgt <- function(j) if (j < 0) 0 else hc$height[j];
  seg <- vector("list", n - 1);
  for (i in seq_len(n - 1)) {
    a <- hc$merge[i, 1];
    b <- hc$merge[i, 2];
    h <- hc$height[i];
    pa <- getPos(a);
    pb <- getPos(b);
    nodePos[i] <- (pa + pb) / 2;
    seg[[i]] <- data.table(pos = c(pa, pb, pa), h = c(getHgt(a), getHgt(b), h),
      pos_end = c(pa, pb, pb), h_end = c(h, h, h));
  }
  rbindlist(seg);
};

panel_grob <- function(p) {
  g <- ggplotGrob(p);
  g$grobs[[which(g$layout$name == "panel")]];
};

# 行聚类树：根在左，叶节点贴近热图；第 1 个叶节点在最上方
row_dendro_grob <- function(hc) {
  n <- length(hc$order);
  s <- dendro_seg(hc);
  p <- ggplot(s) +
    geom_segment(aes(x = h, y = n + 1 - pos, xend = h_end, yend = n + 1 - pos_end),
      linewidth = lw, lineend = "square") +
    scale_x_reverse(expand = expansion(mult = c(0.02, 0))) +
    scale_y_continuous(limits = c(0.5, n + 0.5), expand = c(0, 0)) +
    theme_void();
  panel_grob(p);
};

col_dendro_grob <- function(hc) {
  m <- length(hc$order);
  s <- dendro_seg(hc);
  p <- ggplot(s) +
    geom_segment(aes(x = pos, y = h, xend = pos_end, yend = h_end),
      linewidth = lw, lineend = "square") +
    scale_x_continuous(limits = c(0.5, m + 0.5), expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    theme_void();
  panel_grob(p);
};

group_strip_grob <- function(grp, cols) {
  n <- length(grp);
  sd <- data.table(x = 1, y = n:1, f = cols[grp]);
  p <- ggplot(sd, aes(x = x, y = y, fill = f)) +
    geom_raster() +
    scale_fill_identity() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.5, n + 0.5), expand = c(0, 0)) +
    theme_void();
  panel_grob(p);
};

group_colors <- function(k) {
  if (k <= 8) brewer.pal(max(3, k), "Set2")[seq_len(k)] else hcl.colors(k, "Dark 3");
};

legend_title <- function(mode) {
  if (endsWith(mode, "_row")) return("Row z-score");
  if (mode == "log2_none") return(bquote(log[2] * "(Exp + " * .(format(args$pseudo)) * ")"));
  "Exp";
};

mode_desc <- function(mode) {
  tr <- if (startsWith(mode, "log2")) sprintf("log2(Exp + %s)", format(args$pseudo)) else "Exp";
  sc <- if (endsWith(mode, "_row")) "row scaled" else "no scale";
  sprintf("%s [%s]", tr, sc);
};

# R 的 pdf 设备把 "-" 画成减号：pdf 输出前把文字中紧跟在非空白字符后的 "-"（如 z-score、基因名中的 -）
# 换成软连字符 (U+00AD) 以显示为连字符；开头的负号（如刻度 -1.0）仍画成减号
soft_hyphen <- function(g) {
  if (inherits(g, "text") && is.character(g$label)) g$label <- gsub("(?<=\\S)-", "\u00ad", g$label, perl = TRUE);
  if (inherits(g, "gtable")) g$grobs <- lapply(g$grobs, soft_hyphen);
  if (inherits(g, "gTree")) for (nm in names(g$children)) g$children[[nm]] <- soft_hyphen(g$children[[nm]]);
  g;
};

save_plot <- function(g, file, fmt, W, H, dpi) {
  if (fmt == "pdf") {
    g <- soft_hyphen(g);
    pdf(file, width = W / 25.4, height = H / 25.4, useDingbats = FALSE, family = args$font);
    grid.draw(g);
    invisible(dev.off());
    if (nzchar(Sys.which("gs"))) embed_pdf_fonts(file, fnt$dir);
  } else if (fmt == "png") {
    ragg::agg_png(file, width = W, height = H, units = "mm", res = dpi);
    grid.draw(g);
    invisible(dev.off());
  } else {
    ragg::agg_tiff(file, width = W, height = H, units = "mm", res = dpi, compression = "lzw");
    grid.draw(g);
    invisible(dev.off());
  }
};

##############################
# 逐个处理方式作图
##############################
# 文字宽度等尺寸的计算需要打开一个图形设备
pdf(NULL);

for (mode in modes) {
  x <- transform_mat(expr, mode);
  n <- nrow(x);
  m <- ncol(x);
  if (n == 0) {
    message(sprintf("[%s] 没有可作图的基因，跳过", mode));
    next;
  }

  # 聚类；不聚类时保持输入顺序
  doRow <- !args$no_cluster_rows && n >= 2;
  doCol <- !args$no_cluster_cols && m >= 2;
  hcRow <- if (doRow) hclust(calc_dist(x, args$dist), method = args$hclust_method) else NULL;
  hcCol <- if (doCol) hclust(calc_dist(t(x), args$dist), method = args$hclust_method) else NULL;
  rOrd <- if (doRow) hcRow$order else seq_len(n);
  cOrd <- if (doCol) hcCol$order else seq_len(m);
  xo <- x[rOrd, cOrd, drop = FALSE];

  # 基因分组：组号按热图从上到下的出现顺序重新编号
  grp <- NULL;
  if (!is.null(args$k)) {
    if (!doRow) {
      message(sprintf("  [%s] 未对基因聚类，忽略 --k", mode));
    } else if (args$k < 1 || args$k > n) {
      message(sprintf("  [%s] --k 须在 1..%d 之间，忽略", mode, n));
    } else {
      g0 <- cutree(hcRow, k = args$k)[rOrd];
      grp <- match(g0, unique(g0));
      fn_grp <- sprintf("%s_%s.grp", args$out_prefix, mode);
      fwrite(data.table(Gene_ID = rownames(xo), Group_ID = grp), fn_grp, sep = "\t", quote = FALSE);
      message(sprintf("  [%s] 分组写入 %s", mode, fn_grp));
    }
  }

  # 涂色数值：默认与聚类数值相同；--fill_log2 时统一用 log2(x + pseudo)，行列顺序沿用聚类结果
  if (args$fill_log2) {
    fo <- log2(expr[rownames(xo), colnames(xo), drop = FALSE] + args$pseudo);
    fo[is.nan(fo)] <- NA;
    fillTitle <- legend_title("log2_none");
  } else {
    fo <- xo;
    fillTitle <- legend_title(mode);
  }

  # 色标范围：按行缩放的数值以 0 为中心对称
  if (!args$fill_log2 && endsWith(mode, "_row")) {
    lim <- c(-1, 1) * max(abs(fo), na.rm = TRUE);
  } else {
    lim <- range(fo, na.rm = TRUE);
  }
  plotTitle <- sprintf("%s, %d genes%s", mode_desc(mode), n,
    if (args$fill_log2) sprintf("; fill: log2(Exp + %s)", format(args$pseudo)) else "");

  # 基因名是否显示、字号
  showRow <- switch(args$show_rownames, yes = TRUE, no = FALSE, auto = n <= preset$rowname_max);
  cellPt <- if (is.na(preset$height)) preset$cell_h * .pt else (preset$height - 25) / n * .pt;
  rowFs <- min(afs, cellPt * 0.9);
  if (showRow && rowFs < 5) message(sprintf("  [%s] 基因名字号仅 %.1f pt，可增大 --height 或用 --show_rownames no", mode, rowFs));

  colAngle <- if (!is.null(args$col_angle)) args$col_angle else
    if (m <= 8 && max(nchar(colnames(xo))) <= 8) 0 else 90;

  dt <- data.table(col = rep(seq_len(m), each = n), y = rep(n:1, times = m), value = as.vector(fo));
  p <- ggplot(dt, aes(x = col, y = y, fill = value)) +
    geom_raster() +
    scale_x_continuous(breaks = seq_len(m), labels = colnames(xo),
      limits = c(0.5, m + 0.5), expand = c(0, 0)) +
    scale_y_continuous(breaks = if (showRow) seq_len(n) else NULL, labels = if (showRow) rev(rownames(xo)) else NULL,
      limits = c(0.5, n + 0.5), expand = c(0, 0), position = "right") +
    scale_fill_gradientn(colours = pal, limits = lim, na.value = args$na_col, name = fillTitle,
      guide = guide_colourbar(order = 1, theme = theme(
        legend.key.width = unit(preset$bar_w, "mm"), legend.key.height = unit(preset$bar_h, "mm"),
        legend.ticks = element_blank()))) +
    labs(x = NULL, y = NULL, title = if (preset$title) plotTitle else NULL) +
    theme_minimal(base_size = fs, base_family = args$font) +
    theme(
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_text(size = afs, colour = "black", angle = colAngle,
        hjust = if (colAngle == 0) 0.5 else 1, vjust = if (colAngle == 0) 1 else 0.5,
        margin = margin(t = 1, unit = "mm")),
      axis.text.y.right = element_text(size = rowFs, colour = "black", hjust = 0,
        margin = margin(l = 1, unit = "mm")),
      legend.title = element_text(size = fs),
      legend.text = element_text(size = afs),
      legend.justification = c(0, 1),
      legend.box.just = "left",
      legend.box.spacing = unit(2, "mm"),
      legend.spacing.y = unit(3, "mm"),
      plot.title = element_text(size = fs, hjust = 0.5, margin = margin(b = 2, unit = "mm")),
      plot.title.position = "plot",
      plot.margin = margin(2, 2, 2, 2, unit = "mm"));

  # 分组图例：用不可见的点生成方块图例，与色标放在一起
  if (!is.null(grp)) {
    gcols <- group_colors(args$k);
    gd <- data.table(col = 1, y = 1, grp = factor(seq_len(args$k)));
    p <- p +
      geom_point(data = gd, aes(x = col, y = y, colour = grp), inherit.aes = FALSE, alpha = 0, shape = 15) +
      scale_colour_manual(values = gcols, name = "Group",
        guide = guide_legend(order = 2, override.aes = list(alpha = 1, size = preset$grp_key),
          theme = theme(legend.key.size = unit(preset$grp_key, "mm"),
            legend.key.spacing.y = unit(0.5, "mm"))));
  }

  # 拼接：热图 panel 左侧从外到内为行聚类树、分组色带，上方为列聚类树
  # （每次都紧贴 panel 左侧插入，所以先插外层的聚类树）
  g <- ggplotGrob(p);
  pl <- g$layout$l[g$layout$name == "panel"];
  pt <- g$layout$t[g$layout$name == "panel"];
  if (doRow) {
    g <- gtable_add_cols(g, unit(gap, "mm"), pos = pl - 1);
    g <- gtable_add_cols(g, unit(preset$row_dendro, "mm"), pos = pl - 1);
    g <- gtable_add_grob(g, row_dendro_grob(hcRow), t = pt, l = pl, clip = "off", name = "row-dendro");
    pl <- pl + 2;
  }
  if (!is.null(grp)) {
    g <- gtable_add_cols(g, unit(gap, "mm"), pos = pl - 1);
    g <- gtable_add_cols(g, unit(preset$strip, "mm"), pos = pl - 1);
    g <- gtable_add_grob(g, group_strip_grob(grp, gcols), t = pt, l = pl, clip = "off", name = "row-group");
    pl <- pl + 2;
  }
  if (doCol) {
    g <- gtable_add_rows(g, unit(gap, "mm"), pos = pt - 1);
    g <- gtable_add_rows(g, unit(preset$col_dendro, "mm"), pos = pt - 1);
    g <- gtable_add_grob(g, col_dendro_grob(hcCol), t = pt, l = pl, clip = "off", name = "col-dendro");
    pt <- pt + 2;
  }

  # 画布大小：给定值优先；否则按格子大小固定 panel，再累加各部件尺寸
  if (is.na(preset$width)) {
    g$widths[pl] <- unit(m * preset$cell_w, "mm");
    W <- convertWidth(sum(g$widths), "mm", valueOnly = TRUE);
  } else {
    W <- preset$width;
  }
  if (is.na(preset$height)) {
    g$heights[pt] <- unit(n * preset$cell_h, "mm");
    H <- convertHeight(sum(g$heights), "mm", valueOnly = TRUE);
  } else {
    H <- preset$height;
  }

  for (fmt in fmts) {
    fn <- sprintf("%s_%s.%s", args$out_prefix, mode, fmt);
    save_plot(g, fn, fmt, W, H, preset$dpi);
    message(sprintf("  [%s] %d genes x %d samples -> %s (%.1f x %.1f mm)", mode, n, m, fn, W, H));
  }
}

invisible(dev.off());
