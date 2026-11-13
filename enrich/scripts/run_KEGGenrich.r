#!/usr/bin/env Rscript
# KEGG 通路富集分析（clusterProfiler::enricher）并画点图
#   - KO -> pathway 对应关系、通路名称与分类（br08901）在线下载（--kegg_source rest 或 KEGGREST），缓存到本地
#   - 富集结果表 <prefix>.tsv 的列与旧版脚本一致，末尾增加 KEGG_class_A / KEGG_class_B 两列
#   - 点图：x = GeneRatio_in_term 等（--x_col 指定列），点大小 = Count，颜色 = -log10(p.adjust)，按 KEGG 一级分类分面（--no_facet 不分面）
#   - preset：pub（Nature Genetics 规格，宽 180 mm、高不超过 185 mm）与 screen（放大浏览）
#   - --from_tsv：直接读已有结果表作图，不重跑富集
#
# 用法：
#   Rscript run_KEGGenrich.r --gene_list CPMMV_up.list --kegg_annot_2col MEAM1_v1.2-ko_2col.tsv --output_prefix CPMMV_up-KEGGen
#   Rscript run_KEGGenrich.r --from_tsv CPMMV_up-KEGGen.tsv --preset screen --format pdf,png

suppressPackageStartupMessages({
  library(argparse);
  library(data.table);
  library(ggplot2);
  library(grid);
});

#### >>> 共用代码开始：以下至「共用代码结束」在 run_GOenrich.r、run_KEGGenrich.r、run_IPRenrich.r 中相同 ####

pt2mm <- 25.4 / 72;
lw <- 0.5 / .pt;

# 画图相关的命令行参数
add_plot_args <- function(parser) {
  parser$add_argument("--padj_cutoff", type = "double", default = 0.05,
    help = "p.adjust 阈值：结果表 Tag 列标 OVER、以及上图的 term 都按此阈值 [default: %(default)s]");
  parser$add_argument("--x_col", default = "GeneRatio_in_term",
    help = "点图 x 轴取结果表中的哪一列，x 轴名称即该列名；a/b 形式的列（如 GeneRatio、GeneRatio_in_term）按分数计算。常用：GeneRatio_in_term (= Count / term 背景基因数，本脚本输出，总是存在)、RichFactor (clusterProfiler 新版才有的列，数值与 GeneRatio_in_term 相同)、FoldEnrichment、GeneRatio (= Count / 输入基因数) [default: %(default)s]");
  parser$add_argument("--order_by", default = "x", choices = c("x", "padj"),
    help = "分面内 term 的排列：x = 按 x 值从大到小（自上而下）；padj = 按 p.adjust 从小到大 [default: %(default)s]");
  parser$add_argument("--grid", action = "store_true", default = FALSE,
    help = "画出每个 term 的横向网格线 [default: 不画]");
  parser$add_argument("--no_facet", action = "store_true", default = FALSE,
    help = "不按类别分面（GO 的 BP / MF / CC，KEGG 的 class A；IPR 本来就不分面），所有 term 画在同一个面板里，按 --order_by 统一排序；--top_n、--simplify 也改为在全部 term 中进行 [default: 分面]");
  parser$add_argument("--shape_by_class", action = "store_true", default = FALSE,
    help = "用点的形状区分类别（GO 的 BP / MF / CC，KEGG 的 class A），并加形状图例；常与 --no_facet 一起用。IPR 没有类别，忽略此参数 [default: 不区分]");
  parser$add_argument("--shapes", default = "21,22,24,23,25",
    help = "--shape_by_class 按类别顺序（GO 为 BP, MF, CC）依次使用的形状，逗号分隔的 R pch 编号，须为可填色的 21~25：21 圆、22 方、23 菱形、24 上三角、25 下三角 [default: %(default)s]");
  parser$add_argument("--top_n", type = "integer", default = NULL,
    help = "每个分面（--no_facet 时为全部 term 中）只画 p.adjust 最小的前 N 个 term [default: 全部]");
  parser$add_argument("--min_count", type = "integer", default = 1,
    help = "只画 Count >= 此值的 term [default: %(default)s]");
  parser$add_argument("--simplify", default = "none",
    help = "去冗余：none 不去；0~1 的数值 = 同一分面内两 term 的富集基因 Jaccard 相似度 >= 此值时只保留 p.adjust 最小（其次 term 较小）的一个；1 = 只合并基因完全相同的 term。只影响作图，不改动结果表 [default: %(default)s]");
  parser$add_argument("--wrap", type = "integer", default = NULL,
    help = "term 名称超过此字符数时折行，0 为不折行 [default: pub 70；screen 0]");
  parser$add_argument("--format", default = "pdf",
    help = "输出图片格式，逗号分隔，可选 pdf, png, tiff [default: %(default)s]");
  parser$add_argument("--preset", default = "pub", choices = c("pub", "screen"),
    help = "pub: Nature Genetics 规格（宽 180 mm，高按 term 数自动、不超过 185 mm，7 pt）；screen: 按 term 数自动放大画布，便于浏览 [default: %(default)s]");
  parser$add_argument("--palette", default = "enrichplot",
    help = "颜色：enrichplot = clusterProfiler dotplot 默认配色（低显著蓝 #327EBA -> 高显著红 #E06663）；viridis 系列（viridis, magma, plasma, inferno, cividis, mako, rocket, turbo，自动反转为越显著越深）；或 RColorBrewer 色板名 [default: %(default)s]");
  parser$add_argument("--colors", default = NULL,
    help = "自定义渐变色（低显著 -> 高显著），逗号分隔，如 \"#FDE725,#21908C,#440154\"；指定后忽略 --palette");
  parser$add_argument("--width", type = "double", default = NULL,
    help = "画布宽度 (mm) [default: pub 180；screen 按内容自动计算]");
  parser$add_argument("--panel_width", type = "double", default = NULL,
    help = "面板（图框）宽度 (mm)；指定后画布宽度 = 面板宽度 + term 名称 + 图例等，--width 不再起作用。多张图拼合时各图用相同的值，图框即等长 [default: pub 由 --width 减去其余部分；screen 90]");
  parser$add_argument("--label_width", type = "double", default = NULL,
    help = "画布左缘到面板左缘的宽度 (mm)，即 term 名称区宽度；指定后 term 名称按实际字宽折行以放进此宽度（代替 --wrap），不足部分在左侧留白。与 --panel_width 一起用，各图的图框左右位置完全一致，可直接上下对齐拼合 [default: 按最长的 term 名称自动]");
  parser$add_argument("--legend_width", type = "double", default = NULL,
    help = "面板右缘到画布右缘的宽度 (mm)，即分面名称 + 图例区宽度；不足部分在右侧留白。与 --label_width、--panel_width 一起用，各图名称区、面板、图例区宽度都相同，画布总宽 = 三者之和 [default: 按图例自动]");
  parser$add_argument("--height", type = "double", default = NULL,
    help = "画布高度 (mm)；指定后不再自动计算 [default: 按 term 数自动计算]");
  parser$add_argument("--max_height", type = "double", default = NULL,
    help = "pub 自动高度的上限 (mm)；放不下时 term 名称字号自动缩小（最小 5 pt），仍放不下则超出上限并给出警告 [default: 185]");
  parser$add_argument("--font_size", type = "double", default = NULL,
    help = "基础字号 (pt)，刻度与图例文字比它小 1 pt [default: pub 7；screen 9]");
  parser$add_argument("--dpi", type = "integer", default = NULL,
    help = "png/tiff 分辨率 [default: pub 600；screen 150]");
  parser$add_argument("--font", default = "Arial", choices = c("Arial", "Helvetica"),
    help = "图中文字字体：Arial = pdf 中字体名为 ArialMT，并用 ghostscript 嵌入系统中的 Arial 字体（需装有 gs 和 Arial 字体，缺少时自动改用 Helvetica）；Helvetica = R 自带 pdf 设备，嵌入 Nimbus Sans 作为 Helvetica 的替代字体 [default: %(default)s]");
  invisible(parser);
};

# preset（长度单位均为 mm）
get_preset <- function(args) {
  pr <- switch(args$preset,
    pub = list(width = 180, max_height = 185, font_size = 7, min_font = 5, dpi = 600,
      wrap = 70, panel_w = NA, row_pad = 1.4, dot_max = 3.5, bar_w = 2.5, bar_h = 15, title = FALSE),
    screen = list(width = NA, max_height = Inf, font_size = 9, min_font = 9, dpi = 150,
      wrap = 0, panel_w = 90, row_pad = 1.6, dot_max = 5, bar_w = 4, bar_h = 30, title = TRUE));
  if (!is.null(args$width)) pr$width <- args$width;
  if (!is.null(args$panel_width)) {
    pr$panel_w <- args$panel_width;
    pr$width <- NA;
  }
  if (!is.null(args$max_height)) pr$max_height <- args$max_height;
  if (!is.null(args$font_size)) {
    pr$font_size <- args$font_size;
    pr$min_font <- min(pr$min_font, args$font_size);
  }
  if (!is.null(args$dpi)) pr$dpi <- args$dpi;
  if (!is.null(args$wrap)) pr$wrap <- args$wrap;
  pr$height <- if (is.null(args$height)) NA else args$height;
  pr;
};

check_plot_args <- function(args) {
  fmts <- tolower(trimws(strsplit(args$format, ",")[[1]]));
  fmts[fmts == "tif"] <- "tiff";
  bad <- setdiff(fmts, c("pdf", "png", "tiff"));
  if (length(bad) > 0) stop("未知的 --format: ", paste(bad, collapse = ", "));
  args$format <- unique(fmts);
  if (args$simplify != "none") {
    thr <- suppressWarnings(as.numeric(args$simplify));
    if (is.na(thr) || thr <= 0 || thr > 1) stop("--simplify 须为 none 或 (0, 1] 之间的数值");
  }
  if (!is.null(args$top_n) && args$top_n < 1) stop("--top_n 须 >= 1");
  shp <- suppressWarnings(as.integer(trimws(strsplit(args$shapes, ",")[[1]])));
  if (anyNA(shp) || any(!shp %in% 21:25) || anyDuplicated(shp)) stop("--shapes 须为 21~25 中不重复的整数，逗号分隔");
  args$shapes <- shp;
  if (!is.null(args$panel_width) && args$panel_width <= 0) stop("--panel_width 须 > 0");
  if (!is.null(args$label_width) && args$label_width <= 5) stop("--label_width 须 > 5 mm");
  if (!is.null(args$legend_width) && args$legend_width < 0) stop("--legend_width 须 >= 0");
  if (!is.null(args$panel_width) && !is.null(args$width)) message("[警告] 指定了 --panel_width，忽略 --width");
  if (!is.null(args$label_width) && !is.null(args$wrap)) message("[警告] 指定了 --label_width，term 名称按其宽度折行，忽略 --wrap");
  args$font_dir <- NULL;
  if (args$font == "Arial") {
    args$font_dir <- find_font_dir("Arial");
    if (!nzchar(Sys.which("gs")) || is.null(args$font_dir)) {
      message("[警告] 没有 ghostscript (gs) 或系统中没有 Arial 字体，--font 改用 Helvetica");
      args$font <- "Helvetica";
    } else {
      register_arial();
    }
  }
  args;
};

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

# 渐变色（低显著 -> 高显著）
get_palette <- function(args) {
  if (!is.null(args$colors)) {
    pal <- trimws(strsplit(args$colors, ",")[[1]]);
    ok <- vapply(pal, function(cc) !inherits(try(col2rgb(cc), silent = TRUE), "try-error"), logical(1));
    if (!all(ok)) stop("无法识别的颜色: ", paste(pal[!ok], collapse = ", "));
    if (length(pal) < 2) stop("--colors 至少需要 2 种颜色");
    return(pal);
  }
  if (tolower(args$palette) == "enrichplot") return(c("#327EBA", "#E06663"));
  vir <- c("viridis", "magma", "plasma", "inferno", "cividis", "mako", "rocket", "turbo");
  if (tolower(args$palette) %in% vir) return(rev(viridisLite::viridis(9, option = tolower(args$palette))));
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) stop("未知的色板: ", args$palette);
  info <- RColorBrewer::brewer.pal.info;
  if (!args$palette %in% rownames(info)) stop("未知的色板: ", args$palette);
  RColorBrewer::brewer.pal(min(9, info[args$palette, "maxcolors"]), args$palette);
};

# 基因列表：取每行第 1 列，忽略空行和 # 开头的行
read_id_list <- function(f) {
  if (!file.exists(f)) stop("文件不存在: ", f);
  x <- trimws(readLines(f, warn = FALSE));
  x <- x[nzchar(x) & !startsWith(x, "#")];
  unique(sub("[[:space:]].*$", "", x));
};

# 背景基因集
#   annotated: 至少有一个 term 的基因（clusterProfiler 默认）
#   all      : 注释文件中出现的全部基因（含无 term 的基因）
#   <文件>   : 文件中的基因（clusterProfiler 会再与有 term 的基因取交集）
set_universe <- function(t2g, allGenes, mode, minGS) {
  if (mode == "annotated") return(list(t2g = t2g, universe = NULL));
  if (mode == "all") {
    if (minGS < 2) stop("--universe all 需要 --min_gs >= 2");
    noTerm <- setdiff(allGenes, t2g[[2]]);
    if (length(noTerm) > 0) {
      # 无 term 的基因各放进一个单基因的占位 term：计入背景，但小于 min_gs 不参与检验
      ph <- data.table(sprintf("__NO_TERM_%d__", seq_along(noTerm)), noTerm);
      setnames(ph, names(t2g));
      t2g <- rbind(t2g, ph);
    }
    return(list(t2g = t2g, universe = NULL));
  }
  list(t2g = t2g, universe = read_id_list(mode));
};

# 运行 enricher 并整理结果（不含 Tag 列）
run_enricher <- function(genes, t2g, t2n, args) {
  uni <- set_universe(t2g, args$all_genes, args$universe, args$min_gs);
  ego <- clusterProfiler::enricher(genes, TERM2GENE = as.data.frame(uni$t2g),
    TERM2NAME = as.data.frame(t2n), universe = uni$universe, pAdjustMethod = "BH",
    pvalueCutoff = args$padj_cutoff, qvalueCutoff = args$padj_cutoff,
    minGSSize = args$min_gs, maxGSSize = args$max_gs);
  if (is.null(ego) || nrow(ego@result) == 0) return(NULL);
  res <- as.data.table(ego@result);
  res[!startsWith(ID, "__NO_TERM_")];
};

write_result <- function(res, prefix, cutoff) {
  res <- cbind(Tag = ifelse(!is.na(res$p.adjust) & res$p.adjust < cutoff, "OVER", ""), res);
  fn <- paste0(prefix, ".tsv");
  fwrite(res, fn, sep = "\t", quote = FALSE, na = "NA");
  message(sprintf("结果表 -> %s（检验 %d 个 term，p.adjust < %s 的 %d 个）",
    fn, nrow(res), format(cutoff), sum(res$Tag == "OVER")));
};

# 富集基因 Jaccard 相似度去冗余，返回保留与否
simplify_keep <- function(geneID, padj, setSize, thr) {
  gs <- strsplit(as.character(geneID), "/", fixed = TRUE);
  keep <- logical(length(gs));
  kept <- list();
  for (i in order(padj, setSize)) {
    redundant <- FALSE;
    for (g in kept) {
      if (length(intersect(gs[[i]], g)) / length(union(gs[[i]], g)) >= thr) {
        redundant <- TRUE;
        break;
      }
    }
    if (!redundant) {
      keep[i] <- TRUE;
      kept[[length(kept) + 1]] <- gs[[i]];
    }
  }
  keep;
};

# 图例中 Count 的整数刻度
count_breaks <- function(lim) {
  b <- unique(round(pretty(lim, n = 4)));
  b <- b[b >= lim[1] & b <= lim[2]];
  if (length(b) == 0) b <- round(mean(lim));
  b;
};

save_plot <- function(g, file, fmt, W, H, dpi, font = "Helvetica", fontDir = NULL) {
  if (fmt == "pdf") {
    pdf(file, width = W / 25.4, height = H / 25.4, useDingbats = FALSE, family = font,
      encoding = "WinAnsi.enc");
    grid.draw(g);
    invisible(dev.off());
    if (nzchar(Sys.which("gs"))) embed_pdf_fonts(file, fontDir);
  } else if (requireNamespace("ragg", quietly = TRUE)) {
    if (fmt == "png") {
      ragg::agg_png(file, width = W, height = H, units = "mm", res = dpi);
    } else {
      ragg::agg_tiff(file, width = W, height = H, units = "mm", res = dpi, compression = "lzw");
    }
    grid.draw(g);
    invisible(dev.off());
  } else if (nzchar(Sys.which("gs"))) {
    # 没有 ragg 时先画 pdf，再用 ghostscript 转为 png / tiff
    tmp <- tempfile(fileext = ".pdf");
    save_plot(g, tmp, "pdf", W, H, dpi, font, fontDir);
    gsArgs <- c("-q", "-dNOPAUSE", "-dBATCH", "-dSAFER", "-dTextAlphaBits=4", "-dGraphicsAlphaBits=4",
      if (fmt == "png") "-sDEVICE=png16m" else c("-sDEVICE=tiff24nc", "-sCompression=lzw"),
      paste0("-r", dpi), paste0("-sOutputFile=", shQuote(file)), shQuote(tmp));
    st <- system2("gs", gsArgs);
    unlink(tmp);
    if (st != 0) stop("ghostscript 转换 ", fmt, " 失败");
  } else {
    stop("输出 png/tiff 需要 ragg 包或 ghostscript (gs)");
  }
  if (!file.exists(file)) stop("未能生成 ", file);
};

# 按实际字宽折行，使每行宽度不超过 maxW (mm)；单个单词超长时单独成行
wrap_to_width <- function(s, maxW, gp) {
  tw <- function(x) convertWidth(grobWidth(textGrob(x, gp = gp)), "mm", valueOnly = TRUE);
  lines <- character(0);
  cur <- "";
  for (w in strsplit(s, " ", fixed = TRUE)[[1]]) {
    cand <- if (nzchar(cur)) paste(cur, w) else w;
    if (nzchar(cur) && tw(cand) > maxW) {
      lines <- c(lines, cur);
      cur <- w;
    } else {
      cur <- cand;
    }
  }
  paste(c(lines, cur), collapse = "\n");
};

# 点图：res 为富集结果表，facetCol 为分面列名（不存在时不分面）
plot_enrich <- function(res, facetCol, facetLevels, kind, prefix, args) {
  pr <- get_preset(args);
  fs <- pr$font_size;
  afs <- fs - 1;
  pal <- get_palette(args);

  d <- as.data.table(res);
  need <- c("ID", "Description", "BgRatio", "p.adjust", "Count");
  miss <- setdiff(need, names(d));
  if (length(miss) > 0) stop("结果表缺少列: ", paste(miss, collapse = ", "));
  d[, `:=`(Count = as.numeric(Count), p.adjust = as.numeric(p.adjust),
    Description = as.character(Description), M_ = as.numeric(sub("/.*", "", BgRatio)))];
  # x 轴取 --x_col 列；"a/b" 形式（如 GeneRatio）按分数计算
  if (!args$x_col %in% names(d)) {
    stop(sprintf("结果表中没有 --x_col 指定的列 %s；可选的列: %s", args$x_col, paste(names(d), collapse = ", ")));
  }
  xr <- as.character(d[[args$x_col]]);
  xv <- suppressWarnings(if (any(grepl("/", xr, fixed = TRUE))) {
    vapply(strsplit(xr, "/", fixed = TRUE), function(z) as.numeric(z[1]) / as.numeric(z[2]), 0)
  } else as.numeric(xr));
  if (all(is.na(xv))) stop(sprintf("--x_col 列 %s 不是数值（或 a/b 形式的比例）", args$x_col));
  d[, xval := xv];
  hasFacet <- !isTRUE(args$no_facet) && !is.null(facetCol) && facetCol %in% names(d);
  if (hasFacet) {
    f <- as.character(d[[facetCol]]);
    f[is.na(f) | !nzchar(f)] <- "Other";
    d[, Facet := f];
  } else {
    d[, Facet := "All"];
  }
  useShape <- isTRUE(args$shape_by_class) && !is.null(facetCol) && facetCol %in% names(d);
  if (isTRUE(args$shape_by_class) && !useShape) message(sprintf("[%s] 没有类别列，忽略 --shape_by_class", kind));
  if (useShape) {
    s <- as.character(d[[facetCol]]);
    s[is.na(s) | !nzchar(s)] <- "Other";
    d[, Class_ := s];
  }
  shapeTitle <- switch(kind, GO = "Ontology", KEGG = "KEGG class", "Class");

  # 选择上图的 term
  sel <- d[!is.na(p.adjust) & p.adjust < args$padj_cutoff & Count >= args$min_count & !is.na(xval)];
  if (nrow(sel) == 0) {
    message(sprintf("[%s] 没有 p.adjust < %s 且 Count >= %d 的 term，不作图",
      kind, format(args$padj_cutoff), args$min_count));
    return(invisible(NULL));
  }
  lev <- c(intersect(facetLevels, sel$Facet), setdiff(unique(sel$Facet), facetLevels));
  sel[, Facet := factor(Facet, levels = lev)];
  if (args$simplify != "none") {
    if (!"geneID" %in% names(sel)) stop("--simplify 需要结果表中的 geneID 列");
    n0 <- nrow(sel);
    sel[, keep_ := simplify_keep(geneID, p.adjust, M_, as.numeric(args$simplify)), by = Facet];
    sel <- sel[keep_ == TRUE];
    sel[, keep_ := NULL];
    message(sprintf("[%s] 去冗余 (Jaccard >= %s)：%d -> %d 个 term", kind, args$simplify, n0, nrow(sel)));
  }
  if (!is.null(args$top_n)) {
    setorder(sel, Facet, p.adjust, -xval);
    sel <- sel[, head(.SD, args$top_n), by = Facet];
  }
  if (args$order_by == "x") setorder(sel, Facet, -xval, p.adjust) else setorder(sel, Facet, p.adjust, -xval);
  sel[, Facet := droplevels(Facet)];
  if (useShape) {
    lev <- c(intersect(facetLevels, sel$Class_), setdiff(unique(sel$Class_), facetLevels));
    if (length(lev) > length(args$shapes)) {
      stop(sprintf("上图的类别有 %d 种，--shapes 只有 %d 个形状；可增加 --shapes 或用 --drop_class 排除类别", length(lev), length(args$shapes)));
    }
    sel[, Class_ := factor(Class_, levels = lev)];
  }
  nF <- nlevels(sel$Facet);

  # y 位置：每个 term 占的高度随其名称行数变化；各分面用互不重叠的数值区间，便于共用刻度
  # --label_width：按实际字宽折行，可用宽度扣除左边距 2 mm、刻度线 0.8 mm 与刻度文字间距
  labGp <- gpar(fontsize = fs, fontfamily = args$font);
  sel[, Label := if (!is.null(args$label_width)) {
    vapply(Description, wrap_to_width, "", maxW = args$label_width - 4, gp = labGp)
  } else if (pr$wrap > 0) {
    vapply(Description, function(s) paste(strwrap(s, width = pr$wrap), collapse = "\n"), "")
  } else Description];
  sel[, unit_ := (lengths(regmatches(Label, gregexpr("\n", Label))) + pr$row_pad) / pr$row_pad];
  sel[, yEnd_ := cumsum(unit_), by = Facet];
  facU <- sel[, .(U = sum(unit_)), by = Facet];
  sel[, off_ := (as.integer(Facet) - 1) * 1e5];
  sel[, y := off_ + facU$U[match(Facet, facU$Facet)] - (yEnd_ - unit_ / 2)];
  sel[, mlog := -log10(pmax(p.adjust, .Machine$double.xmin))];

  # padU：分面名称比面板长时，面板上下各留出的空白（以行为单位）
  # softHy：R 的 pdf 设备把 "-" 画成减号，经 pdf 输出时改用软连字符 (U+00AD) 才显示为连字符
  build <- function(fy, rowMM, padU, stripLab, softHy = FALSE, legBox = "vertical") {
    hy <- function(x) if (softHy) gsub("-", "\u00ad", x, fixed = TRUE) else x;
    dMax <- min(pr$dot_max, 0.85 * rowMM);
    off <- (seq_len(nF) - 1) * 1e5;
    lims <- data.table(Facet = factor(rep(levels(sel$Facet), 2), levels = levels(sel$Facet)),
      y = c(off - padU / 2, off + facU$U + padU / 2));
    p <- ggplot(sel, aes(x = xval, y = y)) +
      geom_blank(data = lims, aes(y = y), inherit.aes = FALSE) +
      (if (useShape) {
        geom_point(aes(size = Count, fill = mlog, shape = Class_), colour = "black", stroke = lw)
      } else geom_point(aes(size = Count, fill = mlog), shape = 21, colour = "black", stroke = lw)) +
      facet_grid(Facet ~ ., scales = "free_y",
        labeller = as_labeller(setNames(hy(stripLab), levels(sel$Facet)))) +
      scale_y_continuous(breaks = sel$y, expand = c(0, 0),
        labels = function(b) hy(sel$Label[match(round(b, 6), round(sel$y, 6))])) +
      scale_x_continuous(expand = expansion(mult = 0.08)) +
      scale_fill_gradientn(colours = pal, name = expression("-" * log[10] * "(" * italic(P)[adj] * ")"),
        guide = guide_colourbar(order = 1, theme = theme(
          legend.key.width = unit(pr$bar_w, "mm"), legend.key.height = unit(pr$bar_h, "mm"),
          legend.ticks = element_blank()))) +
      scale_size(range = c(max(0.5, dMax * 0.35), dMax) / 0.75, breaks = count_breaks,
        name = "Count", guide = guide_legend(order = 2, override.aes = list(fill = NA, shape = 21),
          theme = theme(legend.key.height = unit(dMax + 0.5, "mm"),
            legend.key.width = unit(dMax + 0.5, "mm"), legend.key.spacing.y = unit(0.8, "mm")))) +
      (if (useShape) {
        scale_shape_manual(values = args$shapes[seq_len(nlevels(sel$Class_))], name = shapeTitle, labels = hy,
          guide = guide_legend(order = 3, override.aes = list(fill = NA, size = 0.8 * dMax / 0.75),
            theme = theme(legend.key.height = unit(dMax + 0.5, "mm"),
              legend.key.width = unit(dMax + 0.5, "mm"), legend.key.spacing.y = unit(0.8, "mm"))))
      }) +
      labs(x = hy(args$x_col), y = NULL) +
      theme_bw(base_size = fs, base_family = args$font) +
      theme(
        panel.border = element_rect(colour = "black", linewidth = lw, fill = NA),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = if (args$grid) element_line(colour = "grey90", linewidth = lw) else element_blank(),
        panel.spacing = unit(1.5, "mm"),
        axis.ticks = element_line(colour = "black", linewidth = lw),
        axis.ticks.length = unit(0.8, "mm"),
        axis.text.x = element_text(size = afs, colour = "black"),
        axis.text.y = element_text(size = fy, colour = "black", lineheight = 0.9),
        axis.title.x = element_text(size = fs),
        strip.background = element_rect(fill = "grey90", colour = "black", linewidth = lw),
        strip.text.y.right = if (hasFacet) {
          element_text(size = fs, angle = -90, colour = "black", lineheight = 0.9,
            margin = margin(0.8, 0.8, 0.8, 0.8, "mm"))
        } else element_blank(),
        legend.title = element_text(size = fs),
        legend.text = element_text(size = afs),
        legend.spacing.y = unit(3, "mm"),
        legend.spacing.x = unit(3, "mm"),
        legend.box = legBox,
        legend.box.just = if (legBox == "horizontal") "top" else "left",
        plot.title = element_text(size = fs + 1),
        plot.margin = margin(2, 2, 2, 2, "mm"));
    if (pr$title) {
      p <- p + ggtitle(hy(sprintf("%s enrichment: %s (%d terms, p.adjust < %s)",
        kind, basename(prefix), nrow(sel), format(args$padj_cutoff))));
    }
    g <- ggplotGrob(p);
    pIdx <- grep("^panel", g$layout$name);
    pRows <- sort(g$layout$t[pIdx]);
    pCol <- unique(g$layout$l[pIdx]);
    # 面板自上而下即分面水平顺序
    for (j in seq_along(pRows)) g$heights[pRows[j]] <- unit((facU$U[j] + padU[j]) * rowMM, "mm");
    if (is.na(pr$width)) g$widths[pCol] <- unit(pr$panel_w, "mm");
    # --label_width：面板左侧不足的宽度补在左边距上，使面板左缘位于固定位置
    leftW <- sum(convertWidth(g$widths[seq_len(pCol - 1)], "mm", valueOnly = TRUE));
    if (!is.null(args$label_width) && leftW < args$label_width) {
      g$widths[1] <- unit(convertWidth(g$widths[1], "mm", valueOnly = TRUE) + args$label_width - leftW, "mm");
    }
    # --legend_width：面板右侧不足的宽度补在右边距上
    nC <- length(g$widths);
    rightW <- sum(convertWidth(g$widths[seq(pCol + 1, nC)], "mm", valueOnly = TRUE));
    if (!is.null(args$legend_width) && rightW < args$legend_width) {
      g$widths[nC] <- unit(convertWidth(g$widths[nC], "mm", valueOnly = TRUE) + args$legend_width - rightW, "mm");
    }
    list(g = g, pRows = pRows, pCol = pCol, leftW = leftW, rightW = rightW);
  };

  # 先按基础字号算出面板以外的固定高度，再决定 term 名称字号与行高
  rowMM_of <- function(f) f * pt2mm * pr$row_pad;
  b0 <- build(fs, rowMM_of(fs), rep(0, nF), levels(sel$Facet));
  fixedH <- sum(convertHeight(b0$g$heights[-b0$pRows], "mm", valueOnly = TRUE));
  fit_rows <- function(Ut) {
    fy <- fs;
    warn <- "";
    if (!is.na(pr$height)) {
      rowMM <- (pr$height - fixedH) / Ut;
      fy <- min(fs, rowMM / (pt2mm * pr$row_pad));
      if (fy < pr$min_font) {
        fy <- pr$min_font;
        warn <- "height";
      }
    } else {
      rowMM <- rowMM_of(fs);
      if (fixedH + Ut * rowMM > pr$max_height) {
        rowMM <- (pr$max_height - fixedH) / Ut;
        fy <- rowMM / (pt2mm * pr$row_pad);
        if (fy < pr$min_font) {
          # 字号已到下限：先把行距压缩到 1.2 倍字高，仍放不下才超出高度上限
          fy <- pr$min_font;
          if (rowMM < fy * pt2mm * 1.2) {
            rowMM <- fy * pt2mm * 1.2;
            warn <- "max";
          }
        }
      }
    }
    list(rowMM = rowMM, fy = fy, warn = warn);
  };
  # 分面名称竖排在面板右侧：面板高度至少要容纳名称中最长的单词
  stripGp <- gpar(fontsize = fs, fontfamily = args$font);
  needMM <- if (hasFacet) {
    vapply(levels(sel$Facet), function(s) {
      max(vapply(strsplit(s, " ", fixed = TRUE)[[1]], function(w) {
        convertWidth(grobWidth(textGrob(w, gp = stripGp)), "mm", valueOnly = TRUE);
      }, 0)) + 2;
    }, 0);
  } else rep(0, nF);
  padU <- rep(0, nF);
  for (it in 1:10) {
    fr <- fit_rows(sum(facU$U + padU));
    newPad <- pmax(0, needMM / fr$rowMM - facU$U);
    if (max(abs(newPad - padU)) < 1e-6) break;
    padU <- newPad;
  }
  fr <- fit_rows(sum(facU$U + padU));
  rowMM <- fr$rowMM;
  fy <- fr$fy;
  if (fr$warn == "height") {
    message(sprintf("[警告] --height %.0f mm 下 term 名称即使用 %.1f pt 也会重叠，请增大高度或减少 term", pr$height, fy));
  } else if (fr$warn == "max") {
    message(sprintf("[警告] %d 个 term 在 %.0f mm 高度内放不下（最小字号 %.1f pt），画布超出上限；可用 --top_n 或 --min_count 减少 term",
      nrow(sel), pr$max_height, fy));
  }
  if (fy < fs) message(sprintf("[%s] term 名称字号缩小为 %.2f pt", kind, fy));
  stripLab <- levels(sel$Facet);
  if (hasFacet) {
    stripLab <- vapply(seq_len(nF), function(j) {
      wrap_to_width(stripLab[j], (facU$U[j] + padU[j]) * rowMM - 2, stripGp);
    }, "");
  }
  legH_of <- function(g) {
    gb <- grep("^guide-box", g$layout$name);
    if (length(gb) == 0) return(0);
    max(vapply(gb, function(i) convertHeight(grobHeight(g$grobs[[i]]), "mm", valueOnly = TRUE), 0));
  };
  H <- if (is.na(pr$height)) fixedH + sum(facU$U + padU) * rowMM else pr$height;
  legBox <- "vertical";
  b <- build(fy, rowMM, padU, stripLab);
  legH <- legH_of(b$g);
  # term 少、上下叠放的图例比图还高时，图例改为左右并排（指定 --legend_width 时须放得下）
  if (H < legH + 6) {
    b2 <- build(fy, rowMM, padU, stripLab, legBox = "horizontal");
    legH2 <- legH_of(b2$g);
    if (legH2 < legH && (is.null(args$legend_width) || b2$rightW <= args$legend_width + 0.01)) {
      legBox <- "horizontal";
      b <- b2;
      legH <- legH2;
    }
  }
  g <- b$g;
  if (!is.null(args$label_width) && b$leftW > args$label_width + 0.01) {
    message(sprintf("[警告] term 名称中有无法折行的长单词，名称区需要 %.1f mm，超出 --label_width %.1f mm，面板位置与其他图不一致",
      b$leftW, args$label_width));
  }
  if (!is.null(args$legend_width) && b$rightW > args$legend_width + 0.01) {
    message(sprintf("[警告] 分面名称与图例需要 %.1f mm，超出 --legend_width %.1f mm，画布右侧与其他图不一致",
      b$rightW, args$legend_width));
  }
  if (is.na(pr$height) && H < legH + 6) H <- legH + 6;
  fixedW <- sum(convertWidth(g$widths[-b$pCol], "mm", valueOnly = TRUE));
  if (is.na(pr$width)) {
    W <- fixedW + pr$panel_w;
  } else {
    W <- pr$width;
    if (W - fixedW < 25) {
      message(sprintf("[警告] 面板宽度只剩 %.1f mm：term 名称过长，可减小 --wrap 或增大 --width", W - fixedW));
    }
  }

  useRagg <- requireNamespace("ragg", quietly = TRUE);
  gSoft <- NULL;
  for (fmt in args$format) {
    fn <- paste0(prefix, ".", fmt);
    # pdf 以及没有 ragg 时经 pdf 转换的 png / tiff 用软连字符版本
    if (fmt == "pdf" || !useRagg) {
      if (is.null(gSoft)) gSoft <- build(fy, rowMM, padU, stripLab, softHy = TRUE, legBox = legBox)$g;
      save_plot(gSoft, fn, fmt, W, H, pr$dpi, args$font, args$font_dir);
    } else {
      save_plot(g, fn, fmt, W, H, pr$dpi, args$font, args$font_dir);
    }
    message(sprintf("[%s] %d 个 term -> %s (%.1f x %.1f mm)", kind, nrow(sel), fn, W, H));
  }

  # 上图 term 列表，顺序同图中自上而下
  out <- copy(sel);
  out[, Plot_order := seq_len(.N)];
  drop <- intersect(c("xval", "Label", "unit_", "yEnd_", "off_", "y", "mlog", "M_", "Class_"), names(out));
  out[, (drop) := NULL];
  setcolorder(out, c("Plot_order", "Facet"));
  fn <- paste0(prefix, "_plotted.tsv");
  fwrite(out, fn, sep = "\t", quote = FALSE, na = "NA");
  message(sprintf("[%s] 上图 term 列表 -> %s", kind, fn));
};

#### <<< 共用代码结束 ####

##############################
# 命令行参数
##############################
parser <- ArgumentParser(description = "KEGG pathway enrichment (clusterProfiler::enricher) and dot plot.");
parser$add_argument("-g", "--gene_list", default = NULL,
  help = "【富集模式必需】基因列表文件，取每行第 1 列（忽略空行和 # 开头的行）");
parser$add_argument("-a", "--kegg_annot_2col", default = NULL,
  help = "【富集模式必需】KO 注释文件，2 列无标题：geneID, KO ID（如 K03654，无 KO 时留空）");
parser$add_argument("-p", "--output_prefix", default = NULL,
  help = "输出前缀：<prefix>.tsv 结果表、<prefix>.<format> 点图、<prefix>_plotted.tsv 上图 term [default: 富集模式 out；--from_tsv 模式为该 tsv 去掉扩展名]");
parser$add_argument("--from_tsv", default = NULL,
  help = "直接读已有的富集结果表（本脚本或旧版脚本输出的 <prefix>.tsv）作图，不重跑富集；表中没有 KEGG_class_A 列时从 KEGG 缓存补充分类");
parser$add_argument("--universe", default = "all",
  help = "背景基因集：all = 注释文件中的全部基因（含无 KO、或 KO 不属于任何通路的基因）；annotated = 至少注释到一个通路的基因（旧版脚本的做法）；或一个基因列表文件 [default: %(default)s]");
parser$add_argument("--min_gs", type = "integer", default = 3,
  help = "参与检验的通路最少基因数 (minGSSize) [default: %(default)s]");
parser$add_argument("--max_gs", type = "integer", default = 500,
  help = "参与检验的通路最多基因数 (maxGSSize) [default: %(default)s]");
parser$add_argument("--drop_class", default = NULL,
  help = "检验前排除这些 KEGG 分类（一级或二级分类名）下的通路，分号分隔，如 \"Human Diseases;Drug Development\"；--from_tsv 模式下只在作图时排除");
parser$add_argument("--kegg_cache", default = "kegg_ko_pathway.rds",
  help = "KEGG 数据缓存文件：存在时直接读取，否则联网下载并保存；设为 none 则每次下载、不保存 [default: %(default)s]");
parser$add_argument("--kegg_source", default = "rest", choices = c("rest", "KEGGREST"),
  help = "KEGG 数据下载方式：rest = 用 R 自带的 readLines 直接访问 KEGG REST API（https://rest.kegg.jp）；KEGGREST = 用 KEGGREST 包（keggLink / keggList / keggGet）。两者访问的是同一个接口，得到的数据相同 [default: %(default)s]");
parser$add_argument("--kegg_refresh", action = "store_true", default = FALSE,
  help = "忽略已有缓存，重新下载 KEGG 数据并覆盖缓存");
add_plot_args(parser);
args <- parser$parse_args();
if (is.null(args$from_tsv) && (is.null(args$gene_list) || is.null(args$kegg_annot_2col))) {
  parser$print_help();
  message("\n[错误] 必须指定 --gene_list 与 --kegg_annot_2col，或用 --from_tsv 指定已有结果表");
  quit(save = "no", status = 1);
}
args <- check_plot_args(args);
if (is.null(args$output_prefix)) {
  args$output_prefix <- if (is.null(args$from_tsv)) "out" else sub("\\.[^./]*$", "", args$from_tsv);
}
dropCls <- if (is.null(args$drop_class)) character(0) else trimws(strsplit(args$drop_class, ";")[[1]]);
dropCls <- dropCls[nzchar(dropCls)];

kegg_levels <- c("Metabolism", "Genetic Information Processing", "Environmental Information Processing",
  "Cellular Processes", "Organismal Systems", "Human Diseases", "Drug Development");

##############################
# KEGG 数据：KO -> pathway、通路名称、通路分类
##############################
kegg_rest <- function(op) {
  u <- paste0("https://rest.kegg.jp/", op);
  x <- tryCatch(readLines(u, warn = FALSE), error = function(e) stop("无法下载 ", u, "：", conditionMessage(e)));
  if (length(x) == 0) stop("下载内容为空: ", u);
  x;
};

get_kegg_db <- function(cacheFile, refresh, source) {
  useCache <- !is.null(cacheFile) && cacheFile != "none";
  if (useCache && !refresh && file.exists(cacheFile)) {
    db <- readRDS(cacheFile);
    message(sprintf("读取 KEGG 缓存 %s（下载于 %s，来源 %s）", cacheFile, db$date,
      if (is.null(db$source)) "rest" else db$source));
    return(db);
  }
  if (source == "KEGGREST") {
    if (!requireNamespace("KEGGREST", quietly = TRUE)) stop("--kegg_source KEGGREST 需要先安装 KEGGREST 包");
    message("用 KEGGREST 包下载 KO-pathway 对应关系、通路名称与分类 ...");
    lk <- KEGGREST::keggLink("pathway", "ko");
    link <- data.table(ko = sub("^ko:", "", names(lk)), pathway = sub("^path:", "", unname(lk)));
    ds <- KEGGREST::keggList("pathway", "ko");
    desc <- data.table(pathway = sub("^path:", "", names(ds)), desc = unname(ds));
    br <- strsplit(KEGGREST::keggGet("br:br08901"), "\n", fixed = TRUE)[[1]];
  } else {
    message("从 KEGG REST API 下载 KO-pathway 对应关系、通路名称与分类 ...");
    lk <- tstrsplit(kegg_rest("link/pathway/ko"), "\t", fixed = TRUE);
    link <- data.table(ko = sub("^ko:", "", lk[[1]]), pathway = sub("^path:", "", lk[[2]]));
    ds <- tstrsplit(kegg_rest("list/pathway/ko"), "\t", fixed = TRUE);
    desc <- data.table(pathway = sub("^path:", "", ds[[1]]), desc = ds[[2]]);
    br <- kegg_rest("get/br:br08901");
  }
  link <- unique(link[startsWith(pathway, "ko")]);
  # br08901 层级文本：A 一级分类，B 二级分类，C 通路（5 位编号 + 名称）
  tp <- substr(br, 1, 1);
  txt <- trimws(gsub("<[^>]+>", "", substring(br, 2)));
  lastIdx <- function(is) {
    i <- cummax(ifelse(is, seq_along(is), 0));
    i[i == 0] <- NA;
    i;
  };
  aI <- lastIdx(tp == "A");
  bI <- lastIdx(tp == "B");
  cls <- data.table(map = sub("\\s.*$", "", txt), class_A = txt[aI], class_B = txt[bI])[tp == "C"];
  cls <- cls[!duplicated(map)];
  db <- list(link = link, desc = desc, class = cls, date = format(Sys.time(), "%Y-%m-%d"), source = source);
  message(sprintf("  KO-pathway %d 对，通路 %d 条，有分类的通路 %d 条",
    nrow(link), nrow(desc), nrow(cls)));
  if (useCache) {
    saveRDS(db, cacheFile);
    message("  已保存缓存 ", cacheFile);
  }
  db;
};

# 通路 ID（ko05410）对应的分类
pathway_class <- function(ids, db) {
  i <- match(sub("^(ko|map)", "", ids), db$class$map);
  list(A = ifelse(is.na(i), "", db$class$class_A[i]), B = ifelse(is.na(i), "", db$class$class_B[i]));
};

##############################
# 富集分析
##############################
if (is.null(args$from_tsv)) {
  suppressPackageStartupMessages(library(clusterProfiler));
  genes <- read_id_list(args$gene_list);
  kt <- fread(args$kegg_annot_2col, header = FALSE, sep = "\t", quote = "", fill = TRUE,
    colClasses = "character");
  if (ncol(kt) < 2) kt[, V2 := ""];
  kt <- kt[, 1:2];
  setnames(kt, c("gene", "ko"));
  for (cc in names(kt)) set(kt, which(is.na(kt[[cc]])), cc, "");
  kt <- kt[nzchar(gene)];
  kt[, ko := sub("^ko:", "", trimws(ko))];
  args$all_genes <- unique(kt$gene);
  nMiss <- sum(!genes %in% args$all_genes);
  message(sprintf("输入基因 %d 个，其中 %d 个不在注释文件中", length(genes), nMiss));

  kdb <- get_kegg_db(args$kegg_cache, args$kegg_refresh, args$kegg_source);
  t2g <- unique(merge(kt[nzchar(ko)], kdb$link, by = "ko")[, .(pathway, gene)]);
  if (length(dropCls) > 0) {
    pc <- pathway_class(t2g$pathway, kdb);
    bad <- setdiff(dropCls, c(kdb$class$class_A, kdb$class$class_B));
    if (length(bad) > 0) message("[警告] KEGG 中没有这些分类名，已忽略: ", paste(bad, collapse = "; "));
    isDrop <- pc$A %in% dropCls | pc$B %in% dropCls;
    message(sprintf("--drop_class：排除 %d 条通路", uniqueN(t2g$pathway[isDrop])));
    t2g <- t2g[!isDrop];
  }
  t2n <- kdb$desc[pathway %in% t2g$pathway];
  message(sprintf("有通路注释的基因 %d 个，通路 %d 条", uniqueN(t2g$gene), uniqueN(t2g$pathway)));

  res <- run_enricher(genes, t2g, t2n, args);
  if (is.null(res)) {
    message("[警告] 没有可检验的通路（输入基因均无通路注释？），不输出结果");
    quit(save = "no", status = 0);
  }
  res[, GeneRatio_in_term := paste0(Count, "/", as.numeric(sub("/.*", "", BgRatio)))];
  pc <- pathway_class(res$ID, kdb);
  res[, `:=`(KEGG_class_A = pc$A, KEGG_class_B = pc$B)];
  write_result(res, args$output_prefix, args$padj_cutoff);
} else {
  res <- fread(args$from_tsv, sep = "\t", quote = "", header = TRUE);
  message(sprintf("读入结果表 %s（%d 个 term）", args$from_tsv, nrow(res)));
  if (!"KEGG_class_A" %in% names(res) || (length(dropCls) > 0 && !"KEGG_class_B" %in% names(res))) {
    kdb <- tryCatch(get_kegg_db(args$kegg_cache, args$kegg_refresh, args$kegg_source), error = function(e) {
      message("[警告] ", conditionMessage(e), "；不按 KEGG 分类分面");
      NULL;
    });
    if (!is.null(kdb)) {
      pc <- pathway_class(res$ID, kdb);
      res[, `:=`(KEGG_class_A = pc$A, KEGG_class_B = pc$B)];
    }
  }
  if (length(dropCls) > 0 && "KEGG_class_A" %in% names(res)) {
    isDrop <- res$KEGG_class_A %in% dropCls | res$KEGG_class_B %in% dropCls;
    message(sprintf("--drop_class：作图时排除 %d 条通路（p.adjust 仍为原结果表中的值）", sum(isDrop)));
    res <- res[!isDrop];
  }
}

##############################
# 点图
##############################
pdf(NULL);
plot_enrich(res, "KEGG_class_A", kegg_levels, "KEGG", args$output_prefix, args);
invisible(dev.off());
