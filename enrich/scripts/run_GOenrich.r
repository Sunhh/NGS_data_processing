#!/usr/bin/env Rscript
# GO 富集分析（clusterProfiler::enricher）并画点图
#   - 富集结果表 <prefix>.tsv 的列与旧版脚本一致
#   - 点图：x = GeneRatio_in_term 等（--x_col 指定列），点大小 = Count，颜色 = -log10(p.adjust)，按 BP / MF / CC 分面（--no_facet 不分面）
#   - preset：pub（Nature Genetics 规格，宽 180 mm、高不超过 185 mm）与 screen（放大浏览）
#   - --from_tsv：直接读已有结果表作图，不重跑富集
#
# 用法：
#   Rscript run_GOenrich.r --gene_list CPMMV_up.list --go_annot_4col MEAM1_v1.2-b2g.annot-GOinEnrich --output_prefix CPMMV_up-GOen
#   Rscript run_GOenrich.r --from_tsv CPMMV_up-GOen.tsv --preset screen --format pdf,png

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
  parser$add_argument("--height", type = "double", default = NULL,
    help = "画布高度 (mm)；指定后不再自动计算 [default: 按 term 数自动计算]");
  parser$add_argument("--max_height", type = "double", default = NULL,
    help = "pub 自动高度的上限 (mm)；放不下时 term 名称字号自动缩小（最小 5 pt），仍放不下则超出上限并给出警告 [default: 185]");
  parser$add_argument("--font_size", type = "double", default = NULL,
    help = "基础字号 (pt)，刻度与图例文字比它小 1 pt [default: pub 7；screen 9]");
  parser$add_argument("--dpi", type = "integer", default = NULL,
    help = "png/tiff 分辨率 [default: pub 600；screen 150]");
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
  args;
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

save_plot <- function(g, file, fmt, W, H, dpi) {
  if (fmt == "pdf") {
    pdf(file, width = W / 25.4, height = H / 25.4, useDingbats = FALSE, family = "Helvetica",
      encoding = "WinAnsi.enc");
    grid.draw(g);
    invisible(dev.off());
    if (nzchar(Sys.which("gs"))) {
      tmp <- paste0(file, ".tmp.pdf");
      embedFonts(file, outfile = tmp);
      invisible(file.rename(tmp, file));
    }
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
    save_plot(g, tmp, "pdf", W, H, dpi);
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
  nF <- nlevels(sel$Facet);

  # y 位置：每个 term 占的高度随其名称行数变化；各分面用互不重叠的数值区间，便于共用刻度
  sel[, Label := if (pr$wrap > 0) {
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
      geom_point(aes(size = Count, fill = mlog), shape = 21, colour = "black", stroke = lw) +
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
        name = "Count", guide = guide_legend(order = 2, override.aes = list(fill = NA),
          theme = theme(legend.key.height = unit(dMax + 0.5, "mm"),
            legend.key.width = unit(dMax + 0.5, "mm"), legend.key.spacing.y = unit(0.8, "mm")))) +
      labs(x = hy(args$x_col), y = NULL) +
      theme_bw(base_size = fs, base_family = "Helvetica") +
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
    list(g = g, pRows = pRows, pCol = pCol);
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
  stripGp <- gpar(fontsize = fs, fontfamily = "Helvetica");
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
  # term 少、上下叠放的图例比图还高时，两个图例改为左右并排
  if (H < legH + 6) {
    b2 <- build(fy, rowMM, padU, stripLab, legBox = "horizontal");
    legH2 <- legH_of(b2$g);
    if (legH2 < legH) {
      legBox <- "horizontal";
      b <- b2;
      legH <- legH2;
    }
  }
  g <- b$g;
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
      save_plot(gSoft, fn, fmt, W, H, pr$dpi);
    } else {
      save_plot(g, fn, fmt, W, H, pr$dpi);
    }
    message(sprintf("[%s] %d 个 term -> %s (%.1f x %.1f mm)", kind, nrow(sel), fn, W, H));
  }

  # 上图 term 列表，顺序同图中自上而下
  out <- copy(sel);
  out[, Plot_order := seq_len(.N)];
  drop <- intersect(c("xval", "Label", "unit_", "yEnd_", "off_", "y", "mlog", "M_"), names(out));
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
parser <- ArgumentParser(description = "GO enrichment (clusterProfiler::enricher) and dot plot.");
parser$add_argument("-g", "--gene_list", default = NULL,
  help = "【富集模式必需】基因列表文件，取每行第 1 列（忽略空行和 # 开头的行）");
parser$add_argument("-a", "--go_annot_4col", default = NULL,
  help = "【富集模式必需】GO 注释文件，4 列无标题：geneID, GO ID, GO root (Biological Process 等), GO name；无 GO 的基因可只有 geneID");
parser$add_argument("-p", "--output_prefix", default = NULL,
  help = "输出前缀：<prefix>.tsv 结果表、<prefix>.<format> 点图、<prefix>_plotted.tsv 上图 term [default: 富集模式 out；--from_tsv 模式为该 tsv 去掉扩展名]");
parser$add_argument("--from_tsv", default = NULL,
  help = "直接读已有的富集结果表（本脚本或旧版脚本输出的 <prefix>.tsv）作图，不重跑富集");
parser$add_argument("--universe", default = "all",
  help = "背景基因集：all = 注释文件中的全部基因（含无 GO 的基因，同旧版脚本）；annotated = 至少有一个 GO 的基因；或一个基因列表文件 [default: %(default)s]");
parser$add_argument("--min_gs", type = "integer", default = 3,
  help = "参与检验的 term 最少基因数 (minGSSize) [default: %(default)s]");
parser$add_argument("--max_gs", type = "integer", default = 500,
  help = "参与检验的 term 最多基因数 (maxGSSize) [default: %(default)s]");
add_plot_args(parser);
args <- parser$parse_args();
if (is.null(args$from_tsv) && (is.null(args$gene_list) || is.null(args$go_annot_4col))) {
  parser$print_help();
  message("\n[错误] 必须指定 --gene_list 与 --go_annot_4col，或用 --from_tsv 指定已有结果表");
  quit(save = "no", status = 1);
}
args <- check_plot_args(args);
if (is.null(args$output_prefix)) {
  args$output_prefix <- if (is.null(args$from_tsv)) "out" else sub("\\.[^./]*$", "", args$from_tsv);
}

go_levels <- c("Biological Process", "Molecular Function", "Cellular Component");

##############################
# 富集分析
##############################
if (is.null(args$from_tsv)) {
  suppressPackageStartupMessages(library(clusterProfiler));
  genes <- read_id_list(args$gene_list);
  ann <- fread(args$go_annot_4col, header = FALSE, sep = "\t", quote = "", fill = TRUE,
    colClasses = "character");
  if (ncol(ann) < 4) stop("GO 注释文件至少需要 4 列");
  ann <- ann[, 1:4];
  setnames(ann, c("Gene", "GO_ID", "GO_root", "GO_name"));
  for (cc in names(ann)) set(ann, which(is.na(ann[[cc]])), cc, "");
  ann <- ann[nzchar(Gene)];
  t2g <- unique(ann[nzchar(GO_ID), .(GO_ID, Gene)]);
  info <- ann[nzchar(GO_ID)][!duplicated(GO_ID), .(GO_ID, GO_root, GO_name)];
  args$all_genes <- unique(ann$Gene);

  nMiss <- sum(!genes %in% args$all_genes);
  message(sprintf("输入基因 %d 个，其中 %d 个不在注释文件中", length(genes), nMiss));

  res <- run_enricher(genes, t2g, info[, .(GO_ID, GO_name)], args);
  if (is.null(res)) {
    message("[警告] 没有可检验的 GO term（输入基因均无注释？），不输出结果");
    quit(save = "no", status = 0);
  }
  ii <- match(res$ID, info$GO_ID);
  res[, `:=`(GO_root = info$GO_root[ii], GO_name = info$GO_name[ii])];
  res[, GeneRatio_in_term := paste0(Count, "/", as.numeric(sub("/.*", "", BgRatio)))];
  write_result(res, args$output_prefix, args$padj_cutoff);
} else {
  res <- fread(args$from_tsv, sep = "\t", quote = "", header = TRUE);
  message(sprintf("读入结果表 %s（%d 个 term）", args$from_tsv, nrow(res)));
}

##############################
# 点图
##############################
pdf(NULL);
plot_enrich(res, "GO_root", go_levels, "GO", args$output_prefix, args);
invisible(dev.off());
