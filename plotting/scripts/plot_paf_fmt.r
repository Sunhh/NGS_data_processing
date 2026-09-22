#!/usr/bin/env Rscript
# =============================================================================
# plot_paf_fmt.r -- 由 *.paf.fmt 绘制染色体级共线性图（点阵 + 覆盖）
#
# 重构自 example/plot_PUB.r，去除 pafr 与 ggpubr 依赖。
# 设计依据：REFACTOR_DESIGN.md   输入格式：PAF_FMT_SPEC.md   审计：AUDIT_REPORT.md
#
# 产出（视 --plots 而定）：
#   <prefix>.allChr.pdf        全基因组点阵
#   <prefix>.sepChr.<target>.pdf  逐 target 染色体点阵，每个 target 一个文件
#                                 （--sep-onefile on 合成单个多页 PDF）
#   <prefix>.coverage.pdf      target 覆盖条带
#   <prefix>.log.txt           运行日志 + QC 摘要存档（--log off 关闭）
#
# 默认按 Nature 系列的版面尺寸出图（--preset pub），X/Y 锁定同一比例尺。
#
# 用法：Rscript plot_paf_fmt.r --input <in.paf.fmt> --prefix <out_prefix> [选项]
#       Rscript plot_paf_fmt.r --help
#
# 依赖：R >= 4.0，data.table，ggplot2（grid 随 R 分发）
# =============================================================================

VERSION <- "1.1.0"

# grid 随 R 分发，无需安装；data.table 与 ggplot2 需自行安装
local({
  need <- c("data.table", "ggplot2")
  miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) {
    message("ERROR: 缺少 R 包：", paste(miss, collapse = ", "),
            "\n  安装：Rscript -e 'install.packages(c(",
            paste0('"', miss, '"', collapse = ", "), "))'")
    quit(save = "no", status = 1L)
  }
})
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
})

# ggplot2 的 text size / linewidth 以 mm 计，.pt = 72.27/25.4 用于与 pt 互换
PT <- .pt
pt_to_size <- function(pt) pt / PT   # 字号 pt -> ggplot size
pt_to_lw   <- function(pt) pt / PT   # 线宽 pt -> ggplot linewidth

# =============================================================================
# [1] CLI
# =============================================================================

# 参数表：name, default, type, help
#   type: chr / num / int / lgl
ARG_SPEC <- list(
  list("input",            NA,          "chr", "【必需】输入 .paf.fmt 路径"),
  list("prefix",           NA,          "chr", "【必需】输出文件前缀"),

  list("preset",           "pub",       "chr", "尺寸预设：pub（出版，默认）| screen（屏幕探索）"),

  list("keep-tp",          "P",         "chr", "保留的 PrimaryTag 值，逗号分隔。P=仅主比对（等价 pafr::filter_secondary_alignments）"),
  list("min-alen",         50000,       "num", "比对块长下限（BlkLen >）"),
  list("min-mapq",         20,          "num", "比对质量下限（MapQ >）"),
  list("min-qlen",         100000,      "num", "query 序列长度下限（QLen >=）"),
  list("min-pair-nmatch",  300000,      "num", "(query,target) 配对累计 Match 下限"),
  list("min-ident",        0,           "num", "Ident 下限（%），0 = 不过滤"),

  list("plots",            "allChr",    "chr", "要出的图，逗号分隔：allChr,sepChr,coverage。默认只出 allChr"),
  list("chroms",           "",          "chr", "sepChr 要画的 target，逗号分隔。留空=数据中全部"),
  list("x-seqs",           "",          "chr", "自定义 X 轴(query)序列及顺序，逗号分隔"),
  list("y-seqs",           "",          "chr", "自定义 Y 轴(subject)序列及顺序，逗号分隔"),
  list("x-last-seqs",      "chr00",     "chr", "X 轴排序后挪到末尾的序列名，逗号分隔"),
  list("y-last-seqs",      "chr00",     "chr", "Y 轴排序后挪到末尾的序列名，逗号分隔"),

  list("xlab",             "Query",     "chr", "X 轴标题"),
  list("ylab",             "Subject",   "chr", "Y 轴标题"),
  list("show-label",       TRUE,        "lgl", "是否在坐标轴上标序列名"),
  list("label-size",       NA,          "num", "序列名字号(pt)。留空=按 preset"),
  list("axis-title-size",  NA,          "num", "轴标题字号(pt)。留空=按 preset"),
  list("line-size",        NA,          "num", "点阵线宽(pt)。留空=按 preset"),
  list("x-tick-angle",     45,          "num", "X 轴刻度标签旋转角度(度)"),
  list("x-label-strip",    "",          "chr", "从 X 轴刻度标签中去掉的字面前缀"),
  list("y-label-strip",    "",          "chr", "从 Y 轴刻度标签中去掉的字面前缀"),
  list("show-axis-coord",  FALSE,       "lgl", "是否显示每条序列的局部坐标刻度"),
  list("coord-unit",       "Mb",        "chr", "坐标单位：bp | kb | Mb"),
  list("coord-step",       "10Mb",      "chr", "局部坐标刻度间隔，如 10Mb / 500kb；不带单位时按 --coord-unit"),

  list("col-forward",      "#1F4E79",   "chr", "正向比对(+)颜色"),
  list("col-reverse",      "#C00000",   "chr", "反向比对(-)颜色"),
  list("cov-fill",         "#2E7D32",   "chr", "覆盖图中已覆盖区块的填充色"),

  list("units",            "mm",        "chr", "画布尺寸单位：mm | cm | in"),
  list("dot-width",        NA,          "num", "allChr 画布宽。留空=按 preset"),
  list("dot-height",       NA,          "num", "allChr 画布高。留空=按 preset"),
  list("sep-width",        NA,          "num", "sepChr 每页宽。留空=按 preset"),
  list("sep-height",       NA,          "num", "sepChr 每页高。留空=按 preset"),
  list("cov-width",        NA,          "num", "coverage 画布宽。留空=按 preset"),
  list("cov-height",       NA,          "num", "coverage 画布高。留空=按 preset"),
  list("equal-scale",      "on",        "chr", "点阵图 X/Y 锁定同一比例尺：on（默认，长度可比）| off（面板铺满画布）"),
  list("sep-onefile",      "off",       "chr", "sepChr 合成单个多页 PDF：off（默认，每个 target 一个文件）| on"),
  list("fit-canvas",       "on",        "chr", "按字号反推画布尺寸：on（默认）| off。显式给出的 *-width/*-height 始终优先"),
  list("max-canvas",       NA,          "num", "反推画布的尺寸上限（--units 单位）。留空=按 preset"),
  list("min-panel",        NA,          "num", "反推时绘图面板的最小边长（--units 单位）。留空=按 preset"),

  list("font-family",      "Helvetica", "chr", "字体族。Helvetica 为 PDF base-14，与 Arial 度量一致"),
  list("font-file",        "",          "chr", "Arial 等字体的 .ttf 路径，给了就用 showtext 真正嵌入"),

  list("drop-invalid",     FALSE,       "lgl", "遇非法行丢弃并告警，而非报错"),
  list("quiet",            FALSE,       "lgl", "不在屏幕打印 QC 摘要（日志文件里照写）"),
  list("log",              "on",        "chr", "运行日志与 QC 摘要存档到 <prefix>.log.txt：on（默认）| off"),
  list("version",          FALSE,       "lgl", "打印版本号后退出"),
  list("help",             FALSE,       "lgl", "打印本帮助")
)

# preset 决定的默认值（单项参数显式给出时覆盖之）
PRESETS <- list(
  pub = list(`label-size` = 6, `axis-title-size` = 8, `line-size` = 0.3,
             `dot-width` = 183, `dot-height` = 183,
             `sep-width` = 89,  `sep-height` = 89,
             `cov-width` = 183, `cov-height` = 100, units = "mm",
             `max-canvas` = 183, `min-panel` = 80),
  screen = list(`label-size` = 10, `axis-title-size` = 16, `line-size` = 0.6,
                `dot-width` = 21, `dot-height` = 21,
                `sep-width` = 10, `sep-height` = 10,
                `cov-width` = 14, `cov-height` = 7, units = "in",
                `max-canvas` = 20, `min-panel` = 6)
)

# 日志：设备打开前的消息先缓存，log_open() 时补写；之后逐行写入并 flush，
# 中途报错退出也能在日志里看到停在哪一步
LOG <- new.env()
LOG$con <- NULL
LOG$buf <- character(0)

log_write <- function(...) {
  txt <- paste0(...)
  if (is.null(LOG$con)) {
    LOG$buf <- c(LOG$buf, txt)
  } else {
    writeLines(txt, LOG$con)
    flush(LOG$con)
  }
}

log_open <- function(path, header) {
  LOG$con <- file(path, open = "wt", encoding = "UTF-8")
  writeLines(c(header, LOG$buf), LOG$con)
  flush(LOG$con)
  LOG$buf <- character(0)
}

log_close <- function() {
  if (!is.null(LOG$con)) { close(LOG$con); LOG$con <- NULL }
}

stop_cli <- function(...) {
  msg <- paste0("ERROR: ", ...)
  message(msg)
  log_write(msg)
  log_close()
  quit(save = "no", status = 1L)
}

warn_cli <- function(...) {
  msg <- paste0("WARNING: ", ...)
  message(msg)
  log_write(msg)
}

print_help <- function() {
  cat("\nplot_paf_fmt.r -- 由 *.paf.fmt 绘制染色体级共线性图\n\n")
  cat("用法:\n  Rscript plot_paf_fmt.r --input <file.paf.fmt> --prefix <out_prefix> [options]\n\n")
  cat("参数:\n")
  for (a in ARG_SPEC) {
    dflt <- a[[2]]
    dtxt <- if (is.na(dflt[1])) "(按 preset/必需)" else as.character(dflt)
    cat(sprintf("  --%-18s %-14s %s\n", a[[1]], dtxt, a[[4]]))
  }
  cat("\npreset 取值:\n")
  cat("  pub    183x183mm(allChr) 89x89mm(sepChr) 183x100mm(coverage), 标签6pt 标题8pt 线宽0.3pt\n")
  cat("  screen 21x21in           10x10in         14x7in,               标签10pt 标题16pt 线宽0.6pt\n")
  cat("\n示例:\n")
  cat("  # 默认只出 allChr\n")
  cat("  Rscript plot_paf_fmt.r --input a.paf.fmt --prefix out\n\n")
  cat("  # 三图全出\n")
  cat("  Rscript plot_paf_fmt.r --input a.paf.fmt --prefix out --plots allChr,sepChr,coverage\n\n")
  cat("  # 自定义染色体组合（X/Y 各自有序，两组之间无配对关系）\n")
  cat("  Rscript plot_paf_fmt.r --input a.paf.fmt --prefix out.custom \\\n")
  cat("      --x-seqs LA1593Chr01,LA1593Chr03,LA1593Chr02 --y-seqs chr02,chr04\n\n")
  cat("  # 屏幕探索大图 + 显示局部坐标\n")
  cat("  Rscript plot_paf_fmt.r --input a.paf.fmt --prefix out.big --preset screen --show-axis-coord true\n\n")
  quit(save = "no", status = 0L)
}

as_lgl <- function(x, nm) {
  if (is.logical(x)) return(x)
  v <- tolower(as.character(x))
  if (v %in% c("true", "t", "yes", "y", "1"))  return(TRUE)
  if (v %in% c("false", "f", "no", "n", "0")) return(FALSE)
  stop_cli("--", nm, " 只接受 true/false，得到 '", x, "'")
}

as_num <- function(x, nm) {
  v <- suppressWarnings(as.numeric(x))
  if (is.na(v)) stop_cli("--", nm, " 需要数值，得到 '", x, "'")
  v
}

split_csv <- function(x) {
  if (is.null(x) || is.na(x[1]) || !nzchar(x[1])) return(character(0))
  trimws(strsplit(as.character(x), ",", fixed = TRUE)[[1]])
}

parse_args <- function(argv = commandArgs(trailingOnly = TRUE)) {
  spec_names <- vapply(ARG_SPEC, `[[`, "", 1L)
  types <- setNames(vapply(ARG_SPEC, `[[`, "", 3L), spec_names)
  opt <- setNames(lapply(ARG_SPEC, `[[`, 2L), spec_names)
  given <- character(0)

  i <- 1L
  while (i <= length(argv)) {
    tok <- argv[i]
    if (!startsWith(tok, "--")) stop_cli("无法识别的参数 '", tok, "'（参数须以 -- 开头）")
    key <- sub("^--", "", tok); val <- NULL
    if (grepl("=", key, fixed = TRUE)) {
      kv <- strsplit(key, "=", fixed = TRUE)[[1]]
      key <- kv[1]; val <- paste(kv[-1], collapse = "=")
    }
    if (!key %in% spec_names)
      stop_cli("未知参数 --", key, "。用 --help 查看全部参数")
    if (is.null(val)) {
      if (types[[key]] == "lgl" &&
          (i == length(argv) || startsWith(argv[i + 1L], "--"))) {
        val <- "true"                      # 裸开关
      } else {
        if (i == length(argv)) stop_cli("--", key, " 缺少取值")
        i <- i + 1L; val <- argv[i]
      }
    }
    opt[[key]] <- switch(types[[key]],
                         lgl = as_lgl(val, key),
                         num = as_num(val, key),
                         as.character(val))
    given <- c(given, key)
    i <- i + 1L
  }
  attr(opt, "given") <- given
  opt
}

apply_preset <- function(opt) {
  if (!opt$preset %in% names(PRESETS))
    stop_cli("--preset 只接受 ", paste(names(PRESETS), collapse = "/"), "，得到 '", opt$preset, "'")
  given <- attr(opt, "given")
  # --units 在下面就要用来换算，必须先校验；否则非法值会让 switch() 返回 NULL，
  # 尺寸参数变成 NULL，到 validate_args 里报成 "missing value where TRUE/FALSE needed"
  if (!opt$units %in% c("mm", "cm", "in"))
    stop_cli("--units 只接受 mm/cm/in，得到 '", opt$units, "'")
  pre <- PRESETS[[opt$preset]]
  # preset 里的尺寸是按 preset 自己的单位写的；用户显式换了 --units 就得换算，
  # 否则 "--preset pub --units in" 会把 183 mm 当成 183 in
  if ("units" %in% given && !identical(opt$units, pre$units)) {
    for (k in c("dot-width", "dot-height", "sep-width", "sep-height",
                "cov-width", "cov-height", "max-canvas", "min-panel"))
      pre[[k]] <- from_mm(to_mm(pre[[k]], pre$units), opt$units)
  }
  for (k in names(pre)) if (!(k %in% given)) opt[[k]] <- pre[[k]]
  attr(opt, "given") <- given
  opt
}

validate_args <- function(opt) {
  if (is.na(opt$input[1]))  stop_cli("缺少 --input。用 --help 查看用法")
  if (is.na(opt$prefix[1])) stop_cli("缺少 --prefix。用 --help 查看用法")
  if (!file.exists(opt$input)) stop_cli("输入文件不存在：", opt$input)
  if (!nzchar(opt$prefix))  stop_cli("--prefix 不能为空字符串")

  # 输出目录必须已存在且可写——否则要等到第一次 pdf() 才炸，报的是 R 的
  # "cannot open file" 回溯，看不出是路径问题
  odir <- dirname(opt$prefix)
  if (!dir.exists(odir))
    stop_cli("--prefix 的目录不存在：", odir, "\n  先建目录，或换一个 --prefix")
  if (file.access(odir, 2L) != 0L)
    stop_cli("--prefix 的目录不可写：", odir)

  for (k in c("min-alen", "min-mapq", "min-qlen", "min-pair-nmatch", "min-ident",
              "label-size", "axis-title-size", "line-size",
              "dot-width", "dot-height", "sep-width", "sep-height",
              "cov-width", "cov-height", "max-canvas", "min-panel"))
    if (opt[[k]] < 0) stop_cli("--", k, " 不能为负，得到 ", opt[[k]])

  if (opt$`min-mapq` > 255) stop_cli("--min-mapq 应在 0-255，得到 ", opt$`min-mapq`)
  if (opt$`min-ident` > 100) stop_cli("--min-ident 应在 0-100，得到 ", opt$`min-ident`)
  if (!opt$`coord-unit` %in% c("bp", "kb", "Mb"))
    stop_cli("--coord-unit 只接受 bp/kb/Mb，得到 '", opt$`coord-unit`, "'")
  opt$.coord_step_bp <- parse_coord_step(opt$`coord-step`, opt$`coord-unit`)
  if (!opt$units %in% c("mm", "cm", "in"))
    stop_cli("--units 只接受 mm/cm/in，得到 '", opt$units, "'")
  if (!opt$`fit-canvas` %in% c("on", "off"))
    stop_cli("--fit-canvas 只接受 on/off，得到 '", opt$`fit-canvas`, "'")
  if (!opt$`equal-scale` %in% c("on", "off"))
    stop_cli("--equal-scale 只接受 on/off，得到 '", opt$`equal-scale`, "'")
  if (!opt$`sep-onefile` %in% c("on", "off"))
    stop_cli("--sep-onefile 只接受 on/off，得到 '", opt$`sep-onefile`, "'")
  if (!opt$log %in% c("on", "off"))
    stop_cli("--log 只接受 on/off，得到 '", opt$log, "'")
  # 不锁定等比时各页形状一致，没有拆分的理由
  if (identical(opt$`equal-scale`, "off")) opt$`sep-onefile` <- "on"
  if (opt$`min-panel` > opt$`max-canvas`)
    stop_cli("--min-panel (", opt$`min-panel`, ") 不能大于 --max-canvas (",
             opt$`max-canvas`, ")。调小 --min-panel 或调大 --max-canvas")
  if (opt$`label-size` <= 0) stop_cli("--label-size 必须为正，得到 ", opt$`label-size`)

  tp <- split_csv(opt$`keep-tp`)
  if (length(tp) && !all(tp %in% c("P", "S", "I", "i")))
    stop_cli("--keep-tp 只接受 P/S/I/i 的组合，得到 '", opt$`keep-tp`, "'")

  custom <- length(split_csv(opt$`x-seqs`)) > 0 || length(split_csv(opt$`y-seqs`)) > 0
  if (!nzchar(opt$plots)) opt$plots <- "allChr"
  pl <- split_csv(opt$plots)
  bad <- setdiff(pl, c("allChr", "sepChr", "coverage"))
  if (length(bad)) stop_cli("--plots 含未知取值：", paste(bad, collapse = ","),
                            "（合法：allChr,sepChr,coverage）")
  opt$.plots <- pl
  opt$.custom <- custom

  if (nzchar(opt$`font-file`) && !file.exists(opt$`font-file`))
    stop_cli("--font-file 指定的文件不存在：", opt$`font-file`)
  opt
}

# 画布尺寸统一换算成英寸（pdf() 要求）
to_inch <- function(v, units) switch(units, mm = v / 25.4, cm = v / 2.54, `in` = v)
to_mm   <- function(v, units) switch(units, mm = v, cm = v * 10, `in` = v * 25.4)
from_mm <- function(v, units) switch(units, mm = v, cm = v / 10, `in` = v / 25.4)

# =============================================================================
# [2] IO -- 读取与校验 .paf.fmt
# =============================================================================

REQUIRED_COLS <- c("QID", "QLen", "QS", "QE", "Str", "SID", "SLen", "SS", "SE",
                   "Match", "BlkLen", "MapQ")
NUMERIC_COLS  <- c("QLen", "QS", "QE", "SLen", "SS", "SE", "Match", "BlkLen", "MapQ")
IDENT_COLS    <- c("QIdent", "SIdent", "Ident")

#' 补算三个一致性列（输入缺列时用）。定义见 PAF_FMT_SPEC.md §3
derive_ident <- function(dt) {
  if (!"QIdent" %in% names(dt)) dt[, QIdent := round(Match / (QE - QS) * 100, 1)]
  if (!"SIdent" %in% names(dt)) dt[, SIdent := round(Match / (SE - SS) * 100, 1)]
  if (!"Ident"  %in% names(dt)) dt[, Ident  := round(Match / BlkLen   * 100, 1)]
  dt[]
}

#' 读 .paf.fmt 并做全部输入校验
read_paf_fmt <- function(path, drop_invalid = FALSE) {
  # 序列名必须按字符读：纯数字名（如 SID 为 1,2,...）被猜成 integer 后，
  # 下游 offset[SID] 这类具名向量查找会退化成按位置取值，线条错位或变 NA
  hdr <- tryCatch(names(fread(path, sep = "\t", header = TRUE, nrows = 0L)),
                  error = function(e) stop_cli("读取失败：", conditionMessage(e)))
  id_cols <- intersect(c("QID", "SID"), hdr)
  dt <- tryCatch(
    fread(path, sep = "\t", header = TRUE, showProgress = FALSE,
          na.strings = c("", "NA"),
          colClasses = if (length(id_cols)) list(character = id_cols)),
    error = function(e) stop_cli("读取失败：", conditionMessage(e))
  )
  if (nrow(dt) == 0L) stop_cli("输入文件没有数据行：", path)

  miss <- setdiff(REQUIRED_COLS, names(dt))
  if (length(miss))
    stop_cli("输入缺少必需列：", paste(miss, collapse = ", "),
             "\n  实际列：", paste(names(dt), collapse = ", "),
             "\n  期望列见 PAF_FMT_SPEC.md")

  # 数值列可转性
  for (cl in NUMERIC_COLS) {
    if (!is.numeric(dt[[cl]])) {
      v <- suppressWarnings(as.numeric(dt[[cl]]))
      bad <- which(is.na(v) & !is.na(dt[[cl]]))
      if (length(bad))
        stop_cli("列 ", cl, " 含非数值，首个位于数据行 ", bad[1],
                 "（值 '", dt[[cl]][bad[1]], "'）")
      set(dt, j = cl, value = v)
    }
  }

  # 取值域
  bad_str <- which(!dt$Str %in% c("+", "-"))
  bad_q   <- which(dt$QS >= dt$QE)
  bad_s   <- which(dt$SS >= dt$SE)
  bad_all <- sort(unique(c(bad_str, bad_q, bad_s)))
  if (length(bad_all)) {
    first <- bad_all[1]
    reason <- if (first %in% bad_str) "Str 非 +/-" else if (first %in% bad_q) "QS >= QE" else "SS >= SE"
    if (!drop_invalid)
      stop_cli(length(bad_all), " 行非法（首个位于数据行 ", first, "：", reason,
               "）。加 --drop-invalid 可丢弃这些行继续")
    warn_cli("丢弃 ", length(bad_all), " 行非法记录（首个位于数据行 ", first, "：", reason, "）")
    dt <- dt[-bad_all]
    if (nrow(dt) == 0L) stop_cli("丢弃非法行后没有数据了")
  }

  if (!"PrimaryTag" %in% names(dt))
    warn_cli("输入无 PrimaryTag 列（旧的 15 列格式）：跳过 secondary 过滤，",
             "此时次优比对只能靠 --min-mapq 排除")

  derive_ident(dt)

  # 同名序列长度一致性（data.table 的 by= 需要显式 eval，见 ?data.table）
  for (p in list(c("QID", "QLen"), c("SID", "SLen"))) {
    idc <- p[1]; lenc <- p[2]
    n <- dt[, .(nu = uniqueN(.SD[[1]])), by = eval(idc), .SDcols = lenc][nu > 1L]
    if (nrow(n)) warn_cli(idc, " 中有 ", nrow(n), " 个序列的 ", lenc,
                          " 不唯一，取首次出现值：", paste(n[[1]], collapse = ", "))
  }
  dt[]
}

# =============================================================================
# [3] FILTER -- 复现旧脚本的过滤链（AUDIT_REPORT.md §2）
# =============================================================================

#' 记录级过滤。返回过滤后的表，attr("steps") 记录逐级计数供 QC
filter_alignments <- function(dt, opt) {
  steps <- c(input = nrow(dt))

  keep_tp <- split_csv(opt$`keep-tp`)
  if ("PrimaryTag" %in% names(dt) && length(keep_tp)) {
    dt <- dt[PrimaryTag %in% keep_tp]
    steps <- c(steps, setNames(nrow(dt), paste0("keep-tp(", opt$`keep-tp`, ")")))
  }
  dt <- dt[BlkLen > opt$`min-alen` & MapQ > opt$`min-mapq`]
  steps <- c(steps, `min-alen & min-mapq` = nrow(dt))

  dt <- dt[QLen >= opt$`min-qlen`]
  steps <- c(steps, `min-qlen` = nrow(dt))

  if (opt$`min-ident` > 0) {
    dt <- dt[Ident >= opt$`min-ident`]
    steps <- c(steps, `min-ident` = nrow(dt))
  }
  if (nrow(dt) == 0L)
    stop_cli("过滤后没有记录了。当前阈值：--min-alen ", opt$`min-alen`,
             " --min-mapq ", opt$`min-mapq`, " --min-qlen ", opt$`min-qlen`,
             " --min-ident ", opt$`min-ident`, " --keep-tp ", opt$`keep-tp`,
             "\n  放宽阈值后重试")
  setattr(dt, "steps", steps)
  dt[]
}

#' (query,target) 配对的累计 Match
pair_summary <- function(dt) dt[, .(sum_nmatch = sum(Match)), by = .(QID, SID)]

#' 取通过配对阈值的 (QID,SID)
select_pairs <- function(pairs, min_nmatch) pairs[sum_nmatch >= min_nmatch]

# =============================================================================
# [4] LAYOUT -- 序列长度、排序、累积偏移
# =============================================================================

#' 从数据取每条序列的长度（同名取首次出现值）
seq_lengths <- function(dt, id_col, len_col) {
  x <- dt[, .(len = as.numeric(get(len_col))[1]), by = c(id_col)]
  setnames(x, id_col, "id")
  x
}

#' 排序：字典序（名字全为纯数字时按数值），再把 last_seqs 里的名字挪到末尾
order_seqs <- function(ids, last_seqs = character(0)) {
  ids <- unique(as.character(ids))
  ids <- if (length(ids) && all(grepl("^[0-9]+$", ids))) ids[order(as.numeric(ids))] else sort(ids)
  if (!length(last_seqs)) return(ids)
  hit <- ids[ids %in% last_seqs]
  c(ids[!ids %in% last_seqs], hit)
}

#' 累积布局：给定顺序与长度，算 offset / 区间中点 / 边界
cumulative_layout <- function(len_dt, ord) {
  x <- len_dt[match(ord, len_dt$id)]
  x <- x[!is.na(x$id)]
  if (nrow(x) == 0L) return(NULL)
  x[, offset := cumsum(c(0, head(len, -1)))]
  x[, mid := offset + len / 2]
  x[, end := offset + len]
  x[]
}

#' 刻度标签：去掉字面前缀
strip_prefix <- function(labels, prefix) {
  if (!nzchar(prefix)) return(labels)
  hit <- startsWith(labels, prefix)            # 字面匹配，不过正则，省去转义
  labels[hit] <- substring(labels[hit], nchar(prefix) + 1L)
  labels
}

#' 坐标单位换算
unit_div <- function(u) switch(u, bp = 1, kb = 1e3, Mb = 1e6)

#' 解析 --coord-step（"10Mb" / "500kb" / "2000000bp" / 不带单位按 coord_unit），返回 bp
parse_coord_step <- function(x, coord_unit) {
  m <- regmatches(x, regexec("^\\s*([0-9]*\\.?[0-9]+(?:[eE][+-]?[0-9]+)?)\\s*(bp|kb|mb)?\\s*$",
                             x, ignore.case = TRUE, perl = TRUE))[[1]]
  if (!length(m)) stop_cli("--coord-step 格式应为 数值[bp|kb|Mb]，如 10Mb、500kb，得到 '", x, "'")
  u <- if (nzchar(m[3])) c(bp = "bp", kb = "kb", mb = "Mb")[[tolower(m[3])]] else coord_unit
  v <- as.numeric(m[2]) * unit_div(u)
  if (!is.finite(v) || v <= 0) stop_cli("--coord-step 必须为正，得到 '", x, "'")
  v
}

# =============================================================================
# [4b] FIT -- 按字号反推画布尺寸（REFACTOR_DESIGN.md §8.4）
# =============================================================================
# 刻度标签的物理尺寸只由字号(pt)决定，跟画布多大无关。所以"画布定死 183 mm、字缩
# 小去迁就"是反的：画布一缩，序列区块跟着变窄，标签反而更容易撞在一起。正确做法
# 是倒过来——先按可读性定字号，再由"最挤的那对相邻标签必须排得开"反推面板至少要
# 多宽，加上页边距得到画布。字号调大画布跟着长，标签短画布跟着缩，这样任何情况下
# 字号都不必为了塞进固定画布而被压到看不清。

MM_PER_PT <- 25.4 / 72.27

#' 用真实字体度量量一组文字的宽度(mm)。在临时设备上测，不影响输出设备
text_width_mm <- function(txt, pt, family) {
  if (!length(txt)) return(numeric(0))
  tf <- tempfile(fileext = ".pdf")
  ok <- tryCatch({ grDevices::pdf(tf, width = 8, height = 8, family = family); TRUE },
                 error = function(e) FALSE)
  # 量不了就按 0.55 em 估，宁可略宽也不要算窄
  if (!ok) return(nchar(txt) * pt * MM_PER_PT * 0.55)
  on.exit({ grDevices::dev.off(); unlink(tf) }, add = TRUE)
  graphics::par(ps = pt, family = family)
  graphics::strwidth(txt, units = "inches") * 25.4
}

#' 标签沿"轴的走向"占多长(mm)——决定两个相邻刻度会不会撞上。
#' 注意这跟标签的横向宽度不是一回事：X 轴上水平标签沿轴占的是字宽，而旋转之后，
#' 一堆平行的标签只要锚点间距在垂直于文字方向上的投影 >= 行高就不会重叠，与标签
#' 有多长无关（这正是旋转刻度能省地方的原因）；Y 轴上水平标签沿轴占的则是行高。
label_extent <- function(labs, ang, pt, family, axis = c("x", "y")) {
  axis <- match.arg(axis)
  if (!length(labs)) return(0)
  hline <- pt * MM_PER_PT * 1.2
  a <- (ang %% 180) * pi / 180
  w <- text_width_mm(labs, pt, family)
  if (axis == "x") {
    if (abs(sin(a)) < 1e-6) w else rep(hline / abs(sin(a)), length(w))
  } else {
    if (abs(cos(a)) < 1e-6) w else rep(hline / abs(cos(a)), length(w))
  }
}

#' 一组刻度要互不相撞，面板在该轴方向至少需要多长(mm)
#' ext 为各标签沿轴的占位长度，pos 为锚点的数据坐标，total 为该轴数据全长
axis_span_need <- function(ext, pos, total, tol = 0.004) {
  if (length(pos) < 2L || !is.finite(total) || total <= 0) return(0)
  o <- order(pos); ext <- rep_len(ext, length(pos))[o]; pos <- pos[o]
  # 间距小于 tol 的刻度画布再大也会叠在一起（show-axis-coord 下，序列交界处上一条
  # 的末刻度和下一条的 0 刻度几乎同位），先并成一簇再算间距，免得需求发散
  grp <- cumsum(c(TRUE, diff(pos) / total > tol))
  pos <- as.numeric(tapply(pos, grp, mean))
  ext <- as.numeric(tapply(ext, grp, max))
  if (length(pos) < 2L) return(0)
  gapfrac <- pmax(diff(pos) / total, 1e-9)
  need <- (head(ext, -1) + tail(ext, -1)) / 2 * 1.08   # 相邻两个各伸出一半
  max(need / gapfrac)
}

#' 旋转后的刻度标签在垂直于轴方向上吃掉的深度(mm)
tick_depth <- function(labs, ang, pt, family) {
  if (!length(labs)) return(0)
  hline <- pt * MM_PER_PT * 1.2
  a <- (ang %% 180) * pi / 180
  max(text_width_mm(labs, pt, family)) * abs(sin(a)) + hline * abs(cos(a))
}

#' show-axis-coord 模式下每条序列的局部坐标刻度（axis_breaks() 与画布反推共用）
num_ticks <- function(lay, opt) {
  div  <- unit_div(opt$`coord-unit`)
  step <- opt$.coord_step_bp
  if (max(lay$len) / step > 1000)
    stop_cli("--coord-step ", opt$`coord-step`, " 太小：最长序列上会有 ",
             floor(max(lay$len) / step), " 个刻度。调大 --coord-step")
  # 固定间隔、全部保留：序列末端的刻度可能与下一条序列的 0 刻度挨得很近，接受轻微重叠
  pos <- numeric(0); lab <- character(0); sq <- integer(0)
  for (i in seq_len(nrow(lay))) {
    loc <- seq(0, lay$len[i], by = step)
    pos <- c(pos, lay$offset[i] + loc)
    lab <- c(lab, format(loc / div, trim = TRUE))
    sq  <- c(sq, rep(i, length(loc)))
  }
  list(pos = pos, lab = lab, seq = sq)
}

#' 局部坐标刻度的面板需求：只看同一条序列内相邻刻度的间距。
#' 序列交界处（上一条末刻度与下一条 0 刻度）的重叠已接受，不为它放大画布
coord_ticks_need <- function(nt, ext, total) {
  ext <- rep_len(ext, length(nt$pos))
  max(0, vapply(split(seq_along(nt$pos), nt$seq),
                function(ix) axis_span_need(ext[ix], nt$pos[ix], total), 0))
}

#' 面板要多大才不让刻度标签撞上：返回 c(w, h)，单位 mm
panel_need_dot <- function(qlay, slay, opt) {
  pt  <- opt$`label-size`; fam <- opt$`font-family`
  qt  <- max(qlay$end); st <- max(slay$end)
  ang <- opt$`x-tick-angle`
  xn  <- strip_prefix(qlay$id, opt$`x-label-strip`)
  yn  <- strip_prefix(slay$id, opt$`y-label-strip`)
  if (!opt$`show-label`) return(c(0, 0))

  if (opt$`show-axis-coord`) {
    # 主轴是局部坐标数字，序列名挪到次轴（不旋转）
    nx <- num_ticks(qlay, opt); ny <- num_ticks(slay, opt)
    c(max(coord_ticks_need(nx, label_extent(nx$lab, ang, pt, fam, "x"), qt),
          axis_span_need(label_extent(xn, 0, pt, fam, "x"), qlay$mid, qt)),
      max(coord_ticks_need(ny, label_extent(ny$lab, 0, pt, fam, "y"), st),
          axis_span_need(label_extent(yn, 0, pt, fam, "y"), slay$mid, st)))
  } else {
    c(axis_span_need(label_extent(xn, ang, pt, fam, "x"), qlay$mid, qt),
      axis_span_need(label_extent(yn, 0, pt, fam, "y"), slay$mid, st))
  }
}

#' 量 ggplot 自己的布局：面板之外的固定边距吃掉多少(mm)，返回 c(w, h)
#' 从 gtable 里按 layout 索引定位面板所在的行列，其余行列之和即边距。
#' 不能改用"挑出含 null 的行列丢掉"——axis-l / axis-b 这些格子的单位是
#' sum(定长, 1null, ...) 的混合体，整格丢掉会少算好几毫米。
#' 量之前先摘掉 coord_fixed：等比约束会让面板在画布里"信箱式"留白，
#' 那样量到的就不是边距了。边距只由字号与标签文本决定，与画布大小无关。
fixed_margins_mm <- function(p, opt) {
  fallback <- c(14, 20)
  p$coordinates <- ggplot2::coord_cartesian()
  tf <- tempfile(fileext = ".pdf")
  ok <- tryCatch({ grDevices::pdf(tf, width = 8, height = 8,
                                  family = opt$`font-family`); TRUE },
                 error = function(e) FALSE)
  if (!ok) return(fallback)
  on.exit({ grDevices::dev.off(); unlink(tf) }, add = TRUE)
  res <- tryCatch({
    g <- ggplot2::ggplotGrob(p)
    lay <- g$layout[g$layout$name == "panel", ][1, ]
    if (is.na(lay$l)) stop("no panel cell")
    cols <- setdiff(seq_along(g$widths),  lay$l:lay$r)
    rows <- setdiff(seq_along(g$heights), lay$t:lay$b)
    c(grid::convertWidth (sum(g$widths[cols]),  "mm", valueOnly = TRUE),
      grid::convertHeight(sum(g$heights[rows]), "mm", valueOnly = TRUE))
  }, error = function(e) fallback)
  if (!all(is.finite(res)) || any(res < 0)) return(fallback)
  res
}

#' 点阵图：反推画布 c(width, height)，单位 mm
#' 面板尺寸由"标签排得开"和"--equal-scale 要求的数据纵横比"共同决定，
#' 边距取自 ggplot 自己的布局，两者相加即画布。超过 --max-canvas 时等比缩面板。
fit_dot_canvas <- function(qlay, slay, opt, plot_obj) {
  need <- panel_need_dot(qlay, slay, opt)
  minp <- to_mm(opt$`min-panel`, opt$units)
  eq   <- identical(opt$`equal-scale`, "on")

  if (eq) {
    a  <- max(slay$end) / max(qlay$end)            # 数据纵横比 Y/X
    # 面板必须同时满足两轴的标签需求，且 --min-panel 约束较长的一边
    pw <- max(need[[1]], need[[2]] / a, minp / max(1, a))
    ph <- a * pw
  } else {
    pw <- ph <- max(need[[1]], need[[2]], minp)    # 不等比时面板取正方形
  }

  fx <- fixed_margins_mm(plot_obj, opt)
  hi <- to_mm(opt$`max-canvas`, opt$units)
  k  <- min(1, (hi - fx[1]) / pw, (hi - fx[2]) / ph)
  if (k < 1) {
    warn_cli(attr(plot_obj, "what"), "：按 --label-size ", opt$`label-size`, " pt 需要 ",
             sprintf("%.0f x %.0f mm", pw + fx[1], ph + fx[2]),
             "，已按 --max-canvas ", sprintf("%.0f mm", hi), " 等比缩到 ",
             sprintf("%.0f x %.0f mm", pw * k + fx[1], ph * k + fx[2]),
             "，刻度标签可能变挤。可选：调大 --max-canvas、调小 --label-size、",
             "用 --x-label-strip 缩短标签，或加大 --x-tick-angle")
    pw <- pw * k; ph <- ph * k
  }
  c(width = pw + fx[1], height = ph + fx[2])
}

#' 覆盖图：由字号反推画布 c(width, height)，单位 mm。高度由序列条数决定
fit_cov_canvas <- function(lay, opt, plot_obj) {
  hl <- opt$`label-size` * MM_PER_PT * 1.2
  panel_h <- nrow(lay) * hl * 1.9                    # 每条序列一行，行距 1.9 倍行高
  fx <- fixed_margins_mm(plot_obj, opt)
  # 覆盖图横向没有标签碰撞约束（且不受 --equal-scale 影响：两轴量纲本来就不同），
  # 按条带图惯用的扁长比例给宽度，并就地卡在可用宽度内——这不算"挤"，不告警
  room    <- to_mm(opt$`max-canvas`, opt$units) - fx[1]
  panel_w <- min(max(panel_h * 2.5, to_mm(opt$`min-panel`, opt$units)), room)
  c(width = panel_w + fx[1], height = panel_h + fx[2])
}

# =============================================================================
# [5] PLOT-DOT -- 点阵图
# =============================================================================

#' 基础主题：白底（去灰色背景）、无网格、指定字体
theme_paf <- function(opt) {
  theme_bw(base_family = opt$`font-family`) +
    theme(
      panel.background  = element_rect(fill = "white", colour = NA),
      plot.background   = element_rect(fill = "white", colour = NA),
      panel.grid        = element_blank(),
      panel.border      = element_rect(colour = "grey30", fill = NA,
                                       linewidth = pt_to_lw(0.4)),
      axis.title        = element_text(size = pt_to_size(opt$`axis-title-size`),
                                       face = "bold"),
      axis.text         = element_text(size = pt_to_size(opt$`label-size`),
                                       colour = "black"),
      axis.ticks        = element_line(colour = "grey30", linewidth = pt_to_lw(0.3)),
      legend.title      = element_text(size = pt_to_size(opt$`label-size`)),
      legend.text       = element_text(size = pt_to_size(opt$`label-size`)),
      legend.key        = element_blank(),
      legend.background = element_blank(),
      plot.title        = element_text(size = pt_to_size(opt$`axis-title-size`),
                                       face = "bold", hjust = 0.5)
    )
}

#' 把 alignment 转成线段坐标。负链 y 端点对调 -> 倒位表现为负斜率
build_dot_data <- function(dt, qlay, slay) {
  qo <- setNames(qlay$offset, qlay$id)
  so <- setNames(slay$offset, slay$id)
  d <- dt[QID %in% qlay$id & SID %in% slay$id]
  if (nrow(d) == 0L) return(NULL)
  d[, `:=`(x = QS + qo[QID], xend = QE + qo[QID])]
  d[, `:=`(y = fifelse(Str == "-", SE, SS) + so[SID],
           yend = fifelse(Str == "-", SS, SE) + so[SID])]
  d[, Strand := fifelse(Str == "-", "Reverse (-)", "Forward (+)")]
  d[]
}

#' 生成某一轴的刻度：show_coord=FALSE 时刻度即序列名；TRUE 时为局部坐标 + 次轴名
axis_breaks <- function(lay, opt, strip) {
  labs <- strip_prefix(lay$id, strip)
  if (!opt$`show-axis-coord`)
    return(list(breaks = lay$mid, labels = labs, sec_breaks = NULL))
  nt <- num_ticks(lay, opt)
  list(breaks = nt$pos, labels = nt$lab, sec_breaks = lay$mid, sec_labels = labs)
}

#' 点阵图主函数
plot_dotplot <- function(dt, qlay, slay, opt, title = NULL) {
  seg <- build_dot_data(dt, qlay, slay)
  xa <- axis_breaks(qlay, opt, opt$`x-label-strip`)
  ya <- axis_breaks(slay, opt, opt$`y-label-strip`)

  p <- ggplot()
  # 序列边界线
  p <- p +
    geom_vline(xintercept = c(qlay$offset, max(qlay$end)),
               linetype = "dotted", colour = "grey55", linewidth = pt_to_lw(0.3)) +
    geom_hline(yintercept = c(slay$offset, max(slay$end)),
               linetype = "dotted", colour = "grey55", linewidth = pt_to_lw(0.3))
  if (!is.null(seg))
    p <- p + geom_segment(data = seg,
                          aes(x = x, xend = xend, y = y, yend = yend, colour = Strand),
                          linewidth = pt_to_lw(opt$`line-size`), lineend = "butt")

  xsec <- if (!is.null(xa$sec_breaks))
    sec_axis(~ ., breaks = xa$sec_breaks, labels = xa$sec_labels) else waiver()
  ysec <- if (!is.null(ya$sec_breaks))
    sec_axis(~ ., breaks = ya$sec_breaks, labels = ya$sec_labels) else waiver()

  p <- p +
    scale_x_continuous(breaks = if (opt$`show-label`) xa$breaks else numeric(0),
                       labels = if (opt$`show-label`) xa$labels else waiver(),
                       limits = c(0, max(qlay$end)), expand = expansion(mult = 0.01),
                       sec.axis = xsec) +
    scale_y_continuous(breaks = if (opt$`show-label`) ya$breaks else numeric(0),
                       labels = if (opt$`show-label`) ya$labels else waiver(),
                       limits = c(0, max(slay$end)), expand = expansion(mult = 0.01),
                       sec.axis = ysec) +
    scale_colour_manual(values = c("Forward (+)" = opt$`col-forward`,
                                   "Reverse (-)" = opt$`col-reverse`),
                        drop = FALSE, name = NULL) +
    labs(x = opt$xlab, y = opt$ylab, title = title) +
    theme_paf(opt) +
    theme(axis.text.x = element_text(angle = opt$`x-tick-angle`,
                                     hjust = if (opt$`x-tick-angle` != 0) 1 else 0.5,
                                     vjust = if (opt$`x-tick-angle` != 0) 1 else 1),
          legend.position = "top",
          legend.margin = margin(0, 0, 0, 0),
          legend.box.spacing = unit(2, "pt"))

  # X/Y 锁定同一比例尺（bp/mm 相同），两轴上的长度、斜率才可以互相比较。
  # 旧版 pafr::dotplot() 也是等比的（实测旧 PDF 面板宽高比 = 数据跨度比，
  # allChr 1.055、sepChr 第 1 页 2.564，吻合到 0.1%）。
  if (identical(opt$`equal-scale`, "on")) p <- p + coord_fixed(ratio = 1)
  p
}

# =============================================================================
# [6] PLOT-COV -- 覆盖条带图
# =============================================================================

#' 求每条序列被覆盖区间的并集。返回 (id, start, end)
coverage_intervals <- function(dt, id_col = "SID", s_col = "SS", e_col = "SE") {
  x <- dt[, .(id = get(id_col), s = get(s_col), e = get(e_col))]
  setorder(x, id, s, e)
  x[, grp := cumsum(c(TRUE, s[-1] > cummax(e)[-.N])), by = id]
  x[, .(start = min(s), end = max(e)), by = .(id, grp)][, grp := NULL][]
}

#' 覆盖统计：覆盖碱基数、覆盖率、缺口数
coverage_stats <- function(cov, len_dt) {
  s <- cov[, .(covered = sum(end - start), n_block = .N), by = id]
  s <- merge(s, len_dt, by = "id", all.y = TRUE)
  s[is.na(covered), `:=`(covered = 0, n_block = 0L)]
  s[, `:=`(pct = covered / len * 100, n_gap = pmax(n_block - 1L, 0L))]
  s[]
}

#' 覆盖图：每条 target 一行，已覆盖填 --cov-fill，未覆盖留白，全长画细外框
plot_coverage_track <- function(cov, len_dt, ord, opt) {
  lay <- len_dt[match(ord, len_dt$id)]
  lay <- lay[!is.na(lay$id)]
  lay[, row := .I]
  h <- 0.38
  cv <- merge(cov, lay[, .(id, row)], by = "id")
  div <- unit_div(opt$`coord-unit`)

  ggplot() +
    geom_rect(data = lay, aes(xmin = 0, xmax = len, ymin = row - h, ymax = row + h),
              fill = NA, colour = "grey30", linewidth = pt_to_lw(0.3)) +
    geom_rect(data = cv, aes(xmin = start, xmax = end, ymin = row - h, ymax = row + h),
              fill = opt$`cov-fill`, colour = NA) +
    scale_x_continuous(labels = function(v) format(v / div, trim = TRUE),
                       expand = expansion(mult = c(0.005, 0.02))) +
    scale_y_reverse(breaks = lay$row,
                    labels = if (opt$`show-label`)
                      strip_prefix(lay$id, opt$`y-label-strip`) else NULL,
                    expand = expansion(add = 0.6)) +
    labs(x = paste0("Position in sequence (", opt$`coord-unit`, ")"), y = NULL) +
    theme_paf(opt) +
    theme(axis.ticks.y = element_blank())
}

# =============================================================================
# [7] MAIN
# =============================================================================

#' 注册并启用自定义字体（--font-file）
setup_font <- function(opt) {
  if (!nzchar(opt$`font-file`)) return(invisible(FALSE))
  if (!requireNamespace("showtext", quietly = TRUE) ||
      !requireNamespace("sysfonts", quietly = TRUE))
    stop_cli("--font-file 需要 showtext 与 sysfonts 包，但未安装")
  sysfonts::font_add(family = opt$`font-family`, regular = opt$`font-file`)
  showtext::showtext_auto()
  invisible(TRUE)
}

#' 决定某张图最终的画布尺寸（返回 --units 单位）。
#' 优先级：命令行显式给出的 *-width/*-height > 按字号反推 > preset 默认值
resolve_canvas <- function(fit_mm, wkey, hkey, opt) {
  given <- attr(opt, "given")
  w <- opt[[wkey]]; h <- opt[[hkey]]
  if (identical(opt$`fit-canvas`, "on") && !is.null(fit_mm)) {
    if (!(wkey %in% given)) w <- from_mm(fit_mm[[1]], opt$units)
    if (!(hkey %in% given)) h <- from_mm(fit_mm[[2]], opt$units)
  }
  c(w, h)
}

open_pdf <- function(path, w, h, opt) {
  grDevices::pdf(file = path,
                 width = to_inch(w, opt$units), height = to_inch(h, opt$units),
                 family = opt$`font-family`, onefile = TRUE, useDingbats = FALSE)
}

# 注：单个 PDF 内所有页的画布尺寸必须一致（pdf(onefile=TRUE) 的限制），
# 因此 sepChr 各页统一用 --sep-width/--sep-height。旧版"又扁又宽"的问题来自
# pafr 的等比坐标，本实现不锁定长宽比，面板会自动填满画布。

#' 日志文件头：时间、版本、完整命令行、全部生效参数（* = 命令行显式给出）
log_header <- function(opt) {
  ca <- commandArgs(trailingOnly = FALSE)
  script <- sub("^--file=", "", grep("^--file=", ca, value = TRUE)[1])
  argv <- commandArgs(trailingOnly = TRUE)
  q <- vapply(argv, function(a)
    if (grepl("^[A-Za-z0-9_./,=:+@%-]+$", a)) a else shQuote(a), "")
  given <- attr(opt, "given")
  spec <- Filter(function(a) !a[[1]] %in% c("version", "help"), ARG_SPEC)
  pars <- vapply(spec, function(a) {
    k <- a[[1]]
    sprintf("  %s --%-18s %s", if (k %in% given) "*" else " ", k,
            format(opt[[k]], scientific = FALSE, trim = TRUE))
  }, "")
  c("== plot_paf_fmt.r 运行日志 ==",
    sprintf("时间     : %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    sprintf("版本     : plot_paf_fmt.r %s  (%s)", VERSION, R.version.string),
    sprintf("工作目录 : %s", getwd()),
    sprintf("命令     : %s", paste(c("Rscript", script, q), collapse = " ")),
    "",
    "== 生效参数（* = 命令行给出，其余为默认/preset）==",
    pars,
    "",
    "== 运行消息 ==")
}

main <- function() {
  t0 <- proc.time()
  opt <- parse_args()
  if (isTRUE(opt$version)) {
    cat("plot_paf_fmt.r ", VERSION, "\n", sep = "")
    quit(save = "no", status = 0L)
  }
  if (isTRUE(opt$help)) print_help()
  opt <- apply_preset(opt)
  opt <- validate_args(opt)
  log_file <- NULL
  if (identical(opt$log, "on")) {
    log_file <- paste0(opt$prefix, ".log.txt")
    log_open(log_file, log_header(opt))
  }
  setup_font(opt)

  # ---- 读入与过滤 ----
  dt <- read_paf_fmt(opt$input, drop_invalid = opt$`drop-invalid`)
  n_raw <- nrow(dt)
  flt <- filter_alignments(dt, opt)
  steps <- attr(flt, "steps")

  pairs_all <- pair_summary(flt)
  pairs_keep <- select_pairs(pairs_all, opt$`min-pair-nmatch`)

  qlen <- seq_lengths(flt, "QID", "QLen")
  slen <- seq_lengths(flt, "SID", "SLen")

  x_seqs <- split_csv(opt$`x-seqs`); y_seqs <- split_csv(opt$`y-seqs`)
  if (length(x_seqs)) {
    miss <- setdiff(x_seqs, qlen$id)
    if (length(miss)) warn_cli("--x-seqs 中这些名字不在数据里，已跳过：",
                               paste(miss, collapse = ", "))
    x_seqs <- x_seqs[x_seqs %in% qlen$id]
    if (!length(x_seqs)) stop_cli("--x-seqs 指定的序列全都不在数据中")
  }
  if (length(y_seqs)) {
    miss <- setdiff(y_seqs, slen$id)
    if (length(miss)) warn_cli("--y-seqs 中这些名字不在数据里，已跳过：",
                               paste(miss, collapse = ", "))
    y_seqs <- y_seqs[y_seqs %in% slen$id]
    if (!length(y_seqs)) stop_cli("--y-seqs 指定的序列全都不在数据中")
  }

  outs <- character(0)
  sep_pages <- NA_integer_
  sep_files <- NA_integer_
  canvas_log <- list()

  # ---- allChr（自定义模式下即为自定义组合图）----
  if ("allChr" %in% opt$.plots) {
    if (opt$.custom) {
      # 显式点名：顺序即用户给的，且跳过配对阈值（REFACTOR_DESIGN.md §4.3.1）
      qord <- if (length(x_seqs)) x_seqs else order_seqs(qlen$id, split_csv(opt$`x-last-seqs`))
      sord <- if (length(y_seqs)) y_seqs else order_seqs(slen$id, split_csv(opt$`y-last-seqs`))
      sub <- flt[QID %in% qord & SID %in% sord]
    } else {
      sub <- merge(flt, pairs_keep[, .(QID, SID)], by = c("QID", "SID"))
      qord <- order_seqs(unique(sub$QID), split_csv(opt$`x-last-seqs`))
      sord <- order_seqs(unique(sub$SID), split_csv(opt$`y-last-seqs`))
    }
    qlay <- cumulative_layout(qlen, qord); slay <- cumulative_layout(slen, sord)
    if (is.null(qlay) || is.null(slay)) {
      warn_cli("allChr：没有可画的序列，跳过")
    } else {
      f <- paste0(opt$prefix, ".allChr.pdf")
      pl <- plot_dotplot(sub, qlay, slay, opt)
      setattr(pl, "what", "allChr")
      wh <- resolve_canvas(fit_dot_canvas(qlay, slay, opt, pl),
                           "dot-width", "dot-height", opt)
      open_pdf(f, wh[1], wh[2], opt)
      print(pl)
      invisible(dev.off())
      outs <- c(outs, f)
      canvas_log[["allChr"]] <- wh
    }
  }

  # ---- sepChr：逐 target 一页 ----
  if ("sepChr" %in% opt$.plots) {
    tg <- split_csv(opt$chroms)
    if (!length(tg)) tg <- order_seqs(unique(pairs_keep$SID), split_csv(opt$`y-last-seqs`))
    miss <- setdiff(tg, slen$id)
    if (length(miss)) warn_cli("--chroms 中这些 target 不在数据里，已跳过：",
                               paste(miss, collapse = ", "))
    tg <- tg[tg %in% slen$id]

    # 同一个 PDF 内各页尺寸必须一致，所以先把每页的布局都算出来，
    # 取所有页反推结果的最大值作为统一画布，再开设备出图
    pg <- list()
    for (tn in tg) {
      keep_q <- pairs_keep[SID == tn, unique(QID)]
      if (!length(keep_q)) next
      sub <- flt[SID == tn & QID %in% keep_q]
      if (nrow(sub) == 0L) next
      qord <- order_seqs(unique(sub$QID), split_csv(opt$`x-last-seqs`))
      qlay <- cumulative_layout(qlen, qord)
      slay <- cumulative_layout(slen, tn)
      if (is.null(qlay) || is.null(slay)) next
      pg[[length(pg) + 1L]] <- list(tn = tn, sub = sub, qlay = qlay, slay = slay)
    }

    for (i in seq_along(pg)) {
      x <- pg[[i]]
      pl <- plot_dotplot(x$sub, x$qlay, x$slay, opt)
      setattr(pl, "what", paste0("sepChr/", x$tn))
      pg[[i]]$plot <- pl
      pg[[i]]$wh <- resolve_canvas(fit_dot_canvas(x$qlay, x$slay, opt, pl),
                                   "sep-width", "sep-height", opt)
    }

    pages <- 0L
    if (!length(pg)) {
      warn_cli("sepChr：没有任何一页有内容，未生成文件")
    } else if (identical(opt$`sep-onefile`, "on")) {
      # 单文件多页：各页尺寸必须一致（pdf(onefile=TRUE) 的限制），取各页的最大值。
      # --equal-scale on 时各页仍各自等比，跨度比小的页会在面板内留白
      f <- paste0(opt$prefix, ".sepChr.pdf")
      wh <- apply(vapply(pg, function(x) x$wh, numeric(2)), 1, max)
      open_pdf(f, wh[1], wh[2], opt)
      for (x in pg) { print(x$plot); pages <- pages + 1L }
      invisible(dev.off())
      outs <- c(outs, f); sep_pages <- pages
      canvas_log[["sepChr"]] <- wh
    } else {
      # 每个 target 一个文件：各页画布按自己的数据纵横比裁剪，互不迁就
      for (x in pg) {
        f <- paste0(opt$prefix, ".sepChr.", x$tn, ".pdf")
        open_pdf(f, x$wh[1], x$wh[2], opt)
        print(x$plot)
        invisible(dev.off())
        outs <- c(outs, f); pages <- pages + 1L
        canvas_log[[paste0("sepChr.", x$tn)]] <- x$wh
      }
      sep_files <- pages
    }
  }

  # ---- coverage ----
  cov_stat <- NULL
  if ("coverage" %in% opt$.plots) {
    if (opt$.custom) {
      sord <- if (length(y_seqs)) y_seqs else order_seqs(slen$id, split_csv(opt$`y-last-seqs`))
      qsel <- if (length(x_seqs)) x_seqs else unique(flt$QID)
      sub <- flt[SID %in% sord & QID %in% qsel]
    } else {
      sub <- merge(flt, pairs_keep[, .(QID, SID)], by = c("QID", "SID"))
      sord <- order_seqs(unique(sub$SID), split_csv(opt$`y-last-seqs`))
    }
    if (nrow(sub) == 0L) {
      warn_cli("coverage：没有可用记录，跳过")
    } else {
      cov <- coverage_intervals(sub)
      cov_stat <- coverage_stats(cov, slen[id %in% sord])
      f <- paste0(opt$prefix, ".coverage.pdf")
      clay <- slen[match(sord, slen$id)]; clay <- clay[!is.na(clay$id)]
      pl <- plot_coverage_track(cov, slen, sord, opt)
      wh <- resolve_canvas(fit_cov_canvas(clay, opt, pl),
                           "cov-width", "cov-height", opt)
      open_pdf(f, wh[1], wh[2], opt)
      print(pl)
      invisible(dev.off())
      outs <- c(outs, f)
      canvas_log[["coverage"]] <- wh
    }
  }

  # ---- QC 摘要：屏幕（除非 --quiet）与日志共用同一份文本 ----
  el <- (proc.time() - t0)[["elapsed"]]
  qc <- character(0)
  add <- function(...) qc <<- c(qc, sprintf(...))
  add("== plot_paf_fmt.r QC 摘要 ==")
  add("  输入        : %s (%d 条记录)", opt$input, n_raw)
  add("  过滤链      :")
  for (i in seq_along(steps))
    add("     %-24s %7d", names(steps)[i], steps[i])
  if (!"PrimaryTag" %in% names(dt))
    add("     (无 PrimaryTag 列，secondary 过滤依赖 --min-mapq)")
  add("  正向 / 反向 : %d / %d", sum(flt$Str == "+"), sum(flt$Str == "-"))
  add("  (Q,S) 配对  : %d 个，其中 %d 个过 --min-pair-nmatch %g",
      nrow(pairs_all), nrow(pairs_keep), opt$`min-pair-nmatch`)
  if (opt$.custom)
    add("  模式        : 自定义组合（--x-seqs/--y-seqs），已跳过配对阈值筛选")
  if (!is.null(cov_stat)) {
    add("  target 覆盖 :")
    setorder(cov_stat, id)
    for (i in seq_len(nrow(cov_stat)))
      add("     %-10s %5.1f%%  (%.1f/%.1f Mb, %d 段, %d 缺口)",
          cov_stat$id[i], cov_stat$pct[i], cov_stat$covered[i] / 1e6,
          cov_stat$len[i] / 1e6, cov_stat$n_block[i], cov_stat$n_gap[i])
  }
  add("  输出        :")
  for (f in outs) add("     %-34s %8.1f KB", f, file.size(f) / 1024)
  if (!is.na(sep_pages))
    add("     (sepChr 共 %d 页)", sep_pages)
  if (!is.na(sep_files))
    add("     (sepChr 拆成 %d 个单页文件；--sep-onefile on 可合成一个多页 PDF)", sep_files)
  if (!is.null(log_file))
    add("  日志        : %s", log_file)
  add("  预设        : %s（--fit-canvas %s，--equal-scale %s，字号 %g pt）",
      opt$preset, opt$`fit-canvas`, opt$`equal-scale`, opt$`label-size`)
  if (length(canvas_log)) {
    add("  画布 (%s)   :", opt$units)
    for (k in names(canvas_log))
      add("     %-10s %7.1f x %-7.1f", k, canvas_log[[k]][1], canvas_log[[k]][2])
  }
  add("  字体        : %s%s", opt$`font-family`,
      if (nzchar(opt$`font-file`)) paste0(" (嵌入 ", opt$`font-file`, ")")
      else " (PDF base-14, 与 Arial 度量一致)")
  add("  用时        : %.2f 秒", el)
  if (!opt$quiet) cat("\n", paste(qc, collapse = "\n"), "\n\n", sep = "")
  log_write("")
  log_write(paste(qc, collapse = "\n"))

  if (!length(outs))
    stop_cli("没有生成任何输出文件。常见原因：--min-pair-nmatch ", opt$`min-pair-nmatch`,
             " 把所有 (query,target) 配对都滤掉了，或 --plots/--chroms/--x-seqs 选出的范围为空。",
             "\n  上面的 QC 摘要里 \"(Q,S) 配对\" 一行可以确认")
  log_close()
  invisible(NULL)
}

if (sys.nframe() == 0L || identical(environment(), globalenv()))
  withCallingHandlers(main(), warning = function(w)
    log_write("R WARNING: ", conditionMessage(w)))
