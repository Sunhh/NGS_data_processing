#!/usr/bin/Rscript
# 2014-11-21 Edit basic function - .pattern.rangeSet() to find leftmost and right most alignment boundaries. 
# 2014-11-22 Edit Mate-pair data clean functions to separate read pairs with or without junction sequence matched. 

###########################################################
# Storage of basic sequences: 
#  1. Junction sequence used in Silin' mate-pair libraries. 
# junc_seq    <- DNAString("GGTCGATAACTTCGTATAATGTATGCTATACGAAGTTATACA") # Sequence given by Dr. Fei, which should be from Silin Zhong.
# junc_seq    <- DNAString(  "CGTATAACTTCGTATAATGTATGCTATACGAAGTTATACA") # Should be this one according to the file "P3_ec5k_ATTCCT_time2_R1.cleanPE.paired"
# For P1/3_MP2k_R1/2.paired : GGTCGATAACTTCGTATAATGTATGCTATACGAAGTTATACA
# For P1_MP5k_R1/2.paired : GGTCGATAACTTCGTATAATGTATGCTATACGAAGTTATACA
# For P1/3_cre5k_CACTCA_time1/2_R1/2.paired : CGTATAACTTCGTATAATGTATGCTATACGAAGTTATACA (/Some TCAATAACTTCGTATAATGTATGCTATACGAAGTTATACA)
# For P1/3_ec5k_ATGAGC_time1/2_R1/2.paired : CGTATAACTTCGTATAATGTATGCTATACGAAGTTATACA
#  2. Junction sequence used in outer mate-pair libraries. 
# junc_seq    <- DNAString("CTGTCTCTTATACACATCT")
#
#



############################################################
## Start: Artificial sequences (<forward> / <reverse complemented>) list.
## For the 2nd illumina sequencing technology introduced by Shan.
## sequence primer 1 (33 bp)    :                          ACACTCTTTCCCTACACGACGCTCTTCCGATCT / AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
## pcr primer 1                 : AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT / AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
## sequence primer for barcode  :                                                                       GATCGGAAGAGCACACGTCTGAACTCCAGTCAC / GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC
## sequence primer 2 (34 bp)    :                                 GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT / AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
## pcr primer 2                 : CAAGCAGAAGACGGCATACGAGAT <<<<<< GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT / AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC >>>>>> ATCTCGTATGCCGTCTTCTGCTTG
## barcode                      :                                                                                                         >>>>>>
## adapter A (1)                :  ACACTCTTTCCCTACACGACgctcttccgatc t / a gatcggaagagcGTCGTGTAGGGAAAGAGTGT
## adapter B (2)                : GTGACTGGAGTTCAGACGTGTgctcttccgatc   /   gatcggaagagcACACGTCTGAACTCCAGTCAC
## End  : Artificial sequence list.
## Initial the pattern settings.
############################################################
## For paired-end (NOT mate-pair!) sequencing technology.
## Finally, we should clean "pcr primer 2"-reverse complemented sequence at 3' end for reads from "sequence primer 1" (R1), and
##                    clean "pcr primer 1"-reverse complemented sequence at 3' end for reads from "sequence primer 2" (R2).
## Sometimes there are full "pcr primer 2"-reverse complemented artificial sequence near 3' end of R1 reads, followed by mainly polyA. Its proportion is quiet small, so I can address this case later.
############################################################


############################################################
# 初始化library
############################################################
# FUN: 声明library和变量; 
suppressMessages({
  library(ShortRead)
  library(optparse)
  library(foreach)
  library(iterators)
  library(doMC)
  registerDoMC(core=20)
  .mcoptions <- list(preschedule=TRUE, set.seed=FALSE)
}) # End suppressMessages

################################################################################
################################################################################

############################################################
# 定义默认参数; 
############################################################
# Start
############################################################

##################################################
# .get.align.opts: 获得默认比对相关参数的list, 可以指定 thres.mismatch.ratio 来修改错配比例, 暂时不专门为indel错配设置匹配分值阈值; 
#   此函数也可以用来在后面文件处理之前重新生成匹配参数用; 
#   .get.align.opts( thres.width.min=6, thres.width.up=100, thres.mismatch.ratio=0.2, match=1, mismatch=3, gapOpening=5, gapExtension=2  )
# .default.align.opts: 比对参数
# When using thres.width.min=10bp for Mate-Paired data, there are reads with 6bp adaptor ended in the reads.
.get.align.opts <- function ( 
		thres.width.min=6, thres.width.up=100, thres.mismatch.ratio=0.2, 
		match=1, mismatch=3, gapOpening=5, gapExtension=2, 
		... 
	) {
	back.opts <- lapply(list(...), function (x) x)
	names(back.opts) <- names(list(...))
	
	back.opts$thres.width.min      = thres.width.min 
	back.opts$thres.width.up       = thres.width.up 
	back.opts$thres.mismatch.ratio = thres.mismatch.ratio
		
	back.opts$match        = match
	back.opts$mismatch     = mismatch
	back.opts$gapOpening   = gapOpening
	back.opts$gapExtension = gapExtension 
	
	if ( is.null(back.opts$min.score) ) {
		# 生成与匹配片段长度(width)相关联的$min.score系列值; 
		#   允许10个碱基里面有2个(20%)错配时, score值只有8-6=2, 因此不适合简单用score来区分, 应该配合aln.width数据; 
		#   要求匹配长度至少为6个碱基, 不大于最低匹配长度时要求完全相同, 大于最低长度后允许20% mismatch差异; 
		#   暂不考虑gap的罚分结果(那个只会导致threshold更低); 
		for (i in 1:back.opts$thres.width.up) {
			if (i <= back.opts$thres.width.min) {
				# 其实前面的分数没用, 只有[thres.width.min]处的分数有用; 
				back.opts$min.score[i] <- back.opts$match * back.opts$thres.width.min
			} else {
				maxmisNum <- floor(back.opts$thres.mismatch.ratio * i)
				# 允许1个gap+剩余extension; 这个不好掌握, 先不写这个规定了; 
				back.opts$min.score[i] <- min(
					back.opts$match * (i-maxmisNum) - back.opts$mismatch * maxmisNum
				)
			}
			if (back.opts$min.score[i] < 0) back.opts$min.score[i] = 1
		}# End for
		names(back.opts$min.score) <- 1:back.opts$thres.width.up
	}# if ( is.null(back.opts$min.score) ) {
	
	return(back.opts)
}# .get.align.opts()
##################################################
# 获取全局所需默认align.opts参数; 
.default.align.opts <- .get.align.opts() 
##################################################
.reset.align.score <- function ( 
		back.opts = .default.align.opts
	) {
	# 生成与匹配片段长度(width)相关联的$min.score系列值; 
	#   允许10个碱基里面有2个(20%)错配时, score值只有8-6=2, 因此不适合简单用score来区分, 应该配合aln.width数据; 
	#   要求匹配长度至少为6个碱基, 不大于最低匹配长度时要求完全相同, 大于最低长度后允许20% mismatch差异; 
	#   暂不考虑gap的罚分结果(那个只会导致threshold更低); 
	for (i in 1:back.opts$thres.width.up) {
		if (i <= back.opts$thres.width.min) {
			# 其实前面的分数没用, 只有[thres.width.min]处的分数有用; 
			back.opts$min.score[i] <- back.opts$match * back.opts$thres.width.min
		} else {
			maxmisNum <- floor(back.opts$thres.mismatch.ratio * i)
			# 允许1个gap+剩余extension; 这个不好掌握, 先不写这个规定了; 
			back.opts$min.score[i] <- min(
				back.opts$match * (i-maxmisNum) - back.opts$mismatch * maxmisNum
			)
		}
		if (back.opts$min.score[i] < 0) back.opts$min.score[i] = 1
	}# End for
	names(back.opts$min.score) <- 1:back.opts$thres.width.up

	return(back.opts)
}# .reset.align.score()

# .default.qual.opts: 低质量相关参数
.get.qual.opts <- function(
		min.qual   = 20, min.length = 40, wind.size  = 4, 
		java_cmd_pre = 'java -cp /home/Sunhh/src/trimmomatic/ org.usadellab.trimmomatic.Trimmomatic ', # For WWZ server. 
#		java_cmd_pre = 'java -cp /home/Sunhh/tools/clean_reads/trimmomatic/ org.usadellab.trimmomatic.Trimmomatic ', # For Penguin server. 
		java_cores = -1, 
		... 
	) {
	# 先定义(...)指定参数; 
	back.opts <- lapply(list(...), function (x) x)
	names(back.opts) <- names(list(...))
	# 再用已定义参数覆盖; 
	back.opts$min.qual   <- min.qual
	back.opts$min.length <- min.length
	back.opts$wind.size  <- wind.size
	back.opts$java_cmd_pre <- java_cmd_pre
	back.opts$java_cores   <- java_cores
	# 返回list(); 
	return(back.opts)
}# .get.qual.opts()
.default.qual.opts <- .get.qual.opts() 

############################################################
# 定义默认参数; 
############################################################
# End
############################################################

############################################################
# 逐个定义子函数; (独立函数)
############################################################
# Start
############################################################

##################################################
# .unimplemented: Stop and announce this script has not been implemented. 
#      不计算，直接跳出当前运算
.unimplemented <- function() stop("UNIMPLEMENTED")
##################################################

##################################################
# .show.mem.usage1: 给出当前状态下变量的内存占用; 
#      不计算，直接跳出当前运算
.show.mem.usage <- function(units='',...) {
	# ls() can use (envir=.GlobalEnv) as input. 
	return(
		rev(
			sort(
				sapply(
					ls(...), 
					function (obj.name) {
						ss <- object.size( get(obj.name) )
						switch(
							units,
							'K'=round(ss/1024), 
							'k'=round(ss/1024), 
							'M'=round(ss/1024/1024), 
							'm'=round(ss/1024/1024), 
							'G'=round(ss/1024/1024/1024), 
							'g'=round(ss/1024/1024/1024), 
							# bytes. 
							ss
						)
					}
				)#sapply
			)
		)#rev
	)#return 
 }# .show.mem.usage1()
##################################################


##################################################
# FUN: 输出含有时间标签的记录信息
# Other: 
#   message() 作为记录输出命令比直接用cat要好一些，它能够过滤掉一些没用的信息; 
#   另一个时间标签函数: format(Sys.time(), "%Y-%m-%d %T")
## Timestampped message
.tsmsg <- function(...) {
  message("[", date(), "]: ", ...)
	write( paste( "[", date(), "]: ", ..., sep=""), file="R_log_info", append=TRUE )
}# End .tsmsg 
##################################################

##################################################
# .merge.lists: 将l1/l2合并, 合并根据为其names()的异同, 重复的取第一个; 
# Other: 
#   setdiff(x,y) 查找存在于x但不存在于y中的组分; 类似地有, union(x, y)直接合并成非冗余, intersect(x, y)查找二者共有组分(组内重复算1个); 
## Merge l1 and l2 by names
.merge.lists <- function(l1, l2) {
  new.names <- setdiff(names(l2), names(l1))
  l1[new.names] <- l2[new.names]
  l1
}# End .merge.lists()  
##################################################

##################################################
# .merge.3bands: 将l1/l2 (threebands对象)合并, 合并根据为1个threebands
## Merge l1 and l2 threebands objects. 
.merge.3bands <- function(l1, l2) {
	ll        <- list()
	ll$left   <- append(l1$left, l2$left)
	ll$middle <- append(l1$middle, l2$middle)
	ll$right  <- append(l1$right, l2$right)
	return(ll)
}# End .merge.3bands()  
##################################################

##################################################
# .strip.names: 扒掉输入变量的名字属性; 没明白这有什么价值
## Return an object sans names
.strip.names <- function(x) {
  names(x) <- NULL
  return( x )
}# End .strip.names()
##################################################

############################################################
# 逐个定义子函数; (独立函数)
############################################################
# End
############################################################

############################################################
# 逐个定义子函数; 用于多线程(foreach)的基础 (调用、依赖自定义函数，从而顺序相关)
############################################################
# Start 
############################################################

##################################################
# .get.chunk.size: 指定需要分块整体的总大小, 获取每个分块的适合大小; 可以指定分块的最大数量; 
# Other: 
#   getDoParWorkers() 根据 doPar 后台确定并行运行线程数量; 
#   
.get.chunk.size <- function(vec.length,
                           min.chunk.size=NULL, max.chunk.size=NULL,
                           max.chunks=NULL) {
  if (is.null(max.chunks)) {
    max.chunks <- getDoParWorkers()
  }
  size <- vec.length / max.chunks
  if (!is.null(max.chunk.size)) {
    size <- min(size, max.chunk.size)
  }
  if (!is.null(min.chunk.size)) {
    size <- max(size, min.chunk.size)
  }
  num.chunks <- ceiling(vec.length / size)
  actual.size <- ceiling(vec.length / num.chunks)
  return(actual.size)
}# End .get.chunk.size()
##################################################

##################################################
# .ichunk.vectors: 返回一个list对象, 该list含元素"nextElem"函数, 该函数能返回迭代所用分块内容的list; 
#   这是一种调用方法, 看起来不太直接, 但是效率还可以; 
# Other: 
#  nextElem() 是 iterators 包中的函数, 返回迭代器对象(iterator)的下一个迭代; 
#  idiv() 是 iterators 包中的函数, 返回一个iterator对象, 循环分割divide输入值, 分割方式由chunks/chunkSize控制, 可供nextElem()调用操作; 
#  
.ichunk.vectors <- function(vectors=NULL,
		min.chunk.size=NULL,
		max.chunk.size=NULL,
		max.chunks=NULL
	) {
	## Calculate number of chunks
	recycle.length <- max(sapply(vectors, length))
	actual.chunk.size <- .get.chunk.size(recycle.length, min.chunk.size, max.chunk.size, max.chunks)
	num.chunks <- ceiling(recycle.length / actual.chunk.size)

	## Make the chunk iterator
	i <- 1
	it <- idiv(recycle.length, chunks=num.chunks)
	nextEl <- function() {
		n <- nextElem(it)
		ix <- seq(i, length = n)
		i <<- i + n
		vchunks <- foreach(v=vectors) %do% v[1+ (ix-1) %% length(v)]
		names(vchunks) <- names(vectors)
		vchunks
	}
	obj <- list(nextElem = nextEl)
	class(obj) <- c("ichunk", "abstractiter", "iter")
	obj
}# End .ichunk.vectors()
##################################################

##################################################
# .chunkapply: 以chunk方式循环执行输入函数(FUN), 切割输入变量"VECTOR.ARGS", 传给FUN额外输入参数"SCALAR.ARGS"; 
#   (...)给了 .ichunk.vectors() 函数: 可以指定max.chunks, min.chunk.size, max.chunk.size 
#   "*.ARGS"需要以list类别输入; 
#   .chunkapply <- function(FUN, VECTOR.ARGS, SCALAR.ARGS=list(), MERGE=TRUE, .multicombine=FALSE, ...)
# Other: 
#   foreach 也是一个高效的循环方法, 而且能够按照顺序返回结果; 
#     .options.multicore=.mcoptions 语法来自另一个R包(doMC); 
#   append 能够向变量vector中增补元素(包括list, ShortRead对象); 
.chunkapply <- function(FUN, VECTOR.ARGS, SCALAR.ARGS=list(), MERGE=TRUE, .multicombine=FALSE, ...) {
  ## Check that the arguments make sense
  stopifnot(is.list(VECTOR.ARGS))
  stopifnot(length(VECTOR.ARGS) >= 1)
  stopifnot(is.list(SCALAR.ARGS))
  ## Choose appropriate combine function
  if (is.logical(MERGE) ) {
    if (MERGE) {
      combine.fun <- append
    } else {
      combine.fun <- foreach:::defcombine
    }
  } else {
    combine.fun <- MERGE
  }
  ## Chunk and apply, and maybe merge
  foreach(vchunk=.ichunk.vectors(vectors=VECTOR.ARGS, ...), # Old method. 
          .combine=combine.fun,
					.multicombine=.multicombine, 
          .options.multicore = .mcoptions) %dopar%
  {
    do.call(FUN, args=append(vchunk, SCALAR.ARGS))
  }
}# End .chunkapply()
##################################################

##################################################
# maybe.chunkapply: 根据 getDoParWorkers() 返回结果决定是否调用并行运算; 
#   input : 用(max.chunks=NULL)开启自动判断, 用(max.chunks=1/0)设定为单线程运行; 
#     maybe.chunkapply( FUN=多线程运行函数, VECTOR.ARGS=传给FUN()的需要拆分的参量指定如lis(rd=reads), SCALAR.ARGS=传给FUN()的不拆分参量如lis(min.length=40), ...其它传给chunkapply()的参数 )
# Other: 
#   getDoParWorkers() 根据 doPar 后台确定并行运行线程数量; 
## Only do .chunkapply if it will run in parallel
maybe.chunkapply <- function(FUN, VECTOR.ARGS, SCALAR.ARGS=list(), max.chunks=NULL, ...) {
	if ( is.null(max.chunks) ) {
		max.chunks <- getDoParWorkers()
	} else {
		# max.chunks 不能是0; 
		max.chunks <- max( ceiling(max.chunks), 1 )
		# 不考虑实际情况了; 
		# max.chunks <- min( max.chunks, getDoParWorkers() )
	}
	
	# 开始判断
  if ( max.chunks > 1) {
    .chunkapply(FUN, VECTOR.ARGS, SCALAR.ARGS, max.chunks=max.chunks, ...)
  } else {
    do.call(FUN, append(VECTOR.ARGS, SCALAR.ARGS))
  }
}# End maybe.chunkapply
##################################################

############################################################
# 逐个定义子函数; 用于多线程(foreach)的基础 (调用、依赖自定义函数，从而顺序相关)
############################################################
# End
############################################################


############################################################
# 逐个定义子函数; 数据格式转换 (调用、依赖自定义函数，从而顺序相关)
############################################################
# Start
############################################################

##################################################
# FUN: 将连续成对(pair-end在单个文件内)的reads, 分割成两个成对排列的两个list组分; 
## Takes a set of interleaved reads (or anything else) and
## de-interleaves them
.deinterleave.pairs <- function(rd) {
  stopifnot(length(rd) %% 2 == 0)
  mask <- seq(from=1, to=length(rd), by=2)
  return(list(read1=rd[mask], read2=rd[-mask]))
}# End .deinterleave.pairs()
##################################################

##################################################
## Define some missing type coercions
## 重新定义对象类别转换
#  setAs(from, to, def)定义之后, 可以用as(object, Class)来调用定义改变数据类型; 
# "ShortRead"                ==>> "DNAStringSet"
# "PhredQuality"             ==>> "FastqQuality"
# "SolexaQuality"            ==>> "SFastqQuality"
# "QualityScaledXStringSet"  ==>> "ShortReadQ" 
# "ShortReadQ"               ==>>  "QualityScaledDNAStringSet"

# "ShortRead"     ==>> "DNAStringSet"
setAs(from="ShortRead", to="DNAStringSet", def=function(from) sread(from))

# "PhredQuality"  ==>> "FastqQuality"
setAs(from="PhredQuality", to="FastqQuality", def=function(from) FastqQuality(BStringSet(from)))

# "SolexaQuality" ==>> "SFastqQuality"
setAs(from="SolexaQuality", to="SFastqQuality", def=function(from) SFastqQuality(BStringSet(from)))

# "QualityScaledXStringSet"  ==>> "ShortReadQ" 
setAs(
	from="QualityScaledXStringSet", to="ShortReadQ", 
	def=function(from) {
		q <- quality(from)
		new.quality.class <- switch(
			class(q),
			SolexaQuality="SFastqQuality",
			PhredQuality="FastqQuality",
			stop("Unknown quality type: ", class(q))
		)
		q <- as(q, new.quality.class)
		ShortReadQ(
			sread=as(from, "DNAStringSet"),
			quality=q,
			id=BStringSet(names(from))
		)
	}
)
# End setAs(from="QualityScaledXStringSet", to="ShortReadQ")

# "ShortReadQ"  ==>>  "QualityScaledDNAStringSet"
##  虽然系统有默认转换方法, 但默认方法不能保留序列名字; 
## Override the provided method to keep the sequence names
setAs(from="ShortReadQ", to="QualityScaledDNAStringSet",
	def=function (from, to = "QualityScaledDNAStringSet", strict = TRUE) {
		q <- quality(from)
		new.quality.class <- switch(
			class(q),
			SFastqQuality="SolexaQuality",
			FastqQuality="PhredQuality",
			"XStringQuality"
		)
		q <- as(q, new.quality.class)
		x <- QualityScaledDNAStringSet(sread(from), q)
		names(x) <- as.character(id(from))
		x
	}
)
# End setAs(from="ShortReadQ", to="QualityScaledDNAStringSet")

# 重新定义Quality到matrix的转换; 默认要求等长度; 
# "SolexaQuality"/"PhredQuality" ==>> "matrix"
setAs(from="SolexaQuality", to="matrix", def=function (from) as( as(from, "SFastqQuality"), "matrix") )
setAs(from="PhredQuality",  to="matrix", def=function (from) as( as(from,  "FastqQuality"), "matrix") )
##################################################

##################################################
# .read.QualityScaledDNAStringSet: 读入fastq文件并转换格式成"QualityScaledDNAStringSet"; 
## Define functions for reading fastq into standard Biostrings object
## and writing it back out. The standard functions readFastq and
## writeFastq operate on ShortRead objects. These simply wrap them in
## conversion to/from QualityScaledDNAStringSet.
.read.QualityScaledDNAStringSet <- function(filepath, format = "fastq", ...) {
	switch(
		format,
		fastq=as(readFastq(filepath, withIds=TRUE, ...), "QualityScaledDNAStringSet"),
		## Default
		stop("Unknown quality-scaled sequence format: ", format)
	)
}# End .read.QualityScaledDNAStringSet()
##################################################

##################################################
# .write.QualityScaledDNAStringSet: 将"QualityScaledDNAStringSet"写出到文件; 
#   .write.QualityScaledDNAStringSet <- function (x, filepath, append = FALSE, format = "fastq")
# Other: 
#   sink() 似乎是一个比较容易调控的写文件函数; 
.write.QualityScaledDNAStringSet <- function (x, filepath, append = FALSE, format = "fastq") {
	if(length(x) > 0) {
		sr <- as(x, "ShortReadQ")
		switch(
			format,
			fastq={
				if (!append) { unlink(filepath); }
				writeFastq(
					object=sr,
					file=filepath, mode=ifelse(append, "a", "w"),
					compress=FALSE
				)
			},
			## Default
			stop("Unknown quality-scaled sequence format: ", format)
		)
	} else {
		## Zero-length sequence; just truncate/touch the file
		sink(file=filepath, append=append)
		sink()
	}
}# End .write.QualityScaledDNAStringSet() 
##################################################

############################################################
# 逐个定义子函数; 数据格式转换 (调用、依赖自定义函数，从而顺序相关)
############################################################
# End
############################################################


############################################################
# 逐个定义子函数; 去低质量相关基础 (调用、依赖自定义函数，从而顺序相关)
############################################################
# Start
############################################################

##################################################
# .wapply.sum: 节省函数替代过程的wapply(), step by恒定为1, 不要'...'这个柔性传参; 
# http://rmazing.wordpress.com/2013/04/23/wapply-a-faster-but-less-functional-rollapply-for-vector-setups/
#   input : 输入一维数组, 以及sliding window size(width); 
#   output: 一维数组, 数组长度为 input_len-width+1
.wapply.sum <- function(x, width=1)
{
	SEQ1 <- seq_len(length(x)-width+1)
	SEQ2 <- lapply( SEQ1, function(x) x-1+seq_len(width) )
	
	OUT <- lapply(SEQ2, function(a) sum(x[a]))
	OUT <- base:::simplify2array(OUT, higher = TRUE)
	return(OUT)
}# .wapply.sum() 
##################################################

##################################################
# .longestTWind: 在一维数组中寻找最长的连续TRUE片段; 
#  input : 一个由 T/F 组成的一维数组; 
# output : 一维数组中的最大连T区间; (1,0,0) 表示无此类区间(数组内无T); 
.longestTWind <- function (x) {
	ll <- length(x)
	f.idx <- which( is.na(x) | x != TRUE )
	f.idx.r <- c(f.idx, ll+1)
	f.idx.l <- c(0   , f.idx)
	f.dist <- f.idx.r - f.idx.l - 1
	idx <- which(f.dist == max(f.dist))[1]
	# 返回最长连续TRUE的window的(起始,终止,长度), 长度=-1 表示无T-window. (1,0,0), narrow()识别end-start==-1为空; 
	back_range <- c(1,0,0)
	if ( f.dist[idx] > 0) {
		back_range <- c( f.idx.l[idx]+1, f.idx.r[idx]-1, f.dist[idx] ) 
	}
	names(back_range) <- c("start","end","length")
	return(back_range)
}# .longestTWind()
##################################################

##################################################
# .trimTailF: 在一维数组中根据左右位置坐标向内逐个去除非TRUE元素
#  input : (x=一维数组, ss=起始(左)检查位置, ee=终止(右)检查位置)
# output : 内缩去除非TRUE元素后的新位置区间; 
.trimTailF <- function (x, ss=1, ee=NULL) {
	if (is.null(ee)) ee <- length(x)
	back_range <- c(1,0,0)
	if (ss <= ee) {
		t.range <- range( which(x[ss:ee] == TRUE) )
		back_range <- c(t.range[1]+ss-1, t.range[2]+ss-1, t.range[2]-t.range[1]+1)
	}
	names(back_range) <- c("start", "end", "length")
	return(back_range)
}# .trimTailF() 
##################################################

##################################################
# .high.qual.range: 在一维数值数组中指定数值下限, 指定sliding window size, 查找最长、合格区域坐标; 每次步移1bp; 
#   input : 由质量值数字组成的一维数组, 质量值threshold, sliding window size; 
#           (x=由碱基质量数值组成的一维数组, thresV=最低允许碱基质量数值, windL=计算碱基质量平均值时的window大小)
#  output : 返回与一个array, ("start", "end", "length"), 最佳剩余(高质量)区间, (1,0,0)表示无高质量区间可用; 
#   时间分配: x-thresV:0.0004050732 secs, wapply.sum: 0.002195835 secs; generating .Tag: 0.0006275177 secs; longestTWind: 0.0008149147 secs; .trimTailF: 0.0007812977 secs; Total time: 0.004824638 secs; 
#       可见 wapply.sum() 占用了45%接近一半的时间, 以此时间计算, 25M个reads需要33.5个小时, 20个cpu也需要1.5个小时, 这个太耗费时间了, 
.high.qual.range <- function (x, thresV=0, windL=1) {
	# rollmean() 耗费时间比 rollapply(FUN=sum) 少一半多; 
	# 貌似rollmean()的速度也还可以再次改进; http://rmazing.wordpress.com/2013/04/23/wapply-a-faster-but-less-functional-rollapply-for-vector-setups/
	# 
	# x.windMean <- rollmean( x, k=windL, align="left" ) 
	# sum阈值, sum(sliding window), good(Tag) window位置, 
	x <- x-thresV
	x.windSum <- .wapply.sum(x, width=windL)
	x.windSum.Tag <- !is.na(x.windSum) & x.windSum >= 0
	# 最长连续good window 区间; 
	x.windSum.keptP <- .longestTWind(x.windSum.Tag)
	if (x.windSum.keptP[3] > 0) {
		# 如果x.windSum.keptP[3]==0, 说明没有好的window可用; 从连续window区的末端逐个trim低质量碱基; 
		x.Tag <- !is.na(x) & x >= 0
		x.windSum.keptP <- .trimTailF( x.Tag, ss=x.windSum.keptP[1], ee=x.windSum.keptP[2]+windL-1 )
	}
	names(x.windSum.keptP) <- c("start", "end", "length")
	return( x.windSum.keptP )
}# .high.qual.range()
##################################################

##################################################
# .high.qual.rangeSet: 指定Quality threshold, sliding window size, 
#   input : FastqQuality / SFastqQuality 类对象(质量值集合), 
#           质量值threshold, 
#           sliding window size; 
#           (rd=fastq_reads, thresV=最低允许碱基质量值, windL=计算碱基质量平均值时的window大小)
#   output: 返回reads的高质量区间IRanges对象; (1,0,0) 表示该行reads无高质量碱基区间; 
.high.qual.rangeSet <- function (rd.q, thresV=0, windL=1) {
	rd.q.mat <- as(rd.q, "matrix")
	# rd.q.kept.range <- t( apply( rd.q.mat, MARGIN=1,FUN=function(x) .high.qual.range(x, thresV=thresV, windL=windL)) )
	# 如果输入参数说明没什么好trim的, 那就不trim了; 
	if (thresV <= -10 | windL <= 0 ) return (IRanges(start=1, end=width(rd.q))) 
	# 如果还是需要trim, 那就来吧!!! 
	rd.q.mat <- rd.q.mat - thresV
	# 构建一个加和后的质量值矩阵; 
	n.rc <- dim(rd.q.mat)
	zero_i <- seq_len(n.rc[2]-windL+1)
	rd.q.mat.sum <- rd.q.mat[ ,zero_i ]
	if (windL > 1) {
		add_i <- seq_len(windL-1)
		for (wind.startI in add_i) {
			rd.q.mat.sum <- rd.q.mat.sum + rd.q.mat[ ,(wind.startI+zero_i) ]
		}
	}
	# 寻找最长的高质量连续窗口
	rd.q.mat.sum.Tag <- !is.na(rd.q.mat.sum) & rd.q.mat.sum >= 0
	rd.q.mat.sum.keptP <- apply(rd.q.mat.sum.Tag, MARGIN=1, FUN=.longestTWind) # 横向
	rd.q.mat.sum.keptP.ee <- vapply( rd.q.mat.sum.keptP[2,], FUN=function(x) ifelse(x<=0, x, x+windL-1), FUN.VALUE=1 )
	# 收缩窗口; 暂时不做了，可能更耗费时间; 而直接用 trimEnds() 函数非常的快！
	rd.q.kept.range <- IRanges(start=rd.q.mat.sum.keptP[1,], end=rd.q.mat.sum.keptP.ee)
	return( rd.q.kept.range )
}# .high.qual.rangeSet()
#.high.qual.rangeSet.bak <- function (rd.q, thresV=0, windL=1) {
#	rd.q.mat <- as(rd.q, "matrix")
#	rd.q.kept.range <- t( apply( rd.q.mat, MARGIN=1,FUN=function(x) .high.qual.range(x, thresV=thresV, windL=windL)) )
#	rd.q.kept.range <- IRanges(start=rd.q.kept.range[,1], end=rd.q.kept.range[,2])
#	return( rd.q.kept.range )
#}# .high.qual.rangeSet.bak()
##################################################

##################################################
# .high.qual.fq: 指定Quality threshold, sliding window size, 
#    .high.qual.fq( rd1, rd2=NULL, thresV=0, windL=1, min.length=1, max.chunks=NULL, ... ) # ... for maybe.chunkapply() 
#   本函数返回的是list, 而在foreach中, 我没找出适合的, 合并这个list的内置函数, 自写函数我怕效率低, 因此决定将多线程整合到这个函数内部; 
#   用 max.chunks=NULL 控制多线程开关, NULL为开启多线程, 并不限cpu上限, max.chunks=1表示单线程, max.chunks=0表示傻叉在设置函数; 
#   input : QualityScaledDNAStringSet类对象, 
#           质量值threshold, 
#           sliding window size; 
#   output: 返回高质量的reads, reads存储在list(R1.pair=NULL, R1.single=NULL, R2.pair=NULL, R2.single=NULL)之中; 
.high.qual.fq <- function ( rd1, rd2=NULL, thresV=0, windL=1, min.length=1, max.chunks=NULL, ... ) {
	# 返回值为 back.rd=list(); 
	back.rd <- list(R1.pair=NULL, R1.single=NULL, R2.pair=NULL, R2.single=NULL) 
	
	if (is.null(rd2)) {
		# 多线程运行(foreach)在这里, ... 可以进一步用于指定如 min.chunk.size 等参数; 运行foreach之前先清理内存; 
		#   拔掉名字, 有啥用处? 能够减小内存, 增加分解速度? 
		saved.names    <- BStringSet(names(rd1))
		rd1            <- .strip.names(rd1)
		invisible(gc())
		.tsmsg("[Msg]     Begin to run multicore-foreach. trim.low.quality.se")
		hq.rangeSet    <- maybe.chunkapply(
			FUN         = .high.qual.rangeSet, # 返回 IRanges 对象; 
			VECTOR.ARGS = list( rd=quality(rd1) ), 
			SCALAR.ARGS = list( thresV=thresV, windL=windL ), 
			MERGE       = c, # 用c()即可粘贴, 且不必两两粘贴(.multicombine=TRUE)
			.multicombine = TRUE, 
			max.chunks    = max.chunks, 
			...
		)
		.tsmsg("[Msg]     Finish to run multicore-foreach. trim.low.quality.se")
		back.rd$R1.single <- narrow(rd1, start=start(hq.rangeSet), end=end(hq.rangeSet))
		names(back.rd$R1.single) <- as.character(saved.names)
		if (min.length > 0) back.rd$R1.single <- back.rd$R1.single[ width(hq.rangeSet) >= min.length ]
	} else {
		# 多线程运行(foreach)在这里, ... 可以进一步用于指定如 min.chunk.size 等参数; 运行foreach之前先清理内存; 
		#   合并两个rd1/rd2变量到list--original.rds, 拔掉名字到rd.names内; 
		original.rds    <- list(d1=rd1, rd2=rd2)
		rm(rd1, rd2)
		# 保留每条read序列的名字到rd.names中
		rd.names        <- foreach(r=original.rds) %do% BStringSet(names(r))
		# 对应 rd1/rd2 这两个reads类别名称到list; 至此rd.names含有所有需要的名称; 
		names(rd.names) <- names(original.rds)
		# 除掉 original.rds 各个元素内部的名字(即去掉每条read序列的名字); 
		original.rds    <- Map(.strip.names, original.rds)
		#   开始运行多线程. lapply() 简化代码, 分别计算list中的每个元素; 
		invisible(gc())
		.tsmsg("[Msg]     Begin to run multicore-foreach. trim.low.quality.pe")
		.tsmsg("[Msg]       trim.low.quality.pe.rd1")
		hq.rangeSet.rds <- list()
		hq.rangeSet.rds$rd1 <- maybe.chunkapply(
			FUN = .high.qual.rangeSet, 
			VECTOR.ARGS = list(rd=quality(original.rds$rd1)), 
			SCALAR.ARGS = list(thresV=thresV, windL=windL), 
			MERGE        =c, 
			.multicombine=TRUE, 
			max.chunks   = max.chunks, 
			... 
		)
		.tsmsg("[Msg]       trim.low.quality.pe.rd2")
		invisible(gc())
		hq.rangeSet.rds$rd2 <- maybe.chunkapply(
			FUN = .high.qual.rangeSet, 
			VECTOR.ARGS = list(rd=quality(original.rds$rd2)), 
			SCALAR.ARGS = list(thresV=thresV, windL=windL), 
			MERGE        =c, 
			.multicombine=TRUE, 
			max.chunks   = max.chunks, 
			... 
		)
		.tsmsg("[Msg]     Finish to run multicore-foreach. trim.low.quality.pe")
		
		# 仅含高质量区的reads, 部分reads长度为0; 恢复 read names; 
		hq.rds <- lapply( names(original.rds), 
			function (x) {
				trimmed <- narrow( original.rds[[x]], start=start(hq.rangeSet.rds[[x]]), end=end(hq.rangeSet.rds[[x]]) ) 
				names(trimmed) <- as.character(rd.names[[x]])
				return(trimmed)
			}
		)
		names(hq.rds) <- names(rd.names)
		rm(original.rds)
		
		# 考虑 min.length 的要求; 
		if ( min.length > 0 ) {
			kept.width.r1 <- width(hq.rangeSet.rds$rd1) >= min.length
			kept.width.r2 <- width(hq.rangeSet.rds$rd2) >= min.length
			kept.pair <- kept.width.r1 & kept.width.r2
			kept.s1   <- !kept.pair & kept.width.r1
			kept.s2   <- !kept.pair & kept.width.r2
			
			back.rd$R1.pair   <- hq.rds$rd1[ kept.pair ]
			back.rd$R1.single <- hq.rds$rd1[ kept.s1 ]
			back.rd$R2.pair   <- hq.rds$rd2[ kept.pair ]
			back.rd$R2.single <- hq.rds$rd2[ kept.s2 ]
		} else {
			back.rd$R1.pair   <- hq.rds$rd1
			back.rd$R2.pair   <- hq.rds$rd2
			back.rd$R1.single <- hq.rds$rd1[1][FALSE]
			back.rd$R2.single <- hq.rds$rd2[1][FALSE]
		}
	}
	return(back.rd)
}# .high.qual.fq()
##################################################


############################################################
# 逐个定义子函数; 去低质量相关基础 (调用、依赖自定义函数，从而顺序相关)
############################################################
# End
############################################################


############################################################
# 逐个定义子函数; 去"adaptor"相关基础 (调用、依赖自定义函数，从而顺序相关)
############################################################
# Start
############################################################

##################################################
# .pattern.rangeSet: 需要输入fastq reads ("DNAStringSet"/"XStringSet"类别对象), 需要输入检索用pattern序列(adaptor)1条,
#   input : subj=adaptor_seq (pattern), 
#           rd=fastq_reads, 
#   output: threebands list for mapped region. 无匹配reads的left为reads全长; 
.pattern.rangeSet <- function (
		subj, rd, 
		align.opts=list()
	) {
	# 合并比对参数, 重复以align.opts为准
	align.opts <- .merge.lists(align.opts, .default.align.opts)
	
	rd.width <- width(rd)
	rd.num   <- length(rd)
	aln.ranges <- IRanges( start = rd.width+1, end = rd.width )

	pre.aln.mapped <- rep(FALSE, rd.num)
	for (i in 1:length(subj)) {
		# Added 2014-01-10, used for skipping trimming adaptor. 
		if (is.null(subj[i]) | is.na(subj[i]) | subj[i] == "") {
			next 
		}
		invisible(gc())

		# Use recursive pairwiseAlignment() to get the most left match
		rest.aln.mapped <- rep(TRUE, rd.num)
		tmp.rd <- rd
		while ( sum(rest.aln.mapped) ) {
			# Compare pattern(subj[i]) to reads(tmp.rd). 
			aln.fwd <- pairwiseAlignment(
				pattern = tmp.rd[rest.aln.mapped], 
				subject = subj[i], 
				type    = 'overlap', 
				substitutionMatrix=nucleotideSubstitutionMatrix(
					match    = align.opts$match, 
					mismatch = -align.opts$mismatch
				), 
				gapOpening   = -align.opts$gapOpening, 
				gapExtension = -align.opts$gapExtension
			)
			aln.fwd.pat <- pattern(aln.fwd)
			aln.width   <- width(aln.fwd.pat)
			# Determine which read is currently really matched under the score system. 
			rest.aln.mapped.2 <- apply(
				cbind( aln.width, score(aln.fwd) ), 
				MARGIN=1,
				function (x) {
					if (x[1] > 0) {
						if (x[1] > align.opts$thres.width.up) {
							x[1] = align.opts$thres.width.up
						}
						return( x[2] >= align.opts$min.score[ x[1] ] )
					} else {
						return( FALSE )
					}
				}# End function
			)
			aln.start <- start(aln.fwd.pat[rest.aln.mapped.2])
			aln.end   <- end(  aln.fwd.pat[rest.aln.mapped.2])
			# Make a T/F vector "rest.aln.mapped.3", whose length == rest.aln.mapped, and only currently matched reads are TRUE
			# This is a T/F vector for currently rest reads. 
			rest.aln.mapped.3 <- rest.aln.mapped
			rest.aln.mapped.3[rest.aln.mapped][!rest.aln.mapped.2] <- FALSE
			# Mask read data "tmp.rd" at aligned region with "N" strings. 
			subseq(tmp.rd[rest.aln.mapped.3], start=aln.start, end=aln.end) <- DNAStringSet( sapply( aln.width[rest.aln.mapped.2], function (x) { paste(rep("N", x),collapse="") } ) )

			# renew the leftmost and rightmost boundaries in variable "aln.ranges"
			pre.aln.start <- start(aln.ranges[rest.aln.mapped.3])
			pre.aln.end   <-   end(aln.ranges[rest.aln.mapped.3])
			aln.ranges[rest.aln.mapped.3] <- IRanges(
				start = pmin(pre.aln.start, aln.start), 
				end   = ifelse( pre.aln.mapped[rest.aln.mapped.3], 
					pmax(pre.aln.end, aln.end), 
					aln.end
				)
			)
			pre.aln.mapped[rest.aln.mapped.3] <- TRUE

			# Renew variable "rest.aln.mapped" for next cycle. 
			# Knock out currently not matching reads from the previously rest reads set. 
			rest.aln.mapped <- rest.aln.mapped.3
		}#End while (sum(rest.aln.mapped))

	}#End for i in 1:length(subj)

	aln.threebands <- threebands(
		IRanges(start=1, end=rd.width), 
		start = start(aln.ranges), 
		end   = end(aln.ranges)
	)
	return(aln.threebands)
}# .pattern.rangeSet() 
# 为 .pattern.rangeSet() 提供处理不同类别("ShortRead"/"QualityScaledDNAStringSet"/"QualityScaledXStringSet")数据的方法, 方便使用; 
suppressMessages({
	invisible(
		setMethod(
			".pattern.rangeSet", signature=c(rd="ShortRead"),
			function (subj, rd, align.opts) {
				callGeneric(subj, as(rd, "DNAStringSet"), align.opts)
			}
		)
	)
	invisible(
		setMethod(
			".pattern.rangeSet", signature=c(rd="QualityScaledDNAStringSet"),
			function (subj, rd, align.opts) {
				callGeneric(subj, as(rd, "DNAStringSet"), align.opts)
			}
		)
	)
	invisible(
		setMethod(
			".pattern.rangeSet", signature=c(rd="QualityScaledXStringSet"),
			function (subj, rd, align.opts) {
				callGeneric(subj, as(rd, "XStringSet"), align.opts)
			}
		)
	)
})# End suppressMessages( for setMethod(".pattern.rangeSet") )
##################################################

##################################################
# .adaptor.trimmed.pe.fq: 从fastq格式reads中去除PE建库类型的adaptor序列; 
#   应该没有多少步骤, 只是这样封装一下, 以后处理mate-paired数据时更简单一些; 
#   这个函数单线程的速度与 trimLRPattern() 相若, 只是由于做了双向match, 因此时间消耗翻倍; 
#   输出改成了list, 因此将foreach整合到函数内部; 
#    .adaptor.trimmed.pe.fq( rd1, adaptor1, rd2=NULL, adaptor2=NULL, align.opts=list(), min.length=1, use.right=TRUE, ... ) # ... used for maybe.chunkapply()
#  input : fq_reads, adaptor_seq_string, min.length, align.opts; 
#              (use.right=TRUE) 参数表示保留adaptor右侧序列; 
#          
# output : 一个list, 当输入仅有rd1时, 返回 $R1.single 这个list; 
#                    当输入包括rd2时, 返回 $R1.pair, $R1.single, $R2.pair, $R2.single 这个list; 
.adaptor.trimmed.pe.fq <- function (
		rd1,      adaptor1=NULL, 
		rd2=NULL, adaptor2=NULL, 
		align.opts=list(), 
		min.length=1, 
		max.chunks=NULL, 
		use.right=FALSE, 
		...
	) {
	# 由于后面调用函数的时候又做了一面merge, 所以这一句可有可无; 索性它不会循环很多次，所以无所谓了; 
	align.opts <- .merge.lists(align.opts, .default.align.opts)
	
	# 输入可能是单向(rd1), 也可能是双向(rd1 & rd2); 这里仅根据(rd2)的存在来判断, 这样如果漏掉了 adaptor2, 程序会报错跳出; 
	
	# 定义一个返回的 reads list, list包含4个变量(R1.pair, R1.single, R2.pair, R2.single)
	#   当输入仅有rd1时, 返回 $R1.single 这个list; 
	#   当输入包括rd2时, 返回 $R1.pair, $R1.single, $R2.pair, $R2.single 这个list; 
	back.rd <- list( R1.pair=NULL, R1.single=NULL, R2.pair=NULL, R2.single=NULL ) # 呵呵, 这个有意思, 用NULL赋值是删除, 赋值时指定则添加; 
	if ( is.null(rd2) ) {
		# 仅有 rd1 时
		# 获取 trimmed 后剩余序列位置; 
		.tsmsg("[Msg]     Begin to run multicore-foreach. trim.adaptor.se")
		invisible(gc())
		match.3bands.rd1 <- maybe.chunkapply(
			FUN            = .pattern.rangeSet, 
			VECTOR.ARGS    = list(rd=sread(rd1)), 
			SCALAR.ARGS    = list(subj=adaptor1, align.opts=align.opts), 
			MERGE          = .merge.3bands, 
			.multicombine  = FALSE, 
			max.chunks     = max.chunks, 
			... 
		)
		.tsmsg("[Msg]     Finish to run multicore-foreach. trim.adaptor.se")
		
		# 有时候需要去除位于 reads 5'端的 adaptor (例如454的index)
		if (use.right) {
			trimmed.ranges.rd1 <- IRanges(
				start=start(match.3bands.rd1$right), 
				end  =  end(match.3bands.rd1$right)
			)
		} else {
			trimmed.ranges.rd1 <- IRanges(
				start=start(match.3bands.rd1$left), 
				end  =  end(match.3bands.rd1$left)
			)
		}
		
		# 获取干净的reads; 
		back.rd$R1.single <- narrow( rd1, start=start(trimmed.ranges.rd1), end=end(trimmed.ranges.rd1), use.names=TRUE )
		if ( min.length > 0 ) back.rd$R1.single <- back.rd$R1.single[ width(trimmed.ranges.rd1) >= min.length ]
		.tsmsg( "[Rec]   Kept : ", sum(width(trimmed.ranges.rd1) >= min.length) )
	} else {
		# 同时存在 rd1/rd2 时
		# 获取 trimmed 后剩余序列位置; 

		#   开始运行多线程. lapply() 简化代码, 分别计算list中的每个元素; 
		.tsmsg("[Msg]   Begin to run multicore-foreach. trim.adaptor.pe")
		.tsmsg("[Msg]     trim.adaptor.pe.1 Start")
		invisible(gc())
		match.3bands.rd1 <- maybe.chunkapply(
			FUN            = .pattern.rangeSet, 
			VECTOR.ARGS    = list(rd=sread(rd1)), 
			SCALAR.ARGS    = list(subj=adaptor1, align.opts=align.opts), 
			MERGE          = .merge.3bands, 
			.multicombine  = FALSE, 
			max.chunks     = max.chunks, 
			... 
		)
		.tsmsg("[Msg]     trim.adaptor.pe.2 Start")
		invisible(gc())
		match.3bands.rd2 <- maybe.chunkapply(
			FUN            = .pattern.rangeSet, 
			VECTOR.ARGS    = list(rd=sread(rd2)), 
			SCALAR.ARGS    = list(subj=adaptor2, align.opts=align.opts), 
			MERGE          = .merge.3bands, 
			.multicombine  = FALSE, 
			max.chunks     = max.chunks, 
			... 
		)
		.tsmsg("[Msg]     trim.adaptor.pe.2 End")
		if (use.right) {
			trimmed.ranges.rd1 <- IRanges(
				start=start(match.3bands.rd1$right), 
				end  =  end(match.3bands.rd1$right)
			)
			trimmed.ranges.rd2 <- IRanges(
				start=start(match.3bands.rd2$right), 
				end  =  end(match.3bands.rd2$right)
			)
		} else {
			trimmed.ranges.rd1 <- IRanges(
				start=start(match.3bands.rd1$left), 
				end  =  end(match.3bands.rd1$left)
			)
			trimmed.ranges.rd2 <- IRanges(
				start=start(match.3bands.rd2$left), 
				end  =  end(match.3bands.rd2$left)
			)
		}
		.tsmsg("[Msg]   Finish to run multicore-foreach. trim.adaptor.pe")
		
		# 获取干净的reads; 
		rd1 <- narrow( rd1, start=start(trimmed.ranges.rd1), end=end(trimmed.ranges.rd1), use.names=TRUE )
		rd2 <- narrow( rd2, start=start(trimmed.ranges.rd2), end=end(trimmed.ranges.rd2), use.names=TRUE )
		
		# 根据 min.length 处理reads. 
		if ( min.length > 0 ) {
			# 本来不需要这个判断的, 不过考虑到有时候可能会有需要保留与源文件同样多的reads数量和行数, 这里做个min.length的要求, min.length==0时, 空reads也会保留下来; 
			kept.width.r1 <- width(trimmed.ranges.rd1) >= min.length
			kept.width.r2 <- width(trimmed.ranges.rd2) >= min.length
			kept.pair <- kept.width.r1 & kept.width.r2
			back.rd$R1.pair   <- rd1[ kept.pair ]
			back.rd$R2.pair   <- rd2[ kept.pair ]
			back.rd$R1.single <- rd1[ kept.width.r1 & !kept.width.r2 ]
			back.rd$R2.single <- rd2[ kept.width.r2 & !kept.width.r1 ]
			.tsmsg("[Rec] Kept : ", paste( sum(kept.pair), sum(kept.width.r1 & !kept.width.r2), sum(kept.width.r2 & !kept.width.r1), sep=":") )
		} else {
			back.rd$R1.pair   <- rd1
			back.rd$R2.pair   <- rd2
			back.rd$R1.single <- rd1[1][FALSE]
			back.rd$R2.single <- rd2[1][FALSE]
		}# End if (min.length > 0) else 
	}# End if (is.null(rd2)) else 
	
	# 返回结果reads; 
	return( back.rd )
}# .adaptor.trimmed.pe.fq() 
##################################################

##################################################
# .adaptor.trimmed.mp.fq: 从fastq格式reads中去除PE建库类型的adaptor序列; 
#    .adaptor.trimmed.mp.fq( rd1, junction.seq, rd2=NULL, align.opts=list(), min.length=1, ... ) # ... used for maybe.chunkapply()
#           对于单reads输入来说, 取出junction左右两侧reads序列作为成对序列输出; 
#           对于双reads输入来说, 取出junction左侧reads序列作为成对序列输出; 
#  input : fq_reads (ShortFastq格式), junction_seq_string, min.length, align.opts; 
#
# output : 一个list, 返回 $R1.pair + $R2.pair + R1.single + R2.single 这个list; 
#                    当输入仅有 rd1 时, R1 来自 junction 左端, R2 来自 junction 右端的反向互补序列; 
#                    当输入包括 rd2 时, R1 来自 R1     的左端, R2 来自 R2     的左端; 都不做反向互补; 
.adaptor.trimmed.mp.fq <- function (
		rd1,      rd2=NULL, 
		junction.seq=NULL, 
		align.opts=list(), 
		min.length=1, 
		max.chunks=NULL, 
		...
	) {
	# 由于后面调用函数的时候又做了一面merge, 所以这一句可有可无; 索性它不会循环很多次，所以无所谓了; 
	align.opts <- .merge.lists(align.opts, .default.align.opts)
	# 必须指定 junction.seq, 否则跳出; 
	if (is.null(junction.seq)) {
		.tsmsg("[Err] Cannot find junction sequence from junction.seq\n")
		stop("[Err] Junction seq not assigned.\n")
	}
	
	# 定义一个返回的 reads list, list包含4个变量(R1.pair, R2.pair)
	back.rd <- list( R1.pairA=NULL, R1.pairB=NULL, R1.single=NULL, R2.pairA=NULL, R2.pairB=NULL, R2.single=NULL ) # 呵呵, 这个有意思, 用NULL赋值是删除, 赋值时指定则添加; 
	if ( is.null(rd2) ) {
		# 仅有 rd1 时
		# 获取 trimmed 后剩余序列位置; 我不明白 delox 中为什么要trip掉序列名字, 也没有看出来速度变快, 这里就不效仿了; 
		invisible(gc())
		.tsmsg("[Msg]     Begin to run multicore-foreach. trim.junction.se")
		match.3bands.r1 <- maybe.chunkapply(
			FUN         = .pattern.rangeSet, 
			VECTOR.ARGS = list( rd=sread(rd1) ), 
			SCALAR.ARGS = list( subj=junction.seq, align.opts=align.opts ), 
			MERGE       = .merge.3bands, 
			.multicombine = FALSE, 
			max.chunks    = max.chunks, 
			... 
		)
		.tsmsg("[Msg]     Finish to run multicore-foreach. trim.junction.se")
		
		# When width(match.3bands.r1$left[i]) == width(rd1[i]), there should be no junction sequence found in the current read. 
		rd1.left  <- narrow(rd1, start=start(match.3bands.r1$left), end=end(match.3bands.r1$left), use.names=TRUE)
		rd1.right <- reverseComplement( narrow(rd1, start=start(match.3bands.r1$right), end=end(match.3bands.r1$right), use.names=TRUE) )

		width.left  <- width(match.3bands.r1$left)
		width.right <- width(match.3bands.r1$right)

		# 从5'端(left)起始的junction所在reads对不能要; 
		no.5start <- width.left > 0
		# "no.aln" shows the reads without junction adaptor found. 
		no.aln    <- width.left == width(rd1)
		
		if (min.length > 0) {
			save.left   <- width.left  >= min.length
			save.right  <- width.right >= min.length
			save.pairA  <- save.left & save.right
			back.rd$R1.pairA  <- rd1.left[  save.pairA ]
			back.rd$R2.pairA  <- rd1.right[ save.pairA ]
			back.rd$R1.pairB  <- rd1[1][FALSE]
			back.rd$R2.pairB  <- rd1[1][FALSE]
			back.rd$R1.single <- rd1.left[  save.left & !save.right]
			back.rd$R2.single <- rd1.right[!save.left &  save.right]
			.tsmsg("[Rec]   Kept : pairA:pairB:single1:single2")
			.tsmsg("[Rec]        : ", paste( sum(save.pairA), 0, sum(save.left & !save.right), sum(!save.left &  save.right), sep=":" ))
		} else {
			save.pairA <- no.5start  & !no.aln
			back.rd$R1.pairA  <- rd1.left[save.pairA]
			back.rd$R2.pairA  <- rd1.right[save.pairA]
			back.rd$R1.pairB  <- rd1[1][FALSE]
			back.rd$R2.pairB  <- rd1[1][FALSE]
			back.rd$R1.single <- rd1.left[no.aln]
			back.rd$R2.single <- rd1.right[!no.5start]
			.tsmsg("[Rec]   Kept : pairA:pairB:single1:single2")
			.tsmsg("[Rec]   Kept : ", paste( sum(save.pairA), 0, sum(no.aln), sum(!no.5start), sep=":" ))
		}
	} else {
		# 同时存在 rd1/rd2 时

		#   开始运行多线程. lapply() 简化代码, 分别计算list中的每个元素; 
		invisible(gc())
		.tsmsg("[Msg]   Begin to run multicore-foreach. trim.junction.pe")
		.tsmsg("[Msg]     trim.junction.pe.1 Start")
		match.3bands.r1 <- maybe.chunkapply(
			FUN         = .pattern.rangeSet, 
			VECTOR.ARGS = list( rd=sread(rd1) ), 
			SCALAR.ARGS = list( subj=junction.seq, align.opts=align.opts ), 
			MERGE       = .merge.3bands, 
			.multicombine = FALSE, 
			max.chunks    = max.chunks, 
			... 
		)
		.tsmsg("[Msg]     trim.junction.pe.1 End")
		invisible(gc())
		.tsmsg("[Msg]     trim.junction.pe.2 Start")
		match.3bands.r2 <- maybe.chunkapply(
			FUN         = .pattern.rangeSet, 
			VECTOR.ARGS = list( rd=sread(rd2) ), 
			SCALAR.ARGS = list( subj=junction.seq, align.opts=align.opts ), 
			MERGE       = .merge.3bands, 
			.multicombine = FALSE, 
			max.chunks    = max.chunks, 
			... 
		)
		.tsmsg("[Msg]     trim.junction.pe.2 End")
		.tsmsg("[Msg]   Finish to run multicore-foreach. trim.junction.pe")
		
		width.raw1 <- width(rd1)
		width.raw2 <- width(rd2)
		# 获取干净的reads; 
		rd1 <- narrow( rd1, start=start(match.3bands.r1$left), end=end(match.3bands.r1$left), use.names=TRUE )
		rd2 <- narrow( rd2, start=start(match.3bands.r2$left), end=end(match.3bands.r2$left), use.names=TRUE )
		width.new1 <- width(match.3bands.r1$left)
		width.new2 <- width(match.3bands.r2$left)
		# 不要从5'端(left端)起始junction的; 
		no.5start <- width.new1 > 0 & width.new2 > 0
		no.aln    <- (width.raw1 == width.new1) & (width.raw2 == width.new2)
		if (min.length > 0) {
			save.left  <- width.new1  >= min.length & no.5start
			save.right <- width.new2  >= min.length & no.5start
			save.pairA <- save.left & save.right & !no.aln
			save.pairB <- save.left & save.right &  no.aln
			back.rd$R1.pairA <- rd1[ save.pairA ]
			back.rd$R2.pairA <- rd2[ save.pairA ]
			back.rd$R1.pairB <- rd1[ save.pairB ]
			back.rd$R2.pairB <- rd2[ save.pairB ]
			back.rd$R1.single <- rd1[ save.left & !save.right  ]
			back.rd$R2.single <- rd2[ !save.left & save.right  ]
			.tsmsg("[Rec]   Kept : pairA:pairB:single1:single2")
			.tsmsg("[Rec]        : ", paste(sum(save.pairA), sum(save.pairB), sum(save.left & !save.right), sum(!save.left & save.right), sep=":"))
		} else {
			save.pairA <- no.5start & !no.aln
			save.pairB <- no.5start &  no.aln
			back.rd$R1.pairA  <- rd1[save.pairA]
			back.rd$R2.pairA  <- rd2[save.pairA]
			back.rd$R1.pairB  <- rd1[save.pairB]
			back.rd$R2.pairB  <- rd2[save.pairB]
			back.rd$R1.single <- rd1[!no.5start & width.new1 > 0]
			back.rd$R2.single <- rd2[!no.5start & width.new2 > 0]
			.tsmsg("[Rec]   Kept : pairA:pairB:single1:single2")
			.tsmsg("[Rec]        : ", paste(sum(save.pairA), sum(save.pairB), sum(!no.5start & width.new1 > 0), sum(!no.5start & width.new2 > 0), sep=":"))
		}
	}# End if ( is.null(rd2) ) else 
	# 返回结果reads; 
	return( back.rd )
}# .adaptor.trimmed.mp.fq() 
##################################################

############################################################
# 逐个定义子函数; 去"adaptor"相关基础 (调用、依赖自定义函数，从而顺序相关)
############################################################
# End
############################################################


###****************************************************#####
# 暂时不用的函数, 留作备份; 
###****************************************************#####
# Start
###****************************************************#####
if (FALSE) {

library(zoo)

##################################################
# wapply: 类似rollapply的功能, 但sum速度更快, mean的速度就很慢了; 
# http://rmazing.wordpress.com/2013/04/23/wapply-a-faster-but-less-functional-rollapply-for-vector-setups/
wapply <- function(x, width, by = NULL, FUN = NULL, ...)
{
	FUN <- match.fun(FUN)
	if (is.null(by)) by <- width
	
	lenX <- length(x)
	SEQ1 <- seq(1, lenX - width + 1, by = by)
	SEQ2 <- lapply(SEQ1, function(x) x:(x + width - 1))
	
	OUT <- lapply(SEQ2, function(a) FUN(x[a], ...))
	OUT <- base:::simplify2array(OUT, higher = TRUE)
	return(OUT)
}# wapply() 
##################################################

##################################################
# discard.short.reads: 过滤掉过短的reads; 
discard.short.reads <- function(rd, min.length=0) {
	kept.reads <- rd[width(rd) >= min.length]
	return(kept.rd)
}# End discard.short.reads()
##################################################


##################################################
# windSum.mat: 给定sliding window size，计算matrix的window加和值; 
#   实际
windSum.mat <- function (mat, windL=1) {
	# Need rollapply() from "zoo"
	# return( t( apply(mat, MARGIN=1, FUN=function(x) rollapply( x, width=windL, FUN=sum, by=1, align="left" ) ) ) )
	return( t( apply(mat, MARGIN=1, FUN=function(x) .wapply.sum( x, width=windL) ) ) )
}
##################################################

##################################################
# .trimTailF.se: 在一维数组中根据左右位置坐标向内逐个去除非TRUE元素(每行第1/2个元素指定SE)
.trimTailF.se <- function (x) {
	ss <- x[1]+2
	ee <- x[2]+2
	back_range <- c(1,0,0)
	if (ss <= ee) {
		t.range <- range( which(x[ss:ee] == 1) )
		back_range <- c(t.range[1]+ss-1-2, t.range[2]+ss-1-2, t.range[2]-t.range[1]+1)
	}
	names(back_range) <- c("start", "end", "length")
	return(back_range)
}# .trimTailF.se() 
##################################################

##################################################
# .trimTailF.mat: 在矩阵中，按行根据左右位置坐标向内逐个去除非TRUE元素
.trimTailF.mat <- function (x, ss=1, ee=NULL) {
	n.rc <- dim(x)
	if (is.null(ee)) {
		ee <- n.rc[2]
	}
	new.mat <- cbind(ss, ee, x)
	t(apply( new.mat, MARGIN=1, FUN=.trimTailF.se ))
}# trimTailV.mat()
##################################################

##################################################
# .longestTWind.mat: 给定数值矩阵和数值下限，按行(row)返回低于下限的单值所夹最长window位置(起始-终止-长度); 
.longestTWind.mat <- function (mat, thresV=0) {
	good_pos <- !is.na(mat) & mat-thresV >= 0
	t( apply(good_pos, MARGIN=1, FUN=.longestTWind) )
}# .longestTWind.mat()
##################################################

}# End if (FALSE) 

##################################################
# .adaptor.trimmed.rangeSet.pe: 从fastq格式reads中返回去除PE建库类型的adaptor序列后的剩余IRange位置信息; 
#   应该没有多少步骤, 只是这样封装一下, 以后处理mate-paired数据时更简单一些; 
#   这个函数单线程的速度与 trimLRPattern() 相若, 只是由于做了双向match, 因此时间消耗翻倍; 
# .adaptor.trimmed.rangeSet.pe( rd="DNAStringSet"类型的对象就可以了, adaptor, align.opts=list(), min.length=0, use.right=FALSE )
#  input : fq_reads, adaptor_seq_string, min.length, align.opts; 
#           默认使用 adaptor 左侧(5'端)剩余序列, use.right=TRUE 时, 则将同时低于 min.length & right 剩余序列长度的 left 区间换成 right 区间; 
#           因此如果想要使用 right 区间, 通常需要同时设置( use.right=TRUE, min.length=Inf )
# output : 返回剩余(trimmed)区间的位置(IRanges), 顺序和reads数量总是与输入 fq_reads 相同; 
.adaptor.trimmed.rangeSet.pe <- function (
		rd, adaptor, 
		align.opts=list(), 
		min.length=0, 
		use.right=FALSE 
	) {
	# 由于后面调用函数的时候又做了一面merge, 所以这一句可有可无; 索性它不会循环很多次，所以无所谓了; 
	align.opts <- .merge.lists(align.opts, .default.align.opts)
	
	# 调用函数获取pe.adaptor匹配位置; 
	match.3bands <- .pattern.rangeSet(subj=adaptor, rd=rd, align.opts=align.opts)
	trimmed.ranges <- IRanges(
		start=start(match.3bands$left), 
		end  =  end(match.3bands$left)
	)
	
	if ( use.right & min.length > 0 ) {
		# 加个判断就能少做不少计算; 当指定了 min.length>0 时, 可能需要使用adaptor右侧剩余序列; 
		left.width <- width(match.3bands$left)
		right.width <- width(match.3bands$right)
		kept.right <- left.width < min.length & left.width < right.width
		trimmed.ranges[ kept.right ] <- IRanges( start=start(match.3bands$right[kept.right]), end=end(match.3bands$right[kept.right]) )
	}#End if (min.length > 0)
	
	return( trimmed.ranges )
}# .adaptor.trimmed.rangeSet.pe() 
suppressMessages({
	invisible(
		setMethod(
			".adaptor.trimmed.rangeSet.pe", signature=c(rd="ShortRead"),
			function (rd, adaptor, align.opts, min.length, use.right) {
				callGeneric(as(rd, "DNAStringSet"), adaptor, align.opts, min.length, use.right)
			}
		)
	)
	invisible(
		setMethod(
			".adaptor.trimmed.rangeSet.pe", signature=c(rd="QualityScaledDNAStringSet"),
			function (rd, adaptor, align.opts, min.length, use.right) {
				callGeneric(as(rd, "DNAStringSet"), adaptor, align.opts, min.length, use.right)
			}
		)
	)
	invisible(
		setMethod(
			".adaptor.trimmed.rangeSet.pe", signature=c(rd="QualityScaledXStringSet"),
			function (rd, adaptor, align.opts, min.length, use.right) {
				callGeneric(as(rd, "DNAStringSet"), adaptor, align.opts, min.length, use.right)
			}
		)
	)
})# End suppressMessages( for setMethod(".adaptor.trimmed.rangeSet.pe") )
##################################################

###****************************************************#####
# 暂时不用的函数, 留作备份; 
###****************************************************#####
# End
###****************************************************#####


###****************************************************#####
# 暂留垃圾信息
###****************************************************#####
# Start
###****************************************************#####

##################################################
# FUN: 
# Other: 
#   pairwiseAlignment() 函数参数: 
#     pattern=序列, character vector/XString/XStringSet; 
#     subject=序列, 长度为1的字符型向量, 或XString类对象; 
#     type="global/local/overlap/global-local/local-global", 其中
#       global       - 全局比对, 匹配末端剩余gap有罚分
#       local        - 局部比对; 
#       overlap      - 全局比对, 匹配末端剩余gap无罚分; # 适用于一条序列包含于另一条的情况; 代替 vmatchPattern()
#                      似乎也包含局部匹配, 即pattern与subject末端重叠的情况
#       global-local - pattern 为global, subject 为local; 
#       local-global - pattern 为local, subject 为global; 
#     substitutionMatrix=比对所用的取代打分矩阵; 
#     gapOpening/gapExtension=打开/延伸gap的罚分; 
#   nucleotideSubstitutionMatrix() 生成打分矩阵, 可针对 DNA/RNA 字符, 含简并碱基(可指定); 
#   Map(f, ...) 是类似apply()的函数, 能够为输入向量(或list类型)执行给定函数(f)操作向量(...); 
#   pattern() 这里是处理 PairwiseAlignmentsSingleSubject 对象数据, 从 PairwiseAlignments 对象中提取 AlignedXStringSet 对象; 而非默认 ?pattern 得到的那个功能(x为"XStringPartialMatches"类对象); 
#   threebands() 功能比 narrow 更广泛一些, 返回一个list, 包含三类区域位置信息("left", "middle" 和 "right"), 与narrow对应; 
#     输入irange对象, 以irange对象坐标为参照, 根据start/end信息重新位移; 
#   pairwiseAlignment() 这个函数有些问题, 当使用type="overlap"时, 比较pattern="TTTTATATATATATATATA" vs. subject="AT"会得到匹配位置为17-18, 这显然是不合适的, 这种情况下这个函数不是很准确; 似乎还是vmatchPattern更好用一些, 可惜不支持indel. 
#   pairwiseAlignment() 匹配位置的问题, 可以通过正反两次比对来获得最靠前/靠后的匹配位置; 同时, pairwiseAlignment(type="overlap") 已经可以大约包含tail_trimLR的方法, 因此不必再做一遍trimPatternLR了. 
##################################################

###****************************************************#####
# 暂留垃圾信息
###****************************************************#####
# End
###****************************************************#####

###****************************************************#####
# 主体部分; 
###****************************************************#####
# Start
###****************************************************#####

############################################################
# 逐个定义子函数; 操作文件函数, 可看作终端函数, 适宜直接在此函数库之外(另一个文件中)使用 (调用、依赖自定义函数，从而顺序相关)
############################################################
# Start
############################################################

.add.cur.to.num <- function(cur,num) {
	num$in.R1.num <- num$in.R1.num + cur$in.R1.num
	num$in.R2.num <- num$in.R2.num + cur$in.R2.num
	num$in.R1.bp  <- num$in.R1.bp  + cur$in.R1.bp
	num$in.R2.bp  <- num$in.R2.bp  + cur$in.R2.bp
	num$out.p1.num <- num$out.p1.num + cur$in.p1.num
	num$out.p2.num <- num$out.p2.num + cur$in.p2.num
	num$out.s1.num <- num$out.s1.num + cur$in.s1.num
	num$out.s2.num <- num$out.s2.num + cur$in.s2.num
	num$out.p1.bp  <- num$out.p1.bp  + cur$in.p1.bp
	num$out.p2.bp  <- num$out.p2.bp  + cur$in.p2.bp
	num$out.s1.bp  <- num$out.s1.bp  + cur$in.s1.bp
	num$out.s2.bp  <- num$out.s2.bp  + cur$in.s2.bp
	return(num)
}# .add.cur.to.num() 

##################################################
# clean.pe.fq.file : 封装 .adaptor.trimmed.pe.fq() 和 .high.qual.fq() 两个函数; 
#                    将 .high.qual.fq() 函数替换为调用 java 程序, 这个程序速度更快一些; 
#   clean.pe.fq.file( inFqName1, outFqName1, adaptor1, inFqName2=NULL,..., RdPerYield=Inf, align.opts=list(), qual.opts=list(), ... ) # ... for maybe.chunkapply()
#     qual.opts() 包含 list( min.qual=20, min.length=40, wind.size=4 )
#
#    .high.qual.fq( rd1, rd2=NULL, thresV=0, windL=1, min.length=1, max.chunks=NULL, ... ) # ... for maybe.chunkapply() 
#    .adaptor.trimmed.pe.fq( rd1, adaptor1, rd2=NULL, adaptor2=NULL, align.opts=list(), min.length=1, use.right=TRUE, ... ) # ... used for maybe.chunkapply()
clean.pe.fq.file <- function (
		inFqName1,      outFqName1,      adaptor1, 
		inFqName2=NULL, outFqName2=NULL, adaptor2=NULL, 
		RdPerYield=100e6, 
		min.length=40, 
		align.opts=list(), use.right=FALSE,  
		qual.opts=list(), 
		...
	) {
	invisible(gc())
	
	# 合并match比对参数, 去除了原"errRate=0.2"参数, 改为直接使用"align.opts"来控制一切(可以用.get.align.opts()来重新生成相关参数); 
	# 由于qual.opts相关参数需要在每条reads水平逐个使用, 因此如果每次都用.merge.lists整合, 会非常浪费计算资源, 因此仅在主函数中修改此相关参量; 
	align.opts <- .merge.lists(align.opts, .default.align.opts)
	align.opts <- .reset.align.score( back.opts=align.opts )
	qual.opts  <- .merge.lists(qual.opts, .default.qual.opts)
	
	num.rec <- list(
		in.R1.num=0, in.R2.num=0, in.R1.bp=0, in.R2.bp=0, 
		out.p1.num=0, out.s1.num=0, out.p2.num=0, out.s2.num=0, 
		out.p1.bp=0, out.s1.bp=0, out.p2.bp=0, out.s2.bp=0
	)
	cur.rec <- list(
		in.R1.num=0, in.R2.num=0, in.R1.bp=0, in.R2.bp=0, 
		out.p1.num=0, out.s1.num=0, out.p2.num=0, out.s2.num=0, 
		out.p1.bp=0, out.s1.bp=0, out.p2.bp=0, out.s2.bp=0
	)
	if (is.null(inFqName2)) 
	{
		# single cleaning
		.tsmsg( "[Rec] Begin to clean single end file: ", inFqName1 )
		# 逐步读取文件; 
		# 换用郑轶修改过的 trimmomatic 软件去除低质量, 效率高很多; "java_cmd_pre"
		# 如果输出文件已存在, 将其删除重写; 
		oFqName <- list() 
		# oFqName$p1 <- paste0(outFqName1,'.paired')
		oFqName$s1 <- paste0(outFqName1,'.single')
		lapply( oFqName, function (fn) { if (file.exists(fn)) file.remove(fn) } )
		.tsmsg( "[Rec]   New output file: ", oFqName$s1 )
		
		# 去除低质量; 
		if (qual.opts$min.qual > 0 & qual.opts$wind.size > 0) {
			# 如果没有指定quality format, 就预先读进来1w条猜测一下; 
			.tsmsg("[Msg]   Trimming SE quality Start : ")
			if ( is.null(qual.opts$q.format) ) {
				.tsmsg("[Rec]     Guessing the quality format.")
				inFh1 <- FastqStreamer(inFqName1, n=1e4)
				inRd1 <- yield(inFh1)
				qual.opts$q.format <- switch(
					class( quality(inRd1) ), 
					SFastqQuality="-phred64",
					FastqQuality="-phred33", 
					stop("Unknown quality-scaled sequence format: ", class( quality(inRd1) ))
				)
				close(inFh1)
				.tsmsg("[Rec]     The quality format should be : ", qual.opts$q.format)
			}
			
			oFqName$s1.q <- paste0(outFqName1, '.single.hq')
			java.cores <- ''
			if (qual.opts$java_cores > 1) {
				java.cores <- paste("-threads ", as.integer(qual.opts$java_cores))
			}
			system(
				paste(
					qual.opts$java_cmd_pre,"SE", 
					java.cores, qual.opts$q.format, 
					inFqName1, oFqName$s1.q, 
					paste('SLIDINGWINDOW', qual.opts$wind.size, qual.opts$min.qual, sep=':'), 
					paste('MINLEN', qual.opts$min.length ,sep=':') 
				)
			)
			.tsmsg("[Msg]   Trimming SE quality End: ")
		}# End if 低质量 ; 
		
		# 没有指定去除adaptor; 
		if (is.null(adaptor1)) {
			# 没有指定adaptor1, 不是"", 那就不用去除adaptor了; 
			if (is.null(oFqName$s1.q)) {
				# 啥也没让俺干! 
				.tsmsg("[Rec] Nothing has been done. Check parameter and the original single file: ", inFqName1)
				return (1)
			}else{
				# 只是去除了low quality. 
				file.rename( from=oFqName$s1.q, to=oFqName$s1 )
				return (0)
			}
		}
		
		# 开始去除adaptor
		.tsmsg("[Msg]   Trimming SE adaptors Start : ")
		if (!is.null(oFqName$s1.q)) {
			inFh1 <- FastqStreamer(oFqName$s1.q, n=RdPerYield)
		}else{
			inFh1 <- FastqStreamer(inFqName1, n=RdPerYield)
		}
		
		while ( length( rd1 <- yield(inFh1) ) ) {
			## reads 读取完毕, 自此展开各个多线程运行; 
			.tsmsg("[Msg]     Cycle.")
			#  Trmming adaptors. 
			clean.fq <- .adaptor.trimmed.pe.fq(
				rd1=rd1, adaptor1=adaptor1, 
				align.opts=align.opts, min.length=min.length, use.right=use.right, 
				... 
			)
			# 输出结果; ()
			if (length(clean.fq$R1.single) > 0) writeFastq(clean.fq$R1.single,   file=oFqName$s1, mode="a", compress=FALSE)
		}# End while( length( rd1 <- yield(inFh1) ) ) 
		close(inFh1)
		.tsmsg("[Msg]   Trimming SE adapter End.")
		# Edit here 
	} else {
		# paired cleaning
		.tsmsg( "[Rec] Begin to clean paired end files: ", inFqName1, " and ", inFqName2, " ." )
	
		# 如果输出文件已存在, 将其删除重写; 
		oFqName <- list() 
		oFqName$p1 <- paste0(outFqName1,'.paired')
		oFqName$s1 <- paste0(outFqName1,'.single')
		oFqName$p2 <- paste0(outFqName2,'.paired')
		oFqName$s2 <- paste0(outFqName2,'.single')
		lapply( oFqName, function (fn) { if (file.exists(fn)) file.remove(fn) } )
		.tsmsg( "[Rec]   New output files R1: ", oFqName$p1, " and ", oFqName$s1 )
		.tsmsg( "[Rec]   New output files R2: ", oFqName$p2, " and ", oFqName$s2 )
		
		# 去除低质量; 
		if (qual.opts$min.qual > 0 & qual.opts$wind.size > 0) {
			# 如果没有指定quality format, 就预先读进来1w条猜测一下; 
			.tsmsg("[Msg]   Trimming PE quality Start : ")
			if ( is.null(qual.opts$q.format) ) {
				.tsmsg("[Rec]     Guessing the quality format.")
				inFh1 <- FastqStreamer(inFqName1, n=1e4)
				inRd1 <- yield(inFh1)
				qual.opts$q.format <- switch(
					class( quality(inRd1) ), 
					SFastqQuality="-phred64",
					FastqQuality="-phred33", 
					stop("Unknown quality-scaled sequence format: ", class( quality(inRd1) ))
				)
				close(inFh1)
				.tsmsg("[Rec]     The quality format should be : ", qual.opts$q.format)
			}
			oFqName$s1.q <- paste0(outFqName1, '.single.hq')
			oFqName$p1.q <- paste0(outFqName1, '.paired.hq')
			oFqName$s2.q <- paste0(outFqName2, '.single.hq')
			oFqName$p2.q <- paste0(outFqName2, '.paired.hq')
			
			java.cores <- ''
			if (qual.opts$java_cores > 1) {
				java.cores <- paste("-threads ", as.integer(qual.opts$java_cores))
			}
			system(
				paste(
					qual.opts$java_cmd_pre, "PE", 
					java.cores, qual.opts$q.format, 
					inFqName1, inFqName2, 
					oFqName$p1.q, oFqName$s1.q, 
					oFqName$p2.q, oFqName$s2.q, 
					paste('SLIDINGWINDOW', qual.opts$wind.size, qual.opts$min.qual, sep=':'), 
					paste('MINLEN', qual.opts$min.length ,sep=':')
				)
			)# End system command. 
			.tsmsg("[Msg]   Trimming quality End : ")
		}# End if 低质量 
		
		# 没有指定去除 adaptor; 
		if (is.null(adaptor1)) {
			# 没有指定adaptor1, 不是"", 那就不用去除adaptor了; 
			if (is.null(oFqName$s1.q)) {
				# 啥也没让俺干! 
				.tsmsg("Nothing has been done. Check parameter and the original PE files: ", inFqName1, " & ", inFqName2)
				return (1)
			}else{
				# 只是去除了low quality. 
				file.rename( from=oFqName$p1.q, to=oFqName$p1 )
				file.rename( from=oFqName$p2.q, to=oFqName$p2 )
				file.rename( from=oFqName$s1.q, to=oFqName$s1 )
				file.rename( from=oFqName$s2.q, to=oFqName$s2 )
				return (0)
			}
		}
		
		# 开始去除 adaptor
		if (is.null(oFqName$p1.q)) {
			.tsmsg("[Msg] Trimming adaptor from : ", paste(inFqName1, inFqName2, sep=":"))
			inFh.p1 <- FastqStreamer(inFqName1, n=RdPerYield)
			inFh.p2 <- FastqStreamer(inFqName2, n=RdPerYield)
			inFh.s1 <- NULL
			inFh.s2 <- NULL
		}else{
			.tsmsg("[Msg] Trimming adaptor from : ", paste(oFqName$p1.q, oFqName$p2.q, oFqName$s1.q, oFqName$s2.q, sep=":"))
			inFh.p1 <- FastqStreamer(oFqName$p1.q, n=RdPerYield)
			inFh.p2 <- FastqStreamer(oFqName$p2.q, n=RdPerYield)
			inFh.s1 <- FastqStreamer(oFqName$s1.q, n=RdPerYield)
			inFh.s2 <- FastqStreamer(oFqName$s2.q, n=RdPerYield)
		}

		# 逐步读取文件; 
		# 读取并处理 paired 文件; 
		.tsmsg("[Msg]   Trimming PE adaptor Start : ")
		# while ( length( rd1.p <- as(yield(inFh.p1), "QualityScaledDNAStringSet") ) ) {
		while ( length( rd1.p <- yield(inFh.p1) ) ) {
			.tsmsg("[Msg]     Cycle.")
			rd2.p <- yield(inFh.p2)
			#  Trmming adaptors. 
			# 高质量成对reads的trimming; 
			.tsmsg( "[Msg]     Begin trimming adaptor pairs: " )
			clean.fq.P <- .adaptor.trimmed.pe.fq( 
				rd1=rd1.p, adaptor1=adaptor1, 
				rd2=rd2.p, adaptor2=adaptor2, 
				align.opts=align.opts, min.length=min.length, use.right=use.right, 
				... 
			)
			if (length(clean.fq.P$R1.pair) > 0) writeFastq(clean.fq.P$R1.pair,   file=oFqName$p1, mode="a", compress=FALSE)
			if (length(clean.fq.P$R2.pair) > 0) writeFastq(clean.fq.P$R2.pair,   file=oFqName$p2, mode="a", compress=FALSE)
			if (length(clean.fq.P$R1.single) > 0) writeFastq(clean.fq.P$R1.single,   file=oFqName$s1, mode="a", compress=FALSE)
			if (length(clean.fq.P$R2.single) > 0) writeFastq(clean.fq.P$R2.single,   file=oFqName$s2, mode="a", compress=FALSE)
		}#End paired while 
		close(inFh.p1)
		close(inFh.p2)
		.tsmsg("[Msg]   Trimming PE adaptor End : ")

		# clean.fq$R1.single, $R2.single 是保留下的single序列; 
		
		# 读取并处理 single 文件; 
		if (!is.null(inFh.s1)) {
			.tsmsg("[Msg]   Trimming SE1 adaptor Start : ")
			while ( length( rd1.s <- yield(inFh.s1) ) ) {
				.tsmsg("[Msg]     Cycle.")
				#  Trmming adaptors. 
				# 高质量 single reads 的 trimming; 
				clean.fq.S <- .adaptor.trimmed.pe.fq( 
					rd1=rd1.s, adaptor1=adaptor1, 
					align.opts=align.opts, min.length=min.length, use.right=use.right, 
					... 
				)
				if (length(clean.fq.S$R1.single) > 0) writeFastq(clean.fq.S$R1.single,   file=oFqName$s1, mode="a", compress=FALSE)
			}# End while (rd1.s)
			close(inFh.s1)
			.tsmsg("[Msg]   Trimming SE1 adaptor End : ")
		}#End if Single-1
		if (!is.null(inFh.s2)) {
			.tsmsg("[Msg]   Trimming SE2 adaptor Start : ")
			while ( length( rd2.s <- yield(inFh.s2) ) ) {
				.tsmsg("[Msg]     Cycle.")
				#  Trmming adaptors. 
				# 高质量 single reads 的 trimming; 
				clean.fq.S <- .adaptor.trimmed.pe.fq( 
					rd1=rd2.s, adaptor1=adaptor2, 
					align.opts=align.opts, min.length=min.length, use.right=use.right, 
					... 
				)
				if (length(clean.fq.S$R1.single) > 0) writeFastq(clean.fq.S$R1.single,   file=oFqName$s2, mode="a", compress=FALSE)
			}# End while (rd2.s)
			close(inFh.s2)
			.tsmsg("[Msg]   Trimming SE2 adaptor End : ")
		}#End if Single-2
	}# if (is.null(inFqName2)) else 
	.tsmsg("[Rec] Over.")
	return(0)
}# clean.pe.fq.file() 
##################################################

##################################################
# clean.mp.fq.file : 处理 clean 过的 mate-paired reads 文件; 
#   封装 .adaptor.trimmed.mp.fq( rd1, junction.seq, rd2=NULL, align.opts=list(), min.length=1, ... ) # ... used for maybe.chunkapply()
# For align.opts$thres.width.min selection: I grepped 1M R1 reads with 6bp w/ 0.2 rate, and it is fine. (FP=(270+196)/1M~=4.7e-4). It is acceptable. 

clean.mp.fq.file <- function (
		inFqName1,      outFqName1, inFqName2=NULL, 
		junction.seq=NULL, 
		RdPerYield=100e6, 
		align.opts=list(), 
		min.length=40,
		...
	) {
	invisible(gc())
	
	# 合并match比对参数, 去除了原"errRate=0.2"参数, 改为直接使用"align.opts"来控制一切(可以用.get.align.opts()来重新生成相关参数); 
	align.opts <- .merge.lists(align.opts, .default.align.opts)
	align.opts <- .reset.align.score( back.opts=align.opts )
	if (is.null(junction.seq)) {
		.tsmsg("[Err] Cannot find junction sequence from junction.seq\n")
		stop("[Err] Junction seq not found.\n")
	}
	
	# 预设值输出文件; 如果输出文件已存在, 将其删除重写; 
	# p1a/p2a are reads with junction sequences trimmed, p1b/p2b are reads without junction sequences found. 
	oFqName <- list()
	oFqName$p1a <- paste0(outFqName1,'.p1a')
	oFqName$p2a <- paste0(outFqName1,'.p2a')
	oFqName$p1b <- paste0(outFqName1,'.p1b')
	oFqName$p2b <- paste0(outFqName1,'.p2b')
	oFqName$s1  <- paste0(outFqName1,'.s1')
	oFqName$s2  <- paste0(outFqName1,'.s2')
	lapply( oFqName, function (fn) { if (file.exists(fn)) file.remove(fn) } )
	.tsmsg( "[Rec]   Output clean pairs w/  junction sequences trimmed: ", paste(oFqName$p1a, oFqName$p2a, ":") )
	.tsmsg( "[Rec]   Output       pairs w/o junction sequences found  : ", paste(oFqName$p1b, oFqName$p2b, ":") )
	.tsmsg( "[Rec]   Output clean singles   after trimming            : ", paste(oFqName$s1 , oFqName$s2 , ":") )
	
	if (is.null(inFqName2))
	{
		# single cleaning
		.tsmsg( "[Rec] Begin to clean single end file: ", inFqName1 )
		
		# 开始处理 junction sequence. 
		.tsmsg("[Msg]   Trimming SE-MP junction Start : ")
		inFh1 <- FastqStreamer(inFqName1, n=RdPerYield)
		
		# while ( length( rd1 <- as(yield(inFh1), "QualityScaledDNAStringSet") ) ) {
		while ( length( rd1 <- yield(inFh1) ) ) {
			## reads 读取完毕, 自此展开各个多线程运行; 
			.tsmsg("[Msg]     Cycle.")
			#  Trmming junction sequence. 
			clean.fq <- .adaptor.trimmed.mp.fq( 
				rd1=rd1, junction.seq=junction.seq, 
				align.opts=align.opts, min.length=min.length, 
				... 
			)
			# 输出结果; ()
			if ( length( clean.fq$R1.pairA  ) > 0 ) writeFastq(object=clean.fq$R1.pairA , file=oFqName$p1a, mode="a", compress=FALSE) 
			if ( length( clean.fq$R2.pairA  ) > 0 ) writeFastq(object=clean.fq$R2.pairA , file=oFqName$p2a, mode="a", compress=FALSE)
			if ( length( clean.fq$R1.single ) > 0 ) writeFastq(object=clean.fq$R1.single, file=oFqName$s1 , mode="a", compress=FALSE)
			if ( length( clean.fq$R2.single ) > 0 ) writeFastq(object=clean.fq$R2.single, file=oFqName$s2 , mode="a", compress=FALSE)
		}# End while( length( rd1 <- as(yield(inFh1), "QualityScaledDNAStringSet") )
		close(inFh1)
		.tsmsg("[Msg]   Trimming SE-MP junction End.")
	} else {
		# paired cleaning
		.tsmsg( "[Rec] Begin to clean paired end files: ", paste(inFqName1, inFqName2, sep=":") )
		
		# 开始去除 junction
		.tsmsg("[Msg]   Trimming PE-MP junction Start : ")
		inFh.p1 <- FastqStreamer(inFqName1, n=RdPerYield)
		inFh.p2 <- FastqStreamer(inFqName2, n=RdPerYield)
		while ( length( rd1.p <- yield(inFh.p1) ) ) {
			.tsmsg("[Msg]     Cycle.")
			rd2.p <- yield(inFh.p2)
			#  Trmming adaptors. 
			#  Trmming junction sequence. 
			clean.fq.P <- .adaptor.trimmed.mp.fq(
				rd1=rd1.p, junction.seq=junction.seq, rd2=rd2.p, 
				align.opts=align.opts, min.length=min.length, 
				... 
			)
			if (length(clean.fq.P$R1.pairA)  > 0) writeFastq(object=clean.fq.P$R1.pairA, file=oFqName$p1a, mode="a", compress=FALSE)
			if (length(clean.fq.P$R2.pairA)  > 0) writeFastq(object=clean.fq.P$R2.pairA, file=oFqName$p2a, mode="a", compress=FALSE)
			if (length(clean.fq.P$R1.pairB)  > 0) writeFastq(object=clean.fq.P$R1.pairB, file=oFqName$p1b, mode="a", compress=FALSE)
			if (length(clean.fq.P$R2.pairB)  > 0) writeFastq(object=clean.fq.P$R2.pairB, file=oFqName$p2b, mode="a", compress=FALSE)
			if (length(clean.fq.P$R1.single) > 0) writeFastq(object=clean.fq.P$R1.single, file=oFqName$s1, mode="a", compress=FALSE)
			if (length(clean.fq.P$R1.single) > 0) writeFastq(object=clean.fq.P$R2.single, file=oFqName$s2, mode="a", compress=FALSE)
		}#End paired while 
		close(inFh.p1)
		close(inFh.p2)
		.tsmsg("[Msg]   Trimming PE-MP junction End : ")
	}# if (is.null(inFqName2)) else 
	.tsmsg("[Rec] Over.")
	return(0)
}# clean.mp.fq.file() 
##################################################


##################################################
## 测试使用多线程(foreach)效果
##################################################
## Start
##################################################
## 根据检查比较, 可见当 reads 总数量不高于 3e5 时, 多线程的优势不明显, 但当 reads 总数超过 3e6 时, 多线程优势开始体现; 
## 10次运行均值: readsNum=3e5    time.1.proc= 39.8961 s    time.10.proc.default=34.1348 s    time.10.proc.sep=44.9213 s
## 10次运行均值: readsNum=3e6    time.1.proc=439.6487 s    time.10.proc.default=86.9671 s    time.10.proc.sep=92.3013 s
if (FALSE) {

fq1name <- "A32_CAACTA_L005_R1.ndupB"
fq2name <- "A32_CAACTA_L005_R2.ndupB"
pat1 <- DNAString("AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCAACTAATCTCGTATGCCGTCTTCTGCTTG")
pat2 <- DNAString("AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT")
o1name <- "A32_CAACTA_L005_R1.ndupB.clean"
o2name <- "A32_CAACTA_L005_R2.ndupB.clean"


system.time(clean.pe.fq.file (
		inFqName1=fq1name, outFqName1=o1name, adaptor1=pat1, 
		inFqName2=fq2name, outFqName2=o2name, adaptor2=pat2, 
		RdPerYield=100e6, 
		align.opts=list(), use.right=FALSE, 
		qual.opts=list(), 
		max.chunks=16
))	

system.time(clean.pe.fq.file ( inFqName1=fq1name, outFqName1=o1name, adaptor1=pat1, RdPerYield=5e6, align.opts=list(), use.right=FALSE, qual.opts=list(),  max.chunks=16))

gc()
source("using_subfunc.R")
system.time(clean.pe.fq.file ( inFqName1=fq1name, outFqName1=o1name, inFqName2=fq2name, outFqName2=o2name, adaptor1=pat1, adaptor2=pat2, RdPerYield=5e5, align.opts=list(), use.right=FALSE, qual.opts=list(),  max.chunks=30, min.chunk.size=1000))


}

##################################################
## 测试使用多线程(foreach)效果
##################################################
## End
##################################################

############################################################
# 逐个定义子函数; 操作文件函数, 可看作终端函数, 应该挪移到此函数库之外(另一个文件中)使用 (调用、依赖自定义函数，从而顺序相关)
############################################################
# End
############################################################


###****************************************************#####
# 主体部分; 
###****************************************************#####
# End
###****************************************************#####

