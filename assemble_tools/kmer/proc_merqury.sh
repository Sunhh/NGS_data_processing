#!/usr/bin/bash

rdDb=$1
asmFa=$2
outPref=$3
outDir=$4

# Check that there are exactly 4 arguments
if [ "$#" -ne 4 ]; then
    echo "Error: expected 4 arguments, but got $#."
    echo "Usage: $0 <read-mer.meryl> <asm.fa> <out_prefix> <out_dir>"
    exit 1
fi

# Remove assembly meryl databases left behind from earlier runs.
### hifi_pecatR2.pri_only.wig hifi_pecatR2.pri_only.bed hifi_pecatR2.pri.meryl
baseAsmFa=$(basename "$asmFa")
prefAsmFa="${baseAsmFa%.fa}"
prefAsmFa="${prefAsmFa%.fasta}"
for ext in .meryl _only.bed _only.wig; do
  [ -e "${prefAsmFa}${ext}" ] && rm -rf "${prefAsmFa}${ext}"
done


# ./merqury.sh LA1593_m61.meryl assemblies/hifi_pecatR2.pri.fa eval-hifi_pecatR2.pri.hifim61
$MERQURY/merqury.sh $rdDb $asmFa $outPref

if [ ! -d "$outDir/$outPref" ]; then
  mkdir -p "$outDir/$outPref/"
fi
mv $outPref.* logs/$outPref.* $outDir/$outPref/

