#!/usr/bin/bash

outPref=$1
outDir=$2

if [ ! -d "$outDir/$outPref" ]; then
  mkdir -p "$outDir/$outPref/"
fi
# $MERQURY/merqury.sh $rdDb $asmFa $outPref
mv $outPref.* logs/$outPref.* $outDir/$outPref/
# ./merqury.sh LA1593_m61.meryl assemblies/hifi_pecatR2.pri.fa eval-hifi_pecatR2.pri.hifim61

