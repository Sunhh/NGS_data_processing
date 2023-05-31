rFa=$1
rGff=$2
qFa=$3

cp -p $rFa ./r.fa
cp -p $rGff ./r.p.gff3
cp -p $qFa ./q.fa

export PATH=/data/Sunhh/src/align/anchorwave/anchorwave:$PATH
export PATH=/data/Sunhh/src/align/last/install/v869/bin:$PATH
export PATH=/data/Sunhh/src/align/mummer/install/v4.0.0rc/bin:$PATH

java -jar /home/Sunhh/bin/picard.jar CreateSequenceDictionary  -R r.fa  -O r.dict
deal_fasta.pl -attr key:len q.fa > q.fa.kl
deal_fasta.pl -attr key:len r.fa > r.fa.kl
ln -s /home/Sunhh/tools/github/NGS_data_processing/evolution_tools/SV_detection/ hs_tool

minimap2 -t 20 -x asm20 -c r.fa q.fa | perl hs_tool/fmt_paf.pl > q2r.paf
perl hs_tool/make_fakeCDS_fromPAF.pl s0 q2r.paf r.fa
anchorwave gff2seq -r r.fa -i s0.anchor.gff3 -o s0.anchor.fa
perl hs_tool/find_nonOvlCDS.pl s0.anchor.gff3 r.p.gff3 s0add
cat s0.anchor.gff3 s0add.gff3 > s1.anchor.gff3
anchorwave gff2seq -r r.fa -i s1.anchor.gff3 -o s1.anchor.fa
deal_fasta.pl s1.anchor.fa -drawByList -drawList s0add.list -drawWhole -drawLcol 0 > s0add.c.fa
minimap2 -t 20 -a -x asm20 -p 0.8 -N 20 r.fa s0.anchor.fa > s0.toR.sam
minimap2 -t 20 -a -x asm20 -p 0.8 -N 20 q.fa s0.anchor.fa > s0.toQ.sam
minimap2 -t 20 -a -x splice -G 25k -k 12 -p 0.8 -N 20 r.fa s0add.c.fa > s0add.toR.sam
minimap2 -t 20 -a -x splice -G 25k -k 12 -p 0.8 -N 20 q.fa s0add.c.fa > s0add.toQ.sam
samtools merge -O SAM s1.toR.sam s0.toR.sam s0add.toR.sam
samtools merge -O SAM s1.toQ.sam s0.toQ.sam s0add.toQ.sam
anchorwave proali -t 20 -i s1.anchor.gff3 -r r.fa -a s1.toQ.sam -as s1.anchor.fa -ar s1.toR.sam -s q.fa -n t1.anchors -f t1_anc.maf -w 30000 -R 1 -Q 1 -ns
perl /home/Sunhh/tools/github/NGS_data_processing/software_fix/anchorwave/fix_awMAF.pl t1_anc.maf > t1_anc_fix.maf
mv t1_anc_fix.maf t1_anc.maf
maf-convert blasttab t1_anc.maf > t1_anc.maf.blasttab
perl hs_tool/get_shrt_or_ident_mafTab.pl 30000 99 t1_anc.maf.blasttab > t1_anc.maf.blasttab.good
deal_table.pl t1_anc.maf.blasttab -kSrch_idx t1_anc.maf.blasttab.good -kSrch_drop -kSrch_line > t1_anc.maf.blasttab.todo
perl hs_tool/mmp2Aln_anchorsInMafTab.pl t1_mmp q.fa r.fa t1_anc.maf.blasttab.todo
perl hs_tool/restore_sam_position.pl q.fa.kl r.fa.kl t1_mmp.sam > t1_mmp_res.sam
perl hs_tool/filter_maf_by_tab.pl t1_anc.maf.blasttab.good q.fa.kl t1_anc.maf > t1_anc_good.maf
perl hs_tool/join_maf_blk.pl t1_anc_good.maf > t1_anc_jn.maf
perl hs_tool/rm_0span_maf.pl t1_anc_jn.maf > t1_anc_jn1.maf
maf-convert sam -r "ID:q" -f r.dict t1_anc_jn1.maf > t1_anc_jn1.sam
samtools merge -f -O SAM t2.sam t1_anc_jn1.sam t1_mmp_res.sam
perl hs_tool/fix_sam_cigarID.pl t2.sam > t2_fix.sam
mv t2_fix.sam t2.sam
python hs_tool/sam2delta.py t2.sam --ref_fa $PWD/r.fa  --qry_fa $PWD/q.fa
samtools sort -@ 10 -o t2.bam t2.sam
samtools index t2.bam
rm -f q2r.paf s0.anchor* s0* s1.to*sam t1_anc.maf t1_mmp* t1_anc*maf t1_anc*sam r.dict r.fa.kl q.fa.kl
nucdiff --proc 1 --delta_file t2.sam.delta $PWD/r.fa $PWD/q.fa ./ o
rm -f q.fa r.fa r.p.gff3

