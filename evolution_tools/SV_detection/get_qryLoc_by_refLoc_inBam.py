#!/usr/bin/python
# 2024/8/8: Edit to use 1-based positions as input and output.

import pysam
import sys
import argparse
import os;

def get_aligned_regions(bam_file, seqB, start1, end1):
  start1 -= 1;
  end1 -= 1;
  bam = pysam.AlignmentFile(bam_file, "rb");
  back_alignments = [];
  for read in bam.fetch(seqB, start1, end1):
    # Ensure the alignment overlaps target region.
    if read.is_unmapped:
      continue;
    if not (read.reference_start <= end1 and read.reference_end >= start1):
      continue;
    
    # Locate overlapping region on aligned read part.
    #   x1 = read.get_reference_sequence(); This requires MD tag.
    aligned_pairs = read.get_aligned_pairs(matches_only=False, with_seq=False);
    is_reached = 0;
    qry_name  = read.query_name;
    qry_start = None;
    qry_end   = None;
    ref_start = None;
    ref_end   = None;
    for aln1 in aligned_pairs:
      if is_reached == 0 and aln1[0] is not None and aln1[1] is not None and aln1[1] <= end1:
        qry_start = aln1[0];
        ref_start = aln1[1];
      if is_reached == 0 and aln1[0] is not None and aln1[1] is not None and aln1[1] >= start1:
        is_reached = 1;
      if aln1[1] is not None and aln1[1] <= end1 and aln1[0] is not None:
        qry_end = aln1[0];
        ref_end = aln1[1];
      if aln1[1] is not None and aln1[1] > end1:
        break;

    # Adjust query position by adding hard clipped (5) bases. Soft clipped (4) bases should have already been included.
    cigar_tuples = read.cigartuples;
    left_cigar = cigar_tuples[0];
    right_cigar = cigar_tuples[-1];
    left_clipN = 0; 
    right_clipN = 0;
    if left_cigar[0] == 5:
      left_clipN = left_cigar[1];
    if right_cigar[0] == 5:
      right_clipN = right_cigar[1];
    ### Consider the alignment direction.
    str1 = read.is_reverse;
    qry_str = '+';
    if str1 == False:
      qry_start += left_clipN;
      qry_end   += left_clipN;
    else:
      qry_len = read.infer_query_length(); # 10 means 10-bp length;
      qry_start = qry_len-qry_start-1+right_clipN;
      qry_end   = qry_len-qry_end-1+right_clipN;
      qry_start, qry_end = qry_end, qry_start;
      qry_str = '-';

    back_alignments.append([qry_name, qry_start+1, qry_end+1, qry_end-qry_start+1, qry_str, seqB, ref_start+1, ref_end+1, ref_end-ref_start+1, start1+1, end1+1, end1-start1+1]);

  return(back_alignments);


def main():
  parser = argparse.ArgumentParser(description='A script to find query locations that are mapped to a specific reference region.');
  # parser.add_argument('in_bam', type=str, help='Path to the input BAM file'); # Without '--' means this file is directly assigned without parameter tag.
  parser.add_argument('--in_bam', type=str, help='Path to the input BAM file', default='t2.bam'); # Without '--' means this file is directly assigned without parameter tag.
  parser.add_argument('--refID', type=str, help='reference sequence name', default='CLV01_Chr02'); # '--' before 'refID' means this is optional.
  parser.add_argument('--refS', type=int, help='reference sequence start', default=9744425);
  parser.add_argument('--refE', type=int, help='reference sequence end', default=9750055);
  
  args = parser.parse_args();

  if not os.path.exists(args.in_bam):
    parser.print_help(sys.stderr);
    sys.exit(1);

  alns = get_aligned_regions(args.in_bam, args.refID, args.refS, args.refE);
  for aln1 in alns:
    str_aln1 = [str(ele) for ele in aln1];
    line_aln1 = '\t'.join(str_aln1);
    print(line_aln1);

if __name__ == "__main__":
  main();


