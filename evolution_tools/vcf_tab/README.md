# vcf_tab — VCF ⇄ vcf.tab utilities

Convert between a proper **VCF** and a flat **`vcf.tab`** (`chr, pos, ref, <per-sample
genotype…>`), and pull per-site / per-sample fields. `vcf.tab` is the same table format used
across the repo (`bcftools query` / VariantsToTable style; see also `../../reseq_tools/cnvt_tools/`).

| script | role |
|--------|------|
| `cols2vcfTab.pl` | a `cols` SNP table → `vcf.tab` |
| `tab2vcf.pl` | `vcf.tab` → a proper VCF (`-ref_fas -in_tab`) |
| `add_ref_as_indv_in_vcfTab.pl` | add the reference allele as an extra individual column (`-ref2id REF`) |
| `slct_sites_fromVCF.pl` | keep VCF sites on selected sequences / lengths (`-in_vcf -in_keyLen`) |
| `extract_dp_per_indv.pl` | pull the per-individual read depth (DP) table from a VCF |
| `cnt_allele_withPop.pl` | count alleles per site by population / cross (`-vcf_tab`; parents + offspring) |
| `cvt_snp_to_itayNeed.pl` | reformat a SNP table to a collaborator's layout — **project-specific** |

---

## Common uses
```sh
# cols table <-> vcf.tab <-> VCF
perl cols2vcfTab.pl in_sMao.snp.cols > in_sMao.snp.tab
perl tab2vcf.pl -ref_fas ref.fa -in_tab snps.tab > snps.vcf
perl add_ref_as_indv_in_vcfTab.pl -ref2id REF in.vcfTab > in_wRef.vcfTab

# per-site / per-sample extraction
perl slct_sites_fromVCF.pl -in_vcf in.recode.vcf -in_keyLen db/chr.key_len > selected.vcf
perl extract_dp_per_indv.pl outGATK_filtSNPs_PASS.vcf > outGATK_filtSNPs_PASS.vcf.depBySample
perl cnt_allele_withPop.pl -vcf_tab family_geno.tbl > family_geno.tbl.alleleCnt
```

_Hand-written — `gen_script_index.sh` leaves it alone._
