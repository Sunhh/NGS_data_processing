# `assemble_tools/hifi_hic` — HiFi / Hi-C assembly format conversions

Convert the outputs of HiFi and Hi-C assemblers into repo-standard FASTA / AGP.

| script | role |
|--------|------|
| `cnvt_gfa2fa.pl` | hifiasm `*.p_ctg.gfa` -> contig FASTA |
| `get_HiCanu_ctg.pl` | extract primary contigs from a HiCanu assembly |
| `convert_3ddna_assembly_to_agp.pl` | 3D-DNA `.assembly` -> AGP |
| `cnvt_num2tigID.pl` | map numeric IDs back to tig IDs (3D-DNA/JuiceBox bookkeeping) |

```sh
perl cnvt_gfa2fa.pl hifiasm_asm.p_ctg.gfa > hifiasm_asm.p_ctg.fa
perl convert_3ddna_assembly_to_agp.pl C31_F1.assembly > C31_F1.agp
```
_Hand-written._
