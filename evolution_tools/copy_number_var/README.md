# Gene CNV analysis.

## Test by FET.
- Prepare orthologous group count file. (`synFam.cnt`). First line shows accession names.
- Prepare two-column meta file projecting accessions to populations. (`map.acc_pop`)
- Test presence frequency change between two populations.

```sh
Rscript compare_gene_expansion.r -a synFam.cnt -b map.acc_pop  -f landrace    -t cultivar  -o landrace_to_cultivar.CNV.tbl
Rscript compare_gene_expansion.r -a synFam.cnt -b map.acc_pop  -f cordophanus -t landrace  -o cordophanus_to_landrace.CNV.tbl
```

## Combine test results.
```sh
echo -e "landrace_to_cultivar.CNV.tbl\land_cult" > meta.comparison_name
echo -e "cordophanus_to_landrace.CNV.tbl\tCLC_land" >> meta.comparison_name

perl combine_DEGs.pl meta.comparison_name 1 > combined.CNV.tbl
```

## Identify gene families with truely expansion enriched in one group.
- Those gene families with contraction enriched in the other group should be removed.
- For example, expanded in `CLV` instead of contracted in `CA`.

```sh
perl get_CLV_expansion.pl -final_label CLV_high -expan_label CLV_high CA_to_CLV.tbl > CA_to_CLV-CLV_expanded.tbl
```

