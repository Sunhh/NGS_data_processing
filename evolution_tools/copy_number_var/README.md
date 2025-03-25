# Gene CNV analysis.

## Use FET test
- Prepare orthologous group count file. (`synFam.cnt`). First line shows accession names.
- Prepare two-column meta file projecting accessions to populations. (`map.acc_pop`)
- Test presence frequency change between two populations.

```sh
Rscript compare_gene_expansion.r -a synFam.cnt -b map.acc_pop  -f landrace   -t cultivar  -o landrace_to_cultivar.CNV.tbl
```


