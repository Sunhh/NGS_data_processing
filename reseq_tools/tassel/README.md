# tassel — helper for TASSEL GWAS

`cnvt_col_to_TasselTaxaList.pl` converts a `cols` genotype/sample table into a TASSEL taxa
list, for running association in **TASSEL** as an alternative to the FaST-LMM pipeline in
[`../gwas_tools/`](../gwas_tools/).

```sh
perl cnvt_col_to_TasselTaxaList.pl in.cols > taxa.list
```

_Hand-written — `gen_script_index.sh` leaves it alone._
