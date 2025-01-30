# Backbone code for plotting.

## Plotting overlap of data sets.
- Venn diagram. Accept up to five sets. Result file is `output_prefix.pdf`.
```sh
Rscript plot_venn.r output_prefix  <gene_list_1.txt> <gene_list_2.txt> ...
```
- UpSet plot. Accept any number of sets. Result file is `output_prefix.pdf`.
```sh
Rscript plot_upset.r output_prefix <gene_list_1.txt> <gene_list_2.txt> ...
```

## Plotting boxplot or violin plot.
- Input files should be 1-column files storing values (comma-separated in `--pheno_files`). The file names will be used as data set names.
  - Student `t` test P values can be found in tsv files. Consider to use some other tests for non-normal distributions. Group information (`--test_groups` is `:`-separated column groups (comma-separated)).
```sh
Rscript plot_boxplot.r  --pheno_files FleshBrix_22HN_1,FleshBrix_22HN_2,FleshBrix_19YQ_1,FleshBrix_19YQ_2  --test_groups FleshBrix_22HN_1,FleshBrix_22HN_2:FleshBrix_19YQ_1,FleshBrix_19YQ_2  --output_prefix  test_boxplot
```

## Plotting barplot.
- First input column must by `Group` for each bar. The following columns are color groups within each bar.
```sh
Rscript plot_barplot.r 65K_DEL-3class 65Kb_DEL-AF_3class
```

- Plot barplot with two group layers. Add SD lines too.
```sh
Rscript plot_barplot_wiSD_twoGroups.r
```
