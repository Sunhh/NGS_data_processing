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

