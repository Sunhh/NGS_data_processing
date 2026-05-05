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

## Plot GWAS P values with an array of genes below.
- Required columns in `gene_list.tsv`: geneID  chrID   start   end     strand  highlight
- Required columns in `gwas_Pvalues.tsv`: chr     pos     SV.length       P.value
```sh
Rscript plot_gene_array_with_gwas.r --gene_tsv gene_list.tsv --gwas_tsv gwas_Pvalues.tsv --region CLV01_Chr03:31445138-31699695 --cutP 5.96657624451305 --output CBP60_cluster.pdf
```

## Plot gene body graph with GFF3 input.
A Python script for visualizing gene structure from a GFF3 file, with support for overlaying colored elements (e.g., TE insertions, variants, functional domains). Output is a publication-ready PDF.

### Requirements
- Python 3.6+
- `matplotlib`
- `numpy`

## Usage

```sh
draw_gene.py GENE.gff3 <bp>:<label> [options]
```

### Positional arguments

| Argument | Description |
|---|---|
| `GENE.gff3` | GFF3 file containing exactly **one gene and one mRNA**. Required feature: `CDS`. Optional features: `five_prime_UTR`, `three_prime_UTR` (aliases such as `5'UTR`, `UTR5`, `5UTR` are also recognised). |
| `<bp>:<label>` | Scale bar definition: genomic length and display label, e.g. `1000:1Kb` or `500:500bp`. |

### Options

| Option | Default | Description |
|---|---|---|
| `-e ELEMENTS` | — | Elements TSV file (see format below). |
| `-o OUTPUT` | `output.pdf` | Output PDF path. |
| `--name NAME` | mRNA ID from GFF3 | Gene label shown in the top-left corner. |
| `--cds-color HEX` | `#4CAF50` | Fill colour for CDS boxes. |
| `--fig-width W` | auto | Figure width in inches. |
| `--fig-height H` | `2.2` | Figure height in inches. |
| `--xmin INT` | auto | Force the left boundary of the plot (genomic coordinate). |
| `--xmax INT` | auto | Force the right boundary of the plot (genomic coordinate). |


### Elements TSV format

A tab-separated file with **no header**, four columns:

```
seq_name    start    end    color
```


- **`seq_name`** — either the chromosome/scaffold name, or the mRNA `ID` from the GFF3.
- **`start` / `end`** — 1-based coordinates in the system of `seq_name`:
  - If `seq_name` is the chromosome name, coordinates are genomic.
  - If `seq_name` is the mRNA ID, coordinates are mRNA positions counted from the transcript 5′ end. Elements spanning introns are automatically split and drawn per-exon.
- **`color`** — any colour accepted by matplotlib: hex (`#FFA500`) or name (`orange`).

Elements overlapping CDS regions are drawn at CDS height; elements overlapping UTR regions are drawn at UTR height.

### Drawing conventions

| Feature | Appearance |
|---|---|
| CDS exon | Tall filled rectangle |
| UTR exon | Narrow outlined rectangle (white fill) |
| Intron | Black dashed line |
| Elements | Coloured rectangles overlaid on exons |
| Scale bar | Top-right corner, with end caps and label |
| Gene label | Top-left corner (mRNA ID by default) |
| Strand arrow | Blue `start` arrow below the gene label |

### Examples

```bash
# Minimal: gene structure only
python3 draw_gene.py gene.gff3 1000:1Kb -o gene.pdf

# With coloured elements using genomic coordinates
python3 draw_gene.py gene.gff3 1000:1Kb -e elements.tsv -o gene.pdf

# With coloured elements using mRNA coordinates
# (seq_name column = mRNA ID, e.g. Gene53.1)
python3 draw_gene.py gene.gff3 1000:1Kb -e elements_mrna.tsv -o gene.pdf

# Fix the plot region for multi-gene comparison at the same scale
python3 draw_gene.py geneA.gff3 1000:1Kb --xmin 1000 --xmax 20000 --fig-width 10 -o geneA.pdf
python3 draw_gene.py geneB.gff3 1000:1Kb --xmin 1000 --xmax 20000 --fig-width 10 -o geneB.pdf
```


