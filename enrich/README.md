# How to do GO enrichment from a gene list.

## Input files:
- Subset (test) gene list: in\_geneID.list; One column file without headers.
- Background gene to GO annotation file: in\_whole\_genome\_gene.annot; Two-column GO annotation from Blast2GO.
- OBO file: pub-go.obo; Download from Gene ontology website: https://geneontology.org/docs/download-ontology/

## Extract GO mapping information from OBO file.
- Result file: pub-go.obo.tab
```sh
Rscript cnvt_GOobo_to_tab.r --go_obo pub-go.obo  --output_prefix  pub-go.obo
```

## Generate gene to all GO mapping list for whole-genome genes.
```sh
perl extend_GOannot_for_GOenrich.pl  pub-go.obo.tab  in_whole_genome_gene.annot  in_whole_genome_gene.annot-GOinEnrich
```

## Execute GO enrichment. (Not GSEA)
- Result file: in\_geneID-GOen.tsv, in\_geneID-GOen.svg, in\_geneID-GOen-testLg10.svg;
  - testLg10.svg plots enriched GO terms with at least 10 genes present in test gene set (in\_geneID.list).

```sh
Rscript run_GOenrich.r  --gene_list in_geneID.list  --go_annot_4col in_whole_genome_gene.annot-GOinEnrich  --output_prefix in_geneID-GOen

```
