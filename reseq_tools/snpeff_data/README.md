# snpeff_data — post-process SnpEff SV / InDel effect annotation

Turn **SnpEff** structural-variant / InDel annotation into a tidy "one variant × one gene
per line" table with **simplified effect categories**, then pull out a category of interest.

Upstream (external): annotate the VCF with SnpEff, then flatten the `ANN` fields with SnpSift:
```sh
snpsift extractFields -s ',' -e '.' snpeff-ann.vcf \
  CHROM POS ID SVLEN ALT "ANN[*].ALLELE" "ANN[*].FEATUREID" "ANN[*].EFFECT" "ANN[*].DISTANCE" \
  > allInDel-ann_ud2k.tbl1
```

| file | role |
|------|------|
| `simple_sv_class` | lookup table: each SnpEff effect type (SO term) → simplified class (`simple_1..4`: whole/coding/intron/promoter/…) + priority `order_*` + description |
| `simplify_tbl1.pl` | `.tbl1` → melt to one variant × one gene per line and attach the simplified class (uses `simple_sv_class` + a `geneID↔mrnaID` map) |
| `extract_coding.pl` | keep lines of one merged category: `coding` (whole+coding), `intron` (splice+intron), `promoter` (promoter+upstream), `remnant` (mapgene), `downstream` |

---

## Workflow
```sh
# 1. flatten SnpEff ANN fields  (SnpSift, above) -> allInDel-ann_ud2k.tbl1

# 2. melt + classify: one variant vs one gene per line, with simplified effect class
perl simplify_tbl1.pl map.geneID_mrnaID allInDel-ann_ud2k.tbl1 > allInDel-ann_ud2k.tbl2
#    map.geneID_mrnaID: two columns from the SnpEff build (e.g. snpEff/data/<db>/map.geneID_mrnaID)

# 3. pull one category (coding / intron / promoter / remnant / downstream)
perl extract_coding.pl coding allInDel-ann_ud2k.tbl2 > allInDel-ann_ud2k.coding.tbl3
```
`simple_sv_class` is the single place to adjust how SnpEff effect terms map to the merged
categories (edit its `simple_*` / `order_*` columns).

_Hand-written — `gen_script_index.sh` leaves it alone._
