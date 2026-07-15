# `evolution_tools/mummer_tools` — MUMmer/nucmer alignment helpers

Small helpers around a **MUMmer (nucmer)** whole-genome alignment: chain
neighbouring alignment blocks and plot dot-plots.

External programs: **MUMmer** (`nucmer`, `show-coords`, `delta-filter`). The
`mummerplot` here is MUMmer's own gnuplot-driver script, vendored so it is on
hand; it still needs `gnuplot`.

## Scripts

| script | role |
|--------|------|
| `join_coords.pl` | chain neighbouring `show-coords` alignment blocks that are close in **both** Ref and Qry into longer syntenic segments |
| `mummerplot` | vendored MUMmer dot-plot script (nucmer `.delta`/`.coords` -> gnuplot) |

## Typical use

Run nucmer, tabulate coordinates, then chain nearby blocks:
```sh
nucmer -p ref_qry ref.fa qry.fa
delta-filter -1 ref_qry.delta > ref_qry.1delta          # 1-to-1 optional
show-coords -rclT ref_qry.1delta > ref_qry.coords       # -T = tab, with headers
perl join_coords.pl ref_qry.coords > ref_qry.coords.jn
```

`join_coords.pl` options:
- `-maxDist1 [40e3]` / `-maxDist2 [40e3]` — max gap between two neighbouring
  blocks to still join them, in Ref (`1`) / Qry (`2`) coordinates.
- `-maxOvl1 [0]` / `-maxOvl2 [0]` — max allowed overlap of two neighbouring
  blocks in Ref / Qry.
- `-inType [coords]` — input format; also `coordsTab` or already-`joined`.

Dot-plot the alignment:
```sh
./mummerplot --png -p ref_qry ref_qry.1delta        # needs gnuplot
```

---

_Hand-written. `mummerplot` is a vendored 3rd-party MUMmer script; `join_coords.pl`
is the only original tool here._
