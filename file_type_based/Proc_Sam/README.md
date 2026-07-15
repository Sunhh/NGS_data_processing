# `Proc_Sam` — SAM/BAM record processing

| script | role |
|--------|------|
| `get_required_sam.pl` | select/filter SAM records by criteria; **input must be position-sorted** (streams to keep memory low) |
| `trim_rdEnd_inSam.pl` | trim read ends within a SAM (note: mate-mapping fields become inconsistent afterwards) |

_Hand-written._
