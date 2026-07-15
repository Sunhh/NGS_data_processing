# `file_type_based` — utilities organised by file type

Helpers grouped by the file type they operate on, rather than by biological task.

| dir | operates on | purpose |
|-----|-------------|---------|
| [`Proc_Reads/`](Proc_Reads/) | FASTQ reads | QC, adapter/quality trimming, rRNA & poly-A/T removal, barcode search, and thin mapping wrappers (bwa/bowtie/bowtie2/TopHat2) |
| [`Proc_Sam/`](Proc_Sam/) | SAM/BAM | filter/select records and trim read ends inside SAM |

_Hand-written navigation index._
