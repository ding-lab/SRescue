# SRescue

A tool for validating short-read–specific SVs missed by long-read callers via breakpoint force-genotyping in long-read WGS data.

## Details
Long-read SV callers may miss certain somatic SVs with low VAFs due to
relatively lower sequencing coverage. In contrast, short-read SV detection
often suffers from a high false-positive rate caused by mapping ambiguities and
the lack of spanning-read evidence. To retrieve and validate the
short-read–specific SVs that were missed by long-read callers, this tool was
developed. The workflow is summarized as follows:
1.    Identify SVs detected exclusively by short-read callers,
2.    Convert the detected duplications into insertions, because the long-read
      force-genotyping tools typically have limited sensitivity for
      duplications,
3.    Use cuteSV to force-genotype the breakpoints of short-read-specific SVs
      in the long-read normal and tumor data,
4.    Retrieve short-read–specific SVs that have at least 1 supporting read in
      the long-read tumor sample and no supporting reads in the matching normal
      sample,
5.    Verify whether the retrieved duplications (represented as insertions in
      the callsets) are true duplications by comparing the supporting read
      sequences,
6.    Merge the long-read-specific, common, and validated short-read-specific
      SVs into the final SV set.





