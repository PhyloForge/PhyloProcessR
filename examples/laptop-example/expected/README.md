# Expected results

These reference outputs were generated with two CPU threads and a 4 GB SPAdes
memory limit. Runtime is not an acceptance criterion.

A successful preprocessing and assembly run should:

- recognize seven biological samples across 11 sequencing lanes;
- process 10,000 pairs for each CRH sample, 714 pairs for `Test_dataset_1`,
  and 1,424 pairs for `Test_dataset_2`;
- retain 8,777-8,911 cleaned pairs per CRH sample and 572 and 1,140 pairs for
  the original samples;
- recover at least four samples for 39 of the 40 base candidate loci;
- recover `ranoidea-01644_kif2a-ex2` and
  `ranoidea-05399_ckap5-ex3` from all seven samples.

A successful annotation and alignment run should produce 38 candidate
alignments and copy the ranked best 20 to `alignments/best-20/`. Selected
alignments contain five to seven taxa and span 983-12,146 bp in this reference
run. See `best_20_alignment_summary.tsv` for the exact ranking.

The trimming/QC stage should retain all 20 selected alignments without
removing taxa. Final alignments contain five to seven taxa, span 205-11,589 bp,
and contain no more than 13.13% gaps in this reference run. The five CRH
samples occur in all 20 final alignments; both original samples occur in the
two broadly shared loci. See `final_20_qc_summary.tsv` and
`final_20_sample_occupancy.tsv` for exact values.

Minor sequence and count differences can occur across SPAdes, MAFFT, or other
dependency versions. The tables and FASTA/PHYLIP files here are review
references rather than byte-for-byte portability requirements.
