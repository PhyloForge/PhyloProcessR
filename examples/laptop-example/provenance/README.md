# Fixture provenance

The seven-sample fixture was derived without altering the source reads, the
full `Ranoidea-V1_Markers.fa` reference, or the original two-sample laptop
fixture.

The original two divergent samples were selected by searching their draft
assemblies against all 13,515 Ranoidea markers with TBLASTX. The strongest
shared gene-unique markers were used to extract lane-specific paired reads
through empirical contigs. Those source assemblies and full preparation
outputs are not required to run this distributed example.

The five CRH libraries were mapped directly against the full marker reference.
`screen_full_ranoidea.sh` regenerates the primary mapped-alignment counts; the
1.85 MB full count table is omitted from the repository copy because it is not
needed to run or verify the example. `ranked_gene_unique_markers.tsv` contains
markers present in all five samples, ranked by the minimum sample count and
then total support. One marker per gene was retained.
`selected_40_markers.txt` records the 40-candidate panel; it includes the two
loci already supported by both original samples.

The maintainer-only scripts `screen_full_ranoidea.sh`,
`select_candidate_markers.R`, `extract_candidate_reads.sh`, and
`cap_candidate_reads.sh` preserve the screening, ranking, paired-read
extraction, and read-cap operations. Their required path variables are checked
at startup; optional `PREP_WORK`, `THREADS`, and `READ_PAIR_CAP` values control
temporary output and resource use.

For each CRH sample, both mates were retained whenever either mate mapped to
the candidate panel. Each paired file was then capped together at 10,000 read
pairs. Lane identifiers remain in all filenames. `read_extraction_summary.tsv`
records the distributed fixture counts, and `input_sha256.tsv` records its
checksums.

Candidate choice used read mapping only as a screening statistic. Final locus
selection was based on completed PhyloProcessR annotation and MAFFT alignment:
taxon occupancy was ranked first, followed by lower mean missingness, longer
alignment length, and locus name. All 38 passing alignments are preserved, and
the ranked best 20 are copied into the review set.
