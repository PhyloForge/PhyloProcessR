#!/usr/bin/env bash
set -euo pipefail

: "${PHYLOPROCESSR_BIN:?Set PHYLOPROCESSR_BIN to the Conda bin directory}"
: "${NEW_READS:?Set NEW_READS to the directory containing CRH FASTQ files}"
: "${FULL_TARGETS:?Set FULL_TARGETS to Ranoidea-V1_Markers.fa}"

BIN="${PHYLOPROCESSR_BIN}"
WORK="${PREP_WORK:-/private/tmp/phyloprocessr-full-screen}"
THREADS="${THREADS:-4}"

mkdir -p "${WORK}"
cp "${FULL_TARGETS}" "${WORK}/targets.fa"
if [[ ! -f "${WORK}/targets.fa.bwt" ]]; then
  "${BIN}/bwa" index "${WORK}/targets.fa" >/dev/null 2>&1
fi

printf 'sample\tmarker\tmapped_alignments\n' > "${WORK}/marker_counts.tsv"

for read1 in "${NEW_READS}"/*_R1_001.fastq.gz; do
  read2="${read1/_R1_001.fastq.gz/_R2_001.fastq.gz}"
  sample="$(basename "${read1}" | sed -E 's/_.*$//')"
  "${BIN}/bwa" mem -t "${THREADS}" "${WORK}/targets.fa" "${read1}" "${read2}" \
    2> "${WORK}/${sample}.bwa.log" \
    | "${BIN}/samtools" view -F 2308 - \
    | awk -v sample="${sample}" 'BEGIN {OFS="\t"} {count[$3]++} END {for (marker in count) print sample, marker, count[marker]}' \
    | sort -t $'\t' -k2,2 >> "${WORK}/marker_counts.tsv"
  printf 'screened=%s\n' "${sample}"
done

printf 'screen_complete=%s\n' "${WORK}"
