#!/usr/bin/env bash
set -euo pipefail

: "${PHYLOPROCESSR_BIN:?Set PHYLOPROCESSR_BIN to the Conda bin directory}"
: "${NEW_READS:?Set NEW_READS to the directory containing CRH FASTQ files}"
: "${FULL_TARGETS:?Set FULL_TARGETS to Ranoidea-V1_Markers.fa}"
: "${SELECTED_MARKERS:?Set SELECTED_MARKERS to selected_40_markers.txt}"

BIN="${PHYLOPROCESSR_BIN}"
SELECTED="${SELECTED_MARKERS}"
WORK="${PREP_WORK:-/private/tmp/phyloprocessr-40-candidate}"
THREADS="${THREADS:-4}"

mkdir -p "${WORK}/reduced-reads"
cp "${FULL_TARGETS}" "${WORK}/full_targets.fa"
"${BIN}/samtools" faidx "${WORK}/full_targets.fa"
"${BIN}/samtools" faidx -r "${SELECTED}" "${WORK}/full_targets.fa" > "${WORK}/targets.fa"
"${BIN}/bwa" index "${WORK}/targets.fa" >/dev/null 2>&1

printf 'sample\tlane\tread_pairs_retained\n' > "${WORK}/read_extraction_summary.tsv"

for read1 in "${NEW_READS}"/*_R1_001.fastq.gz; do
  read2="${read1/_R1_001.fastq.gz/_R2_001.fastq.gz}"
  sample="$(basename "${read1}" | sed -E 's/_.*$//')"
  sample_dir="${WORK}/reduced-reads/${sample}"
  mkdir -p "${sample_dir}"

  "${BIN}/bwa" mem -t "${THREADS}" "${WORK}/targets.fa" "${read1}" "${read2}" \
    2> "${WORK}/${sample}.bwa.log" \
    | "${BIN}/samtools" view -b -G 12 -F 2304 -o "${WORK}/${sample}.mapped-pairs.bam" -
  "${BIN}/samtools" sort -n -@ "${THREADS}" \
    -o "${WORK}/${sample}.qname.bam" "${WORK}/${sample}.mapped-pairs.bam"
  "${BIN}/samtools" fastq -@ "${THREADS}" -n \
    -1 "${sample_dir}/${sample}_L001_R1.fastq.gz" \
    -2 "${sample_dir}/${sample}_L001_R2.fastq.gz" \
    -0 /dev/null -s /dev/null "${WORK}/${sample}.qname.bam" >/dev/null 2>&1

  retained_pairs="$(gzip -cd "${sample_dir}/${sample}_L001_R1.fastq.gz" | awk 'END {print NR / 4}')"
  printf '%s\tL001\t%s\n' "${sample}" "${retained_pairs}" >> "${WORK}/read_extraction_summary.tsv"
  rm -f "${WORK}/${sample}.mapped-pairs.bam" "${WORK}/${sample}.qname.bam"
  printf 'extracted=%s pairs=%s\n' "${sample}" "${retained_pairs}"
done

gzip -t "${WORK}"/reduced-reads/*/*.fastq.gz
du -sh "${WORK}/reduced-reads"
