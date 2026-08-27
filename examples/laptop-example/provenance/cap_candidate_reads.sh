#!/usr/bin/env bash
set -euo pipefail

: "${PHYLOPROCESSR_BIN:?Set PHYLOPROCESSR_BIN to the Conda bin directory}"

FASTP="${PHYLOPROCESSR_BIN}/fastp"
WORK="${PREP_WORK:-/private/tmp/phyloprocessr-40-candidate}"
CAP="${READ_PAIR_CAP:-10000}"

mkdir -p "${WORK}/capped-reads" "${WORK}/cap-logs"
printf 'sample\tlane\tread_pairs_retained\n' > "${WORK}/capped_read_summary.tsv"

for read1 in "${WORK}"/reduced-reads/*/*_R1.fastq.gz; do
  read2="${read1/_R1.fastq.gz/_R2.fastq.gz}"
  sample="$(basename "$(dirname "${read1}")")"
  output_dir="${WORK}/capped-reads/${sample}"
  mkdir -p "${output_dir}"
  output1="${output_dir}/${sample}_L001_R1.fastq.gz"
  output2="${output_dir}/${sample}_L001_R2.fastq.gz"

  "${FASTP}" -i "${read1}" -I "${read2}" -o "${output1}" -O "${output2}" \
    --reads_to_process "${CAP}" --disable_adapter_trimming \
    --disable_quality_filtering --disable_length_filtering --thread 2 \
    --json "${WORK}/cap-logs/${sample}.json" \
    --html "${WORK}/cap-logs/${sample}.html" >/dev/null 2>&1

  retained="$(gzip -cd "${output1}" | awk 'END {print NR / 4}')"
  mate2="$(gzip -cd "${output2}" | awk 'END {print NR / 4}')"
  [[ "${retained}" == "${mate2}" ]]
  printf '%s\tL001\t%s\n' "${sample}" "${retained}" >> "${WORK}/capped_read_summary.tsv"
done

gzip -t "${WORK}"/capped-reads/*/*.fastq.gz
du -sh "${WORK}/capped-reads"
