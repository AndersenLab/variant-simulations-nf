#!/bin/bash
# Author: Daniel E. Cook
# Fetches test data for variant simulations
set -e

# Modify this variable to download the appropriate reference genome
REFERENCE="c_elegans/PRJNA13758/WS276"

error() {
    >&2 echo -e "\n\t$(tput setaf 1)$*@*$(tput sgr0)\n"
    exit 1
}

msg() {
    >&2 echo -e "$(tput setaf 3)$*$(tput sgr0)"
}

# Copy reference genome from QUEST
quest_host=$(grep -o 'Host quest' ~/.ssh/config | wc -l)
if [[ ! ${quest_host} -eq 1 ]]; then
    error 'You do not have a quest host configured in your profile.'
fi;

# Fetch data
PARENT_DIR=$(git rev-parse --show-toplevel)

msg "Creating dirs"
VCF_DIR="${PARENT_DIR}/_assets/vcf"
BAM_DIR="${PARENT_DIR}/_assets/bam"
GENOMES="${PARENT_DIR}/_assets/genomes"
mkdir -p "${VCF_DIR}"
mkdir -p "${BAM_DIR}"
mkdir -p "${GENOMES}"

msg "Downloading reference genome: ${REFERENCE}"
mkdir -p "${PARENT_DIR}/_assets/genomes/${REFERENCE}"
rsync -rauL --progress --copy-links "quest:/projects/b1059/data/genomes/${REFERENCE}/*" "${GENOMES}/${REFERENCE}"

msg "Downloading N2 Reference"
N2_BAM="${BAM_DIR}/N2.bam"
N2_05_BAM="${BAM_DIR}/N2.s05.bam"
[ ! -s "$N2_BAM" ]; wget -O "${N2_BAM}" "https://s3.us-east-2.amazonaws.com/elegansvariation.org/bam/strain/N2.bam"

# Subsample N2 bam for testing
[ ! -s "${N2_05_BAM}" ]; samtools view -s 0.10 -b "${N2_BAM}" > "${N2_05_BAM}" && samtools index "${N2_05_BAM}"


msg "Download recent WI release"
WI_FNAME="WI.20180527.soft-filter.vcf.gz"
WI_RELEASE="${VCF_DIR}/${WI_FNAME}"
[ ! -s "$WI_RELEASE" ]; wget -O "${WI_RELEASE}" "https://storage.googleapis.com/elegansvariation.org/releases/20180527/variation/${WI_FNAME}"
bcftools index ${WI_RELEASE}