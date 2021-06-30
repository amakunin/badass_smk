#! /bin/bash

set -e

export WD="/lustre/scratch116/tol/projects/badass/users/am60/"
# export PD="/lustre/scratch116/tol/projects/badass/users/am60/"

snakemake \
    -d ${WD} \
    -s ${WD}/badass_smk/Snakefile \
    --profile lsf \
    --use-singularity \
    --singularity-args "--bind /lustre:/lustre,/lustre/scratch116/tol/resources/busco/v5/lineages/:/lineages/" \
    --use-conda \
    ${@}

    # --shadow-prefix ${PD}/.snakemake \
    # --conda-prefix ${PD}/.snakemake \
    # --singularity-prefix ${PD}/.snakemake \
    