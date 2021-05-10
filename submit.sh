#! /bin/bash

set -e

WD="/lustre/scratch116/tol/teams/lawniczak/data/badass"

snakemake \
    -s ${WD}/badass_smk/Snakefile \
    --profile lsf \
    --use-singularity \
    --use-conda \
    --singularity-args "--bind /lustre:/lustre,/lustre/scratch116/tol/resources/busco/v5/lineages/:/lineages/" \
    -d $WD \
    ${@}