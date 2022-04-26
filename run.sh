#!/usr/bin/env bash

# load conda env where nextflow is installed
eval "$(conda shell.bash hook)"
conda activate nextflow

nextflow run nf_deepvariant/pipelines/main.nf \
	 -c nf_deepvariant/nextflow.config \
	 -c inputs.nf \
	 -profile lsf \
	 --nf_ci_loc $PWD \
	 -resume

# temporary cache dir used by singularity when pulling images from dockerhub
# this is different than the nextflow singularity cache dir when it caches images (defined in nf conf)
# export singularity_cachedir="singularity_cache"
# mkdir -p $singularity_cachedir

# tmp dir used by singularity when pulling images from dockerhub
# export tmpdir="tmpdir"
# mkdir -p $tmpdir
