#!/bin/bash 

docker buildx build \
    --target ttk \
    --cache-from ghcr.io/scivislab/ttk:buildcache \
    -o type=oci,dest=output.tar \
    -f scripts/docker/Dockerfile \
    .

singularity build -F /tmp/ttk-private.sif oci-archive://output.tar && rm output.tar



