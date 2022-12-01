#!/bin/bash 

docker buildx build \
    --target ttk \
    --cache-from ghcr.io/scivislab/ttk:buildcache \
    -f scripts/docker/Dockerfile \
    --load \
    .
