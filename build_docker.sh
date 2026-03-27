#!/usr/bin/env bash
# simple wrapper to build the Docker image for the reference-guided assembly

set -euo pipefail
set -x

IMAGE_NAME="rga:latest"

echo "Building Docker image ${IMAGE_NAME}..."
docker build --network=host --no-cache -t "${IMAGE_NAME}" .

echo "Build complete. You can run the image with:"
echo "  docker run --rm -v \$(pwd):/data ${IMAGE_NAME} --reference /data/ref.fasta --reads /data/reads.fastq --output /data/consensus.fasta"
