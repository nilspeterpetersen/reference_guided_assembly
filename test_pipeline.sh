#!/bin/bash
# Test script to run the assembly pipeline locally

mkdir -p test_data/output

echo "Running with reference AY628201.1.fasta"

docker run \
    --rm \
    -v $(pwd)/test_data:/data rga:latest \
    --reference /data/AY628201.1.fasta \
    --reads /data/reads.AY628201.1.fq \
    --output /data/output/consensus.AY628201.fasta \
    --clair3-model /clair3_models/r1041_e82_400bps_hac_g632

echo "DONE"
echo ""

echo "Running with reference NC_004296.1.fasta"

docker run \
    --rm \
    -v $(pwd)/test_data:/data rga:latest \
    --reference /data/NC_004296.1.fasta \
    --reads /data/reads.AY628201.1.fq \
    --output /data/output/consensus.NC_004296.fasta \
    --clair3-model /clair3_models/r1041_e82_400bps_hac_g632

echo "DONE"
echo ""
