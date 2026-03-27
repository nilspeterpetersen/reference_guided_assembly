FROM mambaorg/micromamba:1.4.5

WORKDIR /app

ENV MAMBA_NO_LOW_SPEED_LIMIT=1 \
    MAMBA_RETRIES=20 \
    MAMBA_RETRY_SLEEP=5 \
    MAMBA_NUM_THREADS=1

# Copy environment specifications
COPY env/env.yaml /tmp/env.yaml

# create environment env from env.yaml
RUN micromamba create -y -f /tmp/env.yaml && micromamba clean --all --yes

# Ensure the created environment's bin is on PATH
ENV PATH=/opt/conda/envs/rga/bin:$PATH

COPY clair_models/rerio/clair3_models/ /clair3_models/

## clone rerio repository
#RUN git clone https://github.com/nanoporetech/rerio.git

## download rerio models
#RUN python3 rerio/download_model.py --clair3

# Copy the script
#COPY assembly.py /app/assembly.py

# Make /data a volume for inputs/outputs
#VOLUME ["/data"]

# Default entrypoint: run the assembly script using the micromamba environment
#ENTRYPOINT ["micromamba", "run", "-n", "rga", "--no-capture-output", "python", "/app/assembly.py"]
