FROM mambaorg/micromamba:2.5.0

WORKDIR /app

ENV MAMBA_NO_LOW_SPEED_LIMIT=1 \
    MAMBA_RETRIES=20 \
    MAMBA_RETRY_SLEEP=5 \
    MAMBA_NUM_THREADS=1

# Copy environment specifications
COPY env/env.yaml /tmp/env.yaml

# create environment env from env.yaml
RUN micromamba create -y -vv -f /tmp/env.yaml && micromamba clean --all --yes

# Ensure the created environment's bin is on PATH
#ENV PATH=/opt/conda/envs/rga/bin:$PATH
ENV ENV_NAME=env

COPY clair_models/clair3_models/ /clair3_models/

# Copy the script
COPY assembly.py /app/assembly.py

# Make /data a volume for inputs/outputs
VOLUME ["/data"]

# Default entrypoint: run the assembly script using the micromamba environment
#ENTRYPOINT ["micromamba", "run", "-n", "env", "--no-capture-output", "python", "/app/assembly.py"]
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "python", "/app/assembly.py"]
