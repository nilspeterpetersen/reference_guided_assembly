FROM mambaorg/micromamba:1.4.5

WORKDIR /app

# Copy environment specification and create micromamba environment
COPY env/env.yaml /tmp/env.yaml
# create environment 'rga' from env.yaml
RUN micromamba create -y -n rga -f /tmp/env.yaml && micromamba clean --all --yes

# Ensure the created environment's bin is on PATH
ENV PATH=/opt/conda/envs/rga/bin:$PATH

# Copy the script
COPY assembly.py /app/assembly.py

# Make /data a volume for inputs/outputs
VOLUME ["/data"]

# Default entrypoint: run the assembly script using the micromamba environment
ENTRYPOINT ["micromamba", "run", "-n", "rga", "--no-capture-output", "python", "/app/assembly.py"]
