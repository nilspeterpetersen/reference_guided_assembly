FROM mambaorg/micromamba:2.5.0

WORKDIR /app

ENV MAMBA_NO_LOW_SPEED_LIMIT=1 \
    MAMBA_RETRIES=20 \
    MAMBA_RETRY_SLEEP=5 \
    MAMBA_NUM_THREADS=1

# Copy environment specifications
COPY env/env.yaml /tmp/env.yaml

# Create Conda environment
RUN micromamba create -y -vv -f /tmp/env.yaml && \
    micromamba clean --all --yes

# Select the environment
ENV ENV_NAME=env

# Copy Python package
COPY --chown=mambauser:mambauser pyproject.toml README.md /app/
COPY --chown=mambauser:mambauser src/ /app/src/

# Install Python package into the Conda environment
RUN micromamba run -n env \
    python -m pip install --no-deps --no-build-isolation .

# Copy Clair3 models
COPY clair_models/clair3_models/ /clair3_models/

# Make /data a volume for inputs/outputs
VOLUME ["/data"]

# Run the installed CLI
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
