# reference_guided_assembly — Docker usage

Build the Docker image (from repository root):

```bash
docker build -t rga:latest .
```

Run the container mounting a folder with your input files (example):

```bash
docker run --rm -v $(pwd):/data rga:latest \
  --reference /data/ref.fasta --reads /data/reads.fastq --output /data/consensus.fasta
```

Notes:
- The Docker image creates a micromamba environment named `rga` from `env/env.yaml`.
- The container entrypoint uses `micromamba run -n rga python /app/assembly.py`, so any
  arguments passed to `docker run` are forwarded to `assembly.py`.
- Building the image downloads and installs conda packages; the build can take several minutes.
- If you prefer a different base image or want the image to be slimmer, I can switch to a multi-stage build or reduce packages.
