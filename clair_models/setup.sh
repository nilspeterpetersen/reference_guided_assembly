#!/usr/bin/env bash

# download rerio models
wget -r -np -nH --cut-dirs=2 \
    -R 'index.html*' \
    -P clair_models/ \
    https://www.bio8.cs.hku.hk/clair3/clair3_models_rerio_pytorch/
