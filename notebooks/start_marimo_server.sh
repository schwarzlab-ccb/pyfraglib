#!/bin/sh
# From e.g. a conda environment with `marimo' and all dependencies installed,
# run this script to start a server. On the remote client, SSH forwarding
# needs to be performed using:
# ``$ ssh -nNT -L localhost:8787:$(hostname -s):8787 $USER@ramses1.itcc.uni-koeln.de''
marimo edit --no-token --headless --port 8787 --host ramses15229 wpr_predictor.py
