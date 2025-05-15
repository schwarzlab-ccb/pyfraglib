#!/bin/sh
# From e.g. a conda environment with `marimo' and all dependencies installed,
# run this script to start a server. On the remote client, SSH forwarding
# needs to be performed as indicated.
SCRIPT="$1"
RAMSES_PORT="${2:-8973}"

echo "Perform SSH forwarding as follows:"
echo "ssh -nNT -L localhost:8787:$(hostname -s):${RAMSES_PORT} ${USER}@ramses1.itcc.uni-koeln.de"
marimo edit --no-token --headless --port "${RAMSES_PORT}" --host ramses15229 "${SCRIPT}"
