#!/bin/bash
# Kills all instances of build_hp.jl, including suspended/stopped background tasks.

echo "Killing all build_hp.jl processes..."
pkill -9 -f "build_hp.jl"
echo "Done."
