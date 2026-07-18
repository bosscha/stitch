#!/bin/bash

echo "Starting N-body simulations..."

echo "Running 1D simulation..."
./nbody-rust-1D

echo "Running 2D simulation..."
./nbody-rust-2D

echo "Running 3D simulation..."
./nbody-rust-3D

echo "All simulations finished."
