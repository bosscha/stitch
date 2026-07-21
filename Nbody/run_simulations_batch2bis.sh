#!/bin/bash

echo "Starting N-body simulations..."

echo "Running 5D simulation..."
./nbody-rust-5D

echo "Running 6D simulation..."
./nbody-rust-6D

echo "All simulations finished."
