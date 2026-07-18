#!/bin/bash

echo "Starting N-body simulations..."

echo "Running 25D simulation..."
./nbody-rust-25D

echo "Running 50D simulation..."
./nbody-rust-50D

echo "Running 100D simulation..."
./nbody-rust-100D

echo "All simulations finished."
