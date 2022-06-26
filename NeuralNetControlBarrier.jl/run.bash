#!/bin/bash

set -e


case "$1" in

    pendulum) echo "Running $1"
            julia --project=/NeuralNetControlBarrier.jl Optimization/verification_pendulum.jl
            ;;

    cartpole) echo "Running $1"
            julia --project=/NeuralNetControlBarrier.jl Optimization/verification_cartpole.jl
            ;;

    husky) echo "Running $1"
            julia --project=/NeuralNetControlBarrier.jl Optimization/verification_husky.jl
            ;;
    *)     echo "Invalid flag"
           ;;
esac
