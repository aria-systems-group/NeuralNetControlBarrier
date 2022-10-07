#!/bin/bash

set -e


case "$1" in

    pendulum) echo "Running $1"
            julia --project=/NeuralNetControlBarrier.jl Optimization/verification_pendulum.jl
            ;;

    cartpole) echo "Running $1"
            julia --project=/NeuralNetControlBarrier.jl Optimization/verification_cartpole.jl
            ;;

    husky4d) echo "Running $1"
            julia --project=/NeuralNetControlBarrier.jl Optimization/verification_husky4d.jl
            ;;

    husky5d) echo "Running $1"
        julia --project=/NeuralNetControlBarrier.jl Optimization/verification_husky5d.jl
        ;;

    acrobot) echo "Running $1"
        julia --project=/NeuralNetControlBarrier.jl Optimization/verification_acrobot.jl
        ;;
    *)     echo "Invalid flag"
           ;;
esac
