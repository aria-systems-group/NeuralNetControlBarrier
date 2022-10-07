# NeuralNetControlBarrier

Neural Networks (NNs) have been successfully employed to represent the state evolution of complex dynamical systems. Such models, referred to as NN dynamic models (NNDMs), use iterative noisy predictions of NN to estimate a distribution of system trajectories over time. Despite their accuracy, safety analysis of NNDMs is known to be a challenging problem and remains largely unexplored. To address this issue, in this paper, we introduce a method of providing safety guarantees for NNDMs. The paper (Neurips 2022), titled: "Safety Guarantees for Neural Network Dynamic Systems via Stochastic Barrier Functions", can be found on [ArXiv](https://arxiv.org/abs/2206.07811).

## Purpose of this code
This code generates stochastic barrier functions and controllers for NNDMs, in accordance with the paper. \
Five case studies are included: 

1. **Pendulum 2D**
2. **Cartpole 4D** 
3. **Husky 4D**
4. **Husky 5D**
5. **Acrobot 6D**

## Repeat Experiments
| **`Linux`** | **`Mac OS X`** | **`Windows`** |
|-----------------|---------------------|-------------------------|

Read the description below for repeatability of all the experiments.

### Docker Image
The Dockerfile is provided in the main folder. Build this docker file to obtain all the required Julia and Mosek packages, as specified in the Project.toml. To build the docker image, navigate into the main folder and run the following command 
```sh
docker build -t neurips_nn_barrier .
```

To start a container 

```sh
docker run -it --name NNBarrierContainer neurips_nn_barrier
```

**DO NOT BUILD THE DOCKERFILE BEFORE OBTAINING THE REQUIRED MOSEK LICENSE!**


### EXTERNAL: Mosek 
Notice, to run the optimizations, a Mosek license is required! Visit https://www.mosek.com to download this license. After downloading, move `mosek.lic` from the mosek folder into the *licenseFile* folder. 

## Run through bash

Use the following commands to run the optimization case studies through bash.

```sh
runOptimization pendulum   # To run Pendulum
runOptimization cartpole   # To run Cartpole
runOptimization husky4d    # To run Husky 4D
runOptimization husky5d    # To run Husky 5D
runOptimization acrobot    # To run Acrobot

```

## Run through Julia
Use the following commands to run the optimization case studies through Julia

Navigate to *```/NeuralNetControlBarrier.jl```* \
In terminal call julia and run the following commands:
1. ```julia 
      using Pkg
   ```
2. ```julia 
      Pkg.activate(".") 
   ```

To run the Pendulum experiment for example, use the following command: 
```julia 
   include("Optimization/verification_pendulum")
```
The same command can be run for the Cartpole and Husky by changing the system's name accordingly after *verification_*

## Change experiment setup for each case study
Navigate to the *src* folder inside NeuralNetControlBarrier.jl and open *Systems.jl*.  \
For each system, change the number of hypercubes by adjusting `number_hypercubes` on line 16, 21 or 26. \
Table 1 in the paper includes the possible number of hypercubes for each system. 

## Contributing
All contributions welcome! All content in this repository is licensed under the MIT license.

## Citing

If the package NeuralNetControlBarrier.jl is useful in your research, and you would like to acknowledge it, please cite this [paper](https://arxiv.org/abs/2206.07811):

```
@article{mazouz2022nncbf,
  author  = {Rayan Mazouz, Karan Muvvala, Akash Ratheesh, Luca Laurenti and Morteza Lahijanian},
  title   = {Safety Guarantees for Neural Network Dynamic Systems via Stochastic Barrier Functions},
  year    = {2022},
  url     = {https://arxiv.org/abs/2206.07811}
}
