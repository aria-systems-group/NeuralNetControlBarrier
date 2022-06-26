FROM julia:1.6

# Install bzip. Required for installing MosekTools package
RUN apt update -y && apt install bzip2 -y

# Copy Julia Control barrier Package to Docker Image
COPY ./NeuralNetControlBarrier.jl /NeuralNetControlBarrier.jl

# Make run.bash executable
RUN chmod +x /NeuralNetControlBarrier.jl/run.bash

# Change the workdir to package root
WORKDIR NeuralNetControlBarrier.jl

# Precompile the Julia Package
RUN julia -e 'using Pkg;Pkg.activate("."); Pkg.instantiate(); Pkg.precompile()'

# Add Alias to run experiment 
RUN echo 'alias runOptimization="/NeuralNetControlBarrier.jl/run.bash"' >> ~/.bashrc

# Add Mosek License file to Mosek Lib Dir
COPY licenseFile/mosek.lic /mosek/mosek.lic

# Add License Env variable
RUN echo 'export MOSEKLM_LICENSE_FILE="/mosek/mosek.lic"' >> ~/.bashrc

# Change Entrypoint to bash (Default: julia)
ENTRYPOINT ["bash"]






