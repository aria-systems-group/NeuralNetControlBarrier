# Dynamical Neural Network Verification using Control Barrier Functions

# Initialization Statement
print("\n Computing Control Barrier Certificate & Controller based on Neural Network Bounds ")

# Import Module 
using NeuralNetControlBarrier
using Optim

# System definition
system_flag =  "pendulum"
input_data = inputs(system_flag)

# Optimize certificate
@time certificate, eta_certificate, certificate_beta_vals, system_dimension = optimization(input_data::Tuple)

# Obtain certificate argmin
h = Meta.parse(string(certificate))
res = optimize((@eval x -> $h),  zeros(system_dimension))
x_star = Optim.minimizer(res)

# Optimize controller
@time certificate, controllers, counts = control_loop(input_data::Tuple, certificate, eta_certificate, x_star, true, certificate_beta_vals)

 