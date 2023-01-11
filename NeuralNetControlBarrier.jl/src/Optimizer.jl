# Sum of squares optimization function
function optimization(input_data::Tuple)

    # Extract input data
    number_hypercubes = input_data[1]
    barrier_degree_input = input_data[2]
    safety_threshold_eta = input_data[3]
    safety_threshold_beta = input_data[4]
    system_flag = input_data[5]
    neural_network_bound = input_data[6]
    layer_flag = input_data[7]
    beta_partition = input_data[8]
    large_range_initial = input_data[9]
    initial_set_radius = input_data[10]
    print_to_txt = input_data[11]
    decision_eta_flag = input_data[12]

    # File reading
    filename = "/models/" * system_flag * "/" * neural_network_bound  * "/" * layer_flag* "_layers/partition_data_"  * string(number_hypercubes) * ".mat"
    file = matopen(pwd()*filename)

    # Extract hypercube data (avoid using float64 for precision issues)
    partitions = read(file, "partitions")
    state_space = read(file, "state_space")

    # Number of hypercubes
    number_hypercubes_constraints = number_hypercubes
    hcube_identifier = 1:number_hypercubes
    
    # Extract Neural Network Bounds
    M_h = read(file, "M_h")
    M_l = read(file, "M_l")
    B_h = read(file, "B_h")
    B_l = read(file, "B_l")

    # Define system and control dimensions
    system_dimension::Int64 = Integer(length(state_space[:,1]))

    # Using Mosek as the SDP solver
    model = SOSModel(optimizer_with_attributes(Mosek.Optimizer,
                                               "MSK_DPAR_INTPNT_TOL_STEP_SIZE" => 1e-6,
                                               "MSK_IPAR_OPTIMIZER" => 0,
                                               "MSK_IPAR_BI_CLEAN_OPTIMIZER" => 0,
                                               "MSK_IPAR_NUM_THREADS" => 16,
                                               "MSK_IPAR_PRESOLVE_USE" => 0))

    # Create state space variables
    @polyvar x[1:system_dimension]

    # Create noise variable
    @polyvar z

    # Create CROWN bounds variables
    @polyvar y[1:system_dimension]

    # Create dummy variable for beta in SOS
    @polyvar w[1:2]

    # Create probability decision variables eta  and beta
    @variable(model, eta)
    if beta_partition == true
        @variable(model, beta_parts_var[1:number_hypercubes_constraints])
        @variable(model, beta)
    else
        @variable(model, beta)
    end

    # Create barrier polynomial, specify degree Lagrangian polynomials
    barrier_degree::Int64 = barrier_degree_input
    alpha::Float64 = 1
    lagrange_degree::Int64 = 2

    # Specify noise element (Gaussian)
    standard_deviation::Float64 = 0.1^2

    # Specify initial state
    x_init::Array{Float64, 2} = zeros(1, system_dimension)

    # Create barrier candidate
    barrier_monomial::MonomialVector{true} = monomials(x, 0:barrier_degree)
    @variable(model, c[1:Integer(length(barrier_monomial))])
    BARRIER::DynamicPolynomials.Polynomial{true, AffExpr} = barrier_polynomial(c, barrier_monomial)

    # Add constraints to model for positive barrier, eta and beta
    add_constraint_to_model(model, BARRIER)
    if decision_eta_flag == true
        @constraint(model, eta >= 1e-6)
        @constraint(model, eta <= safety_threshold_eta)
    else
        eta = safety_threshold_beta
    end
    if beta_partition == true
        for betas = 1:number_hypercubes_constraints
            @constraint(model, beta_parts_var[betas] >= 1e-6)
        end
        @constraint(model, beta >= 1e-6)
    else
        @constraint(model, beta >= 1e-6)
    end
    
    # One initial condition and unsafe conditions
    if system_flag == "cartpole" && large_range_initial == true 
        number_decision_vars = ((system_dimension)^2 + 1 + 1)*length(barrier_monomial)
    elseif system_flag == "husky4d" && large_range_initial == true 
        number_decision_vars = ((system_dimension)^2 + 1 + 1 + 1)*length(barrier_monomial)
    elseif system_flag == "husky5d" && large_range_initial == true 
        number_decision_vars = ((system_dimension)^2 + 1 + 1 + 1)*length(barrier_monomial)
    elseif system_flag == "acrobot" && large_range_initial == true 
        number_decision_vars = ((system_dimension)^2 + 1 + 1 + 1 + 1)*length(barrier_monomial)
    elseif system_flag == "pendulum" && large_range_initial == true 
        number_decision_vars = ((system_dimension)^2 + 1 + 1)*length(barrier_monomial)
    end
    @variable(model, l[1:number_decision_vars])

    barrier_constraints_unsafe_initial = system_dimension + 1

    for ii = 1:barrier_constraints_unsafe_initial

        # Barrier initial condition f(eta)
        if ii == barrier_constraints_unsafe_initial

            # Generate sos polynomial
            count_lag = 0

            # Optimize this code: not all x variables needed in lag_poly_i, lag_poly_theta (only 1 2 4 and 3, respectively)
            lag_poly_i::DynamicPolynomials.Polynomial{true, AffExpr} =  sos_polynomial(l::Vector{VariableRef}, x::Array{PolyVar{true},1}, count_lag::Int64, lagrange_degree::Int64)
            add_constraint_to_model(model, lag_poly_i)
            
            # Change the radius ball of theta
            if system_flag == "cartpole" && large_range_initial == true 
                lag_poly_theta::DynamicPolynomials.Polynomial{true, AffExpr} =  sos_polynomial(l::Vector{VariableRef}, x::Array{PolyVar{true},1}, (count_lag+1)::Int64, lagrange_degree::Int64)
                add_constraint_to_model(model, lag_poly_theta)
            elseif (system_flag == "pendulum" && large_range_initial == true) || (system_flag == "pendulum3d" && large_range_initial == true )
                lag_poly_theta_pen::DynamicPolynomials.Polynomial{true, AffExpr} =  sos_polynomial(l::Vector{VariableRef}, x::Array{PolyVar{true},1}, (count_lag+1)::Int64, lagrange_degree::Int64)
                add_constraint_to_model(model, lag_poly_theta_pen)
            elseif  (system_flag == "husky4d" || system_flag == "husky5d") && large_range_initial == true
                lag_poly_x1 =  sos_polynomial(l::Vector{VariableRef}, x::Array{PolyVar{true},1}, (count_lag+1)::Int64, lagrange_degree::Int64)
                lag_poly_x2 =  sos_polynomial(l::Vector{VariableRef}, x::Array{PolyVar{true},1}, (count_lag+2)::Int64, lagrange_degree::Int64)
            elseif  (system_flag == "acrobot") && large_range_initial == true
                lag_poly_x1 =  sos_polynomial(l::Vector{VariableRef}, x::Array{PolyVar{true},1}, (count_lag+1)::Int64, lagrange_degree::Int64)
                lag_poly_x2 =  sos_polynomial(l::Vector{VariableRef}, x::Array{PolyVar{true},1}, (count_lag+2)::Int64, lagrange_degree::Int64)
            end

            # Initial condition radius and ball
            x_initial_radius = (1e-8)
            x_initial_sums = x_initial_radius

            if system_flag == "cartpole" && large_range_initial == true 
                theta_radius = initial_set_radius^2
                theta_initial_sums = x_initial_radius

                for jj = 1:length(x)
                    if jj != 3
                        x_initial_sums += -(x[jj] - x_init[jj])^2
                    else 
                        theta_initial_sums += -(x[jj] - x_init[jj])^2
                    end
                end
            
            elseif system_flag == "pendulum" && large_range_initial == true 
                theta_radius = initial_set_radius^2
                theta_initial_sums = x_initial_radius

                for jj = 1:length(x)
                    if jj != 2
                        x_initial_sums += -(x[jj] - x_init[jj])^2
                    else 
                        theta_initial_sums += -(x[jj] - x_init[jj])^2
                    end
                end

            elseif (system_flag == "husky4d" || system_flag == "husky5d") && large_range_initial == true 
                x1_initial_radius = 0.1^2
                x2_initial_radius = 0.1^2

                for jj = 1:length(x)
                    if jj == 1
                        x1_initial_radius += -(x[jj] - x_init[jj])^2
                    elseif jj == 2
                        x2_initial_radius += -(x[jj] - x_init[jj])^2    
                    else 
                        x_initial_sums += -(x[jj] - x_init[jj])^2
                    end
                end
                
            elseif (system_flag == "acrobot") && large_range_initial == true 
                x1_initial_radius = 0.1^2
                x2_initial_radius = 0.1^2

                for jj = 1:length(x)
                    if jj == 1
                        x1_initial_radius += -(x[jj] - x_init[jj])^2
                    elseif jj == 2
                        x2_initial_radius += -(x[jj] - x_init[jj])^2    
                    else 
                        x_initial_sums += -(x[jj] - x_init[jj])^2
                    end
                end    

            else
                for jj = 1:length(x)
                    x_initial_sums += -(x[jj] - x_init[jj])^2
                end
            end

            # Barrier constraint eta
            if system_flag == "cartpole" && large_range_initial == true
                _barrier_initial = - BARRIER + eta - lag_poly_i * x_initial_sums - lag_poly_theta*theta_initial_sums
            elseif (system_flag == "husky4d" || system_flag == "husky5d") && large_range_initial == true 
                _barrier_initial = - BARRIER + eta - lag_poly_i * x_initial_sums - lag_poly_x1*x1_initial_radius - lag_poly_x2*x2_initial_radius
            elseif (system_flag == "acrobot") && large_range_initial == true 
                _barrier_initial = - BARRIER + eta - lag_poly_i * x_initial_sums - lag_poly_x1*x1_initial_radius - lag_poly_x2*x2_initial_radius
            elseif system_flag == "pendulum" && large_range_initial == true
                _barrier_initial = - BARRIER + eta - lag_poly_i * x_initial_sums - lag_poly_theta_pen*theta_initial_sums
            end

            # Add constraint to model
            add_constraint_to_model(model, _barrier_initial)

        # Barrier unsafe region conditions (cartpole)
        elseif system_flag == "cartpole" 
            
            if ii == 1 || ii == 3  

                # Generate sos polynomials
                count_lag = 2*ii
                lag_poly_i_lower =  sos_polynomial(l::Vector{VariableRef}, x::Array{PolyVar{true},1}, (count_lag - 1)::Int64, lagrange_degree::Int64)
                lag_poly_i_upper =  sos_polynomial(l::Vector{VariableRef}, x::Array{PolyVar{true},1}, (count_lag)::Int64, lagrange_degree::Int64)

                # State space ranges
                x_i_lower = state_space[ii, 1]
                x_i_upper = state_space[ii, 2]

                # Specify constraints for initial and unsafe set
                _barrier_unsafe_lower = BARRIER - lag_poly_i_lower * (x_i_lower - x[ii]) - 1
                _barrier_unsafe_upper = BARRIER - lag_poly_i_upper * (x[ii] - x_i_upper) - 1

                # Add constraints to model
                add_constraint_to_model(model, lag_poly_i_lower)
                add_constraint_to_model(model, lag_poly_i_upper)
                add_constraint_to_model(model, _barrier_unsafe_lower)
                add_constraint_to_model(model, _barrier_unsafe_upper)
            end

        # Barrier unsafe region conditions (husky)
        elseif (system_flag == "husky4d" || system_flag == "husky5d")

            if ii == 1 || ii == 2

                # Generate sos polynomials
                if large_range_initial == true 
                    count_lag = 2*ii + 1
                else
                    count_lag = 2*ii
                end

                lag_poly_i_lower =  sos_polynomial(l::Vector{VariableRef}, x::Array{PolyVar{true},1}, (count_lag - 1)::Int64, lagrange_degree::Int64)
                lag_poly_i_upper =  sos_polynomial(l::Vector{VariableRef}, x::Array{PolyVar{true},1}, (count_lag)::Int64, lagrange_degree::Int64)

                # State space ranges
                x_i_lower = state_space[ii, 1]
                x_i_upper = state_space[ii, 2]

                # Specify constraints for initial and unsafe set
                _barrier_unsafe_lower = BARRIER - lag_poly_i_lower * (x_i_lower - x[ii]) - 1
                _barrier_unsafe_upper = BARRIER - lag_poly_i_upper * (x[ii] - x_i_upper) - 1

                # Add constraints to model
                add_constraint_to_model(model, lag_poly_i_lower)
                add_constraint_to_model(model, lag_poly_i_upper)
                add_constraint_to_model(model, _barrier_unsafe_lower)
                add_constraint_to_model(model, _barrier_unsafe_upper)

            end

        # Barrier unsafe region conditions (acrobot)
        elseif (system_flag == "acrobot")

            if ii == 2 || ii == 4

                # Generate sos polynomials
                if large_range_initial == true 
                    count_lag = 2*ii + 1
                else
                    count_lag = 2*ii
                end

                lag_poly_i_lower =  sos_polynomial(l::Vector{VariableRef}, x::Array{PolyVar{true},1}, (count_lag - 1)::Int64, lagrange_degree::Int64)
                lag_poly_i_upper =  sos_polynomial(l::Vector{VariableRef}, x::Array{PolyVar{true},1}, (count_lag)::Int64, lagrange_degree::Int64)

                # State space ranges
                x_i_lower = -0.6
                x_i_upper = 0.6

                # Specify constraints for initial and unsafe set
                _barrier_unsafe_lower = BARRIER - lag_poly_i_lower * (x_i_lower - x[ii]) - 1
                _barrier_unsafe_upper = BARRIER - lag_poly_i_upper * (x[ii] - x_i_upper) - 1

                # Add constraints to model
                add_constraint_to_model(model, lag_poly_i_lower)
                add_constraint_to_model(model, lag_poly_i_upper)
                add_constraint_to_model(model, _barrier_unsafe_lower)
                add_constraint_to_model(model, _barrier_unsafe_upper)

            end

        # Barrier unsafe region conditions (pendulum)
        elseif system_flag == "pendulum" 
    
            if ii == 1

                # Generate sos polynomials
                    count_lag = 2*ii
                    lag_poly_i_lower =  sos_polynomial(l::Vector{VariableRef}, x::Array{PolyVar{true},1}, (count_lag - 1)::Int64, lagrange_degree::Int64)
                    lag_poly_i_upper =  sos_polynomial(l::Vector{VariableRef}, x::Array{PolyVar{true},1}, (count_lag)::Int64, lagrange_degree::Int64)
        
                    # State space ranges
                    x_i_lower = state_space[ii, 1]
                    x_i_upper = state_space[ii, 2]
        
                    # Specify constraints for initial and unsafe set
                    _barrier_unsafe_lower = BARRIER - lag_poly_i_lower * (x_i_lower - x[ii]) - 1
                    _barrier_unsafe_upper = BARRIER - lag_poly_i_upper * (x[ii] - x_i_upper) - 1
        
                    # Add constraints to model
                    add_constraint_to_model(model, lag_poly_i_lower)
                    add_constraint_to_model(model, lag_poly_i_upper)
                    add_constraint_to_model(model, _barrier_unsafe_lower)
                    add_constraint_to_model(model, _barrier_unsafe_upper)

            end

        else
            continue
        end

    end

    # Variables g and h for Lagrange multipliers
    lagrange_monomial_length::Int64 = length_polynomial(x::Array{PolyVar{true},1}, lagrange_degree::Int64)
    number_of_variables_exp::Int64 = number_hypercubes_constraints * (system_dimension) * lagrange_monomial_length
    @variable(model, g[1:number_of_variables_exp])
    @variable(model, h[1:number_of_variables_exp])

    # Partition beta to extract ith beta values
    if beta_partition == true

        # Variables for beta in SOS
        num_vars_beta_lagrangian = number_hypercubes_constraints * lagrange_monomial_length
        @variable(model, delta[1:num_vars_beta_lagrangian])

        # Number of constraints
        number_constraints_per_loop = (2*system_dimension) + 1 + 1 + 1
        constraints = Array{DynamicPolynomials.Polynomial{true, AffExpr}}(undef, number_hypercubes_constraints, number_constraints_per_loop)

    else

        number_constraints_per_loop = (2*system_dimension) + 1
        constraints = Array{DynamicPolynomials.Polynomial{true, AffExpr}}(undef, number_hypercubes_constraints, number_constraints_per_loop)

    end

    # Counters
    counter_lag::Int64 = 0
    parts_count::Int64 = 0
    counter_beta::Int64 = 0

    for parts = 1:length(hcube_identifier)

        if hcube_identifier[parts] == 0.0
            continue
        else
            parts_count += 1
            identifier = Integer(hcube_identifier[parts])
        end

        # Create SOS polynomials for X (Partition) and Y (Bounds)
        hCubeSOS_X::DynamicPolynomials.Polynomial{true, AffExpr} = 0
        hCubeSOS_Y::DynamicPolynomials.Polynomial{true, AffExpr} = 0
        
        # Define global or explicit upper and lower bound for kth dimension of partition parts
        M_h_ii = transpose(M_h[identifier, :, :])
        M_l_ii = transpose(M_l[identifier, :, :])
        B_h_ii = B_h[identifier, :]
        B_l_ii = B_l[identifier, :]

        # Loop over hcube higher bound
        hyper_matrix_higher = M_h_ii * x + B_h_ii

        # Loop over hcube lower bound
        hyper_matrix_lower = M_l_ii * x + B_l_ii
   
        # Loop of state space and neural network bounds
        for kk = 1:system_dimension

            # Partition bounds
            x_k_hcube_bound::Vector{Float64} = partitions[identifier, :, kk]
            x_k_lower::Float64 = x_k_hcube_bound[1, 1]
            x_k_upper::Float64 = x_k_hcube_bound[2, 1]

            y_k_upper_explicit = hyper_matrix_higher[kk]
            y_k_lower_explicit = hyper_matrix_lower[kk]
          
            # Generate Lagrange polynomial for kth dimension
            lag_poly_X::DynamicPolynomials.Polynomial{true, AffExpr} = sos_polynomial(g::Vector{VariableRef}, x::Array{PolyVar{true},1}, (counter_lag + kk - 1)::Int64, lagrange_degree::Int64)
            lag_poly_Y::DynamicPolynomials.Polynomial{true, AffExpr} = sos_polynomial(h::Vector{VariableRef}, y::Array{PolyVar{true},1}, (counter_lag + kk - 1)::Int64, lagrange_degree::Int64)

            # Add Lagrange polynomial to constraints vector for the state space
            constraints[parts_count, kk] = lag_poly_X
            constraints[parts_count, kk + system_dimension] = lag_poly_Y

            # Generate SOS polynomials for state space
            hCubeSOS_X::DynamicPolynomials.Polynomial{true, AffExpr} += lag_poly_X*(x_k_upper - x[kk])*(x[kk] - x_k_lower)
            hCubeSOS_Y::DynamicPolynomials.Polynomial{true, AffExpr} += lag_poly_Y*(y_k_upper_explicit - y[kk])*(y[kk] - y_k_lower_explicit)
        end

        # Update system counter
        counter_lag += system_dimension

        # SOS for beta partition
        if beta_partition == true
    
            lag_poly_beta::DynamicPolynomials.Polynomial{true, AffExpr} = sos_polynomial(delta::Vector{VariableRef}, w::Array{PolyVar{true},1}, (counter_beta)::Int64, lagrange_degree::Int64)

            constraints[parts_count, (2*system_dimension) + 1] = lag_poly_beta

            counter_beta += 1

        end

        # Compute expectation
        _e_barrier::DynamicPolynomials.Polynomial{true, AffExpr} = BARRIER
        exp_evaluated::DynamicPolynomials.Polynomial{true, AffExpr} = _e_barrier

        for zz = 1:system_dimension
            exp_evaluated = subs(exp_evaluated, x[zz] => y[zz] + z)
        end

        # Extract noise term
        exp_poly, noise = expectation_noise(exp_evaluated, barrier_degree::Int64, standard_deviation::Float64, z::PolyVar{true})

        # Full expectation term
        exp_current = exp_poly + noise

        # Constraint for hypercube
        if beta_partition == true
            hyper_constraint = - exp_current + BARRIER/alpha + beta_parts_var[parts_count] - hCubeSOS_X - hCubeSOS_Y

            w_min = 0
            w_max = 1

            beta_constraint = (- beta_parts_var[parts_count] + beta - lag_poly_beta)*(w_max - w[1])*(w[1] - w_min)

            constraints[parts_count, number_constraints_per_loop - 1] = beta_constraint

        else
            hyper_constraint = - exp_current + BARRIER/alpha + beta - hCubeSOS_X - hCubeSOS_Y
        end


        # Add to model
        constraints[parts_count, number_constraints_per_loop] = hyper_constraint

    end
    
    # Add constraints to model as a vector of constraints
    @time begin
        @constraint(model, constraints .>= 0)
    end
    print("Constraints made\n")

    # Define optimization objective
    time_horizon = 1
    @objective(model, Min, eta + beta*time_horizon)
    print("Objective made\n")

    # Print number of partitions
    print("\n", "Optimizing for number of partitions = " * string(number_hypercubes), "\n")

    # Optimize model
    optimize!(model)

    # Barrier certificate
    certificate = barrier_certificate(barrier_monomial, c)

    # Return beta values
    eta_val = value(eta)
    if beta_partition == true
        beta_values = value.(beta_parts_var)
        max_beta = maximum(beta_values)

        # Print probability values
        println("Solution: [eta = $(value(eta)), beta = $(value(max_beta)), total = $(value(eta) + value(max_beta)) ]")
  
    else

        # Print probability values
        println("Solution: [eta = $(value(eta)), beta = $(value(beta)), total = $(value(eta) + value(beta)) ]")
        
    end

    # Print beta values to txt file
    if print_to_txt == true
        if isfile("probabilities/probs_cert_"*system_flag*string(number_hypercubes)*".txt") == true
            rm("probabilities/probs_cert_"*system_flag*string(number_hypercubes)*".txt")
        end

        open("probabilities/probs_cert_"*system_flag*string(number_hypercubes)*".txt", "a") do io
        println(io, "eta = $(value(eta)), beta = $(value(beta)), total = $(value(eta) + value(beta)) ")
        end

        # Print beta values to txt file
        if beta_partition == true

            if isfile("probabilities/beta_vals_certificate.txt") == true
                rm("probabilities/beta_vals_certificate.txt")
            end

            open("probabilities/beta_vals_certificate.txt","a") do io
                println(io, beta_values)
            end
        end
    end

    # Return optimization results
    if beta_partition == true
        return certificate, value.(eta), beta_values, system_dimension
    else
        return certificate
    end

end

########################################################################
########################################################################
########################################################################

# Sum of squares control function
function control_loop(input_data::Tuple, certificate, eta_certificate, x_star, minimum_interferance, beta_vals_opt)

    # Extract input data
    number_hypercubes = input_data[1]
    barrier_degree_input = input_data[2]
    safety_threshold_beta = input_data[4]
    system_flag = input_data[5]
    neural_network_bound = input_data[6]
    layer_flag = input_data[7]
    beta_partition = input_data[8]
    print_to_txt = input_data[11]

    # File reading
    filename = "/models/" * system_flag * "/" * neural_network_bound  * "/" * layer_flag* "_layers/partition_data_"  * string(number_hypercubes) * ".mat"
    file = matopen(pwd()*filename)

    # Extract hypercube data (avoid using float64 for precision issues)
    partitions = read(file, "partitions")
    state_space = read(file, "state_space")

    # Number of hypercubes
    number_hypercubes_constraints = number_hypercubes
    hcube_identifier = 1:number_hypercubes
    
    # Extract Neural Network Bounds
    M_h = read(file, "M_h")
    M_l = read(file, "M_l")
    B_h = read(file, "B_h")
    B_l = read(file, "B_l")

    # Define system and control dimensions
    system_dimension::Int64 = Integer(length(state_space[:,1]))

    # Using Mosek as the SDP solver
    model = SOSModel(optimizer_with_attributes(Mosek.Optimizer,
                                               "MSK_DPAR_INTPNT_TOL_STEP_SIZE" => 1e-6,
                                               "MSK_IPAR_OPTIMIZER" => 0,
                                               "MSK_IPAR_BI_CLEAN_OPTIMIZER" => 0,
                                               "MSK_IPAR_NUM_THREADS" => 16,
                                               "MSK_IPAR_PRESOLVE_USE" => 0))

    # Create state space variables
    @polyvar x[1:system_dimension]

    # Create noise variable
    @polyvar z

    # Create global CROWN bounds variables
    @polyvar y[1:system_dimension]

    # Create dummy variable for beta in SOS
    @polyvar w[1:2]

    # Create probability decision variables eta  and beta
    eta = eta_certificate
    if beta_partition == true
        @variable(model, beta_parts_var[1:number_hypercubes_constraints])
        @variable(model, beta)
    else
        @variable(model, beta)
    end

    # Create barrier polynomial, specify degree Lagrangian polynomials
    alpha::Float64 = 1
    lagrange_degree::Int64 = 2

    # Specify noise element (Gaussian)
    standard_deviation::Float64 = 0.1

    # Barrier
    BARRIER = certificate

    # Add constraints to model for positive barrier, eta and beta
    if beta_partition == true
        for betas = 1:number_hypercubes_constraints
            @constraint(model, beta_parts_var[betas] >= 1e-6)
        end
        @constraint(model, beta >= 1e-6)
    else 
        @constraint(model, beta >= 1e-6)
    end

    # Variables g and h for Lagrange multipliers
    lagrange_monomial_length::Int64 = length_polynomial(x::Array{PolyVar{true},1}, lagrange_degree::Int64)
    number_of_variables_exp::Int64 = number_hypercubes_constraints * (system_dimension) * lagrange_monomial_length
    @variable(model, g[1:number_of_variables_exp])
    @variable(model, h[1:number_of_variables_exp])

    # Partition beta to extract ith beta values
    if beta_partition == true

        # Variables for beta in SOS
        num_vars_beta_lagrangian = number_hypercubes_constraints * lagrange_monomial_length
        @variable(model, delta[1:num_vars_beta_lagrangian])

        # Number of constraints
        number_constraints_per_loop = (2*system_dimension) + 1 + 1 + 1
        constraints = Array{DynamicPolynomials.Polynomial{true, AffExpr}}(undef, number_hypercubes_constraints, number_constraints_per_loop)

    else

        number_constraints_per_loop = (2*system_dimension) + 1
        constraints = Array{DynamicPolynomials.Polynomial{true, AffExpr}}(undef, number_hypercubes_constraints, number_constraints_per_loop)

    end

    # Counters
    counter_lag::Int64 = 0
    parts_count::Int64 = 0
    counter_beta::Int64 = 0
    count_bad::Int64 = 0
    feedback_control_store = zeros(system_dimension, number_hypercubes)

    for parts = 1:length(hcube_identifier)

        if hcube_identifier[parts] == 0.0
            continue
        else
            parts_count += 1
            identifier = Integer(hcube_identifier[parts])
        end

        # Create SOS polynomials for X (Partition) and Y (Bounds)
        hCubeSOS_X::DynamicPolynomials.Polynomial{true, AffExpr} = 0
        hCubeSOS_Y::DynamicPolynomials.Polynomial{true, AffExpr} = 0
        
        # Define global or explicit upper and lower bound for kth dimension of partition parts
        M_h_ii = transpose(M_h[identifier, :, :])
        M_l_ii = transpose(M_l[identifier, :, :])
        B_h_ii = B_h[identifier, :]
        B_l_ii = B_l[identifier, :]

        # Loop over hcube higher bound
        hyper_matrix_higher = M_h_ii * x + B_h_ii

        # Loop over hcube lower bound
        hyper_matrix_lower = M_l_ii * x + B_l_ii
   
        # Call on feedback law
        if system_flag == "cartpole"
            control_dimensions = [0, 1, 0, 1]
        elseif system_flag == "husky4d"
            control_dimensions = [0, 0, 1, 1]
        elseif system_flag == "husky5d"
            control_dimensions = [0, 0, 0, 1, 1]
        elseif system_flag == "acrobot"
            control_dimensions = [0, 0, 0, 0, 1, 1]
        elseif system_flag == "pendulum"
            control_dimensions = [0 ,1]
        else
            print("System not defined ...")
        end
        if beta_vals_opt[parts] >= safety_threshold_beta
            count_bad += 1
            feedback_control = controller_convex(system_flag, x_star, system_dimension, identifier, control_dimensions, partitions, M_h_ii, M_l_ii, B_h_ii, B_l_ii)
            feedback_control_store[:, parts] = feedback_control
        end

        # Loop of state space and neural network bounds
        for kk = 1:system_dimension

            # Partition bounds
            x_k_hcube_bound::Vector{Float64} = partitions[identifier, :, kk]
            x_k_lower::Float64 = x_k_hcube_bound[1, 1]
            x_k_upper::Float64 = x_k_hcube_bound[2, 1]

            y_k_upper_explicit = hyper_matrix_higher[kk]
            y_k_lower_explicit = hyper_matrix_lower[kk]

            if beta_vals_opt[parts] >= safety_threshold_beta
                y_k_upper_explicit += feedback_control[kk]
                y_k_lower_explicit += feedback_control[kk]
            end
          
            # Generate Lagrange polynomial for kth dimension
            lag_poly_X::DynamicPolynomials.Polynomial{true, AffExpr} = sos_polynomial(g::Vector{VariableRef}, x::Array{PolyVar{true},1}, (counter_lag + kk - 1)::Int64, lagrange_degree::Int64)
            lag_poly_Y::DynamicPolynomials.Polynomial{true, AffExpr} = sos_polynomial(h::Vector{VariableRef}, y::Array{PolyVar{true},1}, (counter_lag + kk - 1)::Int64, lagrange_degree::Int64)

            # Add Lagrange polynomial to constraints vector for the state space
            constraints[parts_count, kk] = lag_poly_X
            constraints[parts_count, kk + system_dimension] = lag_poly_Y

            # Generate SOS polynomials for state space
            hCubeSOS_X::DynamicPolynomials.Polynomial{true, AffExpr} += lag_poly_X*(x_k_upper - x[kk])*(x[kk] - x_k_lower)
            hCubeSOS_Y::DynamicPolynomials.Polynomial{true, AffExpr} += lag_poly_Y*(y_k_upper_explicit - y[kk])*(y[kk] - y_k_lower_explicit)

        end

        # Update system counter
        counter_lag += system_dimension

        # SOS for beta partition
        if beta_partition == true
    
            lag_poly_beta::DynamicPolynomials.Polynomial{true, AffExpr} = sos_polynomial(delta::Vector{VariableRef}, w::Array{PolyVar{true},1}, (counter_beta)::Int64, lagrange_degree::Int64)

            constraints[parts_count, (2*system_dimension) + 1] = lag_poly_beta

            counter_beta += 1

        end

        # Compute expectation
        _e_barrier::DynamicPolynomials.Polynomial{true, AffExpr} = BARRIER
        exp_evaluated::DynamicPolynomials.Polynomial{true, AffExpr} = _e_barrier

        for zz = 1:system_dimension
            exp_evaluated = subs(exp_evaluated, x[zz] => y[zz] + z)
        end

        # Extract noise term
        barrier_degree = barrier_degree_input
        exp_poly, noise = expectation_noise(exp_evaluated, barrier_degree::Int64, standard_deviation::Float64, z::PolyVar{true})

        # Full expectation term
        exp_current = exp_poly + noise

        # Constraint for hypercube
        if beta_partition == true
            hyper_constraint = - exp_current + BARRIER/alpha + beta_parts_var[parts_count] - hCubeSOS_X - hCubeSOS_Y

            w_min = 0
            w_max = 1

            beta_constraint = (- beta_parts_var[parts_count] + beta - lag_poly_beta)*(w_max - w[1])*(w[1] - w_min)

            constraints[parts_count, number_constraints_per_loop - 1] = beta_constraint

        else
            hyper_constraint = - exp_current + BARRIER/alpha + beta - hCubeSOS_X - hCubeSOS_Y
        end

        # Add to model
        constraints[parts_count, number_constraints_per_loop] = hyper_constraint

    end
    
    # Add constraints to model as a vector of constraints
    @time begin
        @constraint(model, constraints .>= 0)
    end
    print("Constraints made\n")

    # Define optimization objective
    time_horizon = 1
    @objective(model, Min, beta*time_horizon) 
    print("Objective made\n")

    # Print number of partitions
    print("\n", "Optimizing for number of partitions = " * string(number_hypercubes), "\n")

    # Optimize model
    optimize!(model)

    # Return beta values
    if beta_partition == true
        beta_values = value.(beta_parts_var)
        max_beta = maximum(beta_values)

        # Print probability values
        println("Solution: [eta = $(value(eta)), beta = $(value(max_beta)), total = $(value(eta) + value(max_beta)) ]")
        # end

    else

        # Print probability values
        println("Solution: [eta = $(value(eta)), beta = $(value(beta)), total = $(value(eta) + value(beta)) ]")
        
    end

    # Print beta values to txt file
    if print_to_txt == true
        if isfile("probabilities/probs_cont_"*system_flag*string(number_hypercubes)*".txt") == true
            rm("probabilities/probs_cont_"*system_flag*string(number_hypercubes)*".txt")
        end

        open("probabilities/probs_cont_"*system_flag*string(number_hypercubes)*".txt", "a") do io
        println(io, "eta = $(value(eta)), beta = $(value(beta)), total = $(value(eta) + value(beta)) ")
        end

    end

    # Return optimization results
    if beta_partition == true
        return certificate, feedback_control_store, count_bad
    else
        return certificate
    end

end
