# Control hypercube optimization 
function controller_convex(system_flag, x_star, system_dimension, identifier, controller_dimension, partitions, M_h_ii, M_l_ii, B_h_ii, B_l_ii)

    # Optimize code later: only include controlled dimensions 

    # Initialize solver
    lp_control_model = Model()
    set_optimizer(lp_control_model, GLPK.Optimizer)

    # Create feedback law
    @variable(lp_control_model, theta[1:system_dimension])
    @variable(lp_control_model, x[1:system_dimension])
    @variable(lp_control_model, x_prime[1:system_dimension])
    @variable(lp_control_model, y[1:system_dimension])
    @variable(lp_control_model, y_prime[1:system_dimension])
    @variable(lp_control_model, u[1:system_dimension])

    # Argmin of barrier 
    # x_star = zeros(1, system_dimension)
    
    # Free variable bounds
    @constraint(lp_control_model, x .== y - y_prime)
    @constraint(lp_control_model, y .>= 0)
    @constraint(lp_control_model, y_prime .>= 0)

    # Loop over hcube higher bound
    hyper_matrix_higher = M_h_ii * (y-y_prime) + B_h_ii
    
    # Loop over hcube lower bound
    hyper_matrix_lower = M_l_ii * (y-y_prime) + B_l_ii
    
    # Control bounds
    if system_flag == "cartpole"
        u_min = -1
        u_max = 1
    else
        u_min = -1
        u_max = 1
    end

    # Create constraints
    for ctrl = 1:system_dimension

        x_k_hcube_bound_jj = partitions[identifier, :, ctrl]
        x_k_lower_jj = x_k_hcube_bound_jj[1, 1]
        x_k_upper_jj = x_k_hcube_bound_jj[2, 1]
        @constraint(lp_control_model,  y[ctrl] - y_prime[ctrl] >=  x_k_lower_jj)
        @constraint(lp_control_model,  y[ctrl] - y_prime[ctrl] <=  x_k_upper_jj)

        @constraint(lp_control_model, x_prime[ctrl] - x_star[ctrl] <= theta[ctrl])
        @constraint(lp_control_model, x_star[ctrl] - x_prime[ctrl] <= theta[ctrl])

        y_k_upper_explicit = hyper_matrix_higher[ctrl]
        y_k_lower_explicit = hyper_matrix_lower[ctrl]
        @constraint(lp_control_model, x_prime[ctrl] >= y_k_lower_explicit + u[ctrl])
        @constraint(lp_control_model, x_prime[ctrl] <= y_k_upper_explicit + u[ctrl])

        # Controller logic
        if controller_dimension[ctrl] == 1
            @constraint(lp_control_model, u[ctrl] >= -1)
            @constraint(lp_control_model, u[ctrl] <= 1)
        else controller_dimension[ctrl] == 0
            @constraint(lp_control_model, u[ctrl] == 0)
        end

    end

    #Define Objective
    @objective(lp_control_model, Min, sum(theta))

    #Run the opimization
    set_silent(lp_control_model)
    optimize!(lp_control_model)

    # Return feedback control law
    feedback_law = value.(u)
    return feedback_law

end