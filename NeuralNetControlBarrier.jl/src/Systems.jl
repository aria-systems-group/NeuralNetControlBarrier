# Systems for Verification

# Input function: change parameters here to run optimizations
function inputs(system_flag::String)::Tuple

    barrier_degree_input = 4
    safety_threshold_eta = (1 - 1e-6) # norm = (1 - 1e-6)
    safety_threshold_beta = 0.05
    neural_network_bound = "alpha"
    beta_partition = true
    large_range_initial = true
    print_to_txt = true 

    # Pendulum
    if system_flag == "pendulum"
       number_hypercubes = 120
       layer_flag = "1"
       initial_set_radius = deg2rad(5.0)
       decision_eta_flag = true

    # Cartpole
    elseif system_flag ==  "cartpole"
        number_hypercubes = 960
        layer_flag = "1"
        initial_set_radius = deg2rad(5.0)
        decision_eta_flag = true

    # Husky 4D
    elseif system_flag == "husky4d"
        number_hypercubes = 900
        layer_flag = "1"
        initial_set_radius = 0.1
        decision_eta_flag = false

    # Husky 5D
    elseif system_flag == "husky5d"
        number_hypercubes = 432
        layer_flag = "1"
        initial_set_radius = 0.1
        decision_eta_flag = false

    # Acrobot
    elseif system_flag == "acrobot"
        number_hypercubes = 144
        layer_flag = "1"
        initial_set_radius = 0.1
        decision_eta_flag = false

    else
        print("System not defined ...")
        return 0
    end

    # Return system input values
    return  number_hypercubes,
            barrier_degree_input,
            safety_threshold_eta,
            safety_threshold_beta,
            system_flag,
            neural_network_bound,
            layer_flag,
            beta_partition,
            large_range_initial,
            initial_set_radius,
            print_to_txt,
            decision_eta_flag

end

               
