# Functions call for Sum of Squares Optimization

# Create Control Barrier Polynomial
function barrier_polynomial(c::Vector{VariableRef}, barrier_monomial)
    barrier_poly = 0
    for cc in 1:Integer(length(barrier_monomial))
        barrier_poly += c[cc] * barrier_monomial[cc]
    end
    return barrier_poly
end

# Create SOS polynomial function
function sos_polynomial(k::Vector{VariableRef}, var, k_count::Int64, lagrange_degree::Int64)
    sos_polynomial = monomials(var, 0:lagrange_degree)
    sos_poly_t = 0
    for sos in 1:Integer(length(sos_polynomial))
        sos_poly_t += k[sos + k_count*length(sos_polynomial)] * sos_polynomial[sos]
    end
    return sos_poly_t
end

# Function to compute number of decision variables per Lagrange function
function length_polynomial(var, degree::Int64)::Int64
    sos_polynomial = monomials(var, 0:degree)
    length_polynomial::Int64 = length(sos_polynomial)
    return length_polynomial
end

# Function to add constraints to the model
function add_constraint_to_model(model::Model, expression)
    @constraint(model, expression >= 0)
end

# Function to compute the expecation and noise element
function expectation_noise(exp_evaluated, standard_deviations, zs)
    exp_poly = 0
    noise = 0

    for term in terms(exp_evaluated)
        z_degs = [MultivariatePolynomials.degree(term, z) for z in zs]

        z_occurs = sum(z_degs) > 0

        if z_occurs == false
            exp_poly = exp_poly + term
        end

        if z_occurs == true
            all_even = all(iseven, z_degs)

            if all_even
                coeff = subs(term, zs => ones(length(zs)))
                exp_z = prod([expected_univariate_noise(z_deg, standard_deviation) for (z_deg, standard_deviation) in zip(z_degs, standard_deviations)])

                noise_exp = coeff * exp_z
                noise = noise + noise_exp
            end
        end
    end

    return exp_poly, noise
end

function expected_univariate_noise(z_deg, standard_deviation)
    if z_deg == 0
        return 1
    else
        return (doublefactorial(z_deg - 1) * standard_deviation^z_deg)
    end
end

# Compute the final barrier certificate
function barrier_certificate(barrier_monomial, c)

    # Control Barrier Certificate
    barrier_certificate = 0
    for cc in 1:Integer(length(barrier_monomial))
        barrier_certificate += value(c[cc]) * barrier_monomial[cc]
    end

    return barrier_certificate

end