# Functions call for Sum of Squares Optimization

# Create Control Barrier Polynomial
function barrier_polynomial(c::Vector{VariableRef}, barrier_monomial::MonomialVector{true})::DynamicPolynomials.Polynomial{true, AffExpr}
    barrier_poly = 0
    for cc in 1:Integer(length(barrier_monomial))
        barrier_poly += c[cc] * barrier_monomial[cc]
    end
    return barrier_poly
end

# Create SOS polynomial function
function sos_polynomial(k::Vector{VariableRef}, var::Array{PolyVar{true},1}, k_count::Int64, lagrange_degree::Int64)::DynamicPolynomials.Polynomial{true, AffExpr}
    sos_polynomial::MonomialVector{true}  = monomials(var, 0:lagrange_degree)
    sos_poly_t = 0
    for sos in 1:Integer(length(sos_polynomial))
        sos_poly_t += k[sos + k_count*length(sos_polynomial)] * sos_polynomial[sos]
    end
    return sos_poly_t
end

# Function to compute number of decision variables per Lagrange function
function length_polynomial(var::Array{PolyVar{true},1}, degree::Int64)::Int64
    sos_polynomial::MonomialVector{true}  = monomials(var, 0:degree)
    length_polynomial::Int64 = length(sos_polynomial)
    return length_polynomial
end

# Function to add constraints to the model
function add_constraint_to_model(model::Model, expression::DynamicPolynomials.Polynomial{true, AffExpr})
    @constraint(model, expression >= 0)
end

# Function to compute the expecation and noise element
function expectation_noise(exp_evaluated::DynamicPolynomials.Polynomial{true, AffExpr}, barrier_degree::Int64, standard_deviation::Float64, z::PolyVar{true})

    exp_poly::DynamicPolynomials.Polynomial{true, AffExpr} = 0
    noise::DynamicPolynomials.Polynomial{true, AffExpr} = 0

    for zz in 1:length(exp_evaluated)
        z_occurs = occursin('z', string(exp_evaluated[zz]))

        if z_occurs == false
            exp_poly = exp_poly + exp_evaluated[zz]
        end

        if z_occurs == true
            for z_deg = 2:2:barrier_degree
                even_order_z = contains(string(exp_evaluated[zz]), "z^$z_deg")
                if even_order_z == true
                    exp_deep_evaluated::Term{true, AffExpr} = exp_evaluated[zz]
                    z_coefficients::Term{true, AffExpr} = subs(exp_deep_evaluated, z => 1)

                    noise_exp = z_coefficients * (doublefactorial(z_deg - 1) * standard_deviation^z_deg)
                    noise = noise + noise_exp
                end
            end
        end
    end
    return exp_poly, noise
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