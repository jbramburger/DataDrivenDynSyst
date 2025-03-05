using DynamicPolynomials, SumOfSquares, MosekTools
include("utils.jl")

# ODE mu parameter
m = 2

# Bounding method parameters
degV = 20  # Degree of auxiliary function
λ = 10^3

# Bisection Method for Finding Upper Bound
function volume(μ, m, d, λ)
    @polyvar u v

    ϵ = 1e-4;

    model = SOSModel(Mosek.Optimizer)
    set_silent(model)

    #Auxiliary function
    @variable(model, V, createPoly([u,v], d))

    # S procedure
    d2 = d + m 
    @variables(model, begin
        s1, createPoly([u,v], d2)
        s2 ,createPoly([u,v], d2)
        s3 , createPoly(u, d2)
    end)

    # Derivatives
    dVdu = differentiate(V, u)
    dVdv = differentiate(V, v)

    # Replacements
    Vv0 = subs(V, v=>0)

    @constraints(model, begin
        subs(V,u=>0, v=>-ϵ) == 0
        subs(V,u=>1,v=>0) == 0
        λ*(dVdu*v - μ*dVdv*v - dVdv*(1-u)*u^m) -V -u*(1-u)*s1 + v*s2 in SOSCone()
        -Vv0 -ϵ*(1-u) -u*(1-u)*s3 in SOSCone()
        s1 in SOSCone()
        s2 in SOSCone()
        s3 in SOSCone()
    end)

    optimize!(model)

    # Return whether we are able to solve the problem or not
    return primal_status(model)
end

#Bisection method
μ_left = 0
μ_right = 2
μ_mid = 0

while abs((μ_right - μ_left)) >= 1e-5
    μ_mid = 0.5 * (μ_left + μ_right)
    flag = volume(μ_mid, m, degV, λ)

    if flag==FEASIBLE_POINT
        μ_left = μ_mid
    else
        μ_right = μ_mid
    end
end

println("The lower bound on the minimum speed for m = $m is $μ_mid found using degree $degV polynomials.\n")