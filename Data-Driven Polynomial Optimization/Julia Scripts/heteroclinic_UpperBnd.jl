# -------------------------------------------------------------------------
# Upper Bounds on the Existence of a Heteroclinic Orbit
#
# The goal of this script is to use semidefinite programming to provide
# upper bounds on the range of mu values for which a heteroclinic 
# connection exists in the system
#      x1' = x2, 
#      x2' = -μ*x2 - x1*(1-x1).
# The script searches for a desired auxiliary function using SOS 
# relaxations and localizations in phase space.   
#
# To run this notebook one requires paths to the freely available software
# package MOSEK, which can be downloaded at:
#           https://www.mosek.com/downloads/
#
# This script accompanies Section 4.2 of Data-Driven Methods for
# Dynamic Systems. 
#
# This script was adapted from the MATLAB version by Jason J. Bramburger
# 
# Author: Mohit Sahu
# -------------------------------------------------------------------------

using DynamicPolynomials, SumOfSquares, MosekTools
include("utils.jl")

# Bisection Method for Finding Upper Bound

# ODE mu parameter
m = 2

# Bounding method parameters
degV = 20  # Degree of auxiliary function
λ = 2

# SOS relaxations for finding V
function volume(μ, m, d, λ)
    @polyvar u v

    ϵ = 1e-4;

    model = SOSModel(Mosek.Optimizer)
    set_silent(model)

    #Auxiliary function
    @variable(model, V, createPoly([u,v],d))

    # S procedure
    d2 = d + m
    @variables(model, begin
        s1, createPoly([u,v],d2)
        s2, createPoly([u,v],d2)
        s3, createPoly(u,d2)
        s4, createPoly(v,d2)
        s5, createPoly(u,d2)
    end)

    # Derivatives
    dVdu = differentiate(V, u)
    dVdv = differentiate(V, v)

    # Function Replacements
    Vu0 = subs(V, u=>0)
    Vv0 = subs(V, v=>0)
    
    # constraints
    @constraints(model, begin
        subs(V, u=>0, v=>0) == 0
        -λ*(dVdu*v - μ*dVdv*v - dVdv*(1-u)*u^m) -V -u*(1-u)*s1 + v*s2 in SOSCone()
        Vu0 + ϵ*v + v*s4 in SOSCone()
        -Vv0 - ϵ*u*(1-u) - u*(1-u)*s5 in SOSCone()
        s1 in SOSCone()
        s2 in SOSCone()
        s3 in SOSCone()
        s4 in SOSCone()
        s5 in SOSCone()
    end)

    optimize!(model)

    # Return whether we are able to solve the problem or not
    primal_status(model)
end

#Bisection method
μ_left = 0
μ_right = 2
μ_mid = 0

# Check why this is giving false always, there must be something wrong in it's implementation
while abs((μ_right - μ_left)) >= 1e-5
    μ_mid = 0.5 * (μ_left + μ_right)
    flag = volume(μ_mid, m, degV, λ)

    if flag==FEASIBLE_POINT
        μ_right = μ_mid
    else
        μ_left = μ_mid
    end

end

println("The upper bound on the minimum speed for m = $m is $μ_mid found using degree $degV polynomials.\n")