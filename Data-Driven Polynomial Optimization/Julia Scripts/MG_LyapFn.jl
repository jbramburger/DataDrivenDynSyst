# -------------------------------------------------------------------------
# Identifying Lyapunov Functions with Semidefinite Programming
#
# The goal of this script is to use semidefinite programming to identify a
# Lyapunov function for the planar Moore-Greitzer system
#      x1' = -x2 - 1.5x1^2 - 0.5x1^3,
#      x2' = 3x1 - x2.
# The script expands the unknown Lyapunov function as a polynomial up to a
# given degree and tunes the coefficients in order to satisfy both of the
# constraints:
#      V - x1^2 - x2^2 is SOS
#      -LV - x1^2 - x1^2 is SOS
# where LV is the Lie derivative associated to the Moore-Greitzer system.
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


using DynamicPolynomials, SumOfSquares
using MosekTools
using CSDP

@polyvar x[1:2]

f = [-x[2] - 1.5*x[1]^2 - 0.5*x[1]^3;
      3*x[1] - x[2]]

# solver = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true)
model = Model(Mosek.Optimizer)
V_basis = monomials(x, 0:6)  # monomials for polynomial V
@variable(model, c[1:length(V_basis)])

# As julia doesn't handle L1 norm, so t here is a dummy variable to solve this problem
# See Tutorials/Linear programs/Tips and Tricks of Jump.jl documentation to know more about this
@variable(model, t)

V = c'*V_basis

using LinearAlgebra
Dv = differentiate(V, x)
Lie = dot(Dv, f)

# Add SOS constraints
@constraints(model, begin
    V .- sum(x.^2) in SOSCone()
    -Lie - sum(x.^2) in SOSCone()
    [t; c] in MOI.NormOneCone(1+length(c))   # L1 Norm
end)

@objective(model, Min, t)
optimize!(model)

# Check if a solution was found
if termination_status(model) == MOI.OPTIMAL
    println("Lyapunov function identified!")
    V_sol = value(V)
    threshold = 1e-1
    filtered_coeffs = [abs(c) > threshold ? c : zero(c) for c in coefficients(V_sol)]
    println("Filtered V(x1, x2) = ", sum(filtered_coeffs .* monomials(V_sol)))
else
    println("Sorry, unable to identify a Lyapunov function")
end

# Plot the Lyapunov function
using Plots
plotlyjs()
N = 1000;
x1_range = LinRange(-5, 5, N)

x = x1_range' .* ones(length(x1_range))
y = ones(length(x1_range))' .* x1_range
V_discovered = value.(c)'*V_basis;

plot(x, y, V_discovered.(x,y), st=:surface)
xlabel!("x1")
ylabel!("x2")
zlabel!("V(x1, x2)")
