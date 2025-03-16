# -------------------------------------------------------------------------
# Identifying Lyapunov Functions from Data with EDMD
#
# The goal of this script is to use EDMD and semidefinite programming to
# identify a Lyapunov function for the planar discrete-time system
#      x1 --> 0.3*x1 ,
#      x2 --> -x1 + 0.5*x2 + (7/18)*x1^2,
# using only data gathered from the system.
#
# To run this notebook one requires paths to the freely available software
# package MOSEK, which can be downloaded at:
#           https://www.mosek.com/downloads/
#
# This script accompanies Section 4.3 of Data-Driven Methods for
# Dynamic Systems.
#
# This script was adapted from the MATLAB version by Jason J. Bramburger
# 
# Author: Mohit Sahu
# -------------------------------------------------------------------------


using DynamicPolynomials, SumOfSquares, MosekTools, LinearAlgebra
include("utils.jl")

# Method parameters
# Dpₘₐₓ = max degree of Dp dictionary of obserables
# Dqₘₐₓ = max degree of Dq dictionary of obserables
# ϵ = hyperparameter specific to Lyapunov function for sharp bounds
# cleanVal = remove coefficients smaller than this value from Lyap function

Dpₘₐₓ = 4;
Dqₘₐₓ = Dpₘₐₓ + 2;
ϵ = 1;
cleanVal = 1e-4;

# Generate Synthetic Data
# Number of data points
N = 100;

# Random points in [-2,2]x[-2,2]
xdat = 2 .*rand(N,2) .- 1;
ydat = zeros(N,2);

# Images of random points under map dynamics
ydat[:,1] = 0.3*xdat[:,1];
ydat[:,2] = -xdat[:,1] + 0.5*xdat[:,2] + (7/18)*xdat[:,1].^2;

# Create Q matrix
pow = monpowers(2, Dqₘₐₓ)
ell = size(pow)[1]  # number of nontrivial monomials in Q
Q = zeros(ell, N)
for i = 1:ell
    zx = xdat.^permutedims(pow[i]);
    Q[i,:] = prod(zx,dims=2);
end

# Create P matrix
pow = monpowers(2, Dpₘₐₓ)
ell = size(pow)[1]  # number of nontrivial monomials in P
P = zeros(ell, N)
for i = 1:ell
    zy = ydat.^permutedims(pow[i]);
    P[i,:] = prod(zy,dims=2);
end

# Create Koopman and Lie derivative matrix approximations
# Koopman approximations
K = P * pinv(Q)

# Symbolic variables
@polyvar x[1:2] # 2d state variables

z = monomials(x, 0:Dpₘₐₓ)  # monomials that make up span(Dp)
w = monomials(x, 0:Dqₘₐₓ)  # monomials that make up span(Dq)

model = SOSModel(Mosek.Optimizer)

@variable(model, c[1:length(z)])  # coeffs to build Lyapunov function in span(Dp)

# As julia doesn't handle L1 norm, so t here is a dummy variable to solve this problem
# See Tutorials/Linear programs/Tips and Tricks of Jump.jl documentation to know more about this
@variable(model, t) 

# Lie approximation
# ---> Lyapunov function = c.'*z in span(Dp)
L = c'*(K*w) - c'*z;

# Identify Lyapunov function with SOS programming
# Inequality constraints posed relaxed to SOS conditions
@constraints(model, begin
    c'*z - ϵ*dot(x,x) in SOSCone()
    -L - ϵ*dot(x, x) in SOSCone()
    [t; c] in MOI.NormOneCone(1+length(c))  # Minimizing L1 Norm of c
end)

# Objective function
@objective(model, Min, t)

# Solve SOS problem
optimize!(model)

# Remove coefficients smaller than cleanVal
c = [abs(c_)>cleanVal ? c_ : zero(c_) for c_ in value.(c)]

# Check identified Lyapunov function is a real Lyapunov function 
# This time we maximize ϵ and check if it is positive
model = SOSModel(Mosek.Optimizer)
@variable(model, ϵ)  # ϵ is now a variable to be maximized

# Discovered Lyapunov function 
v = c'*z

# Reverse the order of x₁ and x₂ to match the theory, as julia create polynomials in lexical order, so it's necessary to do to get to right answer
v = subs(v, x[1]=>x[2], x[2]=>x[1])

# True RHS of map
f = [0.3*x[1]; -x[1] + 0.5*x[2] + (7/18)*x[1]^2]

# True Lie Derivative
Lie_exact = subs(v, x=>f)

# Find Lyapunov function
@constraints(model, begin
    (v - ϵ * dot(x,x)) in SOSCone()
    (-(Lie_exact - v) - ϵ*dot(x,x)) in SOSCone()
end)

@objective(model, Max, ϵ)
optimize!(model)
solution_summary(model)
println("Maximized value of ϵ is $(value(ϵ))")

if value(ϵ) > 0
    println("We have discovered a true Lyapunov function from the data! \n") 
else
    println("Something went wrong - we cannot verify that this is a Lyapunov function. \n")
end