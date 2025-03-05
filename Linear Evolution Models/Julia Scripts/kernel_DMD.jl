# -------------------------------------------------------------------------
# Kernel Dynamic Mode Decomposition (EDMD)
#
# The goal of this script is to apply the kernel method of dynamic mode
# decomposition to solutions to the Burgers equation, given by the
# nonlinear partial differential equation:
#             uₜ = ν*uₓₓ - u*uₓ
# The method is meant to circumvent exploding dictionary sizes coming from
# high-dimensional data by replacing inner products of dictionary
# evaluations with an evaluation of a kernel function. Here we use the
# polynomial kernel h(x,y) = (1 + yᵀx)ᵈ, where we take the degree d to
# be a parameter.
#
# The last block in the code applies EDMD to the Burgers
# equation using the single observable given by the Cole-Hopf
# transformation. The significance is that this transformation v(u) turns the
# Burgers equation into the linear heat equation
            # vₜ = ν*vₓₓ
# and so the Koopman expansion can be full realized using separation of
# variables.
#
# This script accompanies Section 2.5 of Data-Driven Methods for
# Dynamic Systems.
#
# This script was adapted from the MATLAB version by Jason J. Bramburger
# 
# Author: Mohit Sahu
# -------------------------------------------------------------------------


using SparseArrays, DifferentialEquations, LinearAlgebra

# Viscosity parameter
ν = 0.1
# Spatial and temporal domain parameters
dt = 1e-1
tspan = [0.0, 10.0]
L = 2*π # spatial domain is [-L/2,L/2] - L is the length
# spatial points in each dimension
n = 2^8
x = hcat(range(-L/2, L/2, length=n))
h = x[2] - x[1]

# d_x with Dirichlet BCs
D = sparse(1:n-1, 2:n, ones(n-1) ./ 2 , n, n)
D = (D - D')/h

# d_xx with Dirichlet BCs
D2 = sparse(1:n-1,2:n,ones(n-1),n,n) .- sparse(1:n,1:n,ones(n),n,n)
D2 = (D2 + D2')
D2 = D2/h^2

# Initialize for multiple initial conditions
numICs = 5
xdat = []
ydat = []

function burgers!(du, u, p, t)
    ν, D, D2 = p
    
    # Burger equation in Fourier Space
    du .= ν.*(D2*u) - u.*(D*u)
end

function prob_func(prob, i, repeat)
    # Randomized initial conditions
    numSines = rand(1:10)
    u0 = zeros(n,1)
    for jnd = 1:numSines
        u0 = @. u0 + (2*rand() - 1)*sin.(jnd*x);
    end

    remake(prob, u0=u0)
end

# Simulate PDE using Fourier methods
p = (ν, D, D2)
prob = ODEProblem(burgers!,zeros(n,1), tspan, p)
ensemble_prob = EnsembleProblem(prob, prob_func=prob_func)

# To use multithreading make sure that JULIA_NUM_THREADS is more than 1
sim = solve(ensemble_prob, DP5(), EnsembleThreads(), trajectories = numICs, saveat=dt)

# Aggregate data
xdat = zeros((256, 500))
ydat = zeros((256, 500))

ind = 1;
for i=1:numICs
    xdat[:,ind:ind+99] = sim.u[i][:,1,1:end-1]
    ydat[:,ind:ind+99] = sim.u[i][:,1,2:end]
    ind+=100
end

X = x*ones((1,101))
T = ones((256,1))*sim.u[1].t'

# Plot for visualization
using Plots
plotlyjs()
fig = plot()
surface!(X,T,sim.u[1][:,1,:],c=:gray)
plot!(x,zero(x), sim.u[1][:,1,1], c=:red, lw=5)
plot!(x,10*ones((256,1)), sim.u[1][:,1,100], c=:blue, lw=5)
xlabel!("x")
ylabel!("t")
zlabel!("u(x,t)")

# Build Data matrices with kernel
# Number of time elements
N = length(xdat[1,:])

# Polynomial degree
d = 8

# Kernel matrices
ΨxΨx = (ones(N,N) + xdat'*xdat).^d;
ΨyΨx = (ones(N,N) + ydat'*xdat).^d;

# Diagonalize ΨxΨx using SVD to get V and Σ
F = svd(ΨxΨx)
U = F.U
Σ = Diagonal(sqrt.(F.S))
V = F.Vt

# Create hat(A) matrix and analyze spectrum
Ahat = pinv(Σ)*U'*ΨyΨx*V*pinv(Σ)
F = eigen(Ahat)
μ = F.values

# Make axis line
line = -2:2

# Unit Circle
th = LinRange(0.0, 2π, 1000)
xcos = cos.(th)
ysin = sin.(th)

# Plot kernel DMD eigenvalues
plot(zeros(length(line)), line, lw=2, lc=:black)  # imaginary axis
plot!(line, zeros(length(line)), lw=2, lc=:black)  # real axis
plot!(xcos,ysin,ls=:dash,lw=2) # unit circle
scatter!(μ,c=RGB(1,69/255,79/255))
xlims!(-1.2,1.2)
ylims!(-1.2,1.2)

# DMD on Cole-Hopf transformation
using Trapz

#  Create Cole-Hopf transformed variable v(u)
xi = zeros(length(sim.u[1].t),n)
vCH = zeros(length(sim.u[1].t),n)
for ind = 2:n-1
    xi[:,ind] = (-1/(2*ν)).*trapz(x[1:ind],sim.u[1][1:ind,1,:]')'
end
expxi = exp.(xi)
expxiInt = trapz(x,expxi)
vCH = expxi./expxiInt

# Create DMD matrices
X = vCH[1:end-1,:]'
Y = vCH[2:end,:]'

# DMD matrix
A = Y*pinv(X)
F = eigen(A) # compute eigenvalues + eigenvectors
μ = F.values # extract eigenvalues

# Continuous time eigenvalues for reference
ω = @. log(μ)*L^2/dt/ν/pi^2

# Check which are close to integers
ωSort = sort(sqrt.(abs.(ω)))