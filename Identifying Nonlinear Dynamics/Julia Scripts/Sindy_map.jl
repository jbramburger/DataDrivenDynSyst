# -------------------------------------------------------------------------
# SINDy applied to Poincare and stroboscopic mapping discovery
# -------------------------------------------------------------------------

# This script applies the SINDy method to discovering discrete-time
# iterative schemes that govern coarse-grained data from continuous-time
# systems. Two examples are presented here:

# 1) A stroboscopic map of a periodically driven model for an RC circuit.
# 2) A Poincare map of a planar model of the truncated Hopf normal form.

# This script accompanies Section 3.2 of Data-Driven Methods for
# Dynamic Systems.

# For applications of this method to other systems see the repository:
# https://github.com/jbramburger/Poincare-Maps

# This script was adapted from the MATLAB version by Jason J. Bramburger
# 
# Author: Mohit Sahu
# -------------------------------------------------------------------------


using DifferentialEquations, LinearAlgebra, Plots
include("SINDy_module.jl");

# Model Parameters
A = 1
ω = 2*π;
numTraj = 3; #Number of trajectories

#  Generate Trajectories RC circuit equation
function RC!(du,u, p, t)
   x,y = u
   A,ω = p
   du[1] = A * sin(ω * y) - x
   du[2] = 1
end 

function generate_traj(integrator, x_init, tspan, dt)
   # returns back trajectories for different x0
   xdat = zeros((length(x_init), Int((tspan[2]-tspan[1])/dt + 1),2))
   for j in 1:length(x_init)
       reinit!(integrator, x_init[j])
       n = 1;
       for i in tspan[1]:dt:tspan[2]
           xdat[j, n,:] = integrator.u;
           step!(integrator, dt, true);
           n+=1;
       end
   end

   return xdat
end

dt = 0.01;
tspan = (0.0,20000.0*dt);
p = [A, ω];

x_init = Vector{Float64}[];
x0 = [0.0, 0.0]; # At least one solution
push!(x_init,x0)
prob = ODEProblem(RC!, x0, tspan, p)
integrator = init(prob, reltol = 1e-12, abstol = 1e-12)

if numTraj >= 2
   for k = 2:numTraj
       push!(x_init, [10rand()-5, 0.0]);
   end
end

xdat = generate_traj(integrator, x_init, tspan, dt);

plot(xdat[3,1:10000,1])

# Aggregate stroboscopic section data
#Counting parameter
#Initialize
Psec = Float64[0]
PsecNext = Float64[] 

#Create Poincaré section data
for i in 1:numTraj
   for j in 1:(size(xdat, 2) - 1)  # Loop over all time steps except the last
       if (j == 1) && (i > 1)  # Trajectories start in the section
           push!(Psec, xdat[i, j, 1])  # Append to Psec
       elseif (mod(xdat[i,  j, 2,], 2π / ω) >= 2π / ω - dt && mod(xdat[i,  j+1, 2], 2π / ω) <= dt)
           push!(Psec, xdat[i,  j+1, 1,])  # nth iterate
           push!(PsecNext, xdat[i,  j+1, 1,])  # (n+1)st iterate
       end
   end
   Psec = Psec[1:length(Psec)-1];
end
# Create the recurrence data
xn = Psec;
xnp1 = PsecNext;

# Create Θ matrix from monomial upto degree 5
polyorder = 5; # change maximal polynomial order of monomials in library
Θ = ones((polyorder+1,length(xn))); # constant term

for p = 1:polyorder
   Θ[p+1,:] = xn.^p;
end

xnp1 = reshape(xnp1, (1,length(xnp1)));

total_dim = 1;
Xi_new = sindy(xnp1, Θ,total_dim, 1e-5);
mons = ["" "x" "x^2" "x^3" "x^4" "x^5"];
print_model(Xi_new,mons, total_dim, ["f(x)"])

## Hopf Normal Form ODE
function Hopf!(du,u,p,t)
   x,y = u
   du[1] = x - p * y - x * (x^2 + y^2)
   du[2] = p * x + y - y * (x^2 + y^2)
end

# Model parameters
om = 10π;  # Hopf normal form period (make larger for better map result)
m = 2; # Dimension of ODE
n = m-1; # Dimension of Poincare section
dt = 0.005;
tspan = (0.0,10000.0*dt);

x_init = Vector{Float64}[];
push!(x_init, [0.0001, 0.0]); # Initial condition (Guarantee one starts near origin)
push!(x_init, [0.0, 0.0]); #Trajectory at the fixed point

numTraj  = 10;
if numTraj >= 3
   for k = 3:numTraj
       push!(x_init, [5rand(), 0.0]);
   end
end

prob = ODEProblem(Hopf!, x_init[1], tspan, om)
integrator = init(prob, reltol = 1e-12, abstol = 1e-12)
xdat = generate_traj(integrator, x_init, tspan, dt);
xdat = xdat[:,1:end-1,:];

# Aggregate Poincare section
# Initialize
Psec = Float64[];
PsecNext = Float64[];
push!(Psec, xdat[1,1,1]);
# Create Poincare section data
for i = 1:numTraj
   for j = 1:length(xdat[1,:,1])-1
       if (xdat[i,j,2] < 0) && (xdat[i,j+1,2] >= 0)
           push!(Psec,xdat[i,j+1,1]); #nth iterate
           push!(PsecNext, xdat[i,j+1,1]); #(n+1)st iterate
       end
   end
end
# Create the recurrence data

Psec = Psec[1:length(Psec)-1];
xn = Psec;
xnp1 = PsecNext;

# Create Θ matrix from monomial upto degree 5
polyorder = 2; # change maximal polynomial order of monomials in library
Θ = ones((polyorder+1,length(xn))); # constant term

for p = 1:polyorder
   Θ[p+1,:] = xn.^p;
end

xnp1 = reshape(xnp1, (1,length(xnp1)));

total_dim = 1;
Xi_new = sindy(xnp1, Θ,total_dim, 1e-5);
mons = ["" "x" "x^2" "x^3" "x^4" "x^5"];
print_model(Xi_new,mons, total_dim, ["f(x)"])
