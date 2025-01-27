# -------------------------------------------------------------------------
# Sparse Identification of Nonlinear Dynamics (SINDy)
# -------------------------------------------------------------------------
#
# This script applies the SINDy method of Proctor, Brunton, & Kutz (PNAS,
# 2016) to data collected from the Lorenz or Rossler chaotic systems. The
# script is broken into three parts: using a quadratic library, a cubic
# library, and then using the weak/integral SINDy formulation with a
# quadratic library.
#
# This script accompanies Section 3.1 of Data-Driven Methods for
# Dynamic Systems.
#
# To reproduce the noisy computations in the book load LorenzNoise.mat to
# get the noisy Lorenz data.
#
# This script was adapted from the MATLAB version by Jason J. Bramburger
# 
# Author: Mohit Sahu
# -------------------------------------------------------------------------

using DifferentialEquations, Plots, LinearAlgebra

function lorenz!(du, u, p, t)
    x, y, z = u
    σ, ρ, β = p
    du[1] = dx = σ * (y - x)
    du[2] = dy = x * (ρ - z) - y
    du[3] = dz = x * y - β * z
end

function Rossler!(du, u, p, t)
    x,y,z = u
    a,b,c = p
    du[1] = -y -z;
    du[2] = x + a * y;
    du[3] = b + z * (x - c)
end

# Simulating Lorenz system
dt = 0.0001;
tspan = (0.0, 100.0)
σ = 10.0;
β = 8/3;
ρ = 28.0;
p = [σ, ρ, β]
u0 = [1.0; 2.0; ρ-1]

prob = ODEProblem(lorenz!, u0, tspan,p)

# simulating Rossler system
# a = 0.1;
# b = 0.1;
# c = 18;
# p = [a,b,c]
# u0 = [0 -15 0];

# prob = ODEProblem(Rossler!, u0, tspan,p)

xsol = solve(prob, abstol=1e-12, reltol=1e-12, saveat=dt)

# Add noise
var = 0.0; # noise variance
xsol .= xsol .+ sqrt(var)*randn(size(xsol));
xsol = xsol'; # change dimensions to match theory

# or directly read the file from saved noisy data to reproduce the text
using MAT
file = matopen("D:/Study/Intelligence And Learning/PIM/Learning/My Code/DataDrivenDynSyst-main/Identifying Nonlinear Dynamics/LorenzNoise.mat")
xsol = read(file, "xsol")
close(file)

# plot(xsol[1,:], xsol[2,:], xsol[3,:])  # do this if loaded from LorenzNoise.mat file
plot(xsol)

# Estimated Derivatives
dxdt = (xsol[:,2:end] .- xsol[:,1:end-1]) ./ dt; # Y matrix
x = xsol[:,1:end-1];  #X matrix

Θ = ones(10,length(x[1,:])); # constant term
Θ[2:4,:] .= x; #linear term
Θ[5,:] = x[1,:].^2; # x^2 term
Θ[6,:] = x[1,:].*x[2,:]; # xy term
Θ[7,:] = x[1,:].*x[3,:]; # xz term
Θ[8,:] = x[2,:].^2; # y^2 term
Θ[9,:] = x[2,:].*x[3,:]; # yz term
Θ[10,:] = x[3,:].^2; # z^2 term

size(Θ[2,:])

dxdt = reshape(dxdt, (3,size(dxdt)[2]))  # to change dxdt from vector to matrix

function sindy(dxdt, Θ, lam)
    Xi = dxdt * pinv(Θ);

    k=1;
    Xi_new = Xi;

    while sum(abs.(Xi - Xi_new)) > 0  || k == 1
        
        Xi = deepcopy(Xi_new);

        # loop over all 3 variables
        for ind=1:3
            # find all the big indexes, means they should be above some threshold
            biginds = findall(x-> abs(x) > lam, Xi[ind,:])
            smallinds = setdiff(1:length(Xi[ind,:]), biginds)

            Xi_new[ind,smallinds] .= 0.0;

            # Find coefficients for reduced library
            Xi_new[ind, biginds] = reshape(dxdt[ind,:],(1, length(dxdt[ind,:]))) * pinv(Θ[biginds, :])
        end

        k+=1;

    end

    return Xi_new
end

# Printing the model
Xi_new = sindy(dxdt, Θ, 1e-1)
mons2 = ["" "x" "y" "z" "x^2" "xy" "xz" "y^2" "yz" "z^2"]

function print_model(Xi_new, mons)
    println("Discovered Model using SINDy: ")
    for ind = 1:3
        println("d" * mons[ind+1] *"/dt = ")
        bigcoeffs = abs.(Xi_new[ind,:]) .> 1e-3; # chosen small just to weed out zero coefficients
        for jnd = eachindex(bigcoeffs)
            if bigcoeffs[jnd] == 1
                # Print the model by excluding zeroed out terms
                if Xi_new[ind,jnd] < 0
                    print("-" * string(round(abs.(Xi_new[ind, jnd]), digits=4),mons[jnd]))
                    #print('- %s ', strcat(Xi_new[ind,jnd],mons[jnd]));
                else
                    print("+"*string(round(abs.(Xi_new[ind,jnd]), digits=4),mons[jnd]))
                    # print('+ %s', strcat(Xi_new[ind,jnd],mons[jnd]));
                end
                
            end
        end
        print("\n")
    end
end
print_model(Xi_new, mons2)


# Adding cubic terms to the library
Θ = [Θ;  x[1,:]'.^3]; # x^3 term
Θ = [Θ; x[1,:]'.*x[1,:]'.*x[2,:]']; # xxy term
Θ = [Θ; x[1,:]'.*x[2,:]'.*x[2,:]']; # xyy term
Θ = [Θ; x[1,:]'.*x[2,:]'.*x[3,:]']; # xyz term
Θ = [Θ; x[1,:]'.*x[3,:]'.*x[3,:]']; # xzz term
Θ = [Θ; x[2,:]'.*x[2,:]'.*x[2,:]']; # y^3 term
Θ = [Θ; x[2,:]'.*x[2,:]'.*x[3,:]']; # yyz term
Θ = [Θ; x[2,:]'.*x[3,:]'.*x[3,:]']; # yzz term
Θ = [Θ; x[3,:]'.*x[3,:]'.*x[3,:]']; # z^3 term

Xi_new = sindy(dxdt, Θ, 1e-1)
mons3 = ["" "x" "y" "z" "x^2" "xy" "xz" "y^2" "yz" "z^2" "x^3" "x^2y" "xy^2" "xyz" "xz^2" "y^3" "y^2z" "yz^2" "z^3"];
print_model(Xi_new, mons3)

# Now using Integrals instead of derivatives to counter noise
x = xsol[:,1:end]
Y = xsol[:, 2:end] .- xsol[:,1]
∫Θ = zeros(size(Θ[1:10,:])); # Just keep quadratic terms

# Intialize ThetaInt
∫Θ[:,1] = dt .* Θ[1:10,1];

for n = 2:length(Θ[1,:])
    ∫Θ[:,n] = ∫Θ[:,n-1] + dt * Θ[1:10, n];
end

Xi_new = sindy(Y, ∫Θ,1e-1)
print_model(Xi_new, mons2)