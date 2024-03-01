# -------------------------------------------------------------------------
# Delay Coordinate Dynamic Mode Decomposition 
#
# This code applies the delay-coordinate DMD approach to scalar time series.
# There are two examples included: one synthetic example generated from a
# simple linear span of sines and cosines with differing frequencies and 
# another from the measurements of the nonlinear oscillations of the Van 
# der Pol oscillator. We further take the SVD of the Hankel matrix 
# associated to the Van der Pol measurements to demonstrate the nearly 
# linear oscillation of the observables coming from the left singular
# vectors in U.
#
# This script accompanies Section 2.3 of Data-Driven Methods for
# Dynamic Systems. 
# 
# 
# This script was adapted from the MATLAB version by Jason J. Bramburger
# 
# Author: Daniel Fassler
# -------------------------------------------------------------------------
using LinearAlgebra, Random, Plots, FFTW, DifferentialEquations

## Functions
    # Van der Pol system of ODEs
    function VdP(u,p,t)
        return [u[2], -u[1] + 10*(1 - u[1]^2)*u[2]]
    end

    # Create a hankel matrix with first column c and last row r
    function hankel(c, r)
        N, M = length(c), length(r)
        hankelMatrix = zeros(N, M)
        hankelMatrix[:,1] = c
        hankelMatrix[N,:] = r

        for j = 2:M
            for i = 1:N-1
                hankelMatrix[i,j] = hankelMatrix[i+1,j-1]
            end
        end
        return hankelMatrix
    end

## Simple Synthetic example
    # Generate data
    Δt = 0.1
    t = 0:Δt:40
    numOsc = 8
    x = zeros(1, length(t))

    # randomiwe frequencies and choice of sine and cosine
    ω = rand(1:20, numOsc)
    sin_or_cos = rand(0:1, numOsc)
    coeffs = 10*rand(numOsc) .- 5

    # Create the signal x(t) = c₁f₁(ω₁*t) + c₂f₂(ω₂*t) + ... + cₙfₙ(ωₙ*t)
    # fᵢ either sine or cosine
    for i = 1:numOsc
        global x
        if sin_or_cos[i] == 1
            x += coeffs[i]*sin.(ω[i]*t)'
        else
            x += coeffs[i]*cos.(ω[i]*t)'
        end
    end

    # Hankel Matrix
    delays = 2*numOsc # τ = delays - 1
    xd = hankel(x[1:delays + 1], x[delays+1:end])'

    # xₙ₊₁ = a₁x_{n-d+1} + a₂x_{n-d+2} + ... + a_d*xₙ
    xₙ₊₁ = xd[:, delays+1]
    xₙ = xd[:, 1:delays]
    a = xₙ\xₙ₊₁

    # Create predicted signal
    x̂ = zeros(length(x), 1)
    x̂[1:delays] = x[1:delays]

    for j = delays+1:length(x)
        for k = 1:delays
            x̂[j] += a[k]*x̂[j-delays+k-1]
        end
    end

    # Plot Results
    p1 = plot(t, x', label="True Signal", xlabel="t", ylabel="x", title="Synthetic Example", linestyle = :dash,  color=:blue,)
    plot!(t, x̂, label="Reconstructed signal", color=:red)
    display(p1)

## Nonlinear oscillation example
    # Simulate the ODE
    Δt = 0.05
    tspan = (0, 200)
    x₀ = [2, 2]
    prob = ODEProblem(VdP, x₀, tspan)
    sol = solve(prob, Tsit5(), dt = Δt, saveat= Δt)
    t = sol.t
    x = reduce(hcat, sol.u)'

    # Plotting solution
    l = @layout [a ; b]
    p2 = plot(t, x, layout = l, label = ["x" "ẋ"], xlabel = ["t" "t"], ylabel = ["x(t)" "ẋ(t)"], color = [:blue :red])
    display(p2)
    

## Time-Delay coordinate DMD on VdP
    # Hankel Matrix
    delays = 100
    xd = hankel(x[1:delays + 1, 1], x[delays+1:end, 1])'
    xₙ₊₁ = xd[:, delays+1]
    xₙ = xd[:, 1:delays]
    a = xₙ\xₙ₊₁

    # Create predicted signal
    x̂ = zeros(length(x[:,1]), 1)
    x̂[1:delays] = x[1:delays, 1]
    for j = delays+1:length(x[:,1])
        for k = 1:delays
            x̂[j] += a[k]*x̂[j-delays+k-1]
        end
    end

    # Plot Results
    p3 = plot(t, x[:,1], label="True Signal", xlabel="t", ylabel="x", title="Van der Pol Oscillator with time delay embedding", linestyle = :dash,  color=:blue,)
    plot!(t, x̂, label="Reconstructed signal", color=:red)
    display(p3)

## SVD on Hankel Matrix
    delays = 1000
    xd = hankel(x[1:delays, 1], x[delays:end, 1])

    # Apply SVD
    U, Σ, V = svd(xd, full = true)

    # Plot SVD results
    l = @layout [a ; b ; c]
    p4_1 = scatter(Σ.^2 ./sum(Σ.^2), markersize = 3, label = "Singular Values", xlabel = "Singular Value Index", ylabel = "Energy", title = "Singular Value Spectrum", color = :blue)
    p4_2 = plot(t[1:end-delays+1], V[:,1], color = :red, label = "First Right Singular Vector", xlabel = "t", ylabel = "v(t)", title = "Right Singular Vectors")
    plot!(t[1:end-delays+1], V[:,2], color = :blue, label = "Second Right Singular Vector", linestlye = :dash)
    p4_3 = plot(t[1:end-delays+1], V[:,3], color = :red, label = "Third Right Singular Vector", xlabel = "t", ylabel = "v(t)", title = "Right Singular Vectors")
    plot!(t[1:end-delays+1], V[:,4], color = :blue, label = "Fourth Right Singular Vector", linestlye = :dash)

    p4 = plot(p4_1, p4_2, p4_3, layout = l)
    display(p4)

## Time-delay DMD on Low Rank approximation of Hankel matrix
    # Find rank with 95% of Energy
    energy = 0
    energyTotal = sum(Σ.^2)
    r = 0
    while energy <= 0.95
        global r
        global energy
        r += 1
        energy += Σ[r]^2/energyTotal
    end

    # Apply DMD to first r columns of V
    X1 = V[1:end-1, 1:r]'
    X2 = V[2:end, 1:r]'
    U2, Σ2, V2 = svd(X1)
    A = X2*V2*diagm(1 ./ Σ2)

    # Plot eigenvalues of A
    μ = eigvals(A)
    p5 = scatter(real.(μ), imag.(μ), markersize = 8, color = :red, xlabel = "Re(μ)", ylabel = "Im(μ)", title = "Eigenvalues of A", label = "μ")
    display(p5)

