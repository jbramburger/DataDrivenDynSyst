# -------------------------------------------------------------------------
# Extended Dynamic Mode Decomposition (EDMD)
#
# This code applies the Extended DMD (EDMD) approach to a simple nonlinear
# discrete dynamical system for which the Koopman eigenfunctions are known.
# Below we compare results using a dictionary of monomials up to degree 2
# and up to degree 3. Training data is obtained by drawing N (number of
# snapshots) pairs of points from the uniform distribution on [-1,1]^2 and
# tracking their image under the dynamical system. The system depends on
# two parameters, λ₁ and λ₂, which can be varied to produce similar
# results.
#
# This script accompanies Section 2.4 of Data-Driven Methods for
# Dynamic Systems. 
#
# This code was adapted from the MATLAB version by Jason J. Bramburger
#
#
# Author: Daniel Fassler
# -------------------------------------------------------------------------

using LinearAlgebra, Random, Plots, DifferentialEquations, GenericSchur

## Generate snapshot data

    N = 10000
    X = 2*rand(2, N) .- 1
    Y = zeros(2, N)

    # System parameters
    λ₁ = 1.5
    λ₂ = 0.5

    # Populate Y matrix
    for n = 1:N
        Y[1,n] = λ₁*X[1,n]
        Y[2,n] = λ₂*(X[2,n] .- X[1,n].^2)
    end

    ## Create matrices of ovservables Ψ(X) and Ψ(Y)

    Ψ₁X = zeros(5, N)
    Ψ₁Y = zeros(5, N)
    # Monomials up to degree 2
    # Degree 1 monomials (x₁ and x₂ are X and Y)
    Ψ₁X[1:2, :] = X
    Ψ₁Y[1:2, :] = Y

    # x₁² observable
    Ψ₁X[3, :] = X[1,:].^2
    Ψ₁Y[3, :] = Y[1,:].^2

    # x₂² observable
    Ψ₁X[4, :] = X[2,:].^2
    Ψ₁Y[4, :] = Y[2,:].^2

    # x₁x₂ observable
    Ψ₁X[5, :] = X[1,:].*X[2,:]
    Ψ₁Y[5, :] = Y[1,:].*Y[2,:]

## DMD matrix
    A = Ψ₁Y*pinv(Ψ₁X)
    μ = eigvals(A)
    V = eigvecs(A) # Right eigenvectors are columns of V
    W = eigvecs(A') # Left eigenvectors are rows of W

## Koopman eigenfunctions
    # Print the eigenvalues and associated Koopman eigenfunctions
    for p = 1:5
        println("\nEigenvalue: ", μ[p])
        println("Associated (approximate) Koopman eigenfunction")
        println("$(round(W[1, p], digits = 5))x₁ + $(round(W[2, p], digits = 5))x₂ + $(round(W[3, p], digits = 5))x₁² + $(round(W[4, p], digits = 5))x₂² + $(round(W[5, p], digits = 5))x₁x₂")
    end
    println("\n\n\n\n\n\n")
## Add in cubic observables
    Ψ₂X = zeros(9, N)
    Ψ₂Y = zeros(9, N)

    Ψ₂X[1:5, :] = Ψ₁X
    Ψ₂Y[1:5, :] = Ψ₁Y

    # Monomials up to degree 3
    # x₁³ observable
    Ψ₂X[6, :] = X[1,:].^3
    Ψ₂Y[6, :] = Y[1,:].^3

    # x₂³ observable
    Ψ₂X[7, :] = X[2,:].^3
    Ψ₂Y[7, :] = Y[2,:].^3

    # x₁²x₂ observable
    Ψ₂X[8, :] = X[1,:].^2 .* X[2,:]
    Ψ₂Y[8, :] = Y[1,:].^2 .* Y[2,:]

    # x₁x₂² observable
    Ψ₂X[9, :] = X[1,:] .* X[2,:].^2
    Ψ₂Y[9, :] = Y[1,:] .* Y[2,:].^2

## DMD matrix
    A₂ = Ψ₂Y*pinv(Ψ₂X)
    μ₂ = eigvals(A₂)
    V₂ = eigvecs(A₂) # Right eigenvectors are columns of V
    W₂ = eigvecs(A₂') # Left eigenvectors are rows of W

## Koopman eigenfunctions
    # Print the eigenvalues and associated Koopman eigenfunctions
    for p = 1:9
        println("\nEigenvalue: ", μ₂[p])
        println("Associated (approximate) Koopman eigenfunction")
        println("$(round(W₂[1, p], digits = 5))x₁ + $(round(W₂[2, p], digits = 5))x₂ + $(round(W₂[3, p], digits = 5))x₁² + $(round(W₂[4, p], digits = 5))x₂² + $(round(W₂[5, p], digits = 5))x₁x₂ + $(round(W₂[6, p], digits = 5))x₁³ + $(round(W₂[7, p], digits = 5))x₂³ + $(round(W₂[8, p], digits = 5))x₁²x₂ + $(round(W₂[9, p], digits = 5))x₁x₂²")
    end



