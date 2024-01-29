# -------------------------------------------------------------------------
# Dynamic Mode Decomposition (DMD)
#
# This code generates synthetic snapshot data by numerically integrating
# a time-periodic soliton solution of the nonlinear Schrodinger equation
# (NLS). We then applying DMD to the gathered snapshots and view the
# corresponding DMD eigenvalues to understand the long-time dynamics.
# The standard DMD results are compared to the physics-informed DMD (piDMD)
# approach which seeks identify DMD matrix that is unitary by solving a
# Procrustes problem.
#
# This script accompanies Section 2.1 of Data-Driven Methods for
# Dynamic Systems.
#
# This code was translated from the MATLAB version written by Jason J. Bramburger into Julia.
#
# Author: Daniel Fassler
# -------------------------------------------------------------------------


using FFTW, DifferentialEquations, LinearAlgebra, CairoMakie
## Functions
    function f(û,p,t)
        k = p[1]
        u = ifft(û)
        out = -0.5im * k.^2 .* û .+ 1im * fft(abs.(u).^2 .* u)
        return out
    end

    function meshgrid(x,y)
        nx=length(x)
        ny=length(y)
        xout=zeros(ny,nx)
        yout=zeros(ny,nx)
        for j=1:nx
            for i=1:ny
                xout[i,j]=x[j]
                yout[i,j]=y[i]
            end
        end
        return xout, yout
    end

## Generate Data from NLS using spectral methods
    # Space
    L = 10
    M = 64 # Size of snapshots
    Y = LinRange(-L/2, L/2, M+1)
    x = Y[1:M]
    k = (2π/L)*[Int.(LinRange(0, div(M, 2)-1, div(M, 2))) ; Int.(round.(LinRange(-div(M, 2), -1, div(M,2))))]

    # Time
    slices = 20 # Number of snapshots
    tspan = LinRange(0, 2π, slices+1)
    Δt = tspan[2]-tspan[1]

    # Initial Condition
    u = 2*sech.(x)
    û = fft(u)

    # Solve NLS
    prob = ODEProblem(f, û, (0.0, 2π), [k])
    sol = solve(prob, Tsit5())
    t = sol.t
    usol = zeros(ComplexF64, length(sol.u), length(sol.u[1])) 

    for i = 1:length(sol.u)
        usol[i,:] = ifft(sol.u[i,:][1]) # Back to real space
    end

    # Plot solution
    Z, T = meshgrid(x, t)
    fig = Figure()
    ax = Axis(fig[1,1])
    hm = heatmap!(Z, T, abs.(usol)', colormap = :viridis)
    Colorbar(fig[1,2], hm)
    display(fig)

## Create X and Y data matrices
    # Sample the solution at equidistant time
    tsample = zeros(slices+1)
    tsample_indices = zeros(Int64, slices+1)
    for i = 1:slices+1
        tsample[i] = t[div(length(t), slices+1)*i]
        tsample_indices[i] = div(length(t), slices+1)*i
    end
    usol_sample = usol[tsample_indices,:]
    X = usol_sample[1:end-1,:]' # X data matrix
    Y = usol_sample[2:end,:]' # Y data matrix

## Compute the DMD Matrix Y = AX by A = YX†
    A = Y * pinv(X)
    evectors = eigvecs(A)
    evalues = eigvals(A)
    omega = log.(evalues)/Δt # Continuoues time eigenvalues for reference

## Plotting the eigenvalues
    th = LinRange(0, 2π, 100)
    xcos = cos.(th)
    ysin = sin.(th)
    
    fig2 = Figure()
    ax2 = Axis(fig2[1,1], xlabel = "Re(μ)", ylabel = "Im(μ)")
    evalscat = scatter!(ax2, real.(evalues), imag.(evalues), markersize = 8, color = :red)
    unitdisclines = lines!(ax2, xcos, ysin, color = :black, linestyle = :dash)
    Legend(fig2[1,2], [evalscat, unitdisclines], ["Eigenvalues", "Unit Circle"])
    display(fig2)
## Stepping error for DMD Prediction
    Z = zeros(ComplexF64, M, 201)
    Z[:, 1] = X[:, 1]

    # Forecast using DMD matrix
    for n = 1:200
        Z[:, n+1] = A * Z[:, n]
    end

    # Compare {Y}ₘ to {Z}ₘ₊₁, should be the same if A is perfect
    println("Error of DMD at t = 2π : ", norm(Y[:, 20] - Z[:, 21]))

    # Look at the long term forecast for Z, it is not stable
    println("Norm of DMD Solution after long forecast : ", norm(Z[:, 201]))

## Enforce DMD is unitary
    # Solution to procrustes problem
    M = X*Y'
    U₁, Σ₁, V₁ = svd(M)
    Aᵤ = U₁*V₁' # Unitary DMD matrix

    # Compute and plot the eigenvalues
    evectorsᵤ = eigvecs(Aᵤ)
    evaluesᵤ = eigvals(Aᵤ)

    fig3 = Figure()
    ax3 = Axis(fig3[1,1], xlabel = "Re(μ)", ylabel = "Im(μ)")
    evalscatᵤ = scatter!(ax3, real.(evaluesᵤ), imag.(evaluesᵤ), markersize = 8, color = :red)
    unitdisclinesᵤ = lines!(ax3, xcos, ysin, color = :black, linestyle = :dash)
    Legend(fig3[1,2], [evalscatᵤ, unitdisclinesᵤ], ["Eigenvalues", "Unit Circle"])
    display(fig3)

## Compare Forecasts far into the future
    t_forecast = 0:Δt:20π

    # Regular DMD prediction
    u_dmd = zeros(ComplexF64, length(x), length(t_forecast))
    u_dmd[:, 1] = X[:, 1]
    for n = 2:length(t_forecast)
        u_dmd[:, n] = A * u_dmd[:, n-1]
    end

    # Unitary DMD prediction
    u_dmdᵤ = zeros(ComplexF64, length(x), length(t_forecast))
    u_dmdᵤ[:, 1] = X[:, 1]
    for n = 2:length(t_forecast)
        u_dmdᵤ[:, n] = Aᵤ * u_dmdᵤ[:, n-1]
    end

    # Simulating the RHS far into the future
    problong = ODEProblem(f, û, (0.0, 20π), [k])
    sollong = solve(problong, Tsit5())
    tlong = sollong.t
    usollong = zeros(ComplexF64, length(sollong.u), length(sollong.u[1]))
    for i = 1:length(sollong.u)
        usollong[i,:] = ifft(sollong.u[i,:][1])
    end

    X, Tsim = meshgrid(x, tlong)    
    X, Tdmd = meshgrid(x, t_forecast)
    fig4 = Figure()
    ax4_sim = Axis(fig4[1,1], xlabel = "x", ylabel = "t")
    ax4_dmdu = Axis(fig4[2,1], xlabel = "x", ylabel = "t")
    hm4_sim = heatmap!(ax4_sim, X, Tsim, abs.(usollong)', colormap = :viridis)
    hm4_dmdu = heatmap!(ax4_dmdu, X, Tdmd, abs.(u_dmdᵤ), colormap = :viridis)
    Colorbar(fig4[1,2], limits = (0, 4))
    Colorbar(fig4[2,2], limits = (0, 4))
    display(fig4)

