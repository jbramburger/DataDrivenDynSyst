#Here we will model the system of differential equations:


using DifferentialEquations, LinearAlgebra, Plots

dt = 0.01;
tspan = (0.0,5.0);

function consODE!(dx, x,p, t)
    dx[1] = x[2] * x[3];
    dx[2] = -2 * x[1] * x[3];
    dx[3] = x[1] * x[2]
end

# Simulate conservative system 
numInit = 20; # number of initial conditions (More than one intial condition)

x0 = 10*rand(1,3) .- 5;
prob = ODEProblem(consODE!, x0, tspan);
integrator = init(prob, reltol = 1e-12, abstol = 1e-12);

traj_len = Int(((tspan[2]-tspan[1])/dt+1) * numInit)
dx_dt = zeros((traj_len, 3))
x = zeros((traj_len, 3));

c = 1;
for i in 1:numInit
    
    for i in tspan[1]:dt:tspan[2]
        step!(integrator, dt, true);
        x[c,:] .= integrator.u[:];
        dx_dt[c,:] .= get_du(integrator)[:];
        c+=1;
    end

    # Once again start integrating but with new initial condition
    x0 = 10*rand(1,3) .- 5;
    reinit!(integrator, x0)

end

γ = zeros(size(x)[1], 18);
γ[:,1:3] = dx_dt; # x1, x2, & x3 observables
γ[:,4] = 2*x[:,1].*dx_dt[:,1]; # x1^2 observable
γ[:,5] = x[:,2].*dx_dt[:,1] + x[:,1].*dx_dt[:,2]; # x1*x2 observable
γ[:,6] = x[:,3].*dx_dt[:,1] + x[:,1].*dx_dt[:,3]; # x1*x3 observable
γ[:,7] = 2*x[:,2].*dx_dt[:,2]; # x2^2 observable
γ[:,8] = x[:,3].*dx_dt[:,2] + x[:,2].*dx_dt[:,3]; # x2*x3 observable
γ[:,9] = 2*x[:,3].*dx_dt[:,3]; # x3^2 observable
γ[:,10] = 3*x[:,1].*x[:,1].*dx_dt[:,1]; # x1^3 observable
γ[:,11] = 2*x[:,1].*x[:,2].*dx_dt[:,1] + x[:,1].*x[:,1].*dx_dt[:,2]; # x1*x1*x2 observable
γ[:,12] = x[:,2].*x[:,2].*dx_dt[:,1] + 2*x[:,1].*x[:,2].*dx_dt[:,2]; # x1*x2*x2 observable
γ[:,13] = x[:,2].*x[:,3].*dx_dt[:,1] + x[:,1].*x[:,3].*dx_dt[:,2] + x[:,1].*x[:,2].*dx_dt[:,3]; # x1*x2*x3 observable
γ[:,14] = x[:,3].*x[:,3].*dx_dt[:,1] + 2*x[:,1].*x[:,3].*dx_dt[:,3]; # x1*x3*x3 observable
γ[:,15] = 3*x[:,2].*x[:,2].*dx_dt[:,2]; # x2^3 observable
γ[:,16] = 2*x[:,2].*x[:,3].*dx_dt[:,2] + x[:,2].*x[:,2].*dx_dt[:,3]; # x2*x2*x3 observable
γ[:,17] = x[:,3].*x[:,3].*dx_dt[:,2] + 2*x[:,2].*x[:,3].*dx_dt[:,3]; # x2*x3*x3 observable
γ[:,18] = 3*x[:,3].*x[:,3].*dx_dt[:,3]; # x3^3 observable

U, S, V = svd(γ);
scatter(S)

xi1 = [V[:,end-1] V[:,end]] \ [0;0;0;1;0;0;1;0;1; zeros(9,1)]/2; # angular momentum
xi2 = [V[:,end-1] V[:,end]] \ [0;0;0;1;0;0;2;0;3; zeros(9,1)]/2; # kinetic energy (Hamiltonian)

# Plot trajectories on unit sphere
function sphere(n,r=1.0, C=[0.0,0.0,0.0])   # r: radius; C: center [cx,cy,cz]
    u = range(-π, π; length = n)
    v = range(0, π; length = n)
    x = C[1] .+ r*cos.(u) * sin.(v)'
    y = C[2] .+ r*sin.(u) * sin.(v)'
    z = C[3] .+ r*ones(n) * cos.(v)'
    return x, y, z
end

# Generate initial conditions on the sphere
tspan = [0.0, 25.0];
x_init = [-1 -.5 0;
         -1 0.9 -0.9;
         0.1 -1 0.3;
         0.1 -1 -0.3;
         0.18 -0.9 -0.1;
         0.3 0.3 0.80]

x_init .= x_init ./ norm.(eachrow(x_init))

function generate_traj!(xUnit)
    for n in 1:size(x_init)[1]
        reinit!(integrator, reshape(x_init[n,:], length(x_init[n,:]),1))
        c=1;
        for i in tspan[1]:dt:tspan[2]
            step!(integrator, dt, true);
            xUnit[n,c,:] .= integrator.u[:];
            c+=1;
        end
    end
end

xUnit = zeros(size(x_init)[1], 5001, size(x_init)[2]);
generate_traj!(xUnit);

plotlyjs()
surface(sphere(100))
plot!(xUnit[1,:,1],xUnit[1,:,2],xUnit[1,:,3], lw=6, lc=:blueviolet)
plot!(xUnit[2,:,1],xUnit[2,:,2],xUnit[2,:,3], lw=6, lc=:chocolate)
plot!(xUnit[3,:,1],xUnit[3,:,2],xUnit[3,:,3], lw=6, lc=:cornflowerblue)
plot!(xUnit[4,:,1],xUnit[4,:,2],xUnit[4,:,3], lw=6, lc=:darkgrey)
plot!(xUnit[5,:,1],xUnit[5,:,2],xUnit[5,:,3], lw=6, lc=:darkred)
plot!(xUnit[6,:,1],xUnit[6,:,2],xUnit[6,:,3], lw=6, lc=:lightpink1)


# Plot trajectories on unit ellipsoid
function ellipsoid(cx, cy, cz, rx, ry, rz, n=100)
    u = LinRange(0, 2π, n)
    v = LinRange(0, π, n)

    x = [cx + rx * sin(v[j]) * cos(u[i]) for i in 1:n, j in 1:n]
    y = [cy + ry * sin(v[j]) * sin(u[i]) for i in 1:n, j in 1:n]
    z = [cz + rz * cos(v[j]) for i in 1:n, j in 1:n]

    return x, y, z
end

x_init = [0 -.2 0.1;
          0 -.27 0.48;
          0 -.63 0.01;
          0 -.14 0.55;
          0 -.66 -0.19;
          0 -.65 0.1]

x_init[:,1] .= @. -1*sqrt(1 - 2*x_init[:,2]^2 - 3*x_init[:,3]^2);
generate_traj!(xUnit);

surface(ellipsoid(0,0,0, 1,1/sqrt(2),1/sqrt(3)))
plot!(xUnit[1,:,1],xUnit[1,:,2],xUnit[1,:,3], lw=6, lc=:blueviolet)
plot!(xUnit[2,:,1],xUnit[2,:,2],xUnit[2,:,3], lw=6, lc=:chocolate)
plot!(xUnit[3,:,1],xUnit[3,:,2],xUnit[3,:,3], lw=6, lc=:cornflowerblue)
plot!(xUnit[4,:,1],xUnit[4,:,2],xUnit[4,:,3], lw=6, lc=:darkgrey)
plot!(xUnit[5,:,1],xUnit[5,:,2],xUnit[5,:,3], lw=6, lc=:darkred)
plot!(xUnit[6,:,1],xUnit[6,:,2],xUnit[6,:,3], lw=6, lc=:lightpink1)
