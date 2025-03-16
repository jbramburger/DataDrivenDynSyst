# -------------------------------------------------------------------------
# Plotting upper and lower bounds on heteroclinic existence
#
# This script plots the upper and lower bounds over a range of the
# hyperparameter lambda for the existence of heteroclinic orbit. Data is
# loaded in from the .mat files heteroclinic_lambda_ + 'upper' or 'lower' for
# upper or lower bounds + 'm2' or 'm3' for m = 2 or 3, respectively. Each
# file containts the lambda values, the value of m for reference, and
# arrays c2, c3, c4, c5 representing the values of upper or lower bounds
# on mu for degree 2,3,4, or 5 auxiliary functions, respectively.
#
# This script accompanies Section 4.2 of Data-Driven Methods for
# Dynamic Systems.
#
# This script was adapted from the MATLAB version by Jason J. Bramburger
# 
# Author: Mohit Sahu
# -------------------------------------------------------------------------


using MAT, Plots, Parameters

# Upper Bound Data
# m = 2 data
file_path = "DataDrivenDynSyst/data_driven_poly_optimiazation/"
vars = matread("$(file_path)heteroclinic_lambda_upper_m2.mat")
@unpack c2, c3, c4, c5, lambda, m = vars
vars2 = matread("$(file_path)heteroclinic_lambda_upper_m2_part2.mat")
@unpack c3_2, c4_2, c5_2, lambda_2, m = vars2

# m = 3 data
# vars = matread("$(file_path)heteroclinic_lambda_upper_m3.mat")
# @unpack c2, c3, c4, c5, lambda, m = vars
# vars2 = matread("$(file_path)heteroclinic_lambda_upper_m3_part2.mat")
# @unpack c2_2,c3_2, c4_2, c5_2, lambda_2, m = vars2

# Unpack data and aggregate
λ = [lambda lambda_2]
c3 = [c3 c3_2]
c4 = [c4 c4_2]
c5 = [c5[:,1:181] c5_2]

# Move failures off the figure
if m == 2
    c2 = [c2 2*ones(1,length(lambda_2))];
    c2[1:46] .= 2;
    c2[174:end] .= 2;
    c3[260:end] .= 2;
else
    c2 = [c2 c2_2];
end

# Plot Upper bounds
figure1 = plot()

plot!(λ[1:3:end], c2[1:3:end], c=RGB(1, 69/255, 79/255), lw=4, label="c2")
plot!(λ[1:3:end], c3[1:3:end], ls=:dash, c=RGB(36/255, 122/255, 254/255), lw=4, label="c3")
plot!(λ[1:3:end], c4[1:3:end], ls=:dot, c=RGB(0, 168/255, 0), lw=4, label="c4")
plot!(λ[1:3:end], c5[1:3:end], ls=:dashdot, c=RGB(255/255, 66/255, 161/255), lw=4, label="c5")

display(figure1)

# Load Lower bound data
if m == 2
    vars = matread("$(file_path)heteroclinic_lambda_lower_m2.mat")
    @unpack c3, c4, c5, lambda, m = vars
elseif m == 3
    vars = matread("$(file_path)heteroclinic_lambda_lower_m3.mat")
    @unpack c3, c4, c5, lambda, m = vars
end

# Plot Lower Bounds
figure2 = plot()
if m == 2
    plot!(lambda[7:37], c3[7:37], ls=:dash, c=RGB(36/255, 122/255, 254/255), lw=4, label="c3")
    plot!(lambda[7:37], c4[7:37], ls=:dashdot, c=RGB(0, 168/255, 0), lw=4, label="c4")
    plot!(lambda[7:2:39], c5[7:2:39], ls=:dot, c=RGB(255/255, 66/255, 161/255), lw=4, label="c5")
    xlims!(0, lambda[37])
    ylims!(0.5, 1)
    yticks!(0.5:0.1:1)
elseif m == 3
    plot!(lambda[7:end], c3[7:end], ls=:dash, c=RGB(36/255, 122/255, 254/255), lw=4, label="c3")
    plot!(lambda[7:end], c4[7:end], ls=:dashdot, c=RGB(0, 168/255, 0), lw=4, label="c4")
    plot!(lambda[7:2:end], c5[7:2:end], ls=:dot, c=RGB(255/255, 66/255, 161/255), lw=4, label="c5")
    xlims!(0, lambda[37])
    ylims!(0.30, 0.7)
    yticks!(0.3:0.1:0.765)
end

xlabel!("λ", fontsize=24, fontweight=:bold)
ylabel!("μ", fontsize=24, fontweight=:bold)

display(figure2)
