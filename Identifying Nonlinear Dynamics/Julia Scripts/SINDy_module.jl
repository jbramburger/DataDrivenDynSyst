# -------------------------------------------------------------------------
# Sparse Identification of Nonlinear Dynamics (SINDy)
# -------------------------------------------------------------------------
#
# The function in this module is taken from SINDY.jl, for better Code
# structure
#
# Author: Mohit Sahu
# -------------------------------------------------------------------------

using LinearAlgebra
function sindy(dxdt, Θ, total_dim, lam=0.01)
    Xi = dxdt * pinv(Θ);

    k=1;
    Xi_new = Xi;

    while sum(abs.(Xi - Xi_new)) > 0  || k == 1
        
        Xi = deepcopy(Xi_new);

        # loop over all 3 variables
        for ind=1:total_dim
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

function print_model(Xi_new, mons,total_dim, print_f=["dx/dt", "dy/dt", "dz/dt"], lam=1e-3)
    println("Discovered Model using SINDy: ")
    for ind = 1:total_dim
        println(print_f[ind] * " = ")
        bigcoeffs = abs.(Xi_new[ind,:]) .> lam; # chosen small just to weed out zero coefficients
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
