using Combinatorics

function monpowers(n::Int, max_power::Int)
    result = collect(multiexponents(n, 0))
    for s in 1:max_power
        result = vcat(result, collect(multiexponents(n, s)))
    end
    result
end

function createPoly(u, deg, start=0)
    a = monomials(u, deg)

    deg -=1
    while (deg>=start)
        a = vcat(a, monomials(u, deg))
        deg -=1
    end
    Poly(a)
end