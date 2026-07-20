using LinearAlgebra

function main()
    G = [1   0.2 0.1;
         0.1 2   0.1;
         0.3 0.1 3]

    n = size(G, 1)

    sigma = 0.01
    Pmax = 0.1

    D = Diagonal(G)
    K = G - D

    best_S = nothing
    best_p = nothing

    for S in 2.0:0.1:4.0
        p = (D - S * K) \ (S * sigma * ones(n))

        feasible = all((p .>= 0) .& (p .<= Pmax))

        if feasible
            best_S = S
            best_p = p
        end
    end

    println("Largest feasible target SINR = ", best_S)
    println("Power allocation = ", best_p)
end

main()
