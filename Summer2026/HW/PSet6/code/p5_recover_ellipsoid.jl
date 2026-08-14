using LinearAlgebra
using Statistics
using Printf

include(joinpath(@__DIR__, "readclassjson.jl"))

function symmetric_design(X)
    n, N = size(X)
    pairs = Tuple{Int, Int}[]
    Phi = zeros(N, div(n * (n + 1), 2))

    col = 1
    for i in 1:n, j in i:n
        push!(pairs, (i, j))
        if i == j
            Phi[:, col] = X[i, :] .^ 2
        else
            Phi[:, col] = 2 .* X[i, :] .* X[j, :]
        end
        col += 1
    end

    return Phi, pairs
end

function symmetric_matrix(theta, pairs, n)
    A = zeros(n, n)
    for (col, (i, j)) in enumerate(pairs)
        A[i, j] = theta[col]
        A[j, i] = theta[col]
    end
    return Symmetric(A)
end

function main()
    data = readclassjson(joinpath(@__DIR__, "ellip_bdry_data.json"))
    X = data["X"]
    n, N = size(X)

    Phi, pairs = symmetric_design(X)
    theta = Phi \ ones(N)
    A = symmetric_matrix(theta, pairs, n)

    values = [dot(X[:, k], A * X[:, k]) for k in 1:N]
    errors = values .- 1
    evals = eigvals(A)

    println("n = $n, N = $N")
    println("A = ")
    show(stdout, "text/plain", round.(Matrix(A); digits=6))
    println()
    println("eigenvalues(A) = ")
    show(stdout, "text/plain", round.(evals; digits=6))
    println()
    @printf("min x'Ax = %.6f\n", minimum(values))
    @printf("max x'Ax = %.6f\n", maximum(values))
    @printf("mean x'Ax = %.6f\n", mean(values))
    @printf("RMS boundary error = %.6f\n", norm(errors) / sqrt(N))
    @printf("max abs boundary error = %.6f\n", maximum(abs.(errors)))
end

main()
