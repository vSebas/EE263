using LinearAlgebra
using Printf

include(joinpath(@__DIR__, "readclassjson.jl"))

function main()
    data = readclassjson(joinpath(@__DIR__, "sysid_data.json"))
    m = data["m"]
    n = data["n"]
    N = data["N"]
    X = data["X"]
    Y = data["Y"]

    Ahat = Y * X' * inv(X * X')
    residual = Ahat * X - Y

    println("m = $m, n = $n, N = $N")
    println("rank(X) = ", rank(X))
    println("Ahat = ")
    show(stdout, "text/plain", round.(Ahat; digits=4))
    println()
    @printf("RMS residual = %.6f\n", norm(residual) / sqrt(N))
    @printf("Frobenius residual = %.6f\n", norm(residual))
end

main()
