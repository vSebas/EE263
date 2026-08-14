using LinearAlgebra
using Random
using Printf

Random.seed!(263)

A = randn(3, 3)
M = (I + A) * inv(I - A)

direct = eigvals(M)
mapped = (1 .+ eigvals(A)) ./ (1 .- eigvals(A))

max_eigenvalue_error = max(
    maximum(minimum(abs.(direct .- z)) for z in mapped),
    maximum(minimum(abs.(mapped .- z)) for z in direct),
)

println("A = ")
show(stdout, "text/plain", round.(A; digits=6))
println()
println("eig((I + A)(I - A)^-1) = ")
show(stdout, "text/plain", round.(direct; digits=6))
println()
println("(1 + lambda(A))/(1 - lambda(A)) = ")
show(stdout, "text/plain", round.(mapped; digits=6))
println()
@printf("maximum eigenvalue-set discrepancy = %.3e\n", max_eigenvalue_error)
