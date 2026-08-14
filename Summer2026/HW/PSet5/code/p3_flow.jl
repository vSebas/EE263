using LinearAlgebra
using Printf

function main()
    A = [-1.0  1.0 -1.0 -1.0  0.0  0.0  0.0;
          0.0 -1.0  0.0  0.0 -1.0  0.0  0.0;
          0.0  0.0  0.0  1.0  1.0 -1.0  0.0;
          0.0  0.0  1.0  0.0  0.0  1.0 -1.0]

    s = [1.0, 4.0, 10.0, 10.0]
    fsimple = [5.0, 4.0, 0.0, 0.0, 0.0, 10.0, 20.0]

    fopt = -A' * inv(A * A') * s

    println("A*fopt = ")
    show(stdout, "text/plain", round.(A * fopt; digits=4))
    println()
    println("fopt = ")
    show(stdout, "text/plain", round.(fopt; digits=4))
    println()
    @printf("mean square optimal flow = %.6f\n", sum(fopt.^2) / length(fopt))
    @printf("mean square simple flow = %.6f\n", sum(fsimple.^2) / length(fsimple))
end

main()
