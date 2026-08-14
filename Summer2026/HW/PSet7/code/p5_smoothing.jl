using LinearAlgebra
ENV["GKSwstype"] = "100"
using Plots
using Printf

include(joinpath(@__DIR__, "readclassjson.jl"))

function convolution_matrix(c, n)
    h = div(length(c) - 1, 2)
    A = zeros(n, n)
    for i in 1:n, k in -h:h
        j = i + k
        if 1 <= j <= n
            A[i, j] = c[k + h + 1]
        end
    end
    return A
end

function truncated_estimate(F, y, r)
    Ur = F.U[:, 1:r]
    Vr = F.V[:, 1:r]
    return Vr * ((Ur' * y) ./ F.S[1:r])
end

function tikhonov_estimate(F, y, mu)
    factors = F.S ./ (F.S .^ 2 .+ mu)
    return F.V * (factors .* (F.U' * y))
end

function main()
    data = readclassjson(joinpath(@__DIR__, "regl_data.json"))
    n = data["n"]
    c = data["c"]
    x = vec(data["x"])
    w = vec(data["w"])

    A = convolution_matrix(c, n)
    ymeas = A * x + w
    F = svd(A)

    img_dir = normpath(joinpath(@__DIR__, "..", "latex", "img"))
    mkpath(img_dir)

    p_singular = plot(
        1:n,
        F.S;
        yscale=:log10,
        xlabel="k",
        ylabel="sigma_k",
        label=false,
        linewidth=2,
    )
    savefig(p_singular, joinpath(img_dir, "p5_singular_values.pdf"))

    p_vectors = plot(
        1:n,
        F.V[:, 1:6];
        layout=(3, 2),
        xlabel="i",
        ylabel="v_j(i)",
        label=["v$j" for _ in 1:1, j in 1:6],
        size=(800, 650),
    )
    savefig(p_vectors, joinpath(img_dir, "p5_right_vectors.pdf"))

    xls = A \ ymeas
    p_ls = plot(1:n, x; label="true x", linewidth=2)
    plot!(p_ls, 1:n, xls; label="least-squares estimate", linewidth=1.5)
    xlabel!(p_ls, "i")
    savefig(p_ls, joinpath(img_dir, "p5_least_squares.pdf"))

    ranks = [5, 10, 15, 30, 50]
    rank_panels = Any[]
    for r in ranks
        xr = truncated_estimate(F, ymeas, r)
        p = plot(1:n, x; label="true x", linewidth=2)
        plot!(p, 1:n, xr; label="estimate, r = $r", linewidth=1.5)
        push!(rank_panels, p)
    end
    p_ranks = plot(rank_panels...; layout=(3, 2), size=(900, 750))
    savefig(p_ranks, joinpath(img_dir, "p5_truncated_estimates.pdf"))

    tested_ranks = 1:35
    rank_errors = [
        norm(x - truncated_estimate(F, ymeas, r)) for r in tested_ranks
    ]
    best_index = argmin(rank_errors)
    best_r = tested_ranks[best_index]

    p_errors = plot(
        tested_ranks,
        rank_errors;
        xlabel="r",
        ylabel="norm(x - x_est)",
        label=false,
        linewidth=2,
        marker=:circle,
        markersize=2,
    )
    scatter!(p_errors, [best_r], [rank_errors[best_index]];
             label="best r = $best_r")
    savefig(p_errors, joinpath(img_dir, "p5_truncation_error.pdf"))

    xbest = truncated_estimate(F, ymeas, best_r)
    p_best = plot(1:n, x; label="true x", linewidth=2)
    plot!(p_best, 1:n, xbest;
          label="truncated SVD, r = $best_r", linewidth=2)
    xlabel!(p_best, "i")
    savefig(p_best, joinpath(img_dir, "p5_best_truncated.pdf"))

    mus = 10.0 .^ range(-10, 0; length=301)
    regularization_errors = [
        norm(x - tikhonov_estimate(F, ymeas, mu)) for mu in mus
    ]
    best_mu_index = argmin(regularization_errors)
    best_mu = mus[best_mu_index]
    xreg = tikhonov_estimate(F, ymeas, best_mu)

    p_reg = plot(1:n, x; label="true x", linewidth=2)
    plot!(p_reg, 1:n, xreg;
          label="Tikhonov estimate, mu = $(round(best_mu; sigdigits=4))",
          linewidth=2)
    xlabel!(p_reg, "i")
    savefig(p_reg, joinpath(img_dir, "p5_tikhonov.pdf"))

    gains = F.S ./ (F.S .^ 2 .+ best_mu)
    inverse_error_factors = best_mu ./ (F.S .^ 2 .+ best_mu)

    println("n = $n, h = $(div(length(c) - 1, 2))")
    @printf("sigma_max(A) = %.6e\n", F.S[1])
    @printf("sigma_min(A) = %.6e\n", F.S[end])
    @printf("condition number = %.6e\n", F.S[1] / F.S[end])
    @printf("least-squares estimation error = %.6e\n", norm(x - xls))
    for r in ranks
        @printf("truncated error for r = %2d: %.6f\n", r,
                norm(x - truncated_estimate(F, ymeas, r)))
    end
    @printf("best rank in 1:35 = %d, error = %.6f\n",
            best_r, rank_errors[best_index])
    @printf("selected mu = %.6e, error = %.6f\n",
            best_mu, regularization_errors[best_mu_index])
    @printf("norm(B) = %.6e\n", maximum(gains))
    @printf("worst-case relative inversion error = %.6e\n",
            maximum(inverse_error_factors))
end

main()
